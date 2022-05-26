      program stable
!
!     This program perfoms a stability test on an SCF wavefunction.
!
!     L. M. Thompson, 2022
!
      use mqc_gaussian
      use iso_fortran_env, only: int32, int64, real64
!
!****x* Main/Stability
!*    NAME
!*      SCF Stability Calculator
!*
!*    SYNOPSIS
!*      Computes the stability of an SCF wavefunction.
!
      implicit none
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      character(len=:),allocatable::command,fileName,help_path,wf_string,symInstStr
      character(len=1)::wf_type
      character(len=256)::vecString
      integer(kind=int64)::iOut=6,iPrint=1,iUnit,flag,i,j,nAlpha,nBeta,nBasis,ovDim,oRHessDim,&
        occ1,occ2,virt1,virt2,ind,ind2,neigs2print=5
      real(kind=real64)::vecThresh=0.1
      real(kind=real64),parameter::thresh=1.0e-14
      logical::found,wf_complex
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs),dimension(:),allocatable::eris
      type(mqc_scalar)::Vnn
      type(mqc_scf_integral)::core_hamiltonian,mo_coefficients,fock,density,GMat
      type(mqc_matrix)::OrbRotHess,oRVecs
      type(mqc_vector)::oREigs,vec2process

      TYPE(mqc_scf_eigenvalues)::tmpeigs
      TYPE(mqc_scf_integral)::tmpvecs,lmtmocoreh
!
!*    USAGE
!*      stable [-f <matrix_file>] [--print-level <print_level>] [--help]
!*
!*    OPTIONS
!* 
!
!     Print program information.
!
      write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
          '  SCF Stability Calculator',NEW_LINE('a'), &
          ' ',repeat('*',73),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 22.05.1',NEW_LINE('a'),NEW_LINE('a'),&
          ' L. M. Thompson, Louisville KY, 2022.',NEW_LINE('a')
!
!     Parse input options.
!
!*   1. Input/output
!*
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f') then
!
!*      -f matrix_file                   Input matrix file with set of molecular orbitals and two-
!*                                       electron integrals. 
!*
          call mqc_get_command_argument(i+1,fileName)
          j = i+2
        elseif(command.eq.'--print-level') then
!
!*      --print-level print_level        Verbosity of output. Default print level is 1. Options
!*                                       0-4.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I1)') iPrint
          j = i + 2
        elseif(command.eq.'--neigs') then
!
!*      --neigs neigs                    Number of orbital rotation Hessian eigenvalues to print.
!*                                       Default is 5.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I5)') neigs2print
          j = i + 2
        elseif(command.eq.'--vecThresh') then
!
!*      --vecThresh vecThresh            Lowest limit for orbital rotation Hessian eigenvalue 
!*                                       component to be printed. Default is 0.1.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F10.3)') vecThresh
          j = i + 2
        elseIf(command.eq.'--wf-test') then
!
!*      --wf-test test_string            Specifies the type of instability to be tested. Note that
!*                                       only symmetries of the same type or lower than the input
!*                                       wavefunction symmetry are valid. The default is to test 
!*                                       for an internal instability. Options are: 
!*                                       1) R
!*                                          Test for real instabilities.
!*                                       2) C 
!*                                          Test for complex instabilities.
!*
          call mqc_get_command_argument(i+1,command)
          wf_string = trim(command)
          j = i+2
        elseIf(command.eq.'--help') then
!
!*      --help                           Output help documentation to terminal.
!*
          if(command_argument_count().gt.1) call mqc_error_I('Help output requested with multiple arguments',6, &
            'command_argument_count()',command_argument_count())
          call mqc_get_command_argument(0,help_path)
          help_path = 'less ' // trim(help_path(1:scan(help_path,'/',.true.))) // 'doc/stable.txt'
          call execute_command_line(help_path,exitstat=flag)
          if(flag.ne.0) call mqc_error('Help output command failed')
          stop
        else
          call mqc_error_A('Unrecognised input flag',6,'command',command)
        endIf
      endDo
!
!     Parse input file and extract required data from matrix files.
!
      allocate(eris(1))
      call fileInfo%load(fileName)
      call fileInfo%getMolData(moleculeInfo)
      Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
      call moleculeInfo%print(iOut)
      call Vnn%print(iOut,'Nuclear Repulsion Energy (au)')

      nBasis = fileInfo%getVal('nBasis')
      nAlpha = fileInfo%getVal('nAlpha')
      nBeta = fileInfo%getVal('nBeta')

      call fileInfo%getESTObj('mo coefficients',est_integral=mo_coefficients)
      if(iPrint.ge.2) call mo_coefficients%print(iOut,'MO coefficients')
      call fileInfo%getESTObj('fock',est_integral=fock,foundObj=found)
      if(.not.found) then
        call fileInfo%getESTObj('core hamiltonian',est_integral=core_hamiltonian)
        if(iPrint.ge.2) call core_hamiltonian%print(iOut,'Core Hamiltonian')
        call fileInfo%getESTObj('density',est_integral=density,foundObj=found)
        if(.not.found) then
          density = matmul(mo_coefficients%orbitals('occupied',[nAlpha],[nBeta]), &
            dagger(mo_coefficients%orbitals('occupied',[nAlpha],[nBeta])))
        endIf
        if(iPrint.ge.2) call density%print(iOut,'Density matrix')
        call fileInfo%get2ERIs('regular',eris(1))
        Gmat = contraction(eris,density)
        if(iPrint.ge.3) call Gmat%print(iOut,'G(P)')
        fock = core_hamiltonian + Gmat
      else
        call fileInfo%get2ERIs('regular',eris(1),foundERI=found)
        if(.not.found) call fileInfo%get2ERIs('molecular',eris(1))
      endIf
      if(iPrint.ge.2) call fock%print(iOut,'Fock matrix')
      if(iPrint.ge.4) call eris(1)%print(iOut,'AO 2ERIs') 
      if(eris(1)%type().eq.'regular')  call twoERI_trans(iOut,iPrint,mo_coefficients,eris(1),eris(1))
      if(iPrint.ge.4) call eris(1)%print(iOut,'MO 2ERIs') 
!
!     Establish the symmetry of the input wavefunction and the instability to be tested
!
      call MQC_Gaussian_ICGU(fileInfo%icgu,wf_type,wf_complex)
      if(wf_type.eq.'G'.and.wf_complex) then
        if(mqc_matrix_norm(aimag(fock%getBlock())).lt.thresh) then
          wf_complex = .false.
        endIf
      endIf
      if(.not.allocated(wf_string)) wf_string = ''
      if (len(wf_string).eq.0) then
        if(wf_complex) then
          wf_string = 'C'
        else
          wf_string = 'R'
        endIf
      else
        if((wf_string(1:1).ne.'C'.and.wf_string(1:1).ne.'R').or.&
          len(wf_string).ne.1) &
          call mqc_error_A(' Instability test input format incorrect',6,'wf_string',wf_string)
      endIf
      write(iOut,'(A)') ' Testing for instabilities leading to '//wf_string//' wavefunctions'
      if(wf_complex.and.wf_string(1:1).eq.'R') &
        call mqc_error_A(' Requesting real instability test but wavefunction is already complex', &
        6,'wf_string',wf_string)
!
!     Compute MO Fock matrix (orbital rotation gradient)
!
      fock = matmul(matmul(dagger(mo_coefficients),fock),mo_coefficients)
      if(iPrint.ge.2) call fock%print(iOut,'Orbital Rotation Gradient')

      ovDim = 2*nBasis*(nAlpha+nBeta)-2*nAlpha*nBeta-nAlpha**2-nBeta**2 
!
      if(wf_complex) then
        oRHessDim = 2*ovDim
      else
        oRHessDim = ovDim
      endIf 
      call orbRotHess%init(oRHessDim,oRHessDim)
!
      do i = 1,ovDim 
        do j = 1,i
!
!         Extract orbital indices (related to the index in the MO Fock matrix)
!
          occ1 = mod(i-1,nAlpha)+1
          occ2 = mod(j-1,nAlpha)+1
          virt1 = nAlpha+mod((i-1)/nAlpha,nBasis-nAlpha)+1
          virt2 = nAlpha+mod((j-1)/nAlpha,nBasis-nAlpha)+1
          if((i-1)/(2*nAlpha*(nBasis-nAlpha))+1.eq.2) occ1 = occ1+nBasis
          if((j-1)/(2*nAlpha*(nBasis-nAlpha))+1.eq.2) occ2 = occ2+nBasis
          if(mod((i-1)/(nAlpha*(nBasis-nAlpha)),2)+1.eq.2) virt1 = virt1+nBasis
          if(mod((j-1)/(nAlpha*(nBasis-nAlpha)),2)+1.eq.2) virt2 = virt2+nBasis
!
!         Build A block of orbital rotation Hessian
!
          if(occ1.eq.occ2) call orbRotHess%put(orbRotHess%at(i,j)+fock%at(virt1,virt2),i,j)
          if(virt1.eq.virt2) call orbRotHess%put(orbRotHess%at(i,j)-fock%at(occ1,occ2),i,j)
          call orbRotHess%put(orbRotHess%at(i,j)+eris(1)%at(virt1,occ1,occ2,virt2,interactionIn='doublebar'),i,j)
!
!         Build B block of orbital rotation Hessian 
!
          if(wf_complex) then
            call orbRotHess%put(orbRotHess%at(i,j+ovDim) + &
              eris(1)%at(virt1,occ1,virt2,occ2,interactionIn='doublebar'),i,j+ovDim)
              call orbRotHess%put(orbRotHess%at(i,j+ovDim),j,i+ovDim)
          else
            if(wf_string(1:1).eq.'R') then
              call orbRotHess%put(orbRotHess%at(i,j)+eris(1)%at(virt1,occ1,virt2,occ2,interactionIn='doublebar'),i,j)
            elseIf(wf_string(1:1).eq.'C') then
              call orbRotHess%put(orbRotHess%at(i,j)-eris(1)%at(virt1,occ1,virt2,occ2,interactionIn='doublebar'),i,j)
            endIf
          endIf
          call orbRotHess%put(conjg(orbRotHess%at(i,j)),j,i)
        endDo
      endDo
      if(wf_complex) then
        call orbRotHess%mput(conjg(orbRotHess%mat([1,ovDim],[1,ovDim])),[ovDim+1,-1],[ovDim+1,-1])
        call orbRotHess%mput(conjg(orbRotHess%mat([1,ovDim],[ovDim+1,-1])),[ovDim+1,-1],[1,ovDim])
      endIf
      if(iPrint.ge.3) call orbRotHess%print(6,'Orbital rotation Hessian')
      call orbRotHess%diag(oREigs,oRVecs)
      if(iPrint.ge.2) call oREigs%print(6,'Orbital rotation eigenvalues')
      if(iPrint.ge.3) call oRVecs%print(6,'Orbital rotation eigenvectors')

      
      write(6,'(1x,A)') NEW_LINE('A')//' Lowest '//trim(num2char(neigs2print))//&
        ' orbital rotation Hessian eigenvalues'//NEW_LINE('A')
      do i = 1, min(neigs2print,size(oREigs))
        call mqc_print(oREigs%at(i),6,'Eigenvector '//trim(num2char(i)))
        vec2process = oRVecs%vat([0],[i])
        do while(.true.)
          ind = maxLoc(abs(vec2process))
          if(abs(vec2process%at(ind)).lt.vecThresh) then
            write(6,'(1x,A)') NEW_LINE('A')
            exit
          endIf
          ind2 = mod(ind-1,ovDim)+1
          vecString = '     '//trim(num2char(mod(ind2-1,nAlpha)+1))
          if((ind2-1)/(2*nAlpha*(nBasis-nAlpha))+1.eq.2) then
            vecString = trim(vecString)//'b'
          else
            vecString = trim(vecString)//'a'
          endIf
          vecString = trim(vecString)//' --->'
          vecString = trim(vecString)//' '//trim(num2char(nAlpha+mod((ind2-1)/nAlpha,nBasis-nAlpha)+1))
          if(mod((ind2-1)/(nAlpha*(nBasis-nAlpha)),2)+1.eq.2) then
            vecString = trim(vecString)//'b'
          else
            vecString = trim(vecString)//'a'
          endIf
          if(.not.wf_complex) then
            vecString = trim(vecString)//'     '//trim(num2char(real(vec2process%at(ind))))
            if(wf_string.eq.'C') vecString = trim(vecString)//'i'
          else
            vecString = trim(vecString)//'     '//trim(num2char(vec2process%at(ind)))
          endIf
          write(6,'(A)') trim(vecString)
          call vec2process%put(0.0,ind)
        endDo
      endDo
!    
!
!*    NOTES
!*      Compilation of this program requires the MQC library (https://github.com/MQCPack/mqcPack)
!*      and the gauopen utility (http://gaussian.com/g16/gauopen.zip) and compilation with the
!*      f08 standard.
!*
!*      Compilation tested using: gfortran 9.2.0
!*
!*      Note that subroutine Wr_LCBuf needs modifying in gauopen/qcmatrix.F as follows:
!*        line 58:       LenBX = (LenBuf/(2*abs(NR)))*abs(NR)
!*        line 60:       Call Wr_CBuf(IU,NTot*abs(NR),LenBX,X)
!*
!*      Documentation generated with robodoc. To update documentation edit robodoc.rc to
!*      determine documentation output type and then run robodoc at the command line in the
!*      main directory.
!*
!*    AUTHORS
!*      Lee M. Thompson, University of Louisville, lee.thompson.1@lousiville.edu
!*
!*    COPYRIGHT
!*      (c) 2022 by Lee M. Thompson distributed under terms of the MIT license.
!*
!****
!
 999  End Program stable
