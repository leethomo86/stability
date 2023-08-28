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
      character(len=:),allocatable::command,fileName,help_path,wf_string,outputfileIn,outputFile
      character(len=1)::wf_type
      character(len=256)::vecString,otype='chk',file_tmp,coordinate=''
      integer(kind=int64)::iOut=6,iPrint=1,iUnit,flag,i,j,k,l,nAlpha,nBeta,nBasis,ovDim,oRHessDim,&
        occ1,occ2,virt1,virt2,ind,ind2,neigs2print=5,elem1,elem2,maxIters=5000,degen_start, &
      degen_end,ivec,jvec,iter,maxSteps=1,vpos,jEnd
      real(kind=real64)::vecThresh=0.1,initStep=0.05,step,connectThresh=5.0
      real(kind=real64),parameter::thresh=1.0e-10,etaThresh=1.0e-4,zero=1.0e-12,convThresh=1.0e-8, &
        followThresh=1.0e-3
      logical::found,wf_complex,doUHF,doGHF,doComplex,file_exists
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs),dimension(:),allocatable::AOeris,eris
      type(mqc_scalar)::Vnn,phi,theta,energy,vval,optangle
      type(mqc_scf_integral)::core_hamiltonian,mo_coefficients,fock,AOfock,density,GMat,overlap
      type(mqc_matrix)::OrbRotHess,oRVecs,rotation_matrix,sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,&
        conCoTwo,shCoor,initialORVecs
      type(mqc_vector)::oREigs,vec2process,vecfollow,linearFock,cumvec,vec_overlap
      type(mqc_scf_eigenvalues)::mo_energies
      logical::degen_flag,newtonFlag
!
!*    USAGE
!*      stable [-f <matrix_file>] [--print-level <print_level>] [--neigs <neigs>] [--vecThresh <vecThresh>]
!*        [--otype <extension>] [-o <output_file>] [--wf-test <wf_string>] [--follow <vector>] 
!*        [--coordinate <coordinate>] [--step <step>] [--nsteps <nSteps>] [--connect-angle <angle>] [--help]
!*
!*    OPTIONS
!* 
!
!     Print program information.
!
      write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
          '  SCF Stability Calculator',NEW_LINE('a'), &
          ' ',repeat('*',73),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 23.08.1',NEW_LINE('a'),NEW_LINE('a'),&
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
        elseIf(command.eq.'--otype') then
!
!*      --otype extension                Format of output file. Options are:
!*                                       1) chk (default)
!*                                       2) mat
!*
          call mqc_get_command_argument(i+1,command)
          otype = command
          j = i+2
        elseIf(command.eq.'-o') then
!
!*      -o output_file                   Output file name. Default is input file name.
!*
          call mqc_get_command_argument(i+1,outputFileIn)
          j = i+2
!*   2. Stability test
!*
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
!*                                       3) A 
!*                                          Test for all instabilities.
!*
          call mqc_get_command_argument(i+1,command)
          wf_string = trim(command)
          if (.not.any(['R','C','A'].eq.wf_string)) call mqc_error_a('Unrecognized symmetry input',iOut,'wf_string',wf_string)
          j = i+2
!*   3. Vector following
!*
        elseif(command.eq.'--follow') then
!
!*      --follow vector                  Specifies the vectors to follow when returning a matrix
!*                                       file with a perturbed set of molecular orbitals. Input 
!*                                       vectors are given as a column separated list enclosed in 
!*                                       square brackets. The default is not to return a matrix file.
!
!*
          call mqc_get_command_argument(i+1,command)
          call eigenfollow(command,vecFollow)
          j = i + 2
        elseIf(command.eq.'--coordinate') then
!
!*      --coordinate coordinate          Specifies which coordinate type to follow
!*                                       1) eigenvector (default)
!*                                          Step in the direction of the raw Hessian eigenvectors.
!*                                       2) gradient 
!*                                          Step in the direction of the gradient.
!*                                       3) newton
!*                                          Step in the direction of the gradient projected onto
!*                                          the Hessian eigenvectors.
!*                                       4) connect 
!*                                          Step in the direction of the solution connected by
!*                                          the specified Hessian eigenvectors. This option does
!*                                          eigenvector steps until the projected gradient along 
!*                                          the specified Hessian eigenvectors becomes too small.
!*
          call mqc_get_command_argument(i+1,command)
          coordinate = command
          if(coordinate.ne.'eigenvector'.and.coordinate.ne.'gradient'.and.&
            coordinate.ne.'newton'.and.coordinate.ne.'connect') &
            call mqc_error_a('Argument to --coordinate not recognized',6,'coordinate',coordinate)
          j = i+2
        elseif(command.eq.'--step') then
!
!*      --step step-size                 Size of step when following a vector to return perturbed 
!*                                       set of molecular orbitals. As only first order Zassenhaus 
!*                                       formula is implemented, only small step sizes should be used.
!*                                       Default is 0.05. 
!
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F10.5)') initStep
          j = i + 2
        elseif(command.eq.'--nsteps') then
!
!*      --nsteps steps                   Number of steps when following a vector. Default is 1. 
!
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I5)') maxSteps
          j = i + 2
        elseif(command.eq.'--connect-angle') then
!
!*      --connect-angle angle            The angle from perpendicular of the gradient and step vector
!*                                       to specify when the connect algorithm should switch
!*                                       from following the eigenvector to doing Newton steps.
!
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F10.5)') connectThresh
          j = i + 2
!*   4. Help
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
      call fileInfo%load(fileName)
      call fileInfo%getMolData(moleculeInfo)
      Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
      call moleculeInfo%print(iOut)
      call Vnn%print(iOut,'Nuclear Repulsion Energy (au)')

      nBasis = fileInfo%getVal('nBasis')
      nAlpha = fileInfo%getVal('nAlpha')
      nBeta = fileInfo%getVal('nBeta')
      ovDim = 2*nBasis*(nAlpha+nBeta)-2*nAlpha*nBeta-nAlpha**2-nBeta**2 

      call fileInfo%getArray('SHELL TO ATOM MAP',sh2AtMp)
      call fileInfo%getArray('SHELL TYPES',shlTyp)
      call fileInfo%getArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
      call fileInfo%getArray('PRIMITIVE EXPONENTS',prmExp)
      call fileInfo%getArray('CONTRACTION COEFFICIENTS',conCoef)
      call fileInfo%getArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
      call fileInfo%getArray('COORDINATES OF EACH SHELL',shCoor)
      call fileInfo%getESTObj('mo energies',est_eigenvalues=mo_energies)

      allocate(AOeris(1),eris(1))
      call fileInfo%get2ERIs('regular',AOeris(1))
      call fileInfo%getESTObj('overlap',est_integral=overlap)
      call fileInfo%getESTObj('mo coefficients',est_integral=mo_coefficients)
      if(iPrint.ge.2) call mo_coefficients%print(iOut,'MO coefficients')
      call fileInfo%getESTObj('core hamiltonian',est_integral=core_hamiltonian)
      if(iPrint.ge.2) call core_hamiltonian%print(iOut,'Core Hamiltonian')
!
!     If we are performing eigenvector following routine, we start loop here with updated MO 
!     coefficients (all AO integrals are assumed to be the same). 
!
      call cumvec%init(ovDim)
      newtonFlag = .false.
      do iter = 1, maxSteps
        write(*,*) '--------------'
        write(*,*) 'Step: ',iter
        write(*,*) '--------------'

        density = matmul(mo_coefficients%orbitals('occupied',[nAlpha],[nBeta]), &
          dagger(mo_coefficients%orbitals('occupied',[nAlpha],[nBeta])))
        if(iPrint.ge.2) call density%print(iOut,'Density matrix')
        Gmat = contraction(AOeris,density)
        energy = contraction(core_hamiltonian,density)+0.5*contraction(density%swapODB(),Gmat)+Vnn
        call energy%print(6,'Hartree-Fock Energy',FormatStr='F20.12')
!
!     Compute MO Fock matrix (orbital rotation gradient)
!
      AOfock = core_hamiltonian + Gmat
      fock = matmul(matmul(dagger(mo_coefficients),AOfock),mo_coefficients)
      if(iPrint.ge.2) call fock%print(iOut,'MO Fock matrix')
      call twoERI_trans(iOut,iPrint,mo_coefficients,AOeris(1),eris(1))
      if(iPrint.ge.4) call eris(1)%print(iOut,'MO 2ERIs') 
!
!     Establish the symmetry of the input wavefunction and the instability to be tested.
!     If doing eigenvector following, symmetry may have changed, so ignore input symmetry
!     request. 
!
      wf_complex = .false.
      if(mqc_matrix_norm(aimag(fock%getBlock())).gt.thresh) wf_complex = .true.
      if(.not.allocated(wf_string).or.iter.gt.1) wf_string = ''
      if (len(wf_string).eq.0) then
        if(wf_complex) then
          wf_string = 'C'
        else
          wf_string = 'R'
        endIf
      else
        if((wf_string(1:1).ne.'C'.and.wf_string(1:1).ne.'R'.and.wf_string(1:1).ne.'A').or.&
          len(wf_string).ne.1) &
          call mqc_error_A(' Instability test input format incorrect',6,'wf_string',wf_string)
      endIf
      write(iOut,'(A)') ' Testing for instabilities leading to '//wf_string//' wavefunctions'
      if(wf_complex.and.wf_string(1:1).eq.'R') &
        call mqc_error_A(' Requesting real instability test but wavefunction is already complex', &
        6,'wf_string',wf_string)
      if(wf_string.eq.'A') wf_complex=.true.

      if(wf_complex) then
        oRHessDim = 2*ovDim
      else
        oRHessDim = ovDim
      endIf 
      call linearFock%init(oRHessDim)
      call orbRotHess%init(oRHessDim,oRHessDim)
!
      do i = 1,ovDim 
!
!       Extract orbital indices (related to the index in the MO Fock matrix)
!
        occ1 = mod(i-1,nAlpha)+1
        virt1 = nAlpha+mod((i-1)/nAlpha,nBasis-nAlpha)+1
        if((i-1)/(2*nAlpha*(nBasis-nAlpha))+1.eq.2) occ1 = occ1+nBasis
        if(mod((i-1)/(nAlpha*(nBasis-nAlpha)),2)+1.eq.2) virt1 = virt1+nBasis
        call linearFock%put(fock%at(occ1,virt1),i)
          
        do j = 1,i
!
!         Extract orbital indices (related to the index in the MO Fock matrix)
!
          occ2 = mod(j-1,nAlpha)+1
          virt2 = nAlpha+mod((j-1)/nAlpha,nBasis-nAlpha)+1
          if((j-1)/(2*nAlpha*(nBasis-nAlpha))+1.eq.2) occ2 = occ2+nBasis
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
        call linearFock%vput(conjg(linearFock%vat(1,ovDim)),ovDim+1)
      endIf
      if(iPrint.ge.3) call linearFock%print(iOut,'Orbital rotation gradient')
      call mqc_print(sqrt(dot_product(dagger(linearFock),linearFock)),iOut,'Gradient norm',formatStr='F14.8')
      if(iPrint.ge.3) call orbRotHess%print(iOut,'Orbital rotation Hessian')
      call orbRotHess%diag(oREigs,oRVecs)
      if(iPrint.ge.2) call oREigs%print(iOut,'Orbital rotation Hessian eigenvalues')
      if(iPrint.ge.3) call oRVecs%print(iOut,'Orbital rotation Hessian eigenvectors')
!
!     Ensure eigenvectors have vanishing eta norm 
!
      if(wf_complex) then
        degen_flag = .false.
        do i = 1, size(oREigs)-1
          if(abs(oREigs%at(i)-oREigs%at(i+1)).lt.etathresh.and..not.degen_flag) then
!           We have found the start of a degenerate block. Make record as 
!           necessary.
            degen_start = i
            degen_flag = .true.
            cycle
          elseIf(abs(oREigs%at(i)-oREigs%at(i+1)).lt.etathresh.and.degen_flag) then
!           We are in the middle of a degenerate block, nothing to do for now.
            cycle
          elseIf(abs(oREigs%at(i)-oREigs%at(i+1)).ge.etathresh.and.degen_flag) then
!           We have found the end of a degenerate block, time to cycle through 
!           until vectors are in the correct form.
            degen_end = i
            degen_flag = .false.
            do l = 0, maxIters
              do k = 0, degen_end-degen_start
                ivec = degen_start+mod(k,degen_end-degen_start+1)
                jvec = degen_start+mod(k+1,degen_end-degen_start+1)

                do j = ovDim,1,-1
!                 Pick the element with the largest difference between ind and ind+ovDim and make sure there are no zeros. 
                  ind = mqc_vector_Scalar_at(&
                    mqc_vector_argsort(abs(oRVecs%vat([1,ovDim],[ivec])-conjg(oRVecs%vat([ovDim+1,2*ovDim],[ivec])))),j)
                  if(jvec.ne.degen_start) then
                    if(abs(oRVecs%at(ind,ivec)*oRVecs%at(ind+ovDim,ivec)*oRVecs%at(ind,jvec)*oRVecs%at(ind+ovDim,jvec)).lt.&
                      thresh) then 
                      if(j.ne.1) then
                        cycle 
                      else
!                       If we don't have any choice with a non-zero element, we can only do the phi rotation, so pick the
!                       largest element.
                        ind = mqc_vector_Scalar_at(&
                          mqc_vector_argsort(abs(oRVecs%vat([1,ovDim],[ivec])-conjg(oRVecs%vat([ovDim+1,2*ovDim],[ivec])))),ovDim)
                      endIf
                    endIf
                  else
!                   If we are at the last element in the series, just do phi
                    if(abs(oRVecs%at(ind,ivec)*oRVecs%at(ind+ovDim,ivec)).lt.thresh) cycle
                  endIf
                  if(jvec.ne.degen_start.and.j.ne.1) then
                    theta = real(atan((oRVecs%at(ind,ivec)-conjg(oRVecs%at(ind+ovDim,ivec)))/&
                      (oRVecs%at(ind,jvec)-conjg(oRVecs%at(ind+ovDim,jvec)))))
                    vec2process = oRVecs%vat([0],[ivec])
                    call oRVecs%vput(cos(theta)*vec2process-sin(theta)*oRVecs%vat([0],[jvec]),[0],[ivec])
                    call oRVecs%vput(sin(theta)*vec2process+cos(theta)*oRVecs%vat([0],[jvec]),[0],[jvec])
                  endIf
                  phi = real((-1)*cmplx(0.0,1.0)*log(conjg(oRVecs%at(ind+ovDim,ivec))/oRVecs%at(ind,ivec))/2.0)
                  call oRVecs%vput(exp(cmplx(0.0,1.0)*phi)*oRVecs%vat([0],[ivec]),[0],[ivec])
                  exit
                endDo
              endDo

              if(mqc_matrix_norm(oRVecs%mat([1,ovDim],[degen_start,degen_end])-&
                conjg(oRVecs%mat([ovDim+1,2*ovDim],[degen_start,degen_end]))).lt.convThresh) exit
              if(l.eq.maxIters) then
                call mqc_print(mqc_matrix_norm(oRVecs%mat([1,ovDim],[degen_start,degen_end])-&
                  conjg(oRVecs%mat([ovDim+1,2*ovDim],[degen_start,degen_end]))),iOut,&
                  'Failed to put eigenvectors in eta norm zero format, convergence',formatStr='F20.12')
                call mqc_print(oRVecs%mat([0],[degen_start,degen_end]),iOut,'Final eigenvectors')
                call mqc_error('')
              endIf
         
            endDo

          else
!           Only possible option is we have a nondegenerate eigenvector, so we 
!           should already be in correct vector form. Just cycle.
            cycle
          endif
         
        endDo
        if(iPrint.ge.3) call oRVecs%print(6,'Zero eta-norm eigenvectors')
      endIf
      if(iter.eq.1) initialORVecs = orVecs
!
!     Print the eigenvalues and eigenvectors
!
      write(6,'(1x,A)') NEW_LINE('A')//' Lowest '//trim(num2char(neigs2print))//&
        ' orbital rotation Hessian eigenvalues'//NEW_LINE('A')
      do i = 1, min(neigs2print,size(oREigs))
        call mqc_print(oREigs%at(i),6,'Eigenvector '//trim(num2char(i)))
        vec2process = oRVecs%vat([1,ovDim],[i])
        do while(.true.)
          ind = maxLoc(abs(vec2process))
          if(abs(vec2process%at(ind)).lt.vecThresh) then
            write(6,'(1x,A)') NEW_LINE('A')
            exit
          endIf
!          ind2 = mod(ind-1,ovDim)+1
          ind2 = ind
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
            call vec2process%put(0.0,ind)
          else
            vecString = trim(vecString)//'     '//trim(num2char(vec2process%at(ind)))
            call vec2process%put(0.0,ind)
          endIf
          write(6,'(A)') trim(vecString)
        endDo
      endDo
!
!     If eigenvector following, check for convergence to new stationary point and exit if converged.
!
      if(len(coordinate).gt.0.and.iter.gt.1) then
        if(sqrt(dot_product(dagger(linearFock),linearFock)).lt.followThresh) then
          write(iOut,'(1X,A)') 'Convergence to new stationary point'
          exit
        endIf
      endIf
!
!     Build a rotation matrix to follow the specified root
!
      if(size(vecfollow).ne.0) then
        if(wf_string.eq.'C'.and..not.wf_complex) oRVecs = oRVecs*cmplx(0.0,1.0)

        if(iter.eq.1.or.coordinate.eq.'eigenvector'.or.&
          (coordinate.eq.'connect'.and..not.newtonFlag)) then
          jEnd = size(vecFollow)
        else
          jEnd = 1
        endIf

        do j = 1, jEnd
          if(iter.eq.1.or.coordinate.eq.'eigenvector'.or.&
            (coordinate.eq.'connect'.and..not.newtonflag)) then
            vec_overlap = matmul(dagger(initialORVecs%vat([1,size(ORVecs,1)],[int(vecfollow%at(j))])),oRVecs)
            if(iPrint.ge.3) call mqc_print(vec_overlap,iOut,'Overlap of eigenvectors with previous step vector')
            vpos = maxloc(abs(vec_overlap))
            vval = mqc_vector_scalar_at(vec_overlap,vpos)
            if(abs(vval).lt.thresh) exit
            if(abs(real(vval)).le.zero.and.abs(aimag(vval)).gt.zero) then
              theta = (-1)*asin(aimag(vval)/(sqrt(real(vval)**2+aimag(vval)**2)))
            elseIf(abs(aimag(vval)).le.zero.and.abs(real(vval)).gt.zero) then
              theta = (-1)*acos(real(vval)/(sqrt(real(vval)**2+aimag(vval)**2)))
            elseIf(abs(real(vval)).gt.zero.and.abs(aimag(vval)).gt.zero) then
              theta = (-1)*atan2(vval)
            else
              theta = zero
            endIf
            vval = dot_product(dagger(cmplx(cos(theta),sin(theta))*oRVecs%vat([0],[vpos])),&
              initialORVecs%vat([1,size(ORVecs,1)],[int(vecfollow%at(j))]))
            vec2process = cmplx(cos(theta),sin(theta))*oRVecs%vat([0],[vpos])
            write(iOut,'(1X,A)') 'Eigenvector following information -- Initial:'//&
              trim(num2char(vecFollow%at(j),'I3'))//', Current: '//trim(num2char(vpos,'I3'))//&
              ', Overlap: '//trim(num2char(vval,'F10.6'))
!           update the old vector followed to the previous iteration
            call initialoRVecs%vput(vec2process,[0],[int(vecfollow%at(j))]) 
            step = initStep
          elseIf(coordinate.eq.'gradient') then
            vec2process = (-1)*linearFock
            step = initStep
          elseIf(coordinate.eq.'newton'.or.newtonFlag) then
            if(coordinate.eq.'eigenvector') write(iOut,'(1X,A)') &
              'Angle between eigenvector and gradient too small, swapping to Newton step'
!            vec2process = (-1)*matmul(orbRotHess%inv(),linearFock)
!           Below does the Newton step without inverting the Hessian but gives the same result 
            call vec2process%init(ovDim)
            do i = 1, ovDim
              vec2process = vec2process - &
                (dot_product(dagger(oRVecs%vat([0],[i])),linearFock)/oREigs%at(i))*oRVecs%vat([0],[i])
            endDo
            step = 1.0
          endIf

          if(j.eq.jEnd) then
            optangle = acos(dot_product(dagger((-1)*linearFock),vec2process)/&
              (linearFock%norm()*vec2process%norm()))
            call mqc_print(optangle*180.0/pi,6,'Angle of step vector with gradient vector')
            if(abs(optangle*180.0/pi-90.0).lt.connectThresh) newtonflag = .true.
        
            if(iPrint.ge.2) call vec2process%print(iOut,'Step vector')
            cumvec = cumvec + step*vec2process
            if(iPrint.ge.3) call cumvec%print(iOut,'Cumulative displacement vector')
          endIf

          call rotation_matrix%init(nBasis*2,nBasis*2)
          do i = 1, ovDim
            elem1 = mod(i-1,nAlpha)+1+nBasis*(mod((i-1)/(nAlpha*(nBasis-nAlpha)),2))
            elem2 = nAlpha+mod((i-1)/nAlpha,nBasis-nAlpha)+1+nBasis*((i-1)/(2*nAlpha*(nBasis-nAlpha)))
            call rotation_matrix%put((-1)*step*vec2process%at(i),elem1,elem2)
            call rotation_matrix%put(step*conjg(vec2process%at(i)),elem2,elem1)
          endDo
          rotation_matrix = exp(rotation_matrix)
          if(iPrint.ge.3) call rotation_matrix%print(iOut,'Orbital perturbation matrix '//trim(num2char(vecFollow%at(j))))
          mo_coefficients = matmul(mo_coefficients,rotation_matrix)
        endDo
        if(iPrint.ge.2) call mo_coefficients%print(iOut,'Perturbed MOs')
        if(iPrint.ge.3) call mqc_print(matmul(matmul(dagger(mo_coefficients),overlap),mo_coefficients),iOut,&
          'Molecular orbital orthogonality check')

        if(MQC_Matrix_Norm(mo_coefficients%getBlock('alpha')-mo_coefficients%getBlock('beta')).lt.thresh.and.&
          mqc_matrix_norm(mo_coefficients%getBlock('alpha-beta')).lt.thresh.and. &
          mqc_matrix_norm(mo_coefficients%getBlock('beta-alpha')).lt.thresh) then
          doUHF = .false.
          doGHF = .false.
        elseIf(mqc_matrix_norm(mo_coefficients%getBlock('alpha-beta')).lt.thresh.and. &
          mqc_matrix_norm(mo_coefficients%getBlock('beta-alpha')).lt.thresh) then
          doUHF = .true.
          doGHF = .false.
        else
          doUHF = .false.
          doGHF = .true.
        endIf
        if(mqc_matrix_norm(aimag(mo_coefficients%getBlock('full'))).gt.thresh) then
          doComplex = .true.
        else
          doComplex = .false.
        endIf

        if(.not.allocated(outputFileIn)) then
          outputFileIn = fileName
          outputFileIn = outputFileIn(1:(len(fileName)-4))
        endIf
        outputFile = outputFileIn
        i = 1
        file_tmp = trim(outputfile)
        file_exists = .true.
        do while (file_exists)
          inquire(file=trim(file_tmp)//'.mat',exist=file_exists)
          if(file_exists) then
            i = i+1
            file_tmp = trim(outputfile)//'-'
            call build_string_add_int(i,file_tmp,20)
          endIf
        endDo
        outputfile = trim(file_tmp)

        call write_output_file(iPrint,trim(outputfile)//".mat",fileinfo,doUHF,doGHF,doComplex,mo_coefficients,&
          mo_energies,nBasis,sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor)

        select case (otype)
        case ('chk')
          call EXECUTE_COMMAND_LINE("unfchk -matrix "//trim(outputFile)//".mat "//trim(outputFile)//".chk")
          call EXECUTE_COMMAND_LINE("rm "//trim(outputFile)//".mat")
        case ('mat')
!         nothing to do here
        case default
          call mqc_error('Requested output file type not recognized')
        end select
      endIf
!
!     We have finished updating the MOs so loop back to redo the stability if eigenfollowing routine requested.
!
      endDo
!
      contains

      subroutine write_output_file(iPrint,newMatFile,fileinfo,doUHF,doGHF,doComplex,mo_coefficients,mo_energies, &
          nBasis,sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor)
!
!     Write a matrix file containing molecular orbitals that can be visualized in gaussView.
!
!     There is a bug in Wr_LCBuf in qcmatrix.F of gauopen if you want to write complex.
!     NR should only be negative in Wr_Labl and positive everywhere else.
!     Please recomplile MQC after making these changes.
!
      implicit none
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo
      type(mqc_scf_integral),intent(in)::mo_coefficients
      type(mqc_scf_eigenvalues),intent(in)::mo_energies
      integer,intent(in)::iPrint,nBasis
      character(len=*),intent(in)::newMatFile
      logical,intent(in)::doUHF,doGHF,doComplex
      type(mqc_matrix)::tmpMatrix
      type(mqc_matrix)::sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor
      type(mqc_scf_integral)::densitymatrix


      fileinfo%icgu = 111
      if(doUHF) then
        fileinfo%icgu = fileinfo%icgu + 1
      elseIf(doGHF) then
        fileinfo%icgu = fileinfo%icgu + 100
      endIf
      if(doComplex.or.doGHF) then
        fileinfo%icgu = fileinfo%icgu + 10
      endIf

      call fileinfo%create(newMatFile)
      call fileInfo%writeArray('SHELL TO ATOM MAP',sh2AtMp)
      call fileInfo%writeArray('SHELL TYPES',shlTyp)
      call fileInfo%writeArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
      call fileInfo%writeArray('PRIMITIVE EXPONENTS',prmExp)
      call fileInfo%writeArray('CONTRACTION COEFFICIENTS',conCoef)
      call fileInfo%writeArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
      call fileInfo%writeArray('COORDINATES OF EACH SHELL',shCoor)
      If(.not.doUHF.and..not.doGHF.and..not.doComplex) then
        write(iout,'(A)') ' Writing a real restricted solution'
        call fileinfo%writeESTObj('mo coefficients',est_integral=mo_coefficients,override='space',imagORide='real')
        call fileinfo%writeESTObj('mo energies',est_eigenvalues=mo_energies,override='space')
      elseIf(doUHF.and..not.doGHF.and..not.doComplex) then
        write(iout,'(A)') ' Writing a real unrestricted solution'
        call fileinfo%writeESTObj('mo coefficients',est_integral=mo_coefficients,override='spin',imagORide='real')
        call fileinfo%writeESTObj('mo energies',est_eigenvalues=mo_energies,override='spin')
      elseIf(doGHF.and..not.doUHF.and..not.doComplex) then
        write(iout,'(A)') ' Writing a real general solution'
        call fileinfo%writeESTObj('mo coefficients',est_integral=mo_coefficients,override='general',imagORide='real')
        call fileinfo%writeESTObj('mo energies',est_eigenvalues=mo_energies,override='general')
      elseIf(doComplex.and..not.doUHF.and..not.doGHF) then
        write(iout,'(A)') ' Writing a complex restricted solution'
        call fileinfo%writeESTObj('mo coefficients',est_integral=mo_coefficients,override='space',imagORide='complex')
        call fileinfo%writeESTObj('mo energies',est_eigenvalues=mo_energies,override='space')
      elseIf(doComplex.and.doUHF.and..not.doGHF) then
        write(iout,'(A)') ' Writing a complex unrestricted solution'
        call fileinfo%writeESTObj('mo coefficients',est_integral=mo_coefficients,override='spin',imagORide='complex')
        call fileinfo%writeESTObj('mo energies',est_eigenvalues=mo_energies,override='spin')
      elseIf(doComplex.and.doGHF.and..not.doUHF) then
        write(iout,'(A)') ' Writing a complex general solution'
        call fileinfo%writeESTObj('mo coefficients',est_integral=mo_coefficients,override='general',imagORide='complex')
        call fileinfo%writeESTObj('mo energies',est_eigenvalues=mo_energies,override='general')
      else
        call mqc_error(' Unknown symmetry type in guessGen output')
      endIf
      call Close_MatF(fileinfo%UnitNumber)
!
    end subroutine write_output_file
!    
!
    subroutine eigenfollow(string_in,vector_out)
!
!   Determine list of eigenvectors to follow.
!
!   string_in: string to be parsed
!   vector_out: output vector
!
    implicit none
!
!   input/output variables
    character(len=*),intent(in)::string_in
    type(mqc_vector),intent(inOut)::vector_out

!   text parsing variables
    character(len=80)::processString
    integer::i,evec
    logical::newNum

    do i = 1,len(string_in)
      select case (string_in(i:i))
      case('[','\(')
        if(i.eq.1) then
          newNum = .true.
          cycle
        else
          call mqc_error_A('Eigenvector following selection input format incorrect',6,'string_in', &
            string_in)
        endIf
      case(']','\)')
        if(i.eq.len(string_in).and..not.newNum.and.i.ne.1) then
          read(processString,'(I10)') evec
          if(vector_out%size().ge.1) then
            if(vector_out%at(-1).ge.evec) &
              call mqc_error_I('Trying to specify eigenvectors to follow out of order',&
              6,'last',int(vector_out%at(-1)),'entry',evec)
          endIf
          call vector_out%push(evec)
          exit
        else
          call mqc_error_A('Eigenvector following selection input format incorrect',6,'string_in', &
            string_in)
        endIf
      case('0':'9')
        if(i.eq.1.or.i.eq.len(string_in)) then
          call mqc_error_A('Eigenvector following selection input format incorrect',6,'string_in', &
            string_in)
        else
          if(newNum) then
            processString = string_in(i:i)
          else
            processString = trim(processString)//string_in(i:i)
          endIf
          newNum = .false.
          cycle
        endIf
      case(',',' ')
        if(i.eq.1.or.i.eq.len(string_in).or.i.eq.2.or.i.eq.len(string_in)-1.or.newNum) then
          call mqc_error_A('Eigenvector following selection input format incorrect',6,'string_in', &
            string_in)
        else
          read(processString,'(I10)') evec
          if(vector_out%size().ge.1) then
            if(vector_out%at(-1).ge.evec) &
              call mqc_error_I('Trying to specify eigenvectors to follow out of order', &
              6,'last',int(vector_out%at(-1)),'entry',evec)
          endIf
          call vector_out%push(evec)
          newNum = .true.
          cycle
        endIf
      case default
        call mqc_error_A('Unrecognised character in eigenvector following input format', &
          6,'value',string_in(i:i))
      endselect
    endDo
!  
    end subroutine eigenfollow 
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
