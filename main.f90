program spectrum
use hdf5
use, intrinsic :: iso_c_binding
!use mpi
implicit none
include 'fftw3.f03'
integer, parameter :: startTime=19, nTimes=100, n=256
integer, parameter :: kAlias=int((2.0/3.0)* (n/2))
real, parameter :: pi=3.14159265358979, dx=1.0, l=n*dx
logical :: fileExists, drops=.false.
character(len=200) :: filename, runName, weName='ForIanto/', inDir, fileEnd, str
character(len=200) :: outDir='../we_inf/'
integer :: i,j,k,t,im,jm,km,tm,error,intR,ntask,rank
integer :: ios,dirU,specU,forcU,fPosU,fNegU,fImaU,pSpcU
integer(hid_t) :: file_id, dset_id
integer(hsize_t) :: dims(3)=(/n,n,n/),  dims1d(1)=(/1/) 
real :: Cn,r,time,We,res,weight,wavNum(3),window
real, dimension(n,n,n) :: dxxPhase,cPot,dxCPot,dyCPot,dzCPot
real, dimension(0:kAlias,0:nTimes/2) :: FSpec,ESpec,FSpecPos,FSpecNeg,FImag,PSpec
complex, dimension(n/2+1,n,n) :: cH,cPotH,dxCPotH,dyCPotH,dzCPotH
complex, dimension(n,n,n,nTimes) :: phase,u,v,w,suX,suY,suZ
complex :: surPow
type(C_PTR)  :: plan, plan_inverse, plan4
!call mpi_init(error)
!call mpi_comm_rank(mpi_comm_world,rank,error)
!call mpi_comm_size(mpi_comm_world,ntask,error)
plan        =fftwf_plan_dft_r2c_3d(n, n, n, cPot, cPotH, FFTW_ESTIMATE)
plan_inverse=fftwf_plan_dft_c2r_3d(n, n, n, cPotH, cPot, FFTW_ESTIMATE)
plan4       =fftwf_plan_dft(4, [nTimes, n, n, n], u, u, 1, FFTW_ESTIMATE)
call system('ls /home/alberto.velamartin/drop_time/'//trim(weName)//&
                            '/ > '//trim(outDir)//'dir_list.txt')
open(newunit=dirU, file=trim(outDir)//'dir_list.txt', status='old', action='read')
do
  read(dirU, '(A)', iostat=ios) runName
  if (ios /= 0) exit
  if (trim(runName).ne.'field.000.h5') cycle
  runName=''
  inDir='/home/alberto.velamartin/drop_time/'//trim(weName)//trim(runName)
  write(*,*) trim(inDir)
  write(str,'(i3.3)') startTime+nTimes-1
  inquire(file=trim(inDir)//'/field.'//trim(str)//'.h5', exist=fileExists)
  if (.not.fileExists) then
    write(*,*) trim(inDir)//'/field.'//trim(str)//'.h5', 'NoEexist'
    cycle
  endif
  do t=1,nTimes
    write(str,'(i3.3)') t+startTime
    filename=trim(inDir)//'/field.'//trim(str)//'.h5'
    write(*,*) trim(filename)
    flush(6)
    ! Open the file (read-only)
    !call mpi_barrier(mpi_comm_world,error)
    write(*,*)'ntastk',ntask,'rank',rank
    call h5open_f(error)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
      if(drops)then
        call h5dopen_f(file_id, 'Cn', dset_id, error)
          call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
        call h5dclose_f(dset_id, error)
        call h5dopen_f(file_id, 'We', dset_id, error)
          call h5dread_f(dset_id, H5T_NATIVE_REAL, We, dims1d, error)
        call h5dclose_f(dset_id, error)
        call h5dopen_f(file_id, 'c', dset_id, error)
          call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
        call h5dclose_f(dset_id, error)
        do k=1,n
          do j=1,n
            do i=1,n
              phase(i,j,k,t) = dxxPhase(k,j,i)
            enddo
          enddo
        enddo
      else 
        phase(:,:,:,t)=-1
      endif
      call h5dopen_f(file_id, 'res', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, res, dims1d, error)
      call h5dclose_f(dset_id, error)
      call h5dopen_f(file_id, 'time', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, time, dims1d, error)
      call h5dclose_f(dset_id, error)
      call h5dopen_f(file_id, 'u', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
      call h5dclose_f(dset_id, error)
      write(*,*) 'u2', sum(dxxPhase**2)
      window = 1. - cos( 2.0*pi*(t-1.) / (nTimes-1.) )
      !tanh window
      !if(t.lt.10) then 
      !  window = 0.5+0.5*tanh(t-5.)
      !elseif(t.gt.nTimes-10) then
      !  window = 0.5-0.5*tanh(t-nTimes-5.)
      !else 
      !  window=1
      !endif
      do k=1,n
        do j=1,n
          do i=1,n
            u(i,j,k,t) = dxxPhase(k,j,i)*window
          enddo
        enddo
      enddo
      call h5dopen_f(file_id, 'v', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
      call h5dclose_f(dset_id, error)
      write(*,*) 'v2', sum(dxxPhase**2)
      do k=1,n
        do j=1,n
          do i=1,n
            v(i,j,k,t) = dxxPhase(k,j,i)*window
          enddo
        enddo
      enddo
      call h5dopen_f(file_id, 'w', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
      call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    write(*,*) 'w2', sum(dxxPhase**2)
    do k=1,n
      do j=1,n
        do i=1,n
          w(i,j,k,t) = dxxPhase(k,j,i)*window
        enddo
      enddo
    enddo
    dxxPhase = real(phase(:,:,:,t))
    call fftwf_execute_dft_r2c(plan, dxxPhase, cH)
    do k=1,n
      km=k-1
      if (km.gt.n/2) km=km-n
      wavNum(3)=km*2*pi/l
      do j=1,n
        jm=j-1
        if (jm.gt.n/2) jm=jm-n
        wavNum(2)=jm*2*pi/l
        do i=1,n/2+1
          wavNum(1)=(i-1)*2*pi/l
          cH(i,j,k) = - (wavNum(3)**2 + wavNum(2)**2 + wavNum(1)**2) * cH(i,j,k) /n/n/n
        enddo
      enddo
    enddo
    call fftwf_execute_dft_c2r(plan_inverse, cH, dxxPhase)
    cPot = 1/Cn * (phase(:,:,:,t)*phase(:,:,:,t) - 1)*phase(:,:,:,t) - Cn*dxxPhase
    call fftwf_execute_dft_r2c(plan, cPot, cPotH)
    dxCPotH = cmplx(0.0,0.0)
    dyCPotH = cmplx(0.0,0.0)
    dzCPotH = cmplx(0.0,0.0)
    do k=1,n
      km=k-1
      if (km.gt.n/2) km=km-n
      if (abs(km).gt.kAlias) cycle
      wavNum(3)=km*2*pi/l
      do j=1,n
        jm=j-1
        if (jm.gt.n/2) jm=jm-n
        if (abs(jm).gt.kAlias) cycle
        wavNum(2)=jm*2*pi/l
        do i=1,n/2+1
          if (i-1.gt.kAlias) cycle
          wavNum(1)=(i-1)*2*pi/l
          dxCPotH(i,j,k) = cmplx(0.0,1.0) * wavNum(1) * cPotH(i,j,k) /n/n/n
          dyCPotH(i,j,k) = cmplx(0.0,1.0) * wavNum(2) * cPotH(i,j,k) /n/n/n
          dzCPotH(i,j,k) = cmplx(0.0,1.0) * wavNum(3) * cPotH(i,j,k) /n/n/n
        enddo
      enddo
    enddo
    call fftwf_execute_dft_c2r(plan_inverse, dxCPotH, dxCPot)
    call fftwf_execute_dft_c2r(plan_inverse, dyCPotH, dyCPot)
    call fftwf_execute_dft_c2r(plan_inverse, dzCPotH, dzCPot)
    suX(:,:,:,t) = phase(:,:,:,t) * dxCPot(:,:,:)*window
    suY(:,:,:,t) = phase(:,:,:,t) * dyCPot(:,:,:)*window
    suZ(:,:,:,t) = phase(:,:,:,t) * dzCPot(:,:,:)*window
    phase(:,:,:,t) = phase(:,:,:,t) * window
  enddo
  write(*,*) 'fft u'
  call fftwf_execute_dft(plan4, u, u)
  write(*,*) 'fft v'
  call fftwf_execute_dft(plan4, v, v)
  write(*,*) 'fft w'
  call fftwf_execute_dft(plan4, w, w)
  write(*,*) 'fft suX'
  call fftwf_execute_dft(plan4, suX, suX)
  write(*,*) 'fft suY'
  call fftwf_execute_dft(plan4, suY, suY)
  write(*,*) 'fft suZ'
  call fftwf_execute_dft(plan4, suZ, suZ)
  write(*,*) 'fft phase'
  call fftwf_execute_dft(plan4, phase, phase)
  ESpec=0.0
  FSpec=0.0
  FImag=0.0
  FSpecNeg=0.0
  FSpecPos=0.0
  PSpec=0.0
  weight = 1.0/n/n/n/nTimes
  do t=1,nTimes
    tm=t-1
    if (tm.gt.nTimes/2) tm=abs(tm-nTimes)
    do k=1,n
      km=k-1
      if (km.gt.n/2) km=km-n
      wavNum(3)=km*2*pi/l
      do j=1,n
        jm=j-1
        if (jm.gt.n/2) jm=jm-n
        wavNum(2)=jm*2*pi/l
        do i=1,n
          im=i-1
          if (im.gt.n/2) im=im-n
          wavNum(1)=im*2*pi/l
          r = sqrt( wavNum(3)**2 + wavNum(2)**2 + wavNum(1)**2 )
          !Find the bin number for this radius. Bins have width 1.0/pointsPerWvNum.
          intR = int( r*l/2/pi + 0.5) 
          if(intR.le.kAlias) then
            ESpec(intR,tm) = ESpec(intR,tm) + weight* & 
                                      ( abs(u(i,j,k,t))**2 + &
                                        abs(v(i,j,k,t))**2 + &
                                        abs(w(i,j,k,t))**2 )
            surPow = weight * ( conjg(u(i,j,k,t)) * suX(i,j,k,t) + &
                                conjg(v(i,j,k,t)) * suY(i,j,k,t) + &
                                conjg(w(i,j,k,t)) * suZ(i,j,k,t) )
            
            FSpec(intR,tm) = FSpec(intR,tm) + real(surPow)
            FImag(intR,tm) = FImag(intR,tm) + aimag(surPow)
            PSpec(intR,tm) = PSpec(intR,tm) + weight * r**2 * abs(phase(i,j,k,t))**2 
            if(real(surPow).gt.0.0) then
              FSpecPos(intR,tm) = FSpecPos(intR,tm) + real(surPow) 
            else
              FSpecNeg(intR,tm) = FSpecNeg(intR,tm) + real(surPow)
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  write(str,'(a,i0,a)') '(',nTimes/2+1,'(ES16.7E3))'
  fileEnd=trim(runName)//'.txt'
  open(newunit=specU,file=trim(outDir)//'ESpec_'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=forcU,file=trim(outDir)//'FSpec_'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=fImaU,file=trim(outDir)//'FImag_'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=fPosU,file=trim(outDir)//'FPosi_'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=fNegU,file=trim(outDir)//'FNega_'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=pSpcU,file=trim(outDir)//'PSpec_'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
    do k=0,kAlias
      write(specU,str) ESpec(k,:)
      write(forcU,str) FSpec(k,:) 
      write(fImaU,str) FImag(k,:) 
      write(fPosU,str) FSpecPos(k,:)
      write(fNegU,str) FSpecNeg(k,:)
      write(pSpcU,str) PSpec(k,:)
    enddo
  close(specU,status='keep')
  close(forcU,status='keep')
  close(fImaU,status='keep')
  close(fPosU,status='keep')
  close(fNegU,status='keep')
  close(pSpcU,status='keep')
enddo
call fftw_destroy_plan(plan)
call fftw_destroy_plan(plan_inverse)
call fftw_destroy_plan(plan4)
call fftw_cleanup()
call system('rm '//trim(outDir)//'dir_list.txt')
write(6,*) 'This is the end'
!call mpi_finalize(error)
return
end program spectrum
