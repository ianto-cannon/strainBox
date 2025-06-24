program spectrum
use hdf5
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'
integer, parameter :: startTime=20, nTimes=100, n=256
integer, parameter :: kAlias=int((2.0/3.0)* (n/2))
real, parameter :: pi=3.14159265358979, dx=1.0, l=n*dx
logical :: fileExists
character(len=200) :: filename, runName, weName='we_05/', inDir, outDir='output/', fileEnd, str
integer :: i,j,k,jm,km,t,error,intR,ios,dirU,specU,forcU,fPosU,fNegU,runNum
integer(hid_t) :: file_id, dset_id
integer(hsize_t) :: dims(3)=(/n,n,n/),  dims1d(1)=(/1/) 
real :: Cn,r,time,We,res,weight,wavNum(3)
real, dimension(n,n,n) :: u,v,w,phase,dxxPhase
real, dimension(n,n,n) :: cPot,dxCPot,dyCPot,dzCPot,suX,suY,suZ
real, dimension(0:kAlias) :: FSpec,ESpec,FSpecPos,FSpecNeg
complex, dimension(n/2+1,n,n) :: cH,cPotH,dxCPotH,dyCPotH,dzCPotH,suXH,suYH,suZH,uH,vH,wH
complex, dimension(nTimes,n,n,n) :: u4 
complex :: surPow
type(C_PTR)  :: plan, plan_inverse, plan4
plan        =fftwf_plan_dft_r2c_3d(n, n, n, phase, cH, FFTW_ESTIMATE)
plan_inverse=fftwf_plan_dft_c2r_3d(n, n, n, cH, phase, FFTW_ESTIMATE)
plan4       =fftwf_plan_dft(4, [nTimes, n, n, n], u4, u4, 1, FFTW_ESTIMATE)
call system('ls /home/alberto.velamartin/drop_time/'//trim(weName)//&
                            '/ > ../'//trim(weName)//'dir_list.txt')
open(newunit=dirU, file='../'//trim(weName)//'dir_list.txt', status='old', action='read')
do
  read(dirU, '(A)', iostat=ios) runName
  if (ios /= 0) exit
  str = trim( runName(11:) )
  read( str , *) runNum
  inDir='/home/alberto.velamartin/drop_time/'//trim(weName)//trim(runName)
  write(str,'(i3.3)') startTime+nTimes-1
  inquire(file=trim(inDir)//'/field.'//trim(str)//'.h5', exist=fileExists)
  if (.not.fileExists) cycle
  do t=startTime,startTime+nTimes-1
    write(str,'(i3.3)') startTime+nTimes-1
    filename=trim(inDir)//'/field.'//trim(str)//'.h5'
    write(*,*) trim(filename)
    flush(6)
    ! Open the file (read-only)
    call h5open_f(error)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
      call h5dopen_f(file_id, 'Cn', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, Cn, dims1d, error)
      call h5dclose_f(dset_id, error)
      call h5dopen_f(file_id, 'We', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, We, dims1d, error)
      write(*,*) 'We', We
      return
      call h5dclose_f(dset_id, error)
      call h5dopen_f(file_id, 'c', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
      call h5dclose_f(dset_id, error)
      do k=1,n
        do j=1,n
          do i=1,n
            phase(k,j,i) = dxxPhase(i,j,k)
          enddo
        enddo
      enddo
      call h5dopen_f(file_id, 'res', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, res, dims1d, error)
      call h5dclose_f(dset_id, error)
      call h5dopen_f(file_id, 'time', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, time, dims1d, error)
      call h5dclose_f(dset_id, error)
      call h5dopen_f(file_id, 'u', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
      call h5dclose_f(dset_id, error)
      do k=1,n
        do j=1,n
          do i=1,n
            u(k,j,i) = dxxPhase(i,j,k)
          enddo
        enddo
      enddo
      call h5dopen_f(file_id, 'v', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
      call h5dclose_f(dset_id, error)
      do k=1,n
        do j=1,n
          do i=1,n
            v(k,j,i) = dxxPhase(i,j,k)
          enddo
        enddo
      enddo
      call h5dopen_f(file_id, 'w', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, dxxPhase, dims, error)
      call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    do k=1,n
      do j=1,n
        do i=1,n
          w(k,j,i) = dxxPhase(i,j,k)
        enddo
      enddo
    enddo
    call fftwf_execute_dft_r2c(plan, phase, cH)
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
    cPot = 1/Cn * (phase*phase - 1)*phase - Cn*dxxPhase
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
    do k=1,n
      do j=1,n
        do i=1,n
          suX(i,j,k) = phase(i,j,k) * dxCPot(i,j,k)
          suY(i,j,k) = phase(i,j,k) * dyCPot(i,j,k)
          suZ(i,j,k) = phase(i,j,k) * dzCPot(i,j,k)
        enddo
      enddo
    enddo
  enddo
  !call fftwf_execute_dft_r2c(plan4, suX, suXH)
  !call fftwf_execute_dft_r2c(plan4, suY, suYH)
  !call fftwf_execute_dft_r2c(plan4, suZ, suZH)
  call fftwf_execute_dft(plan4, u4, u4)
  !call fftwf_execute_dft_r2c(plan4, v, vH)
  !call fftwf_execute_dft_r2c(plan4, w, wH)
  ESpec=0.0
  FSpec=0.0
  FSpecNeg=0.0
  FSpecPos=0.0
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
        r = sqrt( wavNum(3)**2 + wavNum(2)**2 + wavNum(1)**2 )
        !Find the bin number for this radius. Bins have width 1.0/pointsPerWvNum.
        intR = int( r*l/2/pi + 0.5) 
        if(intR.le.kAlias) then
          !Half of the kx domain is missing from FFT of real, so we must double
          weight = 2.0/n/n/n
          if(i==1.or.i==n/2+1) weight = 1.0/n/n/n
          ESpec(intR) = ESpec(intR) + weight* & 
                                    ( abs(uH(i,j,k))**2 + &
                                      abs(vH(i,j,k))**2 + &
                                      abs(wH(i,j,k))**2 )
          surPow = weight * ( conjg(uH(i,j,k)) * suXH(i,j,k) + &
                              conjg(vH(i,j,k)) * suYH(i,j,k) + &
                              conjg(wH(i,j,k)) * suZH(i,j,k) )
          
          FSpec(intR) = FSpec(intR) + real(surPow)
          if(real(surPow).gt.0.0) then
            FSpecPos(intR) = FSpecPos(intR) + real(surPow) 
          else
            FSpecNeg(intR) = FSpecNeg(intR) + real(surPow)
          endif
        endif
      enddo
    enddo
  enddo
  write(str,'(a,i0,a)') '(',n/2+2,'(ES16.7E3))'
  fileEnd='.txt'
  open(newunit=specU,file=trim(outDir)//'/ESpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
    write(specU,str) time,ESpec 
  close(specU,status='keep')
  open(newunit=forcU,file=trim(outDir)//'/FSpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
    write(forcU,str) time,FSpec 
  close(forcU,status='keep')
  open(newunit=fPosU,file=trim(outDir)//'/ESpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
    write(fPosU,str) time,FSpecPos
  close(fPosU,status='keep')
  open(newunit=fNegU,file=trim(outDir)//'/ESpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
    write(fNegU,str) time,FSpecNeg
  close(fNegU,status='keep')
enddo
call fftw_destroy_plan(plan)
call fftw_destroy_plan(plan_inverse)
call fftw_cleanup()
call system('rm ../'//trim(weName)//'dir_list.txt')
write(6,*) 'This is the end'
return
end program spectrum
