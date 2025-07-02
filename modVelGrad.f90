module modVelGrad
implicit none
integer, dimension(3), parameter :: nt = (/256,256,256/)
real, dimension(3), parameter :: dx = (/1.0,1.0,1.0/), l = nt*dx
real, parameter :: pi=3.14159265358979
contains

subroutine avgVelGAndReStress(jj,vF,vB,width,nVels,dVeldx,ReStress)
integer, intent(in) :: jj
real, intent(in) :: vF(3), vB(3), width
real, intent(inout), dimension(3,3) :: dVeldx, ReStress
integer, intent(inout) :: nVels
real :: vec(3)
nVels=nVels+1
vec(:) = ( vF(:) - vB(:) ) / width
dVeldx(:,jj) = dVeldx(:,jj) + ( vec(:) - dVeldx(:,jj) ) / nVels
vec = ( vF(:)*vF(jj) - vB(:)*vB(jj) ) / width
ReStress(:,jj) = ReStress(:,jj) + ( vec(:) - ReStress(:,jj) ) / nVels
end subroutine avgVelGAndReStress

subroutine velGBox(realPos,rad,vel,dVeldx,ReStress)
real, intent(in) :: realPos(3), rad
real, intent(in), dimension(3,nt(1),nt(2),nt(3)) :: vel
real, intent(out), dimension(3,3) :: dVeldx, ReStress
integer :: boxWid, pos(3)
integer :: i,j,k,jj,ip,jp,kp,im,jm,km,nVels
boxWid=nint(rad/dx(3))
pos=nint(realPos/dx)
dVeldx=0.0
ReStress=0.0
jj=1
ip = pos(1) + boxWid
ip = modulo(ip-1, nt(jj)) + 1
im = pos(1) - boxWid
im = modulo(im-1, nt(jj)) + 1
nVels = 0
do k = pos(3)-boxWid, pos(3)+boxWid
  kp = modulo(k-1, nt(3)) + 1
  do j = pos(2)-boxWid, pos(2)+boxWid
    jp = modulo(j-1, nt(2)) + 1
    call avgVelGAndReStress(jj, vel(:,ip,jp,kp), vel(:,im,jp,kp), 2.0*boxWid*dx(jj), nVels, dVeldx, ReStress)
  enddo
enddo
jj=2
jp = pos(jj) + boxWid
jp = modulo(jp-1, nt(jj)) + 1
jm = pos(jj) - boxWid
jm = modulo(jm-1, nt(jj)) + 1
nVels = 0
do k = pos(3)-boxWid, pos(3)+boxWid
  kp = modulo(k-1, nt(3)) + 1
  do i = pos(1)-boxWid, pos(1)+boxWid
    ip = modulo(i-1, nt(1)) + 1
    call avgVelGAndReStress(jj, vel(:,ip,jp,kp), vel(:,ip,jm,kp), 2.0*boxWid*dx(jj), nVels, dVeldx, ReStress)
  enddo
enddo
jj=3
kp = pos(jj) + boxWid
kp = modulo(kp-1, nt(jj)) + 1
km = pos(jj) - boxWid
km = modulo(km-1, nt(jj)) + 1
nVels = 0
do j = pos(2)-boxWid, pos(2)+boxWid
  jp = modulo(j-1, nt(2)) + 1
  do i = pos(1)-boxWid, pos(1)+boxWid
    ip = modulo(i-1, nt(1)) + 1
    call avgVelGAndReStress(jj, vel(:,ip,jp,kp), vel(:,ip,jp,km), 2.0*boxWid*dx(jj), nVels, dVeldx, ReStress)
  enddo
enddo
end subroutine velGBox

subroutine velGModes(pos,rad,uHat,vHat,wHat,dVeldx,ReStress)
real, intent(in) :: pos(3), rad
complex, intent(in), dimension(nt(1)/2+1,nt(2),nt(3)) :: uHat,vHat,wHat
real, intent(out), dimension(3,3) :: dVeldx, ReStress
integer :: i,j,k,ii,jm,km
real :: wavNum(3), wav, weight
complex :: eikdotx, vel(3)
dVeldx=0.0
ReStress=0.0
do k=1,nt(3)
  km=k-1
  if (km.gt.nt(3)/2) km=km-nt(3)
  wavNum(3)=km*2*pi/l(3)
  do j=1,nt(2)
    jm=j-1
    if (jm.gt.nt(2)/2) jm=jm-nt(2)
    wavNum(2)=jm*2*pi/l(2)
    do i=1,nt(1)/2+1
      wavNum(1)=(i-1)*2*pi/l(1)
      weight = 2.0/nt(1)/nt(2)/nt(3)
      if(i==1.or.i==nt(1)/2+1) weight = 1.0/nt(1)/nt(2)/nt(3)
      wav = sqrt( wavNum(3)**2 + wavNum(2)**2 + wavNum(1)**2 )
      !Find the bin number for this radius. Bins have width 1.0/pointsPerWvNum.
      if(wav.lt.2*pi/rad) then
        vel = (/uHat(i,j,k), vHat(i,j,k), wHat(i,j,k)/)
        eikdotx = exp( cmplx(0.0, sum( wavNum(:)*pos(:))))
        do ii = 1,3
          dVeldx(:,ii) = dVeldx(:,ii) + weight * real( cmplx(0.0, wavNum(ii) ) * vel(:) * eikdotx)
          ReStress(:,ii) = ReStress(:,ii) + weight**2 * real( cmplx(0.0, wavNum(ii) ) * vel(:) * vel(ii) * eikdotx)
        enddo
      endif
    enddo
  enddo
enddo
end subroutine velGModes

subroutine velGBlob(blob,vel,dVeldx,ReStress)
integer, intent(in), dimension(nt(1),nt(2),nt(3)) :: blob
real, intent(in), dimension(3,nt(1),nt(2),nt(3)) :: vel
real, intent(out), dimension(3,3) :: dVeldx, ReStress
integer :: i,j,k,p,nVels,firs0,last0,jj
dVeldx=0.0
ReStress=0.0
jj=1
nVels = 0
do k=1,nt(3)
  do j=1,nt(2)
    do last0=1,nt(jj)
      p = modulo( last0, nt(jj) ) + 1
      if( blob(last0,j,k).eq.0 .and. blob(p,j,k).ne.0 ) then
        do p=1,nt(1)
          firs0=modulo( last0+p-1, nt(jj)) + 1
          if(blob(firs0,j,k).eq.0) then
            call avgVelGAndReStress(jj, vel(:,firs0,j,k), vel(:,last0,j,k), p*dx(jj), nVels, dVeldx, ReStress)
            exit
          endif
        enddo
      endif
    enddo
  enddo
enddo
jj=2
nVels = 0
do k=1,nt(3)
  do i=1,nt(1)
    do last0=1,nt(jj)
      p = modulo( last0, nt(jj) ) + 1
      if(blob(i,last0,k).eq.0 .and. blob(i,p,k).ne.0) then
        do p=1,nt(jj)
          firs0 = modulo( last0+p-1, nt(jj) ) + 1 
          if(blob(i,firs0,k).eq.0) then
            call avgVelGAndReStress(jj, vel(:,i,firs0,k), vel(:,i,last0,k), p*dx(jj), nVels, dVeldx, ReStress)
            exit
          endif
        enddo
      endif
    enddo
  enddo
enddo
jj=3
nVels = 0
do j=1,nt(2)
  do i=1,nt(1)
    do last0=1,nt(jj)
      p = modulo( last0, nt(jj) ) + 1
      if(blob(i,j,last0).eq.0 .and. blob(i,j,p).ne.0) then
        do p=1,nt(jj)
          firs0 = modulo( last0+p-1, nt(jj) ) + 1 
          if(blob(i,j,firs0).eq.0) then
            call avgVelGAndReStress(jj, vel(:,i,j,firs0), vel(:,i,j,last0), p*dx(jj), nVels, dVeldx, ReStress)
            exit
          endif
        enddo
      endif
    enddo
  enddo
enddo
end subroutine velGBlob

subroutine saveStrain(outDir,domain,time,dVeldx,ReStress)
character(*), intent(in) :: outDir, domain
real, intent(in) :: dVeldx(3,3), ReStress(3,3), time
integer :: i,j,k,ii,jj,straU,dVelU,vortU,error,ReStU
real :: strain(3,3),strainEiVals(3),vort(3),QInva,RInva,work(8)
do ii=1,3
  do jj=1,3
    strain(ii,jj) = 0.5*(dVeldx(ii,jj) + dVeldx(jj,ii))
  enddo
enddo
!vorticity = curl(u) = del x u
vort(1) = dVeldx(3,2) - dVeldx(2,3)
vort(2) = dVeldx(1,3) - dVeldx(3,1)
vort(3) = dVeldx(2,1) - dVeldx(1,2)
do j = 1,3
  do i = 1,3
    !Q citerion as equation (4) from Paul22roleOfBreakup
    QInva = QInva - 0.5*dVeldx(i,j)*dVeldx(j,i)
    do k = 1,3
      !Third invariant of velocty grad tensor as equation (5) from Paul22roleOfBreakup
      !Assumes incompressibility so that 3*det(dVeldx)=tr(dVeldx^3)
      RInva = RInva - (1.0/3.0)*dVeldx(i,j)*dVeldx(j,k)*dVeldx(k,i) 
    enddo
  enddo
enddo
call SSYEV('V','U',3,strain,3,strainEiVals,work,8,error)
call system('mkdir -p '//trim(outDir)//'dVeldx')
open(newunit=dVelU,file=trim(outDir)//'dVeldx/'//trim(domain)//'.txt',access='append',form='formatted')
  write(dVelU,'(12ES16.7E3)') time, QInva, RInva, dVeldx
close(dVelU)
call system('mkdir -p '//trim(outDir)//'strain')
open(newunit=straU,file=trim(outDir)//'strain/'//trim(domain)//'.txt',access='append',form='formatted')
  write(straU,'(13ES16.7E3)') time, strainEiVals, strain
close(straU)
call SSYEV('V','U',3,ReStress,3,strainEiVals,work,8,error)
call system('mkdir -p '//trim(outDir)//'ReStress')
open(newunit=ReStU,file=trim(outDir)//'ReStress/'//trim(domain)//'.txt',access='append',form='formatted')
  write(ReStU,'(13ES16.7E3)') time, strainEiVals, ReStress
close(ReStU)
call system('mkdir -p '//trim(outDir)//'vortic')
open(newunit=vortU,file=trim(outDir)//'vortic/'//trim(domain)//'.txt',access='append',form='formatted')
  write(vortU,'( 4ES16.7E3)') time, vort
close(vortU)
end subroutine saveStrain

subroutine makeSphere(realPos,rad,sphere)
!compute the strain in a sphere with same diameter as drop
real, intent(in) :: realPos(3), rad
integer, intent(out) :: sphere(nt(1),nt(2),nt(3))
integer :: intR,i,j,k,ip,jp,kp,pos(3)
pos=nint(realPos/dx)
intR = nint((rad/dx(3))**2)
do k=1,nt(3)
  kp = abs(pos(3) - k)
  if (kp.gt.nt(3)/2) kp = kp - nt(3) 
  do j=1,nt(2)
    jp = abs(pos(2) - j)
    if (jp.gt.nt(2)/2) jp = jp - nt(2) 
    do i=1,nt(1)
      ip = abs(pos(1) - i)
      if (ip.gt.nt(1)/2) ip = ip - nt(1) 
      if ( ip**2+jp**2+kp**2 .lt. intR ) then
        sphere(i,j,k) = 1
      else
        sphere(i,j,k) = 0
      endif
    enddo
  enddo
enddo
end subroutine makeSphere

subroutine makeEllipse(realPos,MoI,MoIEiVals,dropVol,ellipse)
!make a mask of the ellipse with same moment of inertia as drop
real, intent(in) :: realPos(3), MoI(3,3), MoIEiVals(3), dropVol
integer, intent(out) :: ellipse(nt(1),nt(2),nt(3))
integer :: i,j,k,ii,jj,pos(3),sep(3)
real :: ellMat(3,3)=0.0, rad
do i=1,3
  do j=1,3
    do k=1,3
      ellMat(i,j) = ellMat(i,j) + MoI(i,k) * MoI(j,k) * dropVol/5. / (0.5*sum(MoIEiVals) - MoIEiVals(k))
    enddo
  enddo
enddo
pos=nint(realPos/dx)
do k=1,nt(3)
  sep(3) = pos(3) - k
  if (sep(3).gt. nt(3)/2) sep(3) = sep(3) - nt(3) 
  if (sep(3).lt.-nt(3)/2) sep(3) = sep(3) + nt(3) 
  do j=1,nt(2)
    sep(2) = pos(2) - j
    if (sep(2).gt. nt(2)/2) sep(2) = sep(2) - nt(2) 
    if (sep(2).lt.-nt(2)/2) sep(2) = sep(2) + nt(2) 
    do i=1,nt(1)
      sep(1) = pos(1) - i
      if (sep(1).gt. nt(1)/2) sep(1) = sep(1) - nt(1) 
      if (sep(1).lt.-nt(1)/2) sep(1) = sep(1) + nt(1) 
      rad=0.0
      do ii=1,3
        do jj=1,3
          rad = rad + sep(ii)*ellMat(ii,jj)*sep(jj)
        enddo
      enddo
      if(rad.lt.1) then
        ellipse(i,j,k) = 1
      else
        ellipse(i,j,k) = 0
      endif
    enddo
  enddo
enddo
end subroutine makeEllipse

end module modVelGrad
