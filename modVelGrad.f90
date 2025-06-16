module modVelGrad
implicit none
integer, dimension(3), parameter :: nt = (/256,256,256/)
real, dimension(3), parameter :: dx = (/1.0,1.0,1.0/), l = nt*dx
contains

subroutine velGradBox(pos,vel,dVeldxBox)
integer, intent(in) :: pos(3)
real, intent(in), dimension(3,nt(1),nt(2),nt(3)) :: vel
real, intent(out), dimension(3,3) :: dVeldxBox
integer, parameter :: boxWid=nint(nt(3)/5.)
integer :: i,j,k,jj,ip,jp,kp,im,jm,km,nVels
real :: dVeldx(3)
dVeldxBox=0.0
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
    nVels = nVels + 1
    dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,im,jp,kp)) / dx(1) / boxWid
    dVeldxBox(:,jj) = dVeldxBox(:,jj) + ( dVeldx(:) - dVeldxBox(:,jj) ) / nVels
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
    nVels = nVels + 1
    dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,ip,jm,kp)) / dx(2) / boxWid
    dVeldxBox(:,jj) = dVeldxBox(:,jj) + ( dVeldx(:) - dVeldxBox(:,jj) ) / nVels
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
    nVels = nVels + 1
    dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,ip,jp,km)) / dx(3) / boxWid
    dVeldxBox(:,jj) = dVeldxBox(:,jj) + ( dVeldx(:) - dVeldxBox(:,jj) ) / nVels
  enddo
enddo
end subroutine velGradBox

subroutine velGradBlob(blob,vel,dVeldxTot)
integer, intent(in), dimension(nt(1),nt(2),nt(3)) :: blob
real, intent(in), dimension(3,nt(1),nt(2),nt(3)) :: vel
real, intent(out), dimension(3,3) :: dVeldxTot
integer :: i,j,k,p,nVels,firs0,last0
real :: dVeldx(3), width
dVeldxTot=0.0
nVels = 0
do k=1,nt(3)
  do j=1,nt(2)
    firs0=0
    last0=0
    do i=1,nt(1)
      p=i+1
      if(p.gt.nt(1)) p=p-nt(1)
      if(blob(i,j,k).ne.0.and.blob(p,j,k).eq.0) firs0=p
      if(blob(i,j,k).eq.0.and.blob(p,j,k).ne.0) last0=i
    enddo
    if (last0.ne.0) then
      width = (firs0 - last0)*dx(1)
      if (width.lt.0.0) width = width + l(1)
      nVels = nVels + 1
      dVeldx(:) = (vel(:,firs0,j,k) - vel(:,last0,j,k)) / width
      dVeldxTot(:,1) = dVeldxTot(:,1) + ( dVeldx(:) - dVeldxTot(:,1) ) / nVels
    endif
  enddo
enddo
nVels = 0
do k=1,nt(3)
  do i=1,nt(1)
    firs0=0
    last0=0
    do j=1,nt(2)
      p=j+1
      if(p.gt.nt(2)) p=p-nt(2)
      if(blob(i,j,k).ne.0.and.blob(i,p,k).eq.0) firs0=p
      if(blob(i,j,k).eq.0.and.blob(i,p,k).ne.0) last0=j
    enddo
    if (last0.ne.0) then
      width = (firs0 - last0)*dx(2)
      if (width.lt.0.0) width = width + l(2)
      nVels = nVels + 1
      dVeldx(:) = (vel(:,i,firs0,k) - vel(:,i,last0,k)) / width
      dVeldxTot(:,2) = dVeldxTot(:,2) + ( dVeldx(:) - dVeldxTot(:,2) ) / nVels
    endif
  enddo
enddo
nVels = 0
do j=1,nt(2)
  do i=1,nt(1)
    firs0=0
    last0=0
    do k=1,nt(3)
      p=k+1
      if(p.gt.nt(3)) p=p-nt(3)
      if(blob(i,j,k).ne.0.and.blob(i,j,p).eq.0) firs0=p
      if(blob(i,j,k).eq.0.and.blob(i,j,p).ne.0) last0=k
    enddo
    if (last0.ne.0) then
      width = (firs0 - last0)*dx(3)
      if (width.lt.0.0) width = width + l(3)
      nVels = nVels + 1
      dVeldx(:) = (vel(:,i,j,firs0) - vel(:,i,j,last0)) / width
      dVeldxTot(:,3) = dVeldxTot(:,3) + ( dVeldx(:) - dVeldxTot(:,3) ) / nVels
    endif
  enddo
enddo
end subroutine velGradBlob

subroutine saveStrain(outDir,domain,time,dVeldx)
character(*), intent(in) :: outDir, domain
real, intent(in) :: dVeldx(3,3), time
integer :: i,j,k,ii,jj,straU,dVelU,vortU,error
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
open(newunit=dVelU,file=trim(outDir)//'/dVeldx'//trim(domain)//'.txt',access='append',form='formatted')
  write(dVelU,'(12ES16.7E3)') time, QInva, RInva, dVeldx
close(dVelU)
open(newunit=straU,file=trim(outDir)//'/strain'//trim(domain)//'.txt',access='append',form='formatted')
  write(straU,'(13ES16.7E3)') time, strainEiVals, strain
close(straU)
open(newunit=vortU,file=trim(outDir)//'/vortic'//trim(domain)//'.txt',access='append',form='formatted')
  write(vortU,'( 4ES16.7E3)') time, vort
close(vortU)
end subroutine saveStrain

subroutine makeSphere(pos,sphere)
!compute the strain in a sphere with same diameter as drop
integer, intent(in) :: pos(3)
integer, intent(out) :: sphere(nt(1),nt(2),nt(3))
integer :: intR,i,j,k,ip,jp,kp
intR = nint((nt(3)/5.)**2)
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

end module modVelGrad
