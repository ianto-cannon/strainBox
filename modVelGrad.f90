module modVelGrad
implicit none
integer, dimension(3), parameter :: nt = (/256,256,256/)
real, dimension(3), parameter :: dx = (/1.0,1.0,1.0/), l = nt*dx
contains

subroutine velGradBlob(blob,vel,dVeldxTot)
integer, intent(in), dimension(nt(1),nt(2),nt(3)) :: blob
real, intent(in), dimension(3,nt(1),nt(2),nt(3)) :: vel
real, intent(out), dimension(3,3) :: dVeldxTot
integer, parameter :: boxWid=nint(nt(3)/6.)
integer :: i,j,k,jj,ip,jp,kp,im,jm,km,nVels
real :: dVeldx(3)
dVeldxTot=0.0
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
    dVeldxTot(:,jj) = dVeldxTot(:,jj) + ( dVeldx(:) - dVeldxTot(:,jj) ) / nVels
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
    dVeldxTot(:,jj) = dVeldxTot(:,jj) + ( dVeldx(:) - dVeldxTot(:,jj) ) / nVels
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
    dVeldxTot(:,jj) = dVeldxTot(:,jj) + ( dVeldx(:) - dVeldxTot(:,jj) ) / nVels
  enddo
enddo
end subroutine velGradBox

subroutine velGradSphere(pos,vel,dVeldxTot)
integer, intent(in) :: pos(3)
real, intent(in), dimension(3,nt(1),nt(2),nt(3)) :: vel
real, intent(out), dimension(3,3) :: dVeldxTot
integer, parameter :: sphRad=nint(nt(3)/6.)
integer :: i,j,k,jj,ip,jp,kp,im,jm,km,nVels
real :: dVeldx(3)
dVeldxTot=0.0
jj=1
ip = pos(1) + sphRad
ip = modulo(ip-1, nt(jj)) + 1
im = pos(1) - sphRad
im = modulo(im-1, nt(jj)) + 1
nVels = 0
do k = pos(3)-sphRad, pos(3)+sphRad
  kp = modulo(k-1, nt(3)) + 1
  do j = pos(2)-sphRad, pos(2)+sphRad
    jp = modulo(j-1, nt(2)) + 1
    nVels = nVels + 1
    dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,im,jp,kp)) / dx(1) / sphRad
    dVeldxTot(:,jj) = dVeldxTot(:,jj) + ( dVeldx(:) - dVeldxTot(:,jj) ) / nVels
  enddo
enddo
jj=2
jp = pos(jj) + sphRad
jp = modulo(jp-1, nt(jj)) + 1
jm = pos(jj) - sphRad
jm = modulo(jm-1, nt(jj)) + 1
nVels = 0
do k = pos(3)-sphRad, pos(3)+sphRad
  kp = modulo(k-1, nt(3)) + 1
  do i = pos(1)-sphRad, pos(1)+sphRad
    ip = modulo(i-1, nt(1)) + 1
    nVels = nVels + 1
    dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,ip,jm,kp)) / dx(2) / sphRad
    dVeldxTot(:,jj) = dVeldxTot(:,jj) + ( dVeldx(:) - dVeldxTot(:,jj) ) / nVels
  enddo
enddo
jj=3
kp = pos(jj) + sphRad
kp = modulo(kp-1, nt(jj)) + 1
km = pos(jj) - sphRad
km = modulo(km-1, nt(jj)) + 1
nVels = 0
do j = pos(2)-sphRad, pos(2)+sphRad
  jp = modulo(j-1, nt(2)) + 1
  do i = pos(1)-sphRad, pos(1)+sphRad
    ip = modulo(i-1, nt(1)) + 1
    nVels = nVels + 1
    dVeldx(:) = (vel(:,ip,jp,kp) - vel(:,ip,jp,km)) / dx(3) / sphRad
    dVeldxTot(:,jj) = dVeldxTot(:,jj) + ( dVeldx(:) - dVeldxTot(:,jj) ) / nVels
  enddo
enddo
end subroutine velGradSphere

end module modVelGrad
