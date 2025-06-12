module mod_vel_grad
implicit none
integer, dimension(3), parameter :: nt = (/256,256,256/)
integer, parameter :: boxWid=nint(nt(3)/6.)
real, dimension(3), parameter :: l = nt, dx=l/nt
contains

subroutine velGradBox(pos,vel,dVeldxBox)
integer, intent(in) :: pos(3)
real, intent(in), dimension(3,nt(1),nt(2),nt(3)) :: vel
real, intent(out), dimension(3,3) :: dVeldxBox
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
end module mod_vel_grad
