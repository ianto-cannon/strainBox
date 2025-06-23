program dropStrain
use hdf5
!use, intrinsic :: ISO_FORTRAN_ENV
use, intrinsic :: iso_c_binding
use mpi
use modVelGrad, only : velGradBox,velGradBlob,velGradModes,saveStrain,makeSphere,nt,l,dx,pi
implicit none
include 'fftw3.f03'
type ragged_array
  !2D array containing vectors of different lengths
  real,allocatable::v(:)
  character(len=12) :: indexName
end type ragged_array
character(len=200) :: filename, runName, weName='we_05/', inDir, outDir
real :: modnor
integer,parameter :: maxMom=2
integer,parameter :: kAlias=int((2.0/3.0)* (nt(3)/2))
real, dimension(nt(1),nt(2),nt(3)) :: kur,u,v,w,phase,dxxPhase
real, dimension(nt(1),nt(2),nt(3)) :: chemPot,dxChemPot,dyChemPot,dzChemPot,surfFX,surfFY,surfFZ
real, dimension(3,nt(1),nt(2),nt(3)) :: nor, vel
!use int8 for non shared arrays to save memory
integer, dimension(nt(1),nt(2),nt(3)) :: drop
integer, dimension(0:1,0:1,0:1) :: paint, neigh
integer, dimension(3)  :: last0,last1,pos
integer :: i,j,k,ip,jp,kp,iq,jq,kq,iShifted,mom,ii,jj,im,jm,km
integer :: cols,paintIt,faceOnCorner,genus,onInt,error,intR,ios,rank,ntask
integer :: statU,posiU,veloU,MoInU,topoU,dirU,listU,specU,forcU,fPosU,fNegU
integer :: runNum
real, dimension(maxMom,3) :: dropPos,dropVel
real, dimension(3) :: MoIEiVals,wavNum,farPos
real, dimension(3,3) :: MoI, dVeldx
real :: diag,deformation,dropArea,dA,Cn,r,time,work(8),We,res,weight
real :: maxNor,kurMean,kurStdDev,kurInv,kurInvSq,invSize,surPow,sumSurPow
type(ragged_array) :: hist(3) !histogram of drop mass in x, y and z directions
!integer(kind=int64) :: dropSize, faces, edges, vertices
integer :: dropSize, faces, edges, vertices
character(len=200) :: fileEnd, str
integer(hid_t) :: file_id, dset_id
integer(hsize_t) :: dims(3)=(/nt(1),nt(2),nt(3)/),  dims1d(1)=(/1/) 
complex, dimension(nt(1)/2+1,nt(2),nt(3)) :: cHat,chemPotHat,dxChemPotHat,dyChemPotHat,dzChemPotHat
complex, dimension(nt(1)/2+1,nt(2),nt(3)) :: surfFXHat,surfFYHat,surfFZHat,uHat,vHat,wHat
real, dimension(0:kAlias) :: FSpec,ESpec,FSpecPos,FSpecNeg
type(C_PTR)  :: plan, plan_inverse
call mpi_init(error)
call mpi_comm_rank(mpi_comm_world,rank,error)
call mpi_comm_size(mpi_comm_world,ntask,error)
!call MPI_Get_processor_name(procname, namelen,error)
if (rank.eq.0) write(*,'(1x,a)') 'starting number of drops calculation'
do ii=1,3
  allocate(hist(ii)%v(nt(ii)))
enddo
hist(1)%indexName='1'
hist(2)%indexName='2'
hist(3)%indexName='3'
plan        =fftwf_plan_dft_r2c_3d(nt(3), nt(2), nt(1), phase, cHat, FFTW_ESTIMATE)
plan_inverse=fftwf_plan_dft_c2r_3d(nt(3), nt(2), nt(1), cHat, phase, FFTW_ESTIMATE)
if (rank.eq.0) call system('ls /home/alberto.velamartin/drop_time/'//trim(weName)//&
                            '/ > ../'//trim(weName)//'dir_list.txt')
call mpi_barrier(mpi_comm_world,error)
open(newunit=dirU, file='../'//trim(weName)//'dir_list.txt', status='old', action='read')
do
  read(dirU, '(A)', iostat=ios) runName
  if (ios /= 0) exit
  str = trim( runName(11:) )
  read( str , *) runNum
  !if ( modulo( runNum, ntask ) .ne. rank) cycle
  if ( runNum .ne. 0) cycle
  inDir='/home/alberto.velamartin/drop_time/'//trim(weName)//trim(runName)
  !outDir='../'//trim(weName)//trim(runName)
  outDir='output/'
  call system('mkdir '//trim(outDir))
  !call system('ls '//trim(inDir)//'/field*.h5 > '//trim(outDir)//'/file_list.txt')
  call system('ls '//trim(inDir)//'/field.008.h5 > '//trim(outDir)//'/file_list.txt')
  fileEnd='.txt'
  open(newunit=statU,file=trim(outDir)//'/statDrops'//trim(fileEnd), access='append',form='formatted',status='REPLACE')
  open(newunit=posiU,file=trim(outDir)//'/posiDrops'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=veloU,file=trim(outDir)//'/veloDrops'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=MoInU,file=trim(outDir)//'/MoInDrops'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=topoU,file=trim(outDir)//'/topoDrops'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=specU,file=trim(outDir)//'/ESpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=forcU,file=trim(outDir)//'/FSpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=fPosU,file=trim(outDir)//'/FSpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  open(newunit=fNegU,file=trim(outDir)//'/FSpec'//trim(fileEnd),access='append',form='formatted',status='REPLACE')
  !Open the generated file list
  open(newunit=listU, file=trim(outDir)//'/file_list.txt', status='old', action='read')
  !Loop through file names
  do
    read(listU, '(A)', iostat=ios) filename
    if (ios /= 0) exit
    if (rank.eq.0) write(*,*) trim(filename)
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
        call h5dread_f(dset_id, H5T_NATIVE_REAL, kur, dims, error)
      call h5dclose_f(dset_id, error)
      do k=1,nt(3)
        do j=1,nt(2)
          do i=1,nt(1)
            phase(k,j,i) = kur(i,j,k)
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
        call h5dread_f(dset_id, H5T_NATIVE_REAL, kur, dims, error)
      call h5dclose_f(dset_id, error)
      do k=1,nt(3)
        do j=1,nt(2)
          do i=1,nt(1)
            u(k,j,i) = kur(i,j,k)
          enddo
        enddo
      enddo
      call h5dopen_f(file_id, 'v', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, kur, dims, error)
      call h5dclose_f(dset_id, error)
      do k=1,nt(3)
        do j=1,nt(2)
          do i=1,nt(1)
            v(k,j,i) = kur(i,j,k)
          enddo
        enddo
      enddo
      call h5dopen_f(file_id, 'w', dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_REAL, kur, dims, error)
      call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    do k=1,nt(3)
      do j=1,nt(2)
        do i=1,nt(1)
          w(k,j,i) = kur(i,j,k)
        enddo
      enddo
    enddo
    vel(1,:,:,:)=u
    vel(2,:,:,:)=v
    vel(3,:,:,:)=w
    write(*,*)'rmsVel',sqrt(sum(vel**2)/nt(3)**3/3)
    dropSize=0
    do k=1,nt(3)
      do j=1,nt(2)
        do i=1,nt(1)
          if(phase(i,j,k).ge.0.9)then
            drop(i,j,k)=1
            dropSize=dropSize+1
          else
            drop(i,j,k)=0
          endif
          !compute cell-centred normal
          ip=i+1
          jp=j+1
          kp=k+1
          if(ip.gt.nt(1)) ip=ip-nt(1)
          if(jp.gt.nt(2)) jp=jp-nt(2)
          if(kp.gt.nt(3)) kp=kp-nt(3)
          im=i-1
          jm=j-1
          km=k-1
          if(im.lt.1) im=im+nt(1)
          if(jm.lt.1) jm=jm+nt(2)
          if(km.lt.1) km=km+nt(3)
          nor(1,i,j,k)=(phase(ip,j,k)-phase(im,j,k))*(0.5/dx(1))
          nor(2,i,j,k)=(phase(i,jp,k)-phase(i,jm,k))*(0.5/dx(2))
          nor(3,i,j,k)=(phase(i,j,kp)-phase(i,j,km))*(0.5/dx(3))
          modnor=sqrt(nor(1,i,j,k)**2+nor(2,i,j,k)**2+nor(3,i,j,k)**2)
          !outward pointing normal
          nor(:,i,j,k)=-nor(:,i,j,k)/modnor
        enddo
      enddo
    enddo
    do k=1,nt(3)
      do j=1,nt(2)
        do i=1,nt(1)
          im=i-1
          jm=j-1
          km=k-1
          if(im.lt.1) im=im+nt(1)
          if(jm.lt.1) jm=jm+nt(2)
          if(km.lt.1) km=km+nt(3)
          ip=i+1
          jp=j+1
          kp=k+1
          if(ip.gt.nt(1)) ip=ip-nt(1)
          if(jp.gt.nt(2)) jp=jp-nt(2)
          if(kp.gt.nt(3)) kp=kp-nt(3)
          !compute curvature in backward direction
          kur(i,j,k)=(nor(1,ip,j,k)-nor(1,im,j,k))*(0.5/dx(1))+ &
                     (nor(2,i,jp,k)-nor(2,i,jm,k))*(0.5/dx(2))+ &
                     (nor(3,i,j,kp)-nor(3,i,j,km))*(0.5/dx(3))
        enddo
      enddo
    enddo
    !calculate drop statistics if minimum i j k in the drop lies in the processor's domain
    invSize = 1.0/dropSize
    dropArea=0.0
    kurMean=0.0
    kurStdDev=0.0
    kurInv=0.0
    kurInvSq=0.0
    faces=0
    edges=0
    vertices=0
    do ii=1,3
      hist(ii)%v=0.0
    enddo
    dropVel=0.0
    !make histograms of drop mass in each direction so we can avoid overlap with periodic 
    !boundaries in moment of inertia calculations
    do k=1,nt(3)
      do j=1,nt(2)
        do i=1,nt(1)
          if(drop(i,j,k).ne.0) then
            hist(1)%v(i) = hist(1)%v(i) + invSize
            hist(2)%v(j) = hist(2)%v(j) + invSize
            hist(3)%v(k) = hist(3)%v(k) + invSize
            do ii=1,3
              do mom = 1,maxMom
                dropVel(mom,ii) = dropVel(mom,ii) + (vel(ii,i,j,k)**mom)*invSize
              enddo
            enddo
          endif
        enddo
      enddo
    enddo 
    dropPos=0.0
    last0=0
    last1=0
    do ii=1,3
      !find beginning and end of drop in x,y,z directions
      do i=1,nt(ii)
        ip=i+1
        if(ip.gt.nt(ii)) ip=ip-nt(ii)
        if(hist(ii)%v(i).gt.0.5*invSize.and.hist(ii)%v(ip).lt.0.5*invSize) last1(ii)=i
        if(hist(ii)%v(i).lt.0.5*invSize.and.hist(ii)%v(ip).gt.0.5*invSize) last0(ii)=i
      enddo
      !if part of the drop touches i=1, shift that part to other side of domain so drop is contiguous
      do i=1,nt(ii)
        iShifted = i
        if(last1(ii).lt.last0(ii).and.i.le.last1(ii)) iShifted = i+nt(ii)
        do mom = 1,maxMom
          dropPos(mom,ii) = dropPos(mom,ii) + hist(ii)%v(i)*(iShifted**mom)
        enddo
      enddo 
    enddo
    MoI=0.0
    do k=1,nt(3)
      pos(3)=k
      if(last1(3).lt.last0(3).and.k.le.last1(3)) pos(3) = k+nt(3)
      do j=1,nt(2)
        pos(2)=j
        if(last1(2).lt.last0(2).and.j.le.last1(2)) pos(2) = j+nt(2)
        do i=1,nt(1)
          pos(1)=i
          if(last1(1).lt.last0(1).and.i.le.last1(1)) pos(1) = i+nt(1)
          !use https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
          !Bunner & Tryggvason JFM 2003 misses out on diag part of MoI definition
          !don't bother moving dropPos inside domain until after the moment of inertia calculations
          if(drop(i,j,k).ne.0) then
            !contribution from each cell
            diag=dx(1)**5*(1.0/6.0)
            do ii=1,3
              diag = diag + ( pos(ii)-dropPos(1,ii) )**2 *dx(1)**5
            enddo
            do ii=1,3
              do jj=1,3
                MoI(ii,jj) = MoI(ii,jj) - ( pos(ii)-dropPos(1,ii) )*( pos(jj)-dropPos(1,jj) )*dx(1)**5
                if(ii.eq.jj) MoI(ii,jj) = MoI(ii,jj) + diag
              enddo
            enddo
          endif
          !the interface is here if drop changes in any of the 6 directions
          onInt=0
          do ii=1,3
            ip = i
            jp = j
            kp = k
            iq = i
            jq = j
            kq = k
            if (ii.eq.1) then
              ip=i+1
              if(ip.gt.nt(1))ip=ip-nt(1)
              iq=i-1
              if(iq.lt.1)  iq=iq+nt(1)
            endif
            if (ii.eq.2) then
              jp=j+1
              if(jp.gt.nt(2))jp=jp-nt(2)
              jq=j-1
              if(jq.lt.1)  jq=jq+nt(2)
            endif
            if (ii.eq.3) then
              kp=k+1
              if(kp.gt.nt(3))kp=kp-nt(3)
              kq=k-1
              if(kq.lt.1)  kq=kq+nt(3)
            endif
            if (drop(ip,jp,kp).ne.drop(i,j,k)) faces = faces + 1
            if (drop(ip,jp,kp).ne.drop(i,j,k).or.drop(iq,jq,kq).ne.drop(i,j,k)) onInt=1
          enddo
          if(onInt.eq.1)then
            !find direction which is most aligned with interface normal
            maxNor=0.0
            dA=0.0
            do ii=1,3
              if(abs(nor(ii,i,j,k)).gt.maxNor)then
                maxNor = abs(nor(ii,i,j,k))
                if(ii.eq.1)dA=dx(2)*dx(3)
                if(ii.eq.2)dA=dx(3)*dx(1)
                if(ii.eq.3)dA=dx(1)*dx(2)
              endif
            enddo
            !area of interface is 
            dA=dA/(2.0*maxNor)
            dropArea = dropArea + dA
            kurMean = kurMean + dA*kur(i,j,k)
            kurStdDev = kurStdDev + dA*kur(i,j,k)**2
            kurInv = kurInv + dA/kur(i,j,k)
            kurInvSq = kurInvSq + dA/kur(i,j,k)**2
          endif
          !I Cannon 2023 May
          !calculate the genus of the drop using Mendoza et al. 2006 Acta Mater 10.1016/j.actamat.2005.10.010.
          !genus = tunnels - voids 
          !for a shape made of simple polygons:
          !genus = 1 - (faces - edges + corners)/2
          !use a painter algorithm to count the number of contiguous drop corners at this vertex
          !paint the vof with up to 8 colours surrounding a vertex between cells
          cols=0
          neigh=0
          paint=0
          do kq=0,1
            kp = k + kq
            if(kp.gt.nt(3)) kp=kp-nt(3)
            do jq=0,1
              jp = j + jq
              if(jp.gt.nt(2)) jp=jp-nt(2)
              do iq=0,1
                ip = i + iq
                if(ip.gt.nt(1)) ip=ip-nt(1)
                neigh(iq,jq,kq)=drop(ip,jp,kp)
                if(drop(ip,jp,kp).ne.0) then
                  cols=cols+1
                  paint(iq,jq,kq)=cols
                endif
              enddo
            enddo
          enddo
          faceOnCorner=0
          !no corners if vertex entirely inside drop
          if(cols.eq.8) cols=0
          !trivial case with all drop except for one corner
          if(cols.eq.7) cols=1
          !if drop only at the corner, there are three faces touching it
          if(cols.eq.1) faceOnCorner=3
          !if just two voids on opposite sides of the vertex, we connect the voids using a tunnel with
          !Euler characteristic = v-e+f = 3-9+6 = 0 = 0-6+6. I choose v=0 to skip the painting loop below
          if(cols.eq.6) then
            cols=0
            faceOnCorner=6
            do jq=0,1
              do iq=0,1
                if(neigh(iq,jq,1).ne.neigh(1-iq,1-jq,0)) then
                  cols=6
                  faceOnCorner=0
                endif
              enddo
            enddo
          endif
          if(cols.gt.1) then
            !in the worst case, we need 4 iterations to paint the 2x2x2 cube by neighbours
            !loop 4 times instead of recursion
            do paintIt=1,4
              do kq=0,1
                do jq=0,1
                  do iq=0,1
                    if(neigh(iq,jq,kq).gt.0) then
                      do ii=1,3
                        ip=iq
                        jp=jq
                        kp=kq
                        if(ii.eq.1) ip=1-iq
                        if(ii.eq.2) jp=1-jq
                        if(ii.eq.3) kp=1-kq
                        !paint any touching cells the same colour
                        if(neigh(iq,jq,kq).eq.neigh(ip,jp,kp)) then
                          if(paint(ip,jp,kp).ne.paint(iq,jq,kq)) then
                            paint(iq,jq,kq) = min(paint(iq,jq,kq), paint(ip,jp,kp))
                            paint(ip,jp,kp) = min(paint(iq,jq,kq), paint(ip,jp,kp))
                          endif
                        !if drop has a sign change, there is a face and an edge
                        else if(paintIt.eq.1) then
                          faceOnCorner=faceOnCorner+1
                        endif
                      enddo
                    endif
                  enddo
                enddo
              enddo
            enddo
            !count how many colours remain, this is equal to the number of drop corners at the vertex
            cols=0
            do paintIt=1,8
              colourLoop: do kq=0,1
                do jq=0,1
                  do iq=0,1
                    if(paintIt.eq.paint(iq,jq,kq)) then
                      cols=cols+1
                      exit colourLoop
                    endif
                  enddo
                enddo
              enddo colourLoop
            enddo
          endif
          vertices = vertices + cols
          !the number of edges touching the vertex is equal to the number of faces
          edges=edges+faceOnCorner
        enddo
      enddo
    enddo
    !every edge is between two corners so we double counted the edges
    edges = edges/2
    if(modulo(vertices-edges+faces,2).ne.0)then
      write(*,*)'warning: Euler characteristic is odd, v,e,f=',vertices,edges,faces
    endif
    genus = 1-(vertices-edges+faces)/2
    kurMean = kurMean / dropArea
    !kurStdDev = dsqrt( kurStdDev/dropArea - kurMean**2 )
    kurStdDev = sqrt( kurStdDev/dropArea - kurMean**2 )
    kurInv = kurInv / dropArea
    kurInvSq = kurInvSq / dropArea
    !Calculate eigenvalues of moment of inertia
    if(last0(1).eq.0.or.last0(2).eq.0.or.last0(3).eq.0) then
      write(*,*) 'drop spans all of domain, deformation is undefined'
      deformation=-1.0
    else
      call SSYEV('V','U',3,MoI,3,MoIEiVals,work,8,error)
      if (error.ne.0) write(*,*) 'dropSize', dropSize, 'MoIEiVals', MoIEiVals
      deformation=sqrt(MoIEiVals(3)/MoIEiVals(1))
    endif
    do ii=1,3
      !move drop back inside domain
      if(dropPos(1,ii).ge.nt(ii)+1) then
        do mom=1,maxMom
          dropPos(mom,ii)=dropPos(mom,ii)-nt(ii)**mom
        enddo
      endif
      !integer coordinates of drop centre are used to make sphere and box
      pos(ii) = nint(dropPos(1,ii))
      !put in units of simulation domain size
      do mom=1,maxMom
        dropPos(mom,ii)=dropPos(mom,ii) * dx(ii)**mom
      enddo
      farPos(ii) = dropPos(1,ii) + 0.5*l(ii)
      if (farPos(ii).gt.l(ii)) farPos(ii) = farPos(ii) - l(ii)
    enddo
    call fftwf_execute_dft_r2c(plan, phase, cHat)
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
          cHat(i,j,k) = - (wavNum(3)**2 + wavNum(2)**2 + wavNum(1)**2) * cHat(i,j,k) /nt(1)/nt(2)/nt(3)
        enddo
      enddo
    enddo
    call fftwf_execute_dft_c2r(plan_inverse, cHat, dxxPhase)
    chemPot = 1/Cn * (phase*phase - 1)*phase - Cn*dxxPhase
    call fftwf_execute_dft_r2c(plan, chemPot, chemPotHat)
    dxChemPotHat = cmplx(0.0,0.0)
    dyChemPotHat = cmplx(0.0,0.0)
    dzChemPotHat = cmplx(0.0,0.0)
    do k=1,nt(3)
      km=k-1
      if (km.gt.nt(3)/2) km=km-nt(3)
      if (abs(km).gt.kAlias) cycle
      wavNum(3)=km*2*pi/l(3)
      do j=1,nt(2)
        jm=j-1
        if (jm.gt.nt(2)/2) jm=jm-nt(2)
        if (abs(jm).gt.kAlias) cycle
        wavNum(2)=jm*2*pi/l(2)
        do i=1,nt(1)/2+1
          if (i-1.gt.kAlias) cycle
          wavNum(1)=(i-1)*2*pi/l(1)
          dxChemPotHat(i,j,k) = cmplx(0.0,1.0) * wavNum(1) * chemPotHat(i,j,k) /nt(1)/nt(2)/nt(3)
          dyChemPotHat(i,j,k) = cmplx(0.0,1.0) * wavNum(2) * chemPotHat(i,j,k) /nt(1)/nt(2)/nt(3)
          dzChemPotHat(i,j,k) = cmplx(0.0,1.0) * wavNum(3) * chemPotHat(i,j,k) /nt(1)/nt(2)/nt(3)
        enddo
      enddo
    enddo
    call fftwf_execute_dft_c2r(plan_inverse, dxChemPotHat, dxChemPot)
    call fftwf_execute_dft_c2r(plan_inverse, dyChemPotHat, dyChemPot)
    call fftwf_execute_dft_c2r(plan_inverse, dzChemPotHat, dzChemPot)
    do k=1,nt(3)
      do j=1,nt(2)
        do i=1,nt(1)
          surfFX(i,j,k) = phase(i,j,k) * dxChemPot(i,j,k)
          surfFY(i,j,k) = phase(i,j,k) * dyChemPot(i,j,k)
          surfFZ(i,j,k) = phase(i,j,k) * dzChemPot(i,j,k)
        enddo
      enddo
    enddo
    !do k=1,250
    !  write(vortU,'(256ES16.7E3)') u(:,k,100)
    !enddo
    sumSurPow = sum(surfFX*u + surfFY*v + surfFZ*w)
    call fftwf_execute_dft_r2c(plan, surfFX, surfFXHat)
    call fftwf_execute_dft_r2c(plan, surfFY, surfFYHat)
    call fftwf_execute_dft_r2c(plan, surfFZ, surfFZHat)
    call fftwf_execute_dft_r2c(plan, u, uHat)
    call fftwf_execute_dft_r2c(plan, v, vHat)
    call fftwf_execute_dft_r2c(plan, w, wHat)
    ESpec=0.0
    FSpec=0.0
    FSpecNeg=0.0
    FSpecPos=0.0
    dVeldx=0.0
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
          r = sqrt( wavNum(3)**2 + wavNum(2)**2 + wavNum(1)**2 )
          !Find the bin number for this radius. Bins have width 1.0/pointsPerWvNum.
          intR = int( r*l(3)/2/pi + 0.5) 
          if(intR.le.kAlias) then
            !Half of the kx domain is missing from FFT of real, so we must double
            weight = 2.0/nt(1)/nt(2)/nt(3)
            if(i==1.or.i==nt(1)/2+1) weight = 1.0/nt(1)/nt(2)/nt(3)
            ESpec(intR) = ESpec(intR) + weight* & 
                                      ( abs(uHat(i,j,k))**2 + &
                                        abs(vHat(i,j,k))**2 + &
                                        abs(wHat(i,j,k))**2 )
            surPow = weight * ( real( conjg(uHat(i,j,k)) * surfFXHat(i,j,k) )+ &
                                real( conjg(vHat(i,j,k)) * surfFYHat(i,j,k) )+ &
                                real( conjg(wHat(i,j,k)) * surfFZHat(i,j,k) ))
            
            FSpec(intR) = FSpec(intR) + surPow 
            if(surPow.gt.0.0) then
              FSpecPos(intR) = FSpecPos(intR) + surPow 
            else
              FSpecNeg(intR) = FSpecNeg(intR) + surPow 
            endif
          endif
        enddo
      enddo
    enddo
    !call velGradBlob(drop,vel,dVeldx)
    !call saveStrain(outDir,'Drop',time,dVeldx)
    do i=1,5
      r=l(3)/12.*i
      write(filename,'(i3.3)') nint(r)
      call velGradBox(dropPos(1,:),r,vel,dVeldx)
      call saveStrain(outDir,'BoxR'//trim(filename),time,dVeldx)
      call velGradBox(farPos,r,vel,dVeldx)
      call saveStrain(outDir,'FarBoxR'//trim(filename),time,dVeldx)
      call makeSphere(dropPos(1,:),r,drop)
      !write(*,*) 'sphere',sum(1.*drop)/(1.*nt(1))**3
      call velGradBlob(drop,vel,dVeldx)
      call saveStrain(outDir,'SphereR'//trim(filename),time,dVeldx)
      call makeSphere(farPos,r,drop)
      call velGradBlob(drop,vel,dVeldx)
      call saveStrain(outDir,'FarSphereR'//trim(filename),time,dVeldx)
      call velGradModes(dropPos(1,:),r,uHat,vHat,wHat,dVeldx)
      call saveStrain(outDir,'Modes'//trim(filename),time,dVeldx)
      call velGradModes(farPos,r,uHat,vHat,wHat,dVeldx)
      call saveStrain(outDir,'FarModesR'//trim(filename),time,dVeldx)
    enddo
    write(statU,'(ES16.7E3,i16,4ES16.7E3,6i16,3ES16.7E3)') time,dropSize,dropArea,deformation,kurMean,kurStdDev,&
          last0(:),last1(:),kurInv,kurInvSq,sumSurPow
    flush(statU)
    ! generate format string for writing 
    write(str,'(a,i0,a)') '(',1+maxMom*3,'(ES16.7E3))'
    write(posiU,str) time, (dropPos(mom,:),mom=1,maxMom)
    flush(posiU)
    write(veloU,str) time, (dropVel(mom,:),mom=1,maxMom)
    flush(veloU)
    write(MoInU,'(13ES16.7E3)') time, MoIEiVals, MoI
    flush(MoInU)
    write(topoU,'(ES16.7E3, 4i16)') time,vertices,edges,faces,genus
    flush(topoU)
    write(str,'(a,i0,a)') '(',nt(3)/2+2,'(ES16.7E3))'
    write(specU,str) time,ESpec 
    flush(specU)
    write(forcU,str) time,FSpec 
    flush(forcU)
    write(fPosU,str) time,FSpecPos
    flush(fPosU)
    write(fNegU,str) time,FSpecNeg
    flush(fNegU)
  enddo
  close(posiU,status='keep')
  close(veloU,status='keep')
  close(MoInU,status='keep')
  close(statU,status='keep')
  close(topoU,status='keep')
  close(specU,status='keep')
  close(forcU,status='keep')
  close(fPosU,status='keep')
  close(fNegU,status='keep')
  close(listU)
  call system('rm '//trim(outDir)//'/file_list.txt')
  call system('./transpose.awk.sh '//trim(outDir)//'/ESpec.txt > '//trim(outDir)//'/ESpecTransp.txt')
  call system('./transpose.awk.sh '//trim(outDir)//'/FSpec.txt > '//trim(outDir)//'/FSpecTransp.txt')
enddo
do ii=1,3
  deallocate(hist(ii)%v)
enddo
call fftw_destroy_plan(plan)
call fftw_destroy_plan(plan_inverse)
call fftw_cleanup()
write(*,*) 'rank',rank,'done'
call mpi_barrier(mpi_comm_world,error)
if (rank.eq.0) call system('rm ../'//trim(weName)//'dir_list.txt')
if (rank.eq.0) write(6,*) 'This is the end'
call mpi_finalize(error)
return
end program dropStrain
