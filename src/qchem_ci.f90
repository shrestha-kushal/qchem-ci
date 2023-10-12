program qchem_ci
implicit none
integer::nao,nbasis,natom,nelec,exmax,xdim,&
ical,iter,stopi
double precision,allocatable::chg(:),coorx(:),&
coory(:),coorz(:),bb(:,:),cbb(:,:),vlobe(:,:),&
civec(:),svar(:),start(:),step(:),xmin(:),&
cmo2(:,:),civec2(:),cmo(:,:),ovl(:,:),h1(:,:),&
bbc(:,:),cbbc(:,:),Hci(:,:),fciv(:,:),cieig(:),&
pdens(:,:),fock(:,:),jkmat(:,:,:,:),ovx(:,:),&
ovu(:,:),ovs(:),vq(:,:),h1t(:,:)
integer,allocatable::nref(:),ncor(:),ngau(:),&
nlobe(:),nvdet(:),blokij(:),ntype(:)
double precision::toten,enuc,smo,hfen,orthp,&
pfac,curen,ylo,hfenc,gsen
! get data and init vars
call initinp(chg,coorx,coory,coorz,bb,cbb,&
     nref,ncor,ngau,nao,nbasis,natom,cmo,&
     ovl,h1,nelec,nlobe,vlobe,xdim,nvdet,&
     blokij,civec,svar,exmax,pfac,start,step,&
     xmin,cmo2,civec2,bbc,cbbc,ntype,stopi,&
     Hci,fciv,cieig,pdens,fock,jkmat,ovu,&
     ovx,ovs,vq,h1t)
gsen=molprop(stopi)
! optimize parameters and lower energy
!call simplex(emin,xdim,xdim+1,start,xmin,&
!ylo,step,ical,iter)
! deallocate arrays
deallocate(chg,coorx,coory,coorz,bb,cbb,vlobe,&
civec,svar,start,step,xmin,cmo2,civec2,cmo,ovl,&
h1,bbc,cbbc,nref,ncor,ngau,nlobe,nvdet,blokij,&
ntype,fciv,cieig,pdens,fock,jkmat,ovx,ovu,ovs,&
vq,h1t)

contains

! energy function
! all global variables seen!
double precision function emin(x)
implicit none
double precision,intent(in)::x(xdim)
write(*,'(A10,1X,I8,A1,I8)') 'Iteration:',iter,'-',ical
! apply parameters
call applyx(x,svar,bb,ngau,nao,nbasis,&
     pdens,nvdet,civec,bbc,civec2,cmo2,cbbc)
! only do the following once for optimized
! primitive exponents and coefficients
if (.true.) then
   ! normalize all types of ao's
   call normao(nao,coorx,coory,coorz,&
        bb,cbb,nref,ncor,ngau,nlobe,vlobe)
   ! get overlap
   call matovl(ovl,coorx,coory,coorz,bb,cbb,&
        nref,ncor,ngau,nbasis,nlobe,vlobe)
   ! get nuclear repulsion
   call enrep(enuc,natom,chg,coorx,coory,&
        coorz)
   ! get kinetic mat
   call matkin(h1,coorx,coory,coorz,bb,cbb,&
        nref,ncor,ngau,nbasis,nlobe,vlobe)
   ! add attraction to kinetic matrix
   call attr_h1(chg,coorx,coory,coorz,bb,cbb,&
        nref,ncor,ngau,nbasis,natom,h1,nlobe,&
        vlobe)
   ! get repulsion integrals
   call get_jkmat(jkmat,coorx,coory,coorz,bb,&
        cbb,nref,ncor,ngau,nbasis,nlobe,vlobe)
end if
! normalize mos
!call normmos(cmo,coorx,coory,coorz,bb,cbb,&
!     nref,ncor,ngau,ovl,nbasis)
! orthogonalize mos
!call gramschmidt(cmo,nbasis,ovl)
! normalize orthogonalized mos
!call normmos(cmo,coorx,coory,coorz,bb,cbb,&
!     nref,ncor,ngau,ovl,nbasis)
! run an HF SCF cycle
call scf(jkmat,h1,ovl,fock,ovu,ovx,ovs,pdens,&
     cmo,nbasis,nelec,hfen)
! old way to calculate energy (does restricted
! closed and open shell)
!call deten(cmo,h1,nbasis,coorx,coory,coorz,&
!     bb,cbb,nref,ncor,ngau,nelec,enuc,hfen,&
!     natom,nlobe,vlobe,jkmat)
! store HF 1 det energy
hfenc=hfen+enuc
! get CI result up to exmax excitations if requested
if (exmax/=0) then
   ! normalize GS CI coefficients
   !call normcic(civec,size(civec))
   civec=1.0d0
   ! unscaled contrib to CI energy
   hfen=hfen*(civec(1)**2)
   ! get CI Hamiltonian
   call fullci(coorx,coory,coorz,bb,cbb,nref,&
        ncor,ngau,nelec,natom,nlobe,vlobe,cmo,&
        h1,nbasis,toten,exmax,civec,blokij,Hci)
   ! get All CI vectors
   Hci(1,1)=hfenc-enuc
   call diagci(Hci,fciv,cieig)
   ! final GS energy expression for simplex
   emin=cieig(1)+enuc+orthp
else
   ! final GS energy expression for simplex
   emin=hfenc+orthp
end if
! print stuff
call printdata(iter,ical,hfenc,hfen,toten,enuc,&
     orthp,bb,nao,ngau,cmo,nbasis,civec,1,ntype,&
     nref,ncor,nlobe,Hci,fciv,cieig)
! stop if ical past
if (ical>=stopi) then
   stop
end if
end function emin

! time propagation
double precision function molprop(maxt)
implicit none
integer,intent(in)::maxt
integer::t
double precision::dt,ct
write(*,'(A10,1X,I8,A1,I8)') 'Iteration:',iter,'-',ical
! only do the following once for optimized
! primitive exponents and coefficients
! normalize all types of ao's
call normao(nao,coorx,coory,coorz,&
     bb,cbb,nref,ncor,ngau,nlobe,vlobe)
! get overlap
call matovl(ovl,coorx,coory,coorz,bb,cbb,&
     nref,ncor,ngau,nbasis,nlobe,vlobe)
! get nuclear repulsion
call enrep(enuc,natom,chg,coorx,coory,&
     coorz)
! get kinetic mat
call matkin(h1,coorx,coory,coorz,bb,cbb,&
     nref,ncor,ngau,nbasis,nlobe,vlobe)
! add attraction to kinetic matrix
call attr_h1(chg,coorx,coory,coorz,bb,cbb,&
     nref,ncor,ngau,nbasis,natom,h1,nlobe,&
     vlobe)
! get repulsion integrals
call get_jkmat(jkmat,coorx,coory,coorz,bb,&
     cbb,nref,ncor,ngau,nbasis,nlobe,vlobe)
! get external potential
call makevq(vq,nbasis,coorx,coory,coorz,bb,&
     cbb,nref,ncor,ngau,nlobe,vlobe)
dt=0.01
ct=-dt
do t=1,maxt
   ct=ct+dt
   write(*,'(A,1X,I8,1X,A,1X,ES16.8)') 'Iteration:',t,'time =',ct
   call update_h1(ct,h1,h1t,vq,nbasis)
   ! run an HF SCF cycle (only closed shell)
   call scf(jkmat,h1t,ovl,fock,ovu,ovx,ovs,pdens,&
        cmo,nbasis,nelec,hfen)
   ! old way to calculate energy (does restricted
   ! closed and open shell) no mos with this
   !call deten(cmo,h1,nbasis,coorx,coory,coorz,&
   !     bb,cbb,nref,ncor,ngau,nelec,enuc,hfen,&
   !     natom,nlobe,vlobe,jkmat)
   ! store HF energy
   hfenc=hfen+enuc
   ! get CI result up to exmax excitations if requested
   if (exmax/=0) then
      ! normalize GS CI coefficients
      !call normcic(civec,size(civec))
      civec=1.0d0
      ! unscaled contrib to CI energy
      hfen=hfen*(civec(1)**2)
      ! get CI Hamiltonian
      call fullci(coorx,coory,coorz,bb,cbb,nref,&
           ncor,ngau,nelec,natom,nlobe,vlobe,cmo,&
           h1t,nbasis,toten,exmax,civec,blokij,Hci)
      ! get All CI vectors
      Hci(1,1)=hfenc-enuc
      call diagci(Hci,fciv,cieig)
   end if
   ! print stuff
   call printdata(iter,ical,hfenc,hfen,toten,enuc,&
        orthp,bb,nao,ngau,cmo,nbasis,civec,1,ntype,&
        nref,ncor,nlobe,Hci,fciv,cieig)
end do
molprop=cieig(1)+enuc
end function molprop

! update one electron operator with vq
subroutine update_h1(t,h1,h1t,vq,nbasis)
implicit none
double precision,intent(in)::vq(:,:),h1(:,:)
double precision,intent(inout)::h1t(:,:)
integer,intent(in)::nbasis
double precision,intent(in)::t
double precision::tval
integer::i,j
tval=abs(sin(t))
do i=1,nbasis
   do j=1,nbasis
      h1t(i,j)=h1(i,j)+(tval*vq(i,j))
   end do
end do
end subroutine

! effect of external potential on basis
subroutine makevq(vq,nbasis,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nlobe,vlobe)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::vlobe(:,:)
double precision,intent(out)::vq(:,:)
integer,intent(in)::nlobe(:),nbasis
integer,intent(in)::nref(:),ncor(:),ngau(:)
integer::i,j
double precision::rq(3),vijrq,gexp,gcoef
rq(1)=0.0
rq(2)=0.0
rq(3)=2.0
gexp=1.0d0
gcoef=1.0d0
do i=1,nbasis
   do j=1,nbasis
      vijrq=0.0
      call repl3(i,j,rx,ry,rz,bb,cbb,gexp,&
           gcoef,rq,nref,ncor,ngau,nlobe,vlobe,&
           vijrq)
      vq(i,j)=vijrq
   end do
end do
end subroutine

! calculates repulsion integral between ao's
! with 'lobes' with a spherical gaussian 
! potential at rq
subroutine repl3(ao1,ao2,rx,ry,rz,bb,cbb,bb2,&
           cbb2,rq,nref,ncor,ngau,nlobe,vlobe,&
           vijrq)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::bb2,cbb2,rq(:)
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
integer,intent(in)::ao1,ao2
double precision,intent(out)::vijrq
integer::i,j,aoi,aoj,aci,acj
integer::ngi,ngj,ilob,jlob
integer::nli,nlj
double precision::pi,a,b,ax,ay,az,fact
double precision::bx,by,bz,c,d,cx,cy,cz
double precision::dx,dy,dz,anorm,bnorm,ab2
double precision::alpha,sab,cnorm,dnorm
double precision::cd2,beta,scd,px,py,pz
double precision::qx,qy,qz,gam,pq2,r,vabcd
double precision::signi,signj
double precision::dispa,dispb
pi=3.14159265358979d0
aoi=nref(ao1)
aoj=nref(ao2)
aci=ncor(ao1)
acj=ncor(ao2)
ngi=ngau(aoi)
ngj=ngau(aoj)
nli=nlobe(aoi)
nlj=nlobe(aoj)
c=bb2/2.0d0
cnorm=2.0d0*c/pi
cnorm=cnorm**0.75
cx=rq(1)
cy=rq(2)
cz=rq(3)
d=bb2/2.0d0
dnorm=2.0d0*d/pi
dnorm=dnorm**0.75
beta=c+d
dx=rq(1)
dy=rq(2)
dz=rq(3)
vijrq=0.0d0
do i=1,ngi
   a=bb(aoi,i)
   anorm=2.0d0*a/pi
   anorm=anorm**0.75
   dispa=0.03d0/dsqrt(a)
   do j=1,ngj
      b=bb(aoj,j)
      bnorm=2.0d0*b/pi
      bnorm=bnorm**0.75
      dispb=0.03d0/dsqrt(b)
      alpha=a+b
      signi=1.0d0
      do ilob=1,nli
         ax=rx(aci)+vlobe(ao1,1)*dispa*signi
         ay=ry(aci)+vlobe(ao1,2)*dispa*signi
         az=rz(aci)+vlobe(ao1,3)*dispa*signi
         signj=1.0d0
         do jlob=1,nlj
            bx=rx(acj)+vlobe(ao2,1)*dispb*signj
            by=ry(acj)+vlobe(ao2,2)*dispb*signj
            bz=rz(acj)+vlobe(ao2,3)*dispb*signj
            ab2=(ax-bx)**2+(ay-by)**2+(az-bz)**2
            sab=anorm*bnorm*dexp(-a*b*ab2/alpha)*(pi/alpha)**1.5
            px=(a*ax+b*bx)/alpha
            py=(a*ay+b*by)/alpha
            pz=(a*az+b*bz)/alpha
            cd2=(cx-dx)**2+(cy-dy)**2+(cz-dz)**2
            scd=cnorm*dnorm*dexp(-c*d*cd2/beta)*(pi/beta)**1.5
            qx=(c*cx+d*dx)/beta
            qy=(c*cy+d*dy)/beta
            qz=(c*cz+d*dz)/beta
            pq2=(px-qx)**2+(py-qy)**2+(pz-qz)**2
            gam=(1.0d0/alpha)+(1.0d0/beta)
            gam=1.0d0/gam
            r=gam*pq2
            fact=signi*signj
            vabcd=1.1283791670955d0*sab*scd*dsqrt(gam)*fun(r)
            vijrq=vijrq+cbb(aoi,i)*cbb(aoj,j)*cbb2&
            *cbb2*vabcd*fact
            signj=-signj
         end do
         signi=-signi
      end do
   end do
end do
end subroutine repl3

! scf procedure
subroutine scf(jkmat,h1,ovl,fock,ovu,ovx,ovs,pdens,&
           cmo,nbasis,nelec,hfen)
implicit none
integer,intent(in)::nbasis,nelec
double precision,intent(in)::jkmat(:,:,:,:),h1(:,:),&
ovl(:,:)
double precision,intent(out)::fock(:,:),ovu(:,:),&
ovx(:,:),ovs(:),pdens(:,:),hfen
double precision,intent(inout)::cmo(:,:)
double precision::thresh,olde
thresh=1.0E-7
olde=0.0d0
! make initial density matrix
call makepdens(cmo,pdens,nbasis,nelec)
do
   ! make fock matrix (restricted closed shell)
   call makefock(pdens,jkmat,fock,h1,nbasis)
   ! get HF MOs to define virtual space
   ! ovu, ovx is destroyed, cmo is remade
   call roothan(fock,ovl,ovu,ovx,ovs,cmo,nbasis)
   ! make density matrix
   call makepdens(cmo,pdens,nbasis,nelec)
   ! get HF energy
   call enfock(pdens,fock,h1,nbasis,hfen)
   if ((abs(hfen)-abs(olde))<thresh) then
      exit
   end if
   olde=hfen
end do
end subroutine

! Roothan Equation
subroutine roothan(F,S,U,X,E,C,N)
implicit none
double precision,intent(in)::S(:,:),F(:,:)
double precision,intent(out)::U(:,:),X(:,:),&
E(:),C(:,:)
integer,intent(in)::N
double precision::temp
integer::i,j,k
! diagonalize S
call diagci(S,U,E)
! make transformation matrix
do i=1,N
   if (E(i)<0.000001) then
      E(i)=0.0001
   end if
   X(:,i)=U(:,i)*(E(i)**-0.5d0)
end do
! left multiply by X'
do i=1,N
   do j=1,N
      temp=0.0d0
      do k=1,N
         temp=temp+X(k,i)*F(k,j)
      end do
      U(i,j)=temp
   end do
end do
! right multiply X'F by X
do i=1,N
   do j=1,N
      temp=0.0d0
      do k=1,N
         temp=temp+U(i,k)*X(k,j)
      end do
      C(i,j)=temp
   end do
end do
! diagonalize F'
call diagci(C,U,E)
! left multiply C' by X
do i=1,N
   do j=1,N
      temp=0.0d0
      do k=1,N
         temp=temp+X(i,k)*U(k,j)
      end do
      C(i,j)=temp
   end do
end do
end subroutine

! HF ground state energy from fock
! restricted (closed shell)
subroutine enfock(pdens,fock,h1,nbasis,toten)
implicit none
double precision,intent(in)::pdens(:,:)
double precision,intent(in)::fock(:,:)
double precision,intent(in)::h1(:,:)
double precision,intent(out)::toten
integer,intent(in)::nbasis
integer::u,v
toten=0.0d0
do u=1,nbasis
   do v=1,nbasis
      toten=toten+(pdens(u,v)*(h1(u,v)+fock(u,v)))
   end do
end do
toten=toten*0.5d0
end subroutine enfock

! prepare Fock matrix
subroutine makefock(pdens,jkmat,fock,&
           h1,nbasis)
implicit none
double precision,intent(in)::jkmat(:,:,:,:)
double precision,intent(in)::pdens(:,:),h1(:,:)
double precision,intent(out)::fock(:,:)
integer,intent(in)::nbasis
integer::u,v,r,l
double precision::puv,huv,e2uv,plr,jint,kint
do u=1,nbasis
   ! matrix is real symm
   do v=1,u
      puv=pdens(u,v)
      huv=h1(u,v)
      e2uv=0.0d0
      do l=1,nbasis
         do r=1,nbasis
            plr=pdens(l,r)
            jint=jkmat(u,v,l,r)
            kint=jkmat(u,r,l,v)
            e2uv=e2uv+(plr*(jint-0.50d0*kint))
         end do
      end do
      fock(u,v)=huv+e2uv
      fock(v,u)=fock(u,v)
   end do
end do
end subroutine makefock

! prepare SCF charge density matrix
subroutine makepdens(cmo,pdens,nbasis,nelec)
implicit none
double precision,intent(out)::pdens(:,:)
double precision,intent(in)::cmo(:,:)
integer,intent(in)::nbasis,nelec
integer::i,j,k,kmo
double precision::pij
do i=1,nbasis
   ! matrix is real symm
   do j=1,i
      pij=0.0d0
      do k=1,nelec
         kmo=(k+1)/2
         pij=pij+(cmo(i,kmo)*cmo(j,kmo))
      end do
      pdens(i,j)=pij
      pdens(j,i)=pij
   end do
end do
end subroutine makepdens

! Diagonalize square matrix A
! Fill square matrix C with eigenvectors
! and W with eigenvalues
subroutine diagci(A,C,W)
implicit none
double precision,intent(in)::A(:,:)
double precision,intent(out)::C(:,:)
double precision,intent(out)::W(:)
integer::INFO,LDA,LWORK,N
character(len=1)::UPLO,JOBZ
double precision,allocatable::WORK(:)
double precision::iniw(10)
N=size(A(:,1))
LDA=N
C=A
UPLO='U'
JOBZ='V'
LWORK=-1
call DSYEV( JOBZ, UPLO, N, C, LDA, W, iniw, LWORK, INFO )
LWORK=iniw(1)
allocate(WORK(LWORK))
call DSYEV( JOBZ, UPLO, N, C, LDA, W, WORK, LWORK, INFO )
if (info/=0) then
   write(*,*) 'routine dsyev failed.'
   write(*,*) 'Dont trust Eigenvecs!'
end if
deallocate(WORK)
end subroutine

! print iteration info
subroutine printdata(iter,ical,hfenc,hf,toten,enuc,&
           orthp,bb,nao,ngau,cmo,nbasis,civec,iprint,&
           ntype,nref,ncor,nlobe,Hci,fciv,cieig)
implicit none
integer,intent(in)::iter,ical,nao,nbasis,iprint,ngau(:)
integer,intent(in)::ntype(:),nref(:),ncor(:),nlobe(:)
double precision,intent(in)::hfenc,hf,toten,enuc,orthp
double precision,intent(in)::bb(:,:),cmo(:,:),civec(:)
double precision,intent(in)::Hci(:,:),fciv(:,:),cieig(:)
integer::i,j,bj,ej,ncol,labv(5),k,ncol2,cidim
character(len=100)::fmtstr
if (mod(ical,iprint)/=0) return
! Energy info for current iteration
write(*,*) ' '
write(*,'(A7,1X,ES15.8)') 'HF(GS):',hfenc
write(*,'(A,1X,ES15.8)') 'equal wieght CI(HF):',hfen
write(*,'(A,1X,ES15.8)') 'equal wieght CI(rest):',toten
write(*,'(A4,1X,ES15.8)') 'Vnn:',enuc
write(*,'(A8,1X,ES15.8)') 'E(Corr):',(cieig(1)+enuc)-hfenc
write(*,'(A8,1X,ES15.8)') 'Penalty:',orthp
write(*,*) ' '
! AO exponents and coeffs for current interation
write(*,*) 'AO expansion exponents and coefficients:'
do i=1,nao
   write(*,'(5(I5))') nref(i),ncor(i),nlobe(i),&
   ngau(i),ntype(i)
   fmtstr=''
   write(fmtstr,'(A,I4,A)') '(', ngau(i),'(ES16.8))'
   write(*,fmtstr) bb(i,1:ngau(i))
   write(*,fmtstr) cbb(i,1:ngau(i))
end do
write(*,*) ' '
! MO coefficients for current iteration
write(*,*) 'MO expansion coefficients:'
ncol=nbasis/5
do j=1,ncol
   bj=(j-1)*5+1
   ej=(bj-1)+5
   labv=-1
   do k=1,5
      labv(k)=bj-1+k
   end do
   write(*,'(5X,5(1X,I8,7X))') labv
   do i=1,nbasis
      write(*,'(I5,5(ES16.8))') i,cmo(i,bj:ej)
   end do
   write(*,*) ' '
end do
ncol2=mod(nbasis,5)
if (ncol2/=0) then
   bj=ncol*5+1
   ej=bj+ncol2-1
   labv=-1
   do k=1,ncol2
      labv(k)=bj-1+k
   end do
   write(*,'(5X,5(1X,I8,7X))') labv
   do i=1,nbasis
      write(*,'(I5,5(ES16.8))') i,cmo(i,bj:ej)
   end do
   write(*,*) ' '
end if
! Ground State CI expansion coefficients
!write(*,*) 'GS CI expansion coefficients:'
!ncol=size(civec)/5
!do j=1,ncol
!   bj=(j-1)*5+1
!   ej=(bj-1)+5
!   write(*,'(5(ES16.8))') civec(bj:ej)
!end do
!ncol2=mod(size(civec),5)
!if (ncol2/=0) then
!   bj=ncol*5+1
!   ej=bj+ncol2-1
!   write(*,'(5(ES16.8))') civec(bj:ej)
!end if
!write(*,*) ' '
! CI state energies
write(*,*) 'CI state energies:'
cidim=size(Hci(:,1))
ncol=cidim/5
do j=1,ncol
   bj=(j-1)*5+1
   ej=(bj-1)+5
   write(*,'(5(ES16.8))') cieig(bj:ej)+enuc
end do
ncol2=mod(cidim,5)
if (ncol2/=0) then
   bj=ncol*5+1
   ej=bj+ncol2-1
   write(*,'(5(ES16.8))') cieig(bj:ej)+enuc
end if
write(*,*) ' '
write(*,*) '---------------------------------------'
write(*,*) ' '
return
! CI Hamiltonian for current iteration
write(*,*) 'Full CI Hamiltonian:'
ncol=cidim/5
do j=1,ncol
   bj=(j-1)*5+1
   ej=(bj-1)+5
   labv=-1
   do k=1,5
      labv(k)=bj-1+k
   end do
   write(*,'(5X,5(1X,I8,7X))') labv
   do i=1,cidim
      write(*,'(I5,5(ES16.8))') i,Hci(i,bj:ej)
   end do
   write(*,*) ' '
end do
ncol2=mod(cidim,5)
if (ncol2/=0) then
   bj=ncol*5+1
   ej=bj+ncol2-1
   labv=-1
   do k=1,ncol2
      labv(k)=bj-1+k
   end do
   write(*,'(5X,5(1X,I8,7X))') labv
   do i=1,cidim
      write(*,'(I5,5(ES16.8))') i,Hci(i,bj:ej)
   end do
   write(*,*) ' '
end if
! CI coefficients for current iteration
write(*,*) 'Full CI Coefficients:'
ncol=cidim/5
do j=1,ncol
   bj=(j-1)*5+1
   ej=(bj-1)+5
   labv=-1
   do k=1,5
      labv(k)=bj-1+k
   end do
   write(*,'(5X,5(1X,I8,7X))') labv
   do i=1,cidim
      write(*,'(I5,5(ES16.8))') i,fciv(i,bj:ej)
   end do
   write(*,*) ' '
end do
ncol2=mod(cidim,5)
if (ncol2/=0) then
   bj=ncol*5+1
   ej=bj+ncol2-1
   labv=-1
   do k=1,ncol2
      labv(k)=bj-1+k
   end do
   write(*,'(5X,5(1X,I8,7X))') labv
   do i=1,cidim
      write(*,'(I5,5(ES16.8))') i,fciv(i,bj:ej)
   end do
   write(*,*) ' '
end if
write(*,*) '---------------------------------------'
write(*,*) ' '
end subroutine printdata

! store 2e ints
subroutine get_jkmat(jkmat,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nbasis,nlobe,&
           vlobe)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:),nbasis
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
double precision,intent(out)::jkmat(:,:,:,:)
double precision::vijkl
integer::i,j,k,l
do i=1,nbasis
   do j=1,nbasis
      do k=1,nbasis
         do l=1,nbasis
            call repl2(i,j,k,l,rx,ry,rz,bb,cbb,&
                 nref,ncor,ngau,nlobe,vlobe,vijkl)
            jkmat(i,j,k,l)=vijkl
         end do
      end do
   end do
end do
end subroutine get_jkmat

! gram schmidt orth for mos
subroutine gramschmidt(A,n,ovl)
implicit none
double precision,intent(inout)::A(:,:)
double precision,intent(in)::ovl(:,:)
integer,intent(in)::n
integer::i,j,k,l
double precision::sij,ck,cl
do i=1,n
   do j=1,i-1
      sij=0.0d0
      do k=1,n
         ck=A(k,i)
         do l=1,n
            cl=A(l,j)
            sij=sij+ck*cl*ovl(k,l)
         end do
      end do
      A(:,i)=A(:,i)-sij*A(:,j)
   end do
end do
end subroutine

! normalize CI coefficients
subroutine normcic(civec,n)
implicit none
double precision,intent(inout)::civec(:)
integer,intent(in)::n
integer::i
double precision::fac
fac=0.0d0
do i=1,n
   fac=fac+civec(i)**2.0d0
end do
fac=1.0d0/(fac**0.5d0)
do i=1,n
   civec(i)=civec(i)*fac
end do
end subroutine normcic

! apply new parameters
subroutine applyx(x,svar,bb,ngau,nao,nbasis,&
           cmo,nvdet,civec,bbc,civec2,cmo2,cbbc)
implicit none
double precision,intent(in)::x(:),bbc(:,:),cbbc(:,:)
double precision,intent(in)::civec2(:),cmo2(:,:)
double precision,intent(inout)::svar(:),civec(:)
double precision,intent(inout)::bb(:,:),cmo(:,:)
integer,intent(in)::nvdet(:),ngau(:),nao,nbasis
integer::i,j,xvar
! reset to original
bb(:,:)=bbc(:,:)
cbb(:,:)=cbbc(:,:)
!cmo=cmo2
!civec=civec2
! scaling parameter
xvar=0
do i=1,nao
   xvar=xvar+1
   svar(i)=svar(i)*(1.0d0+x(xvar))
   !call scaleao(i,svar(i),bb,ngau)
end do
! ao exponents and coefficients
do i=1,nao
   ! exponents
   do j=1,ngau(i)
      xvar=xvar+1
     !if (i>2.and.j>0) then
     bb(i,j)=bb(i,j)*abs(1.0d0+x(xvar))+1.0e-12
     !end if
   end do
   ! coefficients
   do j=1,ngau(i)
      xvar=xvar+1
     !if (i>2.and.j>0) then
     cbb(i,j)=cbb(i,j)+x(xvar)+1.0e-12
     !end if
   end do
end do
! mo coefficients or density matrix
do i=1,nbasis
   do j=1,nbasis
      xvar=xvar+1
      !cmo(i,j)=cmo(i,j)+x(xvar)
   end do
end do
! ci coefficients
do i=1,size(civec)
   xvar=xvar+1
   !civec(i)=civec(i)+x(xvar)+1.0e-12
end do
end subroutine applyx

subroutine orthopen(cmo,rx,ry,rz,bb,cbb,nref,&
           ncor,ngau,ovl,nbasis,orthp,fac)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::ovl(:,:),cmo(:,:)
double precision,intent(in)::fac
double precision,intent(out)::orthp
integer,intent(in)::nbasis
double precision::smo
integer::i,j
orthp=0.0d0
do i=1,nbasis
   do j=1,nbasis
      if (j>i) then
         smo=0.0d0
         call mover(i,j,cmo,coorx,coory,coorz,bb,&
         cbb,nref,ncor,ngau,ovl,nbasis,smo)
         orthp=orthp+abs(smo)
      end if
   end do
end do
orthp=orthp*fac
end subroutine orthopen

! Computes CI energy up to any ex level
! needs HF energy
subroutine fullci(rx,ry,rz,bb,cbb,nref,ncor,&
           ngau,nelec,natom,nlobe,vlobe,cmo,&
           h1,nbasis,toten,exmax,civec,&
           blokij,Hci)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::h1(:,:)
double precision,intent(in)::cmo(:,:)
double precision,intent(in)::vlobe(:,:)
double precision,intent(in)::civec(:)
integer,intent(in)::nlobe(:),natom,nbasis
integer,intent(in)::nelec,exmax,blokij(:)
double precision,intent(inout)::Hci(:,:)
double precision,intent(out)::toten
double precision::cisen,fcien,cisen2
toten=0.0d0
cisen=0.0d0
cisen2=0.0d0
fcien=0.0d0
if (exmax>1) then
   ! first the singles block
   call ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
        nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
        cisen,1,1,civec,blokij,Hci,jkmat)
   write(*,'(A2,I3,A1,I3,A2,ES15.8)') 'E(',1,',',1,'):',cisen
   ! singles with ground (should be zero)
   !call ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
   !     nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
   !     cisen2,0,1,civec,blokij)
   ! then all other blocks
   call cidnup(rx,ry,rz,bb,cbb,vlobe,nref,ncor,&
        ngau,nlobe,nbasis,nelec,natom,cmo,h1,&
        fcien,exmax,civec,blokij,Hci,jkmat)
elseif (exmax==1) then
   ! just CIS
   call ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
        nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
        cisen,1,1,civec,blokij,Hci,jkmat)
   write(*,'(A2,I3,A1,I3,A2,ES15.8)') 'E(',1,',',1,'):',cisen
   ! singles with ground (should be zero)
   !call ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
   !     nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
   !     cisen2,0,1,civec,blokij)
else
   cisen=0.0d0
   cisen2=0.0d0
   fcien=0.0d0
end if
toten=cisen+2.0*cisen2+fcien
end subroutine fullci

! Takes care of all CI blocks except singles block
! and HF block. Excitation levels go up to exmax.
subroutine cidnup(rx,ry,rz,bb,cbb,vlobe,nref,ncor,&
           ngau,nlobe,nbasis,nelec,natom,cmo,h1,&
           en,exmax,civec,blokij,Hci,jkmat)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::jkmat(:,:,:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::h1(:,:),cmo(:,:)
double precision,intent(in)::vlobe(:,:),civec(:)
integer,intent(in)::nlobe(:),natom,nbasis,nelec
integer,intent(in)::blokij(:),exmax
double precision,intent(inout)::Hci(:,:)
double precision,intent(out)::en
double precision::en1,en2,en3,eii,eij,eij2
integer::i
en=0.0d0
eii=0.0d0
eij=0.0d0
eij2=0.0d0
do i=2,exmax
   en1=0.0d0
   en2=0.0d0
   en3=0.0d0
   ! diagonal CI block for ex level i
   call ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
        nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
        en1,i,i,civec,blokij,Hci,jkmat)
   ! off diagonal block btw ex levi i and i-1
   call ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
        nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
        en2,i-1,i,civec,blokij,Hci,jkmat)
   ! off diagonal block btw ex levi i and i-2
   call ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
        nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
        en3,i-2,i,civec,blokij,Hci,jkmat)
   eii=eii+en1
   eij=eij+2.0*en2
   eij2=eij2+2.0*en3
   write(*,'(A2,I3,A1,I3,A2,ES15.8)') 'E(',i,',',i,'):',en1
   write(*,'(A2,I3,A1,I3,A2,ES15.8)') 'E(',i-1,',',i,'):',en2
   write(*,'(A2,I3,A1,I3,A2,ES15.8)') 'E(',i-2,',',i,'):',en3
end do
en=eii+eij+eij2
end subroutine cidnup

! initializes the ci coeff vector and other stuff
subroutine initci(v,nv,n,nelec,nbasis,tdet)
implicit none
integer,intent(out)::v(:),nv(:),tdet
integer,intent(in)::n,nelec,nbasis
integer::exocc,exvir,i,nso,nocc,nvir,ndet
double precision::vmag,normci
nso=nbasis*2
nocc=nelec
nvir=nso-nocc
tdet=1
v(1)=tdet
do i=1,n
   v(i+1)=tdet+1
   exocc=factorial(nocc)
   exocc=exocc/(factorial(i)*factorial(nocc-i))
   exvir=factorial(nvir)
   exvir=exvir/(factorial(i)*factorial(nvir-i))
   ndet=exocc*exvir
   nv(i)=ndet
   tdet=tdet+ndet
end do
end subroutine initci

! Calculates energy contribution from <EL1|H|EL2>
! CI block. EL(1 or 2) can be S, D, T, Q, etc.
subroutine ciblock(cmo,h1,nbasis,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nelec,natom,nlobe,vlobe,&
           en,nex1,nex2,civec,startis,Hci,jkmat)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::h1(:,:),cmo(:,:)
double precision,intent(in)::vlobe(:,:),civec(:)
double precision,intent(in)::jkmat(:,:,:,:)
double precision,intent(inout)::Hci(:,:)
integer,intent(in)::nlobe(:),natom,nbasis
integer,intent(in)::nelec,nex1,nex2,startis(:)
double precision,intent(out)::en
double precision::tempe,cij
integer::nso,nocc,nvir,ci,cj
integer,allocatable::oarr1(:),oarr2(:)
integer,allocatable::varr1(:),varr2(:)
integer::ndiff,i,j,k,l,so1,so2,so3,so4
integer::exocc1,exocc2,exvir1,exvir2
integer::sodiff(2,2)
logical::odone1,odone2,vdone1,vdone2,undo
integer,allocatable::sdeti(:),sdetj(:)
! number of spin, occ, and virt orbtials
nso=nbasis*2
nocc=nelec
nvir=nso-nocc
! number of occ combos for ex type 1
exocc1=factorial(nocc)
exocc1=exocc1/(factorial(nex1)*factorial(nocc-nex1))
! number of vir combos for ex type 1
exvir1=factorial(nvir)
exvir1=exvir1/(factorial(nex1)*factorial(nvir-nex1))
! number of occ combos for ex type 2
exocc2=factorial(nocc)
exocc2=exocc2/(factorial(nex2)*factorial(nocc-nex2))
! number of vir combos for ex type 2
exvir2=factorial(nvir)
exvir2=exvir2/(factorial(nex2)*factorial(nvir-nex2))
! CI coeff indeces
allocate(sdeti(nelec),sdetj(nelec))
allocate(oarr1(nex1),varr1(nex1))
allocate(oarr2(nex2),varr2(nex2))
en=0.0d0
! loop thru bra dets for excitation type 1
ci=startis(nex1+1)-1
odone1=.true.
do i=1,exocc1
   call comb_next(nocc,nex1,oarr1,odone1)
   vdone1=.true.
   do j=1,exvir1
      call comb_next(nvir,nex1,varr1,vdone1)
      ci=ci+1
      ! sdeti is bra
      call excitedet(sdeti,nelec,oarr1,varr1,nex1)
      call sort(sdeti,nelec)
      ! loop thru ket dets of ex type 2
      cj=startis(nex2+1)-1
      odone2=.true.
      do k=1,exocc2
         call comb_next(nocc,nex2,oarr2,odone2)
         vdone2=.true.
         do l=1,exvir2
            call comb_next(nvir,nex2,varr2,vdone2)
            cj=cj+1
            cij=civec(ci)*civec(cj)
            ! sdetj is ket
            call excitedet(sdetj,nelec,oarr2,varr2,nex2)
            call sort(sdetj,nelec)
            ! how many orbitals are different?
            call detdiff(sdeti,sdetj,nelec,ndiff,sodiff)
            ! zero SO difference
            if (ndiff==0) then
               call deten_gen(sdetj,cmo,h1,nbasis,rx,&
                    ry,rz,bb,cbb,nref,ncor,ngau,nelec,&
                    natom,nlobe,vlobe,tempe,jkmat)
               en=en+cij*tempe
               Hci(ci,cj)=tempe
               Hci(cj,ci)=tempe
               !write(*,'(A,I5,A,I5,A)') 'Completed CI element H(',ci,',',cj,').'
            ! one SO difference
            elseif (ndiff==1) then
               call deten_1diff(sodiff,sdeti,cmo,h1,rx,ry,rz,&
                    bb,cbb,vlobe,nref,ncor,ngau,nlobe,nbasis,&
                    nelec,natom,tempe,jkmat)
               en=en+cij*tempe
               Hci(ci,cj)=tempe
               Hci(cj,ci)=tempe
               !write(*,'(A,I5,A,I5,A)') 'Completed CI element H(',ci,',',cj,').'
            ! two SO difference
            elseif (ndiff==2) then
               call deten_2diff(sodiff,cmo,rx,ry,rz,bb,cbb,&
                    vlobe,nref,ncor,ngau,nlobe,nelec,nbasis,&
                    natom,tempe,jkmat)
               en=en+cij*tempe
               Hci(ci,cj)=tempe
               Hci(cj,ci)=tempe
               !write(*,'(A,I5,A,I5,A)') 'Completed CI element H(',ci,',',cj,').'
            end if
         end do
      end do
   end do
end do
deallocate(oarr1,varr1)
deallocate(oarr2,varr2)
deallocate(sdeti,sdetj)
end subroutine         

! int energy of 2 dts with 2 diff so's
! (restricted)
subroutine deten_2diff(sodiff,cmo,rx,ry,rz,bb,cbb,&
           vlobe,nref,ncor,ngau,nlobe,nelec,nbasis,&
           natom,en2e,jkmat)
implicit none
integer,intent(in)::sodiff(2,2)
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::jkmat(:,:,:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::vlobe(:,:)
double precision,intent(in)::cmo(:,:)
integer,intent(in)::nlobe(:),nbasis,nelec,natom
double precision,intent(out)::en2e
double precision::j2e,k2e,spik,spjl,spil,spjk
integer::iso,jso,kso,lso,imo,jmo,kmo,lmo
iso=sodiff(1,1)
jso=sodiff(2,1)
kso=sodiff(1,2)
lso=sodiff(2,2)
imo=(iso+1)/2
jmo=(jso+1)/2
kmo=(kso+1)/2
lmo=(lso+1)/2
! get repulsion integrals
call rrepl(imo,kmo,lmo,jmo,rx,ry,rz,bb,&
     cbb,nref,ncor,ngau,j2e,k2e,nbasis,&
     cmo,nlobe,vlobe,jkmat)
! spin integrals
spik=dble(abs(mod((iso+kso),2)-1))
spjl=dble(abs(mod((jso+lso),2)-1))
spil=dble(abs(mod((iso+lso),2)-1))
spjk=dble(abs(mod((jso+kso),2)-1))
en2e=(j2e*spik*spjl)-(k2e*spil*spjk)
end subroutine

! interaction energy of 2 dets with 1 diff
! (restricted)
subroutine deten_1diff(sodiff,so,cmo,h1,rx,ry,rz,&
           bb,cbb,vlobe,nref,ncor,ngau,nlobe,nbasis,&
           nelec,natom,energy,jkmat)
implicit none
integer,intent(in)::so(:),sodiff(2,2)
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::vlobe(:,:),jkmat(:,:,:,:)
double precision,intent(in)::cmo(:,:),h1(:,:)
integer,intent(in)::nlobe(:),nbasis,nelec,natom
double precision,intent(out)::energy
integer::iso,jso,imo,jmo,g,gso,gmo
double precision::en1e,j2e,k2e,en2e
double precision::spij,spig,spjg
energy=0.0d0
! restricted HF indexing:
! odds are alpha, evens beta
iso=sodiff(1,1)
jso=sodiff(1,2)
imo=(iso+1)/2
jmo=(jso+1)/2
!  kinetic and attaction interaction energy
call tv1e_gen(imo,jmo,cmo,h1,en1e,nbasis)
! spin integral for 1e part
spij=dble(abs(mod((iso+jso),2)-1))
en1e=en1e*spij
! repulsion energy
en2e=0.0d0
do g=1,nelec
   ! restricted HF indexing
   gso=so(g)
   gmo=(gso+1)/2
   ! get repulsion integrals
   call rrepl(imo,jmo,gmo,gmo,rx,ry,rz,bb,&
        cbb,nref,ncor,ngau,j2e,k2e,nbasis,&
        cmo,nlobe,vlobe,jkmat)
   ! spin integrals
   spig=dble(abs(mod((iso+gso),2)-1))
   spjg=dble(abs(mod((jso+gso),2)-1))
   ! add to total repulsion
   en2e=en2e+(j2e*spij)-(k2e*spig*spjg)
end do
energy=en1e+en2e
end subroutine

! determine differences between 2 determinants
! warning: subroutine assumes both v1 and v2
! have been sorted in ascending order.
subroutine detdiff(v1,v2,n,ndif,diffs)
implicit none
integer,intent(in)::v1(:),v2(:),n
integer,intent(out)::diffs(2,2),ndif
integer::i
ndif=0
do i=1,n
   if (v1(i)/=v2(i)) then
      ndif=ndif+1
      if (ndif>2) then
         exit
      end if
      diffs(ndif,1)=v1(i)
      diffs(ndif,2)=v2(i)
   end if
end do
end subroutine detdiff

! create excited determinant with m
! excitations
subroutine excitedet(det,n,occ,vir,m)
implicit none
integer,intent(in)::occ(:),vir(:),n,m
integer,intent(out)::det(:)
integer::i
! fill electrons
do i=1,n
   det(i)=i
end do
! excite electrons
do i=1,m
   det(occ(i))=vir(i)+n
end do
end subroutine excitedet

! Total energy of any single determinant
! needs integer vector of spin orbital indeces
! (restricted)
subroutine deten_gen(so,cmo,h1,nbasis,rx,&
           ry,rz,bb,cbb,nref,ncor,ngau,nelec,&
           natom,nlobe,vlobe,dete,jkmat)
implicit none
integer,intent(in)::so(:)
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::vlobe(:,:),jkmat(:,:,:,:)
integer,intent(in)::nlobe(:)
double precision,intent(in)::cmo(:,:),h1(:,:)
integer,intent(in)::nbasis,nelec,natom
double precision,intent(out)::dete
double precision::e1el,e2en,tvmo,j2e,k2e
integer::i,j,imo,jmo,spn
dete=0.0d0
e1el=0.0d0
e2en=0.0d0
do i=1,nelec
   ! restricted HF indexing:
   ! odds are alpha, evens beta
   imo=(so(i)+1)/2
   ! electronic kinetic and attaction energy
   call tv1e(imo,cmo,h1,tvmo,nbasis)
   e1el=e1el+tvmo
   ! electronic repulsion energy
   do j=1,nelec
      ! restricted HF indexing
      jmo=(so(j)+1)/2
      call rrepl(imo,imo,jmo,jmo,rx,ry,rz,bb,&
           cbb,nref,ncor,ngau,j2e,k2e,nbasis,&
           cmo,nlobe,vlobe,jkmat)
      ! get spin part for exchange
      spn=abs(mod((so(i)+so(j)),2)-1)
      e2en=e2en+(j2e-(k2e*dble(spn)*dble(spn)))
   end do
end do
! single determinant energy
dete=e1el+0.5d0*e2en
end subroutine deten_gen

! testing combinations routine
subroutine test()
implicit none
integer::k,n,t,i
logical::done
integer,allocatable::c(:)
n=10
k=0
t=factorial(n)
t=t/(factorial(k)*factorial(n-k))
allocate(c(k))
c=0
done=.true.
do i=1,t+1
 call comb_next(n,k,c,done)
 write(*,*) c(:),done
end do
deallocate(c)
end subroutine test

! get overlap between 2 mos
subroutine mover(imo,jmo,cmo,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,ovl,nbasis,smo)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::ovl(:,:),cmo(:,:)
double precision,intent(out)::smo
integer,intent(in)::imo,jmo,nbasis
double precision::sij,ci,cj
integer::i,j
smo=0.0d0
do i=1,nbasis
   ci=cmo(i,imo)
   do j=1,nbasis
      cj=cmo(j,jmo)
      sij=ovl(i,j)
      smo=smo+(ci*cj*sij)
   end do
end do
end subroutine mover

! normalize mos
subroutine normmos(cmo,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,ovl,nbasis)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::ovl(:,:)
double precision,intent(inout)::cmo(:,:)
integer,intent(in)::nbasis
integer::imo
do imo=1,nbasis
   call nc_mo(imo,cmo,rx,ry,rz,bb,cbb,&
        nref,ncor,ngau,ovl,nbasis)
end do
end subroutine normmos

! get HF energy
subroutine deten(cmo,h1,nbasis,rx,ry,rz,bb,&
           cbb,nref,ncor,ngau,nelec,enuc,een,&
           natom,nlobe,vlobe,jkmat)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::vlobe(:,:)
double precision,intent(in)::jkmat(:,:,:,:)
integer,intent(in)::nlobe(:)
double precision,intent(in)::cmo(:,:)
double precision,intent(in)::enuc,h1(:,:)
double precision,intent(out)::een
integer,intent(in)::nbasis,nelec,natom
double precision::e1el,e2en
! get 1e energy
call eh1(cmo,h1,nbasis,rx,ry,rz,bb,cbb,nref,&
     ncor,ngau,natom,nelec,e1el)
! get 2e energies
call erepl(e2en,rx,ry,rz,bb,cbb,nref,&
     ncor,ngau,nbasis,cmo,nelec,nlobe,&
     vlobe,jkmat)
!een=enuc+e1el+(e2ej-e2ek)
een=enuc+e1el+e2en
end subroutine deten

! calculate HF 1e energy
subroutine eh1(cmo,h1,nbasis,rx,ry,rz,bb,&
           cbb,nref,ncor,ngau,natom,nelec,&
           e1el)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::cmo(:,:)
double precision,intent(in)::h1(:,:)
integer,intent(in)::nbasis,nelec,natom
double precision,intent(out)::e1el
integer::el,imo
double precision::tvmo
e1el=0.0d0
imo=0
do el=1,nelec
   if (mod(el,2)/=0) then
      imo=imo+1
      tvmo=0.0d0
      ! kinetic energy and attaction
      call tv1e(imo,cmo,h1,tvmo,nbasis)
      e1el=e1el+tvmo
   else
      e1el=e1el+tvmo
   endif
end do
end subroutine eh1

! get kinetic and attaction energy btw 2 mos
subroutine tv1e_gen(imo,jmo,cmo,h1,tvmo,nbasis)
implicit none
double precision,intent(in)::cmo(:,:)
double precision,intent(in)::h1(:,:)
integer,intent(in)::nbasis,imo,jmo
double precision,intent(out)::tvmo
integer::i,j
double precision::ci,cj
tvmo=0.0d0
do i=1,nbasis
   ci=cmo(i,imo)
   do j=1,nbasis
      cj=cmo(j,jmo)
      tvmo=tvmo+(ci*cj*h1(i,j))
   end do
end do
end subroutine tv1e_gen

! get kinetic and attaction energy for an mo
subroutine tv1e(imo,cmo,h1,tvmo,nbasis)
implicit none
double precision,intent(in)::cmo(:,:)
double precision,intent(in)::h1(:,:)
double precision,intent(out)::tvmo
integer,intent(in)::nbasis,imo
integer::i,j
double precision::ci,cj
tvmo=0.0d0
do i=1,nbasis
   ci=cmo(i,imo)
   do j=1,nbasis
      cj=cmo(j,imo)
      tvmo=tvmo+(ci*cj*h1(i,j))
   end do
end do
end subroutine tv1e

! normalize mo imo. need overlap and mo coeffs
subroutine nc_mo(imo,cmo,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,ovl,nbasis)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::ovl(:,:)
double precision,intent(inout)::cmo(:,:)
integer,intent(in)::imo,nbasis
double precision::nc,smo,sij,ci,cj
integer::i,j
smo=0.0d0
do i=1,nbasis
   ci=cmo(i,imo)
   do j=1,nbasis
      cj=cmo(j,imo)
      sij=ovl(i,j)
      smo=smo+(ci*cj*sij)
   end do
end do
nc=1.0d0/(smo**0.50d0)
do i=1,nbasis
  cmo(i,imo)=cmo(i,imo)*nc
end do
end subroutine nc_mo

! get kinetic matrix
subroutine matkin(kin,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nbasis,nlobe,vlobe)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
integer,intent(in)::nbasis
double precision,intent(out)::kin(:,:)
double precision::tij
integer::i,j
kin=0.0d0
do i=1,nbasis
   do j=1,nbasis
      call t_chi(i,j,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,tij,nlobe,vlobe)
      kin(i,j)=tij
   end do
end do
end subroutine matkin

! get overlap matrix
subroutine matovl(ovl,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nbasis,nlobe,vlobe)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
integer,intent(in)::nbasis
double precision,intent(out)::ovl(:,:)
double precision::sij
integer::i,j
ovl=0.0d0
do i=1,nbasis
   do j=1,nbasis
      call ovlp_chi(i,j,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nlobe,vlobe,sij)
      ovl(i,j)=sij
   end do
end do
end subroutine matovl

! scales an aotype with s
subroutine scaleao(ao,s,bb,ngau)
implicit none
integer,intent(in)::ao,ngau(:)
double precision,intent(in)::s
double precision,intent(inout)::bb(:,:)
integer::i
do i=1,ngau(ao)
   bb(ao,i)=bb(ao,i)*(s**2.0d0)
end do
end subroutine scaleao

! subroutine normalizes all types of ao's
subroutine normao(nao,coorx,coory,coorz,&
           bb,cbb,nref,ncor,ngau,nlobe,vlobe)
implicit none
double precision,intent(inout)::coorx(:),&
coory(:),coorz(:),bb(:,:),cbb(:,:),vlobe(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:),nao
integer,intent(in)::nlobe(:)
integer::i,j
double precision::sao,nc
do i=1,nao
   sao=0.0
   call ovlp_chi(i,i,coorx,coory,coorz,bb,cbb,&
        nref,ncor,ngau,nlobe,vlobe,sao)
   nc=1.0d0/(sao**0.5d0)
   do j=1,ngau(i)
      cbb(i,j)=nc*cbb(i,j)
   end do
end do
end subroutine normao

! calculates overlap between two ao's
subroutine ovlp_chi(ao1,ao2,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nlobe,vlobe,sij)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
integer,intent(in)::ao1,ao2,nref(:),ncor(:),ngau(:)
double precision,intent(out)::sij
integer::i,j,aoi,aoj,aci,acj,ngi,ngj
integer::ilob,jlob,nli,nlj
double precision::pi,a,b,ax,ay,az
double precision::bx,by,bz,anorm,bnorm
double precision::ab2,alpha,sab
double precision::signi,signj,dispa,dispb 
pi=3.14159265358979d0
aoi=nref(ao1)
aoj=nref(ao2)
aci=ncor(ao1)
acj=ncor(ao2)
ngi=ngau(aoi)
ngj=ngau(aoj)
nli=nlobe(aoi)
nlj=nlobe(aoj)
sij=0.0d0
do i=1,ngi
   a=bb(aoi,i)
   anorm=2.0d0*a/pi
   anorm=anorm**0.75
   dispa=0.03d0/dsqrt(a)
   do j=1,ngj
      b=bb(aoj,j)
      bnorm=2.0d0*b/pi
      bnorm=bnorm**0.75
      dispb=0.03d0/dsqrt(b)
      alpha=a+b
      signi=1.0d0
      do ilob=1,nli
         ax=rx(aci)+vlobe(ao1,1)*dispa*signi
         ay=ry(aci)+vlobe(ao1,2)*dispa*signi
         az=rz(aci)+vlobe(ao1,3)*dispa*signi
         signj=1.0d0
         do jlob=1,nlj
            bx=rx(acj)+vlobe(ao2,1)*dispb*signj
            by=ry(acj)+vlobe(ao2,2)*dispb*signj
            bz=rz(acj)+vlobe(ao2,3)*dispb*signj
            ab2=(ax-bx)**2+(ay-by)**2+(az-bz)**2
            sab=anorm*bnorm*dexp(-a*b*ab2/alpha)*(pi/alpha)**1.5
            sij=sij+cbb(aoi,i)*cbb(aoj,j)*sab*signi*signj
            signj=-signj
         end do
         signi=-signi
      end do
   end do
end do
end subroutine ovlp_chi

! calculates kinetic between 2 ao's
subroutine t_chi(ao1,ao2,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,tij,nlobe,vlobe)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:),vlobe(:,:)
integer,intent(in)::ao1,ao2,nref(:),ncor(:),ngau(:)
integer,intent(in)::nlobe(:)
double precision,intent(out)::tij
integer::i,j,aoi,aoj,aci,acj,ngi,ngj
integer::ilob,jlob,nli,nlj
double precision::pi,a,b,ax,ay,az
double precision::bx,by,bz,anorm,bnorm
double precision::ab2,alpha,sab,tab
double precision::signi,signj,dispa,dispb 
pi=3.14159265358979d0
aoi=nref(ao1)
aoj=nref(ao2)
aci=ncor(ao1)
acj=ncor(ao2)
ngi=ngau(aoi)
ngj=ngau(aoj)
nli=nlobe(aoi)
nlj=nlobe(aoj)
tij=0.0d0
do i=1,ngi
   a=bb(aoi,i)
   anorm=2.0d0*a/pi
   anorm=anorm**0.75
   dispa=0.03d0/dsqrt(a)
   do j=1,ngj
      b=bb(aoj,j)
      bnorm=2.0d0*b/pi
      bnorm=bnorm**0.75
      dispb=0.03d0/dsqrt(b)
      alpha=a+b
      signi=1.0d0
      do ilob=1,nli
         ax=rx(aci)+vlobe(ao1,1)*dispa*signi
         ay=ry(aci)+vlobe(ao1,2)*dispa*signi
         az=rz(aci)+vlobe(ao1,3)*dispa*signi
         signj=1.0d0
         do jlob=1,nlj
            bx=rx(acj)+vlobe(ao2,1)*dispb*signj
            by=ry(acj)+vlobe(ao2,2)*dispb*signj
            bz=rz(acj)+vlobe(ao2,3)*dispb*signj
            ab2=(ax-bx)**2+(ay-by)**2+(az-bz)**2
            sab=anorm*bnorm*dexp(-a*b*ab2/alpha)*(pi/alpha)**1.5
            tab=sab*(a*b/alpha)*(3.0d0-2.0d0*a*b*ab2/alpha)
            tij=tij+cbb(aoi,i)*cbb(aoj,j)*tab*signi*signj
            signj=-signj
         end do
         signi=-signi
      end do
   end do
end do
end subroutine t_chi

! add attraction to 1e atomic matrix
subroutine attr_h1(chg,rx,ry,rz,bb,cbb,nref,ncor,ngau,&
           nbasis,natom,h1,nlobe,vlobe)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:),vlobe(:,:)
double precision,intent(in)::bb(:,:),cbb(:,:),chg(:)
integer,intent(in)::nref(:),ncor(:),ngau(:),nlobe(:)
integer,intent(in)::natom,nbasis
double precision,intent(inout)::h1(:,:)
integer::i,j,k,ii,jj,ao1,ao2,aoi,aoj,aci,acj,ngi,ngj
integer::nuc,ilob,jlob,nli,nlj
double precision::pi,a,b,ax,ay,az,ex,ey,ez,px,py,pz
double precision::bx,by,bz,anorm,bnorm,ep2
double precision::ab2,alpha,sab,atij,akij,r,ichg
double precision::dispa,dispb,signi,signj 
pi=3.14159265358979d0
do i=1,nbasis
   ao1=i
   aoi=nref(ao1)
   aci=ncor(ao1)
   ngi=ngau(aoi)
   nli=nlobe(aoi)
   do j=1,nbasis
      ao2=j
      aoj=nref(ao2)
      acj=ncor(ao2)
      ngj=ngau(aoj)
      nlj=nlobe(aoj)
      akij=0.0d0
      do k=1,natom
         nuc=k
         ichg=chg(nuc)
         ex=rx(nuc)
         ey=ry(nuc)
         ez=rz(nuc)
         atij=0.0d0
         do ii=1,ngi
            a=bb(aoi,ii)
            anorm=2.0d0*a/pi
            anorm=anorm**0.75
            dispa=0.03d0/dsqrt(a)
            do jj=1,ngj
               b=bb(aoj,jj)
               bnorm=2.0d0*b/pi
               bnorm=bnorm**0.75
               alpha=a+b
               dispb=0.03d0/dsqrt(b)
               signi=1.0d0
               do ilob=1,nli
                  ax=rx(aci)+vlobe(ao1,1)*dispa*signi
                  ay=ry(aci)+vlobe(ao1,2)*dispa*signi
                  az=rz(aci)+vlobe(ao1,3)*dispa*signi
                  signj=1.0d0
                  do jlob=1,nlj
                     bx=rx(acj)+vlobe(ao2,1)*dispb*signj
                     by=ry(acj)+vlobe(ao2,2)*dispb*signj
                     bz=rz(acj)+vlobe(ao2,3)*dispb*signj
                     ab2=(ax-bx)**2+(ay-by)**2+(az-bz)**2
                     sab=anorm*bnorm*dexp(-a*b*ab2/alpha)*(pi/alpha)**1.5
                     px=(a*ax+b*bx)/alpha
                     py=(a*ay+b*by)/alpha
                     pz=(a*az+b*bz)/alpha
                     ep2=(ex-px)**2+(ey-py)**2+(ez-pz)**2
                     r=ep2*alpha
                     atij=atij-(cbb(aoi,ii)*cbb(aoj,jj)*&
                     1.1283791670955d0*sab*(alpha**0.5d0)*fun(r))*&
                     signi*signj
                     signj=-signj
                  end do
                  signi=-signi
               end do
            end do
         end do
         akij=akij+atij*ichg
      end do
      h1(i,j)=h1(i,j)+akij
   end do
end do
end subroutine attr_h1

! calculate the 1e expectation value
subroutine vns1e(chg,rx,ry,rz,bb,cbb,nref,ncor,ngau,&
           nbasis,natom,at1e,imo,cmo,nlobe,vlobe)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:),chg(:)
integer,intent(in)::nref(:),ncor(:),ngau(:),imo
integer,intent(in)::natom,nbasis,nlobe(:)
double precision,intent(in)::cmo(:,:),vlobe(:,:)
double precision,intent(out)::at1e
integer::i,j,k,ii,jj,ao1,ao2,aoi,aoj,aci,acj,ngi,ngj
integer::nuc
double precision :: pi,a,b,ax,ay,az,ex,ey,ez,px,py,pz
double precision :: bx,by,bz,anorm,bnorm,ep2,ci,cj
double precision :: ab2,alpha,sab,atij,atmo,r,ichg
pi=3.14159265358979d0
at1e=0.0d0
do k=1,natom
   nuc=k
   ichg=chg(nuc)
   ex=rx(nuc)
   ey=ry(nuc)
   ez=rz(nuc)
   atmo=0.0d0
   do i=1,nbasis
      ci=cmo(i,imo)
      ao1=i
      aoi=nref(ao1)
      aci=ncor(ao1)
      ngi=ngau(aoi)
      ax=rx(aci)
      ay=ry(aci)
      az=rz(aci)
      do j=1,nbasis
         cj=cmo(j,imo)
         ao2=j
         aoj=nref(ao2)
         acj=ncor(ao2)
         ngj=ngau(aoj)
         bx=rx(acj)
         by=ry(acj)
         bz=rz(acj)
         atij=0.0d0
         do ii=1,ngi
            a=bb(aoi,ii)
            anorm=2.0d0*a/pi
            anorm=anorm**0.75
            do jj=1,ngj
               b=bb(aoj,jj)
               bnorm=2.0d0*b/pi
               bnorm=bnorm**0.75
               ab2=(ax-bx)**2+(ay-by)**2+(az-bz)**2
               alpha=a+b
               sab=anorm*bnorm*dexp(-a*b*ab2/alpha)*(pi/alpha)**1.5
               px=(a*ax+b*bx)/alpha
               py=(a*ay+b*by)/alpha
               pz=(a*az+b*bz)/alpha
               ep2=(ex-px)**2+(ey-py)**2+(ez-pz)**2
               r=ep2*alpha
               atij=atij-(cbb(aoi,ii)*cbb(aoj,jj)*&
               &1.1283791670955d0*sab*(alpha**0.5d0)*fun(r))
            end do
         end do
         atmo=atmo+(ci*cj*atij)
      end do
   end do
   at1e=at1e+(ichg*atmo)
end do
end subroutine vns1e

! calculates nuclear repulsion
subroutine enrep(nrep,natom,chg,rx,ry,rz)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::chg(:)
integer,intent(in)::natom
double precision,intent(out)::nrep
integer::i,j
double precision::rijx,rijy,rijz,rij,chgij
nrep=0.0d0
do i=1,natom
   do j=1,natom
      if (j>i) then
         rijx=(rx(i)-rx(j))**2.0d0
         rijy=(ry(i)-ry(j))**2.0d0
         rijz=(rz(i)-rz(j))**2.0d0
         rij=(rijx+rijy+rijz)**0.5d0
         chgij=chg(i)*chg(j)
         nrep=nrep+(chgij/rij)
      end if
   end do
end do
end subroutine enrep

! calculate repulsion energies in system (restricted)
subroutine erepl(e2en,rx,ry,rz,bb,cbb,nref,&
           ncor,ngau,nbasis,cmo,nelec,nlobe,&
           vlobe,jkmat)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::jkmat(:,:,:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
double precision,intent(in)::cmo(:,:)
integer,intent(in)::nbasis,nelec
double precision,intent(out)::e2en
double precision::j2e,k2e
integer::imo,jmo,i,j,spn
imo=0
jmo=0
e2en=0.0d0
do i=1,nelec
   imo=(i+1)/2
   do j=1,nelec
      jmo=(j+1)/2
      call rrepl(imo,imo,jmo,jmo,rx,ry,rz,bb,&
           cbb,nref,ncor,ngau,j2e,k2e,nbasis,&
           cmo,nlobe,vlobe,jkmat)
      ! get spin part for exchange
      spn=abs(mod((i+j),2)-1)
      e2en=e2en+(j2e-k2e*dble(spn)*dble(spn))
    end do
end do
e2en=e2en/2.0d0
end subroutine erepl

! calculates repulsion integrals between mo's
subroutine rrepl(mo1,mo2,mo3,mo4,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,j2e,k2e,nbasis,cmo,nlobe,&
           vlobe,jkmat)
implicit none
integer,intent(in)::mo1,mo2,mo3,mo4
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:),nbasis
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
double precision,intent(in)::cmo(:,:),jkmat(:,:,:,:)
double precision,intent(out)::j2e,k2e
double precision::ci,cj,ck,cl,cijkl,vijkl,vikjl
integer::i,j,k,l
j2e=0.0d0
k2e=0.0d0
do i=1,nbasis
   ci=cmo(i,mo1)
   do j=1,nbasis
      cj=cmo(j,mo2)
      do k=1,nbasis
         ck=cmo(k,mo3)
         do l=1,nbasis
            cl=cmo(l,mo4)
            cijkl=ci*cj*ck*cl
            vijkl=jkmat(i,j,k,l)
            vikjl=jkmat(i,k,j,l)
           !call repl2(i,j,k,l,rx,ry,rz,bb,cbb,&
           !     nref,ncor,ngau,nlobe,vlobe,vijkl)
           !call repl2(i,k,j,l,rx,ry,rz,bb,cbb,&
           !     nref,ncor,ngau,nlobe,vlobe,vikjl)
            vijkl=cijkl*vijkl
            vikjl=cijkl*vikjl
            j2e=j2e+vijkl
            k2e=k2e+vikjl
         end do
      end do
   end do
end do
end subroutine rrepl

! calculates repulsion integral between ao's with 'lobes'
subroutine repl2(ao1,ao2,ao3,ao4,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,nlobe,vlobe,vijkl)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
double precision,intent(in)::vlobe(:,:)
integer,intent(in)::nlobe(:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
integer,intent(in)::ao1,ao2,ao3,ao4
double precision,intent(out)::vijkl
integer::i,j,k,l,aoi,aoj,aok,aol,aci,acj,ack
integer::acl,ngi,ngj,ngk,ngl,ilob,jlob,klob,llob
integer::nli,nlj,nlk,nll
double precision::pi,a,b,ax,ay,az,fact
double precision::bx,by,bz,c,d,cx,cy,cz
double precision::dx,dy,dz,anorm,bnorm,ab2
double precision::alpha,sab,cnorm,dnorm
double precision::cd2,beta,scd,px,py,pz
double precision::qx,qy,qz,gam,pq2,r,vabcd
double precision::signi,signj,signk,signl
double precision::dispa,dispb,dispc,dispd
pi=3.14159265358979d0
aoi=nref(ao1)
aoj=nref(ao2)
aok=nref(ao3)
aol=nref(ao4)
aci=ncor(ao1)
acj=ncor(ao2)
ack=ncor(ao3)
acl=ncor(ao4)
ngi=ngau(aoi)
ngj=ngau(aoj)
ngk=ngau(aok)
ngl=ngau(aol)
nli=nlobe(aoi)
nlj=nlobe(aoj)
nlk=nlobe(aok)
nll=nlobe(aol)
vijkl=0.0d0
do i=1,ngi
   a=bb(aoi,i)
   anorm=2.0d0*a/pi
   anorm=anorm**0.75
   dispa=0.03d0/dsqrt(a)
   do j=1,ngj
      b=bb(aoj,j)
      bnorm=2.0d0*b/pi
      bnorm=bnorm**0.75
      dispb=0.03d0/dsqrt(b)
      alpha=a+b
      do k=1,ngk
         c=bb(aok,k)
         cnorm=2.0d0*c/pi
         cnorm=cnorm**0.75
         dispc=0.03d0/dsqrt(c)
         do l=1,ngl
            d=bb(aol,l)
            dnorm=2.0d0*d/pi
            dnorm=dnorm**0.75
            dispd=0.03d0/dsqrt(d)
            beta=c+d
            signi=1.0d0
            do ilob=1,nli
               ax=rx(aci)+vlobe(ao1,1)*dispa*signi
               ay=ry(aci)+vlobe(ao1,2)*dispa*signi
               az=rz(aci)+vlobe(ao1,3)*dispa*signi
               signj=1.0d0
               do jlob=1,nlj
                  bx=rx(acj)+vlobe(ao2,1)*dispb*signj
                  by=ry(acj)+vlobe(ao2,2)*dispb*signj
                  bz=rz(acj)+vlobe(ao2,3)*dispb*signj
                  signk=1.0d0
                  ab2=(ax-bx)**2+(ay-by)**2+(az-bz)**2
                  sab=anorm*bnorm*dexp(-a*b*ab2/alpha)*(pi/alpha)**1.5
                  px=(a*ax+b*bx)/alpha
                  py=(a*ay+b*by)/alpha
                  pz=(a*az+b*bz)/alpha
                  do klob=1,nlk
                     cx=rx(ack)+vlobe(ao3,1)*dispc*signk
                     cy=ry(ack)+vlobe(ao3,2)*dispc*signk
                     cz=rz(ack)+vlobe(ao3,3)*dispc*signk
                     signl=1.0d0
                     do llob=1,nll
                        dx=rx(acl)+vlobe(ao4,1)*dispd*signl
                        dy=ry(acl)+vlobe(ao4,2)*dispd*signl
                        dz=rz(acl)+vlobe(ao4,3)*dispd*signl
                        cd2=(cx-dx)**2+(cy-dy)**2+(cz-dz)**2
                        scd=cnorm*dnorm*dexp(-c*d*cd2/beta)*(pi/beta)**1.5
                        qx=(c*cx+d*dx)/beta
                        qy=(c*cy+d*dy)/beta
                        qz=(c*cz+d*dz)/beta
                        pq2=(px-qx)**2+(py-qy)**2+(pz-qz)**2
                        gam=(1.0d0/alpha)+(1.0d0/beta)
                        gam=1.0d0/gam
                        r=gam*pq2
                        fact=signi*signj*signk*signl
                        vabcd=1.1283791670955d0*sab*scd*dsqrt(gam)*fun(r)
                        vijkl=vijkl+cbb(aoi,i)*cbb(aoj,j)*cbb(aok,k)&
                        *cbb(aol,l)*vabcd*fact
                        signl=-signl
                      end do
                      signk=-signk
                   end do
                   signj=-signj
               end do
               signi=-signi
            end do
         end do
      end do
   end do
end do
end subroutine repl2

! calculates repulsion integral between ao's
subroutine repl(ao1,ao2,ao3,ao4,rx,ry,rz,bb,cbb,&
           nref,ncor,ngau,vijkl)
implicit none
double precision,intent(in)::rx(:),ry(:),rz(:)
double precision,intent(in)::bb(:,:),cbb(:,:)
integer,intent(in)::nref(:),ncor(:),ngau(:)
integer,intent(in)::ao1,ao2,ao3,ao4
double precision,intent(out)::vijkl
integer::i,j,k,l,aoi,aoj,aok,aol,aci,acj,ack
integer::acl,ngi,ngj,ngk,ngl
double precision::pi,a,b,ax,ay,az
double precision::bx,by,bz,c,d,cx,cy,cz
double precision::dx,dy,dz,anorm,bnorm,ab2
double precision::alpha,sab,cnorm,dnorm
double precision::cd2,beta,scd,px,py,pz
double precision::qx,qy,qz,gam,pq2,r,vabcd
pi=3.14159265358979d0
aoi=nref(ao1)
aoj=nref(ao2)
aok=nref(ao3)
aol=nref(ao4)
aci=ncor(ao1)
acj=ncor(ao2)
ack=ncor(ao3)
acl=ncor(ao4)
ngi=ngau(aoi)
ngj=ngau(aoj)
ngk=ngau(aok)
ngl=ngau(aol)
ax=rx(aci)
ay=ry(aci)
az=rz(aci)
bx=rx(acj)
by=ry(acj)
bz=rz(acj)
cx=rx(ack)
cy=ry(ack)
cz=rz(ack)
dx=rx(acl)
dy=ry(acl)
dz=rz(acl)
vijkl=0.0d0
do i=1,ngi
   a=bb(aoi,i)
   anorm=2.0d0*a/pi
   anorm=anorm**0.75
   do j=1,ngj
      b=bb(aoj,j)
      bnorm=2.0d0*b/pi
      bnorm=bnorm**0.75
      ab2=(ax-bx)**2+(ay-by)**2+(az-bz)**2
      alpha=a+b
      sab=anorm*bnorm*dexp(-a*b*ab2/alpha)*(pi/alpha)**1.5
      px=(a*ax+b*bx)/alpha
      py=(a*ay+b*by)/alpha
      pz=(a*az+b*bz)/alpha
      do k=1,ngk
         c=bb(aok,k)
         cnorm=2.0d0*c/pi
         cnorm=cnorm**0.75
         do l=1,ngl
            d=bb(aol,l)
            dnorm=2.0d0*d/pi
            dnorm=dnorm**0.75
            cd2=(cx-dx)**2+(cy-dy)**2+(cz-dz)**2
            beta=c+d
            scd=cnorm*dnorm*dexp(-c*d*cd2/beta)*(pi/beta)**1.5
            qx=(c*cx+d*dx)/beta
            qy=(c*cy+d*dy)/beta
            qz=(c*cz+d*dz)/beta
            pq2=(px-qx)**2+(py-qy)**2+(pz-qz)**2
            gam=(1.0d0/alpha)+(1.0d0/beta)
            gam=1.0d0/gam
            r=gam*pq2
            vabcd=1.1283791670955d0*sab*scd*dsqrt(gam)*fun(r)
            vijkl=vijkl+cbb(aoi,i)*cbb(aoj,j)*cbb(aok,k)&
            *cbb(aol,l)*vabcd
         end do
      end do
   end do
end do
end subroutine repl

! Read in coords, exps, coeffs, has a lot of
! hard coded stuff...beware.
subroutine initinp(chg,coorx,coory,coorz,bb,cbb,&
           nref,ncor,ngau,nao,nbasis,natom,mo,&
           ovl,kin,nelec,nlobe,vlobe,xdim,nvdet,&
           blokij,civec,svar,exmax,pfac,start,step,&
           xmin,cmo2,civec2,bbc,cbbc,ntype,stopi,&
           Hci,fciv,cieig,pdens,fock,jkmat,ovu,ovx,&
           ovs,vq,h1t)
implicit none
double precision,allocatable,intent(out)::chg(:),&
coorx(:),coory(:),coorz(:),bb(:,:),cbb(:,:),&
vlobe(:,:),bbc(:,:),mo(:,:),ovl(:,:),kin(:,:),&
civec(:),svar(:),start(:),step(:),xmin(:),&
cmo2(:,:),civec2(:),cbbc(:,:),Hci(:,:),fciv(:,:),&
cieig(:),pdens(:,:),fock(:,:),jkmat(:,:,:,:),&
ovu(:,:),ovx(:,:),ovs(:),vq(:,:),h1t(:,:)
integer,allocatable,intent(out)::nref(:),ncor(:),&
ngau(:),nlobe(:),blokij(:),nvdet(:),ntype(:)
integer,intent(out)::nao,nbasis,natom,nelec,exmax,&
xdim,stopi
double precision,intent(out)::pfac
integer::ifile,jw,iw,i,j,maxgau,tdet,ncol,guess
integer::bj,ej,ncol2
character(len=100)::fmtstr
ifile=2
guess=0
open(unit=ifile,file='input.dat',status='old')
! read in number of atoms, number of electrons,
! and excitation level
read(ifile,'(I10)') stopi
read(ifile,'(6(I5))') natom,nelec,nbasis,maxgau,&
exmax,guess
! number of types of ao's
nao=nbasis
! allocate arrays
allocate(ntype(nbasis),nref(nbasis),ncor(nbasis),&
nlobe(nbasis),vlobe(nbasis,3),ngau(nbasis),&
bb(nbasis,maxgau),cbb(nbasis,maxgau),bbc(nbasis,maxgau),&
cbbc(nbasis,maxgau),mo(nbasis,nbasis),&
cmo2(nbasis,nbasis),ovl(nbasis,nbasis),kin(nbasis,nbasis),&
blokij(exmax+1),nvdet(exmax),coorx(natom),coory(natom),&
coorz(natom),chg(natom),svar(nao),pdens(nbasis,nbasis),&
fock(nbasis,nbasis),jkmat(nbasis,nbasis,nbasis,nbasis),&
ovu(nbasis,nbasis),ovx(nbasis,nbasis),ovs(nbasis),&
vq(nbasis,nbasis),h1t(nbasis,nbasis))
! fill everything with zeros
chg=0.0;coorx=0.0;coory=0.0;coorz=0.0;bb=0.0
cbb=0.0;vlobe=0.0;nref=0;ncor=0;ngau=0
nlobe=0;ovl=0.0;kin=0.0;pdens=0.0;fock=0.0
jkmat=0.0;ovs=0.0;ovx=0.0;ovu=0.0;vq=0.0
h1t=0.0
! initialize ci vector
call initci(blokij,nvdet,exmax,nelec,nbasis,tdet)
allocate(civec(tdet),civec2(tdet),Hci(tdet,tdet),&
fciv(tdet,tdet),cieig(tdet))
Hci=0.0d0;fciv=0.0d0;cieig=0.0d0;civec=1.0
do i=1,natom
   ! read in nuclei charge and coordinates 
   read(ifile,'(4(ES16.8))') chg(i),coorx(i),coory(i),coorz(i) 
end do
! read in basis
read(ifile,*)
do i=1,nbasis
   ! read in basis, atom, and lobe refs
   read(ifile,'(5(I5))') nref(i),ncor(i),nlobe(i),&
   ngau(i),ntype(i)
   ! set number of lobes
   if (ntype(i)>0) vlobe(i,ntype(i))=1.0d0
   ! read in expansion (exponents/coefficients)
   fmtstr=''
   write(fmtstr,'(A,I4,A)') '(', ngau(i),'(ES16.8))'
   read(ifile,fmtstr) bb(i,1:ngau(i))
   read(ifile,fmtstr) cbb(i,1:ngau(i))
end do
! init svar
svar=1.20d0;
! get number of parameters
xdim=0
! scaling par.
do i=1,nao
   xdim=xdim+1
end do
! exponents and coefficients
do i=1,nao
   do j=1,ngau(i)
      xdim=xdim+1
   end do
   do j=1,ngau(i)
      xdim=xdim+1
   end do
end do
! mo coefficients
do i=1,nbasis
   do j=1,nbasis
      xdim=xdim+1
   end do
end do
! ci coefficients
do i=1,sum(nvdet)+1
   xdim=xdim+1
end do
! allocate for simplex
allocate(start(xdim),step(xdim),xmin(xdim))
! initialize simplex
start=0.0d0
step=1.0d0
! orthogonality penaly
pfac=1.0d0
! read mos and civec
if (guess==1) then
   ! MO coefficients 
   read(ifile,*)
   ncol=nbasis/5
   do j=1,ncol
      bj=(j-1)*5+1
      ej=(bj-1)+5
      read(ifile,*)
      do i=1,nbasis
         read(ifile,'(I5,5(ES16.8))') iw,cmo(i,bj:ej)
      end do
      read(ifile,*)
   end do
   ncol2=mod(nbasis,5)
   if (ncol2/=0) then
      bj=ncol*5+1
      ej=bj+ncol2-1
      read(ifile,*)
      do i=1,nbasis
         read(ifile,'(I5,5(ES16.8))') iw,cmo(i,bj:ej)
      end do
      read(ifile,*)
   end if
  !! CI expansion coefficients
  !read(ifile,*)
  !ncol=size(civec)/5
  !do j=1,ncol
  !   bj=(j-1)*5+1
  !   ej=(bj-1)+5
  !   !read(ifile,'(5(ES16.8))') civec(bj:ej)
  !   read(ifile,*)
  !end do
  !ncol2=mod(size(civec),5)
  !if (ncol2/=0) then
  !   bj=ncol*5+1
  !   ej=bj+ncol2-1
  !   !read(ifile,'(5(ES16.8))') civec(bj:ej)
  !   read(ifile,*)
  !end if
else
   ! init mo coefficients
   mo=0.0d0
   do i=1,nbasis
      mo(i,i)=1.0d0
   end do
   ! init ci vars
   civec=1.0
end if
! make backups of mos,civecs,and exps
cmo2=mo;civec2=civec;bbc=bb;cbbc=cbb
close(unit=ifile)
end  subroutine initinp


! Below code is taken from the wonderful collection of algorithms
! implemented by John Burkardt:
! https://github.com/johannesgerer/jburkardt-f/blob/master/subset/subset.f90
subroutine comb_next ( n, k, a, done )
!  COMB_NEXT computes combinations of K things out of N.
!  Discussion:
!    The combinations are computed one at a time, in lexicographical order.
!    10 April 1009: Thanks to "edA-qa mort-ora-y" for supplying a
!    correction to this code!
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    10 April 2009
!  Author:
!    John Burkardt
!  Reference:
!    Charles Mifsud,
!    Algorithm 154:
!    Combination in Lexicographic Order,
!    Communications of the ACM,
!    March 1963.
!  Parameters:
!    Input, integer N, the total number of things.
!    Input, integer K, the number of things in each combination.
!    Input/output, integer A(K), contains the list of elements in
!    the current combination.
!    Input/output, logical DONE.  On first call, set DONE to TRUE,
!    and thereafter, its input value should be the output value from
!    the previous call.  The output value will normally be FALSE,
!    indicating that there are further combinations that can be
!    returned.  When DONE is returned TRUE, the sequence is exhausted.
implicit none
integer::k
integer::a(k)
logical done
integer::i,j,n
if ( done ) then
   if ( k <= 0 ) then
      return
   end if
   call i4vec_indicator ( k, a )
   done = .false.
else
   if ( a(k) < n ) then
      a(k) = a(k) + 1
      return
   end if
   do i = k, 2, -1
      if ( a(i-1) < n-k+i-1 ) then
         a(i-1) = a(i-1) + 1
         do j = i, k
            a(j) = a(i-1) + j - ( i-1 )
         end do
         return
      end if
   end do
   done = .true.
end if
end subroutine

subroutine perm_next ( n, p, more, even )
!*****************************************************************************80
!  PERM_NEXT computes all of the permutations of N objects, one at a time.
!  Discussion:
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    30 March 2001
!  Author:
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!  Reference:
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!  Parameters:
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  On the first call, the input value is unimportant.
!    On subsequent calls, the input value should be the same
!    as the output value from the previous call.  In other words, the
!    user should just leave P alone.
!    On output, contains the "next" permutation.
!    Input/output, logical MORE.
!    Set MORE = FALSE before the first call.
!    MORE will be reset to TRUE and a permutation will be returned.
!    Each new call produces a new permutation until
!    MORE is returned FALSE.
!    Input/output, logical EVEN.
!    The input value of EVEN should simply be its output value from the
!    previous call; (the input value on the first call doesn't matter.)
!    On output, EVEN is TRUE if the output permutation is even, that is,
!    involves an even number of transpositions.
  implicit none
  integer::n
  logical::even
  logical::more
  integer::p(n)
  integer::i
  integer::i1
  integer::ia
  integer::id
  integer::is
  integer::j
  integer::l
  integer::m
  if ( .not. more ) then
    call i4vec_indicator ( n, p )
    more = .true.
    even = .true.
    if ( n == 1 ) then
      more = .false.
      return
    end if
    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if
    do i = 1, n - 3
      if ( p(i+1) /= p(i)+1 ) then
        return
      end if
    end do
    more = .false.
  else
    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if
    if ( even ) then
      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.
      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if
      do i = 1, n - 3
        if ( p(i+1) /= p(i)+1 ) then
          return
        end if
      end do
      more = .false.
      return
    else
      more = .false.
      is = 0
      do i1 = 2, n
        ia = p(i1)
        i = i1 - 1
        id = 0
        do j = 1, i
          if ( ia < p(j) ) then
            id = id + 1
          end if
        end do
        is = id + is
        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if
      end do
      if ( .not. more ) then
        p(1) = 0
        return
      end if
    end if
    m = mod ( is + 1, 2 ) * ( n + 1 )
    do j = 1, i
      if ( sign ( 1, p(j)-ia ) /= sign ( 1, p(j)-m ) ) then
        m = p(j)
        l = j
      end if
    end do
    p(l) = ia
    p(i1) = m
    even = .true.
  end if
end subroutine

subroutine i4vec_indicator ( n, a )
! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    09 November 2000
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
  implicit none
  integer::n
  integer::a(n)
  integer::i
  do i = 1, n
    a(i) = i
  end do
  return
end subroutine
! end of the Burkardt extract

double precision function fun(r)
   implicit none
   double precision, intent(inout):: r
   double precision:: rer
   if(r.le.0.0002d0)fun=1.0d0-r/3.0d0 +0.1d0*r*r -r*r*r/42.0d0
   if(r.gt.0.0002d0)then
      r=dsqrt(r)
      rer=derf(r)/r
      fun=rer/1.1283791670955d0
   end if 
   return
end function fun

integer function factorial(n)
implicit none
integer,intent(in)::n
integer::i,ans
ans=1
do i=1,n
   ans=ans*i
end do
factorial=ans
end function factorial

INTEGER FUNCTION  FindMinimum(x, Start, E)
   IMPLICIT  NONE
   INTEGER, DIMENSION(1:), INTENT(IN) :: x
   INTEGER, INTENT(IN)                :: Start, E
   INTEGER                            :: Minimum
   INTEGER                            :: Location
   INTEGER                            :: i
   Minimum  = x(Start)     ! assume the first is the min
   Location = Start        ! record its position
   DO i = Start+1, E     ! start with next elements
      IF (x(i) < Minimum) THEN   !   if x(i) less than the min?
         Minimum  = x(i)      !      Yes, a new minimum found
         Location = i                !      record its position
      END IF
   END DO
   FindMinimum = Location           ! return the position
END FUNCTION  FindMinimum

! This subroutine swaps the values of its two formal arguments.
SUBROUTINE  Swap(a, b)
   IMPLICIT  NONE
   INTEGER, INTENT(INOUT) :: a, b
   INTEGER                :: Temp

   Temp = a
   a    = b
   b    = Temp
END SUBROUTINE  Swap

! This subroutine receives an array x() and sorts it into ascending order.
SUBROUTINE  Sort(x, n)
   IMPLICIT  NONE
   INTEGER, DIMENSION(1:), INTENT(INOUT) :: x
   INTEGER, INTENT(IN)                   :: n
   INTEGER                               :: i
   INTEGER                               :: Location
   DO i = 1, n-1        ! except for the last
      Location = FindMinimum(x, i, n)  ! find min from this to last
      CALL  Swap(x(i), x(Location)) ! swap this and the minimum
   END DO
END SUBROUTINE  Sort

end program
