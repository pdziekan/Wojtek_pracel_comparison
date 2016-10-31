      program moist_convec

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Simple 2d periodic domain flow model
c  with stretched vertical grid
c
c   Large Eddy Model for singlu CU as in Brenguier/Grabowski JAS 1993
c
c  Author: W. Grabowski (grabow@ncar.ucar.edu) 
c
c      modified initialization
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc
cc  nember of time steps:
      parameter(ntime0=10*30)
cc  major dimensions
cc HAVE TO BE CHANGED IN ALL ROUTINES:
      parameter(nx=181,nz=121)
      parameter(nxz=nx*nz)
cc  grid (for regular dzeta isosurfaces)
      data dx,dz,dt /20.,20.,2.0/
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      dimension xx(nx),zz(nz)

      dimension cntrm(ntime0),cntti(ntime0)
      dimension qcav(ntime0)

CC  MODEL VARIABLES
cc thermodynamics
      dimension theta(nx,nz),qv(nx,nz),qc(nx,nz),qr(nx,nz)
cc dynamics
      dimension ux(nx,nz),uy(nx,nz),uz(nx,nz)          ! current time level
      dimension uxp(nx,nz),uyp(nx,nz),uzp(nx,nz)       ! previous time level
cc  forces for model variables
      dimension ft(nx,nz),fx(nx,nz),fy(nx,nz),fz(nx,nz),
     *          fqv(nx,nz),fqc(nx,nz),fqr(nx,nz)
cc  advective velocities:
      dimension uxa(nx+1,nz),uza(nx,nz+1)

cc profiles  
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common /prof_a/ tau(nz)

cc required variables:
      dimension den(nx,nz),p(nx,nz),scr1(nx,nz),scr2(nx,nz)

cc constants
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

cc random noise variables (in pbl)
      dimension thn(nx,nz),qvn(nx,nz)

      call opngks
      call gsclip(0)

cc grid:
      time=0.
      dxi=1./dx
      dzi=1./dz
      dti=1./dt
      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
cc zz is regular (ustretched) dzeta grid
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo

      call zstrtch(zz,nz,dz)

cc initialize moisture parameters:
      call moist_init

cc initialize model profiles:
      call prof_init_dry
c      call prof_init_1layer

cc absorber (tau=1/time scale)
      zab = 3500000.     ! height at which absorber begins
      t_abs=500.      ! time scale
      towi=1./t_abs
      do k=1,nz
      tau(k)=towi*amax1(0.,height(k)-zab)/(height(nz)-zab)
      enddo

       do k=1,nz
       do i=1,nx
       den(i,k)=rho0(k)*gac(k)
       enddo
       enddo

cc initial fields
cc call noise:
c        call noise(thn,nx,nz,nz_noise)
c        call integxz(thn,scr2,nx,nz)
      a=rg/rv
      c=hlatv/cp
      b=hlatv/(rv*tt0)
      d=hlatv/rv
      e=-cp/rg

       zcen=800.
       xcen=(nx-1)*dx/2.
       rad1=200.
       rad2=300.
       pi=4.*atan(1.)

       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       do i=1,nx

        radd=(xx(i)-xcen)**2+(zz(k)-zcen)**2
        radd=sqrt(radd)
         if(radd.gt.rad2) rh=0.
         if(radd.le.rad1) rh=1.
         if(radd.gt.rad1 .and. radd.le.rad2)
     1   rh=1.*(cos(pi/2. * (radd-rad1)/(rad2-rad1)))**2.

        theta(i,k)=th_e(k) + rh*3.

        qc(i,k)=0.
        qr(i,k)=0.

        uy(i,k)=uy_e(k)
        ux(i,k)=ux_e(k)
        uz(i,k)=0.
        uxp(i,k)=ux_e(k)
        uzp(i,k)=0.

        ft(i,k)=0.
        fx(i,k)=0.
        fy(i,k)=0.
        fz(i,k)=0.

        fqv(i,k)=0.
        fqc(i,k)=0.
        fqr(i,k)=0.

        p(i,k)=0.
       enddo
       enddo

cc total water:
       sum=0.
       do i=1,nx-1
       do k=1,nz
       coe=1.
       if(k.eq.1 .or. k.eq.nz) coe=0.5
       sum=sum+coe*den(i,k)*(qv(i,k)+qc(i,k))
       enddo
       enddo
       print*,'------ total water:',sum

cc plot initial fields:
       call diagno_1(ux,uy,uz,theta,nx,nz,scr1,scr2,den)
       call diagno_2(qv,qc,qr,nx,nz,scr1,scr2,den)
       call plot_1(ux,uy,uz,theta)
c       call plot_2(theta,qv,qc,qr)

ccc save initial data:
c       write(17) time,nx,nz
c       write(17) ux,uy,uz,theta,qv,qc,qr


CCCC MARCH FORWARD IN TIME:
              ntime=ntime0
              do itime=1,ntime   ! TIME LOOP
               print*,'*** itime, time: ',itime,time

cc extrapolate in time to get advective momentums:
       call velprd_1(ux,uxp,uxa,uz,uzp,uza,nx,nz,den)
          
cc save previous velocities:
       do i=1,nxz
        uxp(i,1)=ux(i,1)
        uzp(i,1)=uz(i,1)
       enddo

c surface flux:
c       call surfflux(theta,qv,ft,fqv,ux,uy,fx)

cc add half of the force:
       do i=1,nxz
        theta(i,1)=theta(i,1)+.5*dt*ft(i,1)
        ux(i,1)   =   ux(i,1)+.5*dt*fx(i,1)
        uy(i,1)   =   uy(i,1)+.5*dt*fy(i,1)
        uz(i,1)   =   uz(i,1)+.5*dt*fz(i,1)
        qv(i,1)   =   qv(i,1)+.5*dt*fqv(i,1)
        qc(i,1)   =   qc(i,1)+.5*dt*fqc(i,1)
c        qr(i,1)   =   qr(i,1)+.5*dt*fqr(i,1)
       enddo

CC ADVECTION:
c liner: 1-iord=1, 0-iord prescribed inside mpdata
        liner=0
        if(itime/10*10.eq.itime) liner=1

cc advect velocities:
        call mpdat_2d(uxa,uza,   ux,den,1,liner)
        call mpdat_2d(uxa,uza,   uz,den,2,liner)
cc advect thermodynamic variables:
        call mpdat_2d(uxa,uza,theta,den,3,liner)
        call mpdat_2d(uxa,uza,   qv,den,4,liner)
        call mpdat_2d(uxa,uza,   qc,den,5,liner)
c          call rain_fall(qr,tm_e,rho0,uza)
c        call mpdat_2d(uxa,uza,   qr,den,6,liner)

cc save velocities after advection into advective velocities:
cc (used as scratch)
       do k=1,nz
       do i=1,nx
       uxa(i,k)=ux(i,k)
       uza(i,k)=uz(i,k)
       enddo
       enddo

cc finish thermodynamics 
        call thermo(theta,qv,qc,qr,ft,fqv,fqc,fqr)
       
cc add absorber:
       if(zab.lt.zz(nz)) 
     1        call absor(ux,uz,theta,qv,qc,qr,ft,fqv,fqc,fqr)

cc add buoyancy
       epsb=rv/rg-1.
       do k=1,nz
       do i=1,nx
c       scr1(i,k)=gg*( (theta(i,k)-th_e(k))/th0(k)
c     *   + epsb*(qv(i,k)-qv_e(k))-qc(i,k)-qr(i,k) )
       scr1(i,k)=gg*( (theta(i,k)-th_e(k))/th0(k))
       enddo
       enddo

cc filter in vertical
cc       call integz(scr1,scr2,nx,nz)
       call integxz(scr1,scr2,nx,nz)

cc apply
       do k=1,nz
       do i=1,nx
       uz(i,k) = uz(i,k)+.5*dt*scr1(i,k)
       enddo
       enddo

cc calculate pressure gradient force:
      epp=1.e-6
c      epp=1.e-7
      itp=100
      call gcrk_1(p,scr1,scr2,ux,uz,nx,nz,itp,epp)
      call prforc_1(p,scr1,scr2,ux,uz,nx,nz)
      do k=1,nz
      do i=1,nx
      ux(i,k)=scr1(i,k)
      uz(i,k)=scr2(i,k)
      enddo
      enddo

cc calculate velocity forces (using saved velocities after advection):
       do k=1,nz
       do i=1,nx
       fx(i,k)=(ux(i,k)-uxa(i,k))  *2./dt
       fz(i,k)=(uz(i,k)-uza(i,k))  *2./dt
       enddo
       enddo

cc update clock (in minutes...)
       time = float(itime)*dt/60. 

cc output and plot:
       dtout=1. ! in min
       dtape=50.  ! in min
       if(amod(time+.1*dt/60.,dtout).lt.0.5*dt/60.) then
ccc plot selected fields:
       call plot_1(ux,uy,uz,theta)
c       call plot_2(theta,qv,qc,qr)
cc analysis of output:
       print*,'   '
       call diagno_1(ux,uy,uz,theta,nx,nz,scr1,scr2,den)
       call diagno_2(qv,qc,qr,nx,nz,scr1,scr2,den)
cc 
       endif
       if(amod(time+.1*dt/60.,dtape).lt.0.5*dt/60.) then
c       write(17) time,nx,nz
c       write(17) ux,uy,uz,theta,qv,qc,qr
c       print*,'wrote tape for t = ',time
       endif

      cntti(itime)=time
cc center of mass of qc:
      sum1=0.
      sum2=0.
      do i=1,nx
      do k=1,nz
      sum1=sum1 + qc(i,k)
      sum2=sum2 + height(k)*qc(i,k)
      enddo
      enddo
      if(sum1.gt.1.e-3) then
      cntrm(itime)=sum2/sum1  *1.e-3
      else
      cntrm(itime)=0.
      endif

cc conditionally-averaged qc:
      sum1=0.
      sum2=0.
      do i=1,nx
      do k=1,nz
      if(qc(i,k).ge.1.e-5) then
      sum1=sum1 + 1.
      sum2=sum2 + qc(i,k)
      endif
      enddo
      enddo
      if(sum1.gt.1.) then
      qcav(itime)=sum2/sum1  *1.e3
      else
      qcav(itime)=0.
      endif

cc total water:
       sum=0.
       do i=1,nx-1
       do k=1,nz
       coe=1.
       if(k.eq.1 .or. k.eq.nz) coe=0.5
       sum=sum+coe*den(i,k)*(qv(i,k)+qc(i,k))
       enddo
       enddo
       print*,'------ total water:',sum
      
 
           enddo      ! TIME LOOP

cc    finished...

ccc plot cntrm and qcav:
c      call setusv('LW',2000)
c      call set(.15,.95,.1,.5,0.,10.,0.,2.,1)           
c      call labmod('(f3.0)','(f3.0)',4,4,2,2,20,20,0)     
c      call periml(2,5,2,5)                           
c      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)      
c      call curved(cntti,cntrm,ntime)                  
c      call set(.15,.95,.55,.95,0.,10.,0.,1.,1)           
c      call labmod('(f3.0)','(f3.1)',4,4,2,2,20,20,0)     
c      call periml(2,5,2,5)                           
c      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)      
c      call curved(cntti,qcav,ntime)                  
c      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      CALL plchhq(.55,0.03, 'time (min)', 0.016,0.,0)
c      CALL plchhq(.03,0.30, 'qc center of mass (km)', 0.016,90.,0)
c      CALL plchhq(.03,0.75, 'averaged qc (g/kg)', 0.016,90.,0)
c      call frame                                   

        call clsgks
        stop
        end

      subroutine noise(ff,nx,nz,nzn)
      dimension ff(nx,nz)
c      double precision rand

      do i=1,nx
      do k=1,nz
      ff(i,k)=0.
      enddo
      enddo

      do i=1,nx-1
      do k=2,nzn
      ff(i,k)=2.*(rand()-.5)
      enddo
      enddo
      do k=1,nzn
      ff(nx,k)=ff(1,k)
      enddo

      return
      end

      subroutine velprd_1(ux,uxp,uxa,uz,uzp,uza,nx,nz,rho)
      dimension ux(nx,nz),uz(nx,nz)     
      dimension uxp(nx,nz),uzp(nx,nz)   
      dimension uxa(nx+1,nz),uza(nx,nz+1),rho(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

       do k=1,nz
       do i=2,nx
       uxa(i,k) =(0.75*(ux(i-1,k)*rho(i-1,k)+ux(i,k)*rho(i,k)) 
     .          - 0.25*(uxp(i-1,k)*rho(i-1,k)+uxp(i,k)*rho(i,k)) )*dt/dx 
       enddo
cc cyclic in horizontal
       uxa(1,k) = uxa(nx,k)
       uxa(nx+1,k) = uxa(2,k)
       enddo
          
       do i=1,nx
       do k=2,nz
       uza(i,k) =(0.75*(uz(i,k-1)*rho(i,k-1)+uz(i,k)*rho(i,k)) 
     .          - 0.25*(uzp(i,k-1)*rho(i,k-1)+uzp(i,k)*rho(i,k)) )*dt/dz 
       enddo
cc zero flux in vertical
       uza(i,1) = - uza(i,2)
       uza(i,nz+1) = - uza(i,nz)
       enddo

       return
       end

      subroutine mpdat_2d(u1,u2,x,h,iflg,liner)
      parameter(nx=181,nz=121)
      parameter(n1=nx+1,n2=nz+1)
      parameter(n1m=n1-1,n2m=n2-1)
      dimension u1(n1,n2m),u2(n1m,n2),x(n1m,n2m),h(n1m,n2m)
      common// v1(n1,n2m),v2(n1m,n2),f1(n1,n2m),f2(n1m,n2),
     *         cp(n1m,n2m),cn(n1m,n2m),
     *         mx(n1m,n2m),mn(n1m,n2m)
      real mx,mn
      parameter(iord0=2,isor=1,nonos=1,idiv=0)
      data ep/1.e-12/
c
      donor(y1,y2,a)=cvmgm(y2,y1,a)*a
      vdyf(x1,x2,a,r)=(abs(a)-a**2/r)*(abs(x2)-abs(x1))
     1                               /(abs(x2)+abs(x1)+ep)
      vcorr(a,b,y1,y2,r)=-0.125*a*b*y1/(y2*r)
      vcor31(a,x0,x1,x2,x3,r)= -(a -3.*abs(a)*a/r+2.*a**3/r**2)/3.
     1                         *(abs(x0)+abs(x3)-abs(x1)-abs(x2))
     2                         /(abs(x0)+abs(x3)+abs(x1)+abs(x2)+ep)
      vcor32(a,b,y1,y2,r)=0.25*b/r*(abs(a)-2.*a**2/r)*y1/y2
      vdiv1(a1,a2,a3,r)=0.25*a2*(a3-a1)/r
      vdiv2(a,b1,b2,b3,b4,r)=0.25*a*(b1+b2-b3-b4)/r
      pp(y)= amax1(0.,y)
      pn(y)=-amin1(0.,y)
 
      iord=iord0
      if(isor.eq.3) iord=max0(iord,3)
      if(liner.eq.1) iord=1

      do j=1,n2-1
        do i=1,n1
          v1(i,j) = u1(i,j)
        end do 
      end do 
      do i=1,n1-1
        do j=1,n2
          v2(i,j) = u2(i,j)
        end do 
      enddo

      if(nonos.eq.1) then
      do j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do i=1,n1m
      im=(i-1+(n1-i)/n1m*(n1-2))
      ip=(i+1    -i /n1m*(n1-2))
      mx(i,j)=amax1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
      mn(i,j)=amin1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
      end do
      end do
      endif
 
                         do 3 k=1,iord
 
      do 331 j=1,n2-1
      do 331 i=2,n1-1
  331 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      do j=1,n2-1
      f1(1 ,j)=f1(n1-1,j)
      f1(n1,j)=f1(2,j)
      enddo
      do 332 j=2,n2-1
      do 332 i=1,n1-1
  332 f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
      if (iflg.eq.6) then
        do i=1,n1m
          f2(i, 1)=donor(x(i,  1),x(i,  1),v2(i, 1))
          f2(i,n2)=donor(x(i,n2m),x(i,n2m),v2(i,n2))
        end do
      else
        do i=1,n1m
          f2(i, 1)=-f2(i,  2)
          f2(i,n2)=-f2(i,n2m)
        end do
      end if
  
      do 333 j=1,n2-1
      do 333 i=1,n1-1
  333 x(i,j)=x(i,j)-(f1(i+1,j)-f1(i,j)+f2(i,j+1)-f2(i,j))/h(i,j)
 
      if(k.eq.iord) go to 6

      do 49 j=1,n2-1
      do 49 i=1,n1
      f1(i,j)=v1(i,j)
   49 v1(i,j)=0.
      do 50 j=1,n2
      do 50 i=1,n1-1
      f2(i,j)=v2(i,j)
   50 v2(i,j)=0.
      do 51 j=2,n2-2
      do 51 i=2,n1-1
   51 v1(i,j)=vdyf(x(i-1,j),x(i,j),f1(i,j),.5*(h(i-1,j)+h(i,j)))
     *       +vcorr(f1(i,j), f2(i-1,j)+f2(i-1,j+1)+f2(i,j+1)+f2(i,j),
     *   abs(x(i-1,j+1))+abs(x(i,j+1))-abs(x(i-1,j-1))-abs(x(i,j-1)),
     *   abs(x(i-1,j+1))+abs(x(i,j+1))+abs(x(i-1,j-1))+abs(x(i,j-1))+ep,
     *                 .5*(h(i-1,j)+h(i,j)))
      if(idiv.eq.1) then
      do 511 j=2,n2-2
      do 511 i=2,n1-1
  511 v1(i,j)=v1(i,j)
     *    -vdiv1(f1(i-1,j),f1(i,j),f1(i+1,j),.5*(h(i-1,j)+h(i,j)))
     *    -vdiv2(f1(i,j),f2(i-1,j+1),f2(i,j+1),f2(i-1,j),f2(i,j),
     *                 .5*(h(i-1,j)+h(i,j)))
      endif
      do 52 j=2,n2-1
      do 52 i=2,n1-2
   52 v2(i,j)=vdyf(x(i,j-1),x(i,j),f2(i,j),.5*(h(i,j-1)+h(i,j)))
     *       +vcorr(f2(i,j), f1(i,j-1)+f1(i,j)+f1(i+1,j)+f1(i+1,j-1),
     *   abs(x(i+1,j-1))+abs(x(i+1,j))-abs(x(i-1,j-1))-abs(x(i-1,j)),
     *   abs(x(i+1,j-1))+abs(x(i+1,j))+abs(x(i-1,j-1))+abs(x(i-1,j))+ep,
     *                 .5*(h(i,j-1)+h(i,j)))
      i0=n1-2
      do j=2,n2-1
      v2(1,j)=vdyf(x(1,j-1),x(1,j),f2(1,j),.5*(h(1,j-1)+h(1,j)))
     *       +vcorr(f2(1,j), f1(1,j-1)+f1(1,j)+f1(2,j)+f1(2,j-1),
     *   abs(x(2,j-1))+abs(x(2,j))-abs(x(i0,j-1))-abs(x(i0,j)),
     *   abs(x(2,j-1))+abs(x(2,j))+abs(x(i0,j-1))+abs(x(i0,j))+ep,
     *                 .5*(h(1,j-1)+h(1,j)))
      v2(n1-1,j)=v2(1,j)
      enddo

      if(idiv.eq.1) then
      do 521 j=2,n2-1
      do 521 i=1,n1-1
  521 v2(i,j)=v2(i,j)
     *    -vdiv1(f2(i,j-1),f2(i,j),f2(i,j+1),.5*(h(i,j-1)+h(i,j)))
     *    -vdiv2(f2(i,j),f1(i+1,j),f1(i+1,j-1),f1(i,j-1),f1(i,j),
     *                 .5*(h(i,j-1)+h(i,j)))
      endif
      if(isor.eq.3) then
      do 61 j=2,n2-2
      do 61 i=3,n1-2
   61 v1(i,j)=v1(i,j)     +vcor31(f1(i,j),
     1        x(i-2,j),x(i-1,j),x(i,j),x(i+1,j),.5*(h(i-1,j)+h(i,j)))
      do j=2,n2-2
      v1(2,j)=v1(2,j)     +vcor31(f1(2,j),
     1        x(n1-2,j),x(1,j),x(2,j),x(3,j),.5*(h(1,j)+h(2,j)))
      v1(n1-1,j)=v1(n1-1,j) +vcor31(f1(n1-1,j),x(n1-3,j),x(n1-2,j),
     1                  x(n1-1,j),x(2,j),.5*(h(n1-2,j)+h(n1-1,j)))
      enddo
      do 62 j=2,n2-2
      do 62 i=2,n1-1
   62 v1(i,j)=v1(i,j)
     1 +vcor32(f1(i,j),f2(i-1,j)+f2(i-1,j+1)+f2(i,j+1)+f2(i,j),
     *   abs(x(i,j+1))-abs(x(i,j-1))-abs(x(i-1,j+1))+abs(x(i-1,j-1)),
     *   abs(x(i,j+1))+abs(x(i,j-1))+abs(x(i-1,j+1))+abs(x(i-1,j-1))+ep,
     *                   .5*(h(i-1,j)+h(i,j)))
      do 63 j=3,n2-2
      do 63 i=1,n1-1
   63 v2(i,j)=v2(i,j)     +vcor31(f2(i,j),
     1        x(i,j-2),x(i,j-1),x(i,j),x(i,j+1),.5*(h(i,j-1)+h(i,j)))
      do 64 j=3,n2-2
      do 64 i=2,n1-2
   64 v2(i,j)=v2(i,j)
     1 +vcor32(f2(i,j),f1(i,j-1)+f1(i+1,j-1)+f1(i+1,j)+f1(i,j),
     *   abs(x(i+1,j))-abs(x(i-1,j))-abs(x(i+1,j-1))+abs(x(i-1,j-1)),
     *   abs(x(i+1,j))+abs(x(i-1,j))+abs(x(i+1,j-1))+abs(x(i-1,j-1))+ep,
     *                   .5*(h(i,j-1)+h(i,j)))
      do 641 j=3,n2-2
      v2(1,j)=v2(1,j)
     1 +vcor32(f2(1,j),f1(1,j-1)+f1(2,j-1)+f1(2,j)+f1(1,j),
     *   abs(x(2,j))-abs(x(n1-2,j))-abs(x(2,j-1))+abs(x(n1-2,j-1)),
     *   abs(x(2,j))+abs(x(n1-2,j))+abs(x(2,j-1))+abs(x(n1-2,j-1))+ep,
     *                   .5*(h(1,j-1)+h(1,j)))
  641 v2(n1-1,j)=v2(1,j)
      endif
 
        do j=1,n2m
          v1( 1,j)=v1(n1m,j)
          v1(n1,j)=v1(  2,j)
        end do

      if (iflg.ne.6) then
        do i=1,n1m
          v2(i, 1)=-v2(i,  2)
          v2(i,n2)=-v2(i,n2m)
        end do
      end if

                  if(nonos.eq.1) then
c                 non-osscilatory option
      do 401 j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do 401 i=1,n1m
      im=(i-1+(n1-i)/n1m*(n1-2))
      ip=(i+1    -i /n1m*(n1-2))
      mx(i,j)=amax1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mx(i,j))
  401 mn(i,j)=amin1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mn(i,j))

      do 402 j=1,n2m 
      do 4021 i=2,n1-1
 4021 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      f1(1 ,j)=f1(n1m,j)
      f1(n1,j)=f1(2  ,j)
  402 continue
     
      do 403 i=1,n1m
      do 4031 j=2,n2m
 4031 f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
      if(iflg.ne.6) then
      f2(i, 1)=-f2(i,  2)
      f2(i,n2)=-f2(i,n2m)
      else
      f2(i, 1)=0.
      f2(i,n2)=0.
      endif
  403 continue

      do 404 j=1,n2m   
      do 404 i=1,n1m
      cp(i,j)=(mx(i,j)-x(i,j))*h(i,j)/
     1(pn(f1(i+1,j))+pp(f1(i,j))+pn(f2(i,j+1))+pp(f2(i,j))+ep)
      cn(i,j)=(x(i,j)-mn(i,j))*h(i,j)/
     1(pp(f1(i+1,j))+pn(f1(i,j))+pp(f2(i,j+1))+pn(f2(i,j))+ep)
  404 continue
      do 405 j=1,n2m 
      do 4051 i=2,n1m 
 4051 v1(i,j)=pp(v1(i,j))*
     1  ( amin1(1.,cp(i,j),cn(i-1,j))*pp(sign(1., x(i-1,j)))
     1   +amin1(1.,cp(i-1,j),cn(i,j))*pp(sign(1.,-x(i-1,j))) )
     2       -pn(v1(i,j))*
     2  ( amin1(1.,cp(i-1,j),cn(i,j))*pp(sign(1., x(i ,j )))
     2   +amin1(1.,cp(i,j),cn(i-1,j))*pp(sign(1.,-x(i ,j ))) )
      v1( 1,j)=v1(n1m,j)
      v1(n1,j)=v1( 2 ,j)
  405 continue

      do 406 i=1,n1m 
      do 406 j=2,n2m 
  406 v2(i,j)=pp(v2(i,j))*
     1  ( amin1(1.,cp(i,j),cn(i,j-1))*pp(sign(1., x(i,j-1)))
     1   +amin1(1.,cp(i,j-1),cn(i,j))*pp(sign(1.,-x(i,j-1))) )
     1       -pn(v2(i,j))*
     2  ( amin1(1.,cp(i,j-1),cn(i,j))*pp(sign(1., x(i ,j )))
     2   +amin1(1.,cp(i,j),cn(i,j-1))*pp(sign(1.,-x(i ,j ))) )
                  endif
    3                      continue
    6 continue
      return
      end   


      subroutine gcrk_1(p,pfx,pfz,u,w,n1,n3,itr,eps0)
      real p(*),pfx(*),pfz(*),u(*),w(*)
      parameter(nx=181,nz=121)
      parameter(n=nx,l=nz)
      parameter(nn=n*l,nl=n*l)
      common// r(nn),qr(nn),ar(nn)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /strtch/ height(nz),gac(nz)
      common/itero/ niter,nitsm,icount,eer,eem
      parameter (lord=3)
      dimension x(nn,lord),ax(nn,lord),ax2(lord),axar(lord),del(lord)
      dimension rho2d(nx,nz)
convergence test modes **************************************************
      logical ctest                                                     *
      data ctest/.false./                                               *
c     data ctest/.true./                                                *
      parameter (nplt=100)                                              *
      dimension err(0:nplt),xitr(0:nplt)                                *
      if(ctest) then                                                    *
      itr=6000/lord                                                     *
      ner=60                                                            *
      snorm=1./float(n*l)                                               *
      eps0=1.e-15                                                       *
      endif                                                             *
convergence test modes **************************************************
 
      do k=1,l
      do i=1,n
      rho2d(i,k)=rho0(k)*gac(k)
      enddo
      enddo

      eps=eps0*dti
      epa=1.e-30
      nlc=0
      itmn=5
cc iprc 0-no preconditioning, 1-with precon
      iprc=1

      call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,0)

      do k=1,nl
        r(k)=0.
       ar(k)=0.
       qr(k)=0.
      enddo
      do i=1,lord
       do k=1,nl
         x(k,i)=0.
        ax(k,i)=0.
       enddo
      enddo
      call prforc_1(p,pfx,pfz,u,w,n1,n3)
       call rhsdiv_1(pfx,pfz,rho2d,r,n1,n3,-1)
        call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1)
      eer0=0.
      eem0=-1.e15
      rl20=0.
      do k=1,nl
      eer0=eer0+qr(k)**2
      eem0=amax1(eem0,abs(qr(k)))
      rl20=rl20+r(k)**2
      enddo
      eer0=amax1(eer0,epa)
      eem0=amax1(eem0,epa)
      rl20=amax1(rl20,epa)
convergence test modes **************************************************
      if(ctest) then                                                    *
      do ier=0,nplt                                                     *
      err(ier)=eps                                                      *
      enddo                                                             *
      eer=-1.e15                                                        *
      do 3 k=1,nl                                                       *
    3 eer=amax1(eer,abs(r(k)))                                          *
      err(0)=eer                                                        *
      print 300,  err(0)                                                *
  300 format(4x,e11.4,' residual error at it=1')                        *
      endif                                                             *
convergence test modes **************************************************
       do k=1,nl
        x(k,1)=qr(k)
       enddo
      call laplc_1(x(1,1),ax(1,1),pfx,pfz,n1,n3)

      do 100 it=1,itr
       do i=1,lord
        ax2(i)=0.
        rax=0.
         do k=1,nl
          rax=rax+r(k)*ax(k,i)
          ax2(i)=ax2(i)+ax(k,i)*ax(k,i)
         enddo
        ax2(i)=amax1(epa,ax2(i))
        beta=-rax/ax2(i)
        dvmx=-1.e15
        rl2=0.
         do k=1,nl
          p(k)=p(k)+beta* x(k,i)
          r(k)=r(k)+beta*ax(k,i)
          dvmx=amax1(dvmx,abs(r(k)))
          rl2=rl2+r(k)*r(k)
         enddo
       if(dvmx.le.eps.and.it.ge.itmn) go to 200
       if(rl2.ge.rl20.and..not.ctest) go to 200
          rl20=amax1(rl2,epa)
       call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1)
       call laplc_1(qr,ar,pfx,pfz,n1,n3)
        nlc=nlc+1
         do ii=1,i
          axar(ii)=0.
           do k=1,nl
            axar(ii)=axar(ii)+ax(k,ii)*ar(k)
           enddo
          del(ii)=-axar(ii)/ax2(ii)
         enddo
        if(i.lt.lord) then
          do k=1,nl
            x(k,i+1)=qr(k)
           ax(k,i+1)=ar(k)
          enddo
           do ii=1,i
            do k=1,nl
              x(k,i+1)= x(k,i+1)+del(ii)* x(k,ii)
             ax(k,i+1)=ax(k,i+1)+del(ii)*ax(k,ii)
            enddo
           enddo
        else
          do k=1,nl
            x(k,1)=qr(k)+del(1)* x(k,1)
           ax(k,1)=ar(k)+del(1)*ax(k,1)
          enddo
           do ii=2,i
            do k=1,nl
              x(k,1 )= x(k,1)+del(ii)* x(k,ii)
              x(k,ii)=0.
             ax(k,1 )=ax(k,1)+del(ii)*ax(k,ii)
             ax(k,ii)=0.
            enddo
           enddo
        endif
convergence test modes **************************************************
      if(ctest) then                                                    *
      if(nlc/ner*ner.eq.nlc) then                                       *
      ier=nlc/ner                                                       *
      eer=-1.e15                                                        *
      do 50 k=1,nl                                                      *
   50 eer=amax1(eer,abs(r(k)))                                          *
      err(ier)=eer                                                      *
      endif                                                             *
      endif                                                             *
convergence test modes **************************************************
       enddo
  100 continue
  200 continue
      eer=0.
      eem=-1.e15
      do k=1,nl
      eer=eer+qr(k)**2
      eem=amax1(eem,abs(qr(k)))
      enddo
      eer=eer/eer0
      eem=eem/eem0
      niter=nlc
      nitsm=nitsm+niter
      icount=icount+1

convergence test modes **************************************************
      if(ctest) then                                                    *
      print 301, (err(ier),ier=1,nplt,1)                                *
  301 format(4x,5e11.4)                                                 *
      do 400 ier=0,nplt                                                 *
      xitr(ier)=ier*ner                                                 *
  400 err(ier)=alog10(err(ier)*dt )                                     *
      plmx=float(itr*lord)                                              *
      call set(.1,.9,.1,.9,0.,plmx,-10.,0.,1)                           *
      call labmod('(i4)','(f5.0)',4,4,2,2,20,20,0)                      *
      call periml(4,10,5,2)                                             *
      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)                         *
      call curved(xitr,err,nplt+1)                                      *
      i1=int(102.4+409.6)                                               *
      call wtstr(cpux(i1),cpuy(50),'niter',3,0,0)                       *
      call wtstr(cpux(17),cpuy(i1),'log e',3,90,0)                      *
      call frame                                                        *
      endif                                                             *
convergence test modes **************************************************
      return
      end

      subroutine precon_1(rhs,p,r,c11,c33,d,iflg,jfl)
      parameter(nx=181,nz=121)
      parameter(n=nx,l=nz,nl=nx*nz)
      dimension rhs(n,l),p(n,l),r(n,l),
     .          c11(n,l),c33(n,l),d(n,l)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
compute available storage 
      dimension e(n,0:l-1),f(n,0:l-1)
     .       ,px(n,l),dgh(n,l),po(n,l)

      data beta/-1.e15/
      data itr,line/2,1/

      if(iflg.eq.0) then
       do i=1,nl
        p(i,1)=rhs(i,1)
       enddo
      return
      endif
    
      omg=.7
      oms=1.-omg
      dxi2=0.25*dxi*dxi
      dzi2=0.25*dzi*dzi
      do i=1,nl
       c33(i,1)=d(i,1)*dzi2
       c11(i,1)=d(i,1)*dxi2
       dgh(i,1)=0.
        po(i,1)=0.
         p(i,1)=0.
         r(i,1)=0.
      enddo
      if(line.eq.1) then
       do k=1,l
         do i=2,n-1
          dgh(i,k)=c11(i+1,k)+c11(i-1,k)
         enddo
          dgh(1,k)=c11(2,k)+c11(n-1,k)
          dgh(n,k)=c11(2,k)+c11(n-1,k)
       enddo
      endif

      if(jfl.eq.0) then
      if(line.eq.0) then
      beta=-1.e15
      do i=1,nl
      beta=amax1(beta,abs(c11(i,1))/d(i,1))
      enddo
      beta=0.5/beta
      else
      beta=1.
      endif
      return
      endif
      beti=1./beta*(1-line)

      do 100 it=1,itr
      do i=1,nl
       r(i,1)=r(i,1)+d(i,1)*(beti*p(i,1)-rhs(i,1))
     .                  +dgh(i,1)*p(i,1)
      enddo
      do i=1,n
       e(i,0)=1.
       f(i,0)=0.
       dn=d(i,1)*beti+2.*c33(i,2)+dgh(i,1)
       e(i,1)=2.*c33(i,2)/dn
       f(i,1)=     r(i,1)/dn
      enddo
      do k=2,l-1
       do i=1,n
        dn=c33(i,k+1)+c33(i,k-1)*(1.-e(i,k-2))+d(i,k)*beti
     .                                + dgh(i,k)
        e(i,k)=                      c33(i,k+1)/dn
        f(i,k)=(c33(i,k-1)*f(i,k-2)+r(i,k))/dn
       enddo
      enddo
       do i=1,n
        dn=d(i,l)*beti+2.*(1.-e(i,l-2))*c33(i,l-1)
     .                                + dgh(i,l)
        p(i,l)=(r(i,l)+2.*f(i,l-2)*c33(i,l-1))/dn
        p(i,l-1)=f(i,l-1)/(1.-e(i,l-1))
       enddo
      do k=l-2,1,-1
       do i=1,n
        p(i,k)=e(i,k)*p(i,k+2)+f(i,k)
       enddo
      enddo


      if(line.eq.1) then
       do i=1,nl
        p(i,1)=oms*po(i,1)+omg*p(i,1)
       po(i,1)=     p(i,1)
       enddo
      endif

      if(it.eq.itr) go to 101
      do k=1,l
      do i=2,n-1
      px(i,k)=c11(i,k)*(p(i+1,k)-p(i-1,k))
      enddo
      px(1,k)=c11(1,k)*(p(2,k)-p(n-1,k))
      px(n,k)=c11(n,k)*(p(2,k)-p(n-1,k))
      enddo

      do k=1,l
      do 91 i=2,n-1
   91 r(i,k)=px(i+1,k)-px(i-1,k)
      r(1,k)=(px(2,k)-px(n-1,k))
      r(n,k)=(px(2,k)-px(n-1,k))
      enddo

  100 continue
  101 continue

      return
      end


      subroutine laplc_1(p,r,u,w,n1,l1)
      dimension p(n1,l1),r(n1,l1),u(n1,l1),w(n1,l1)
      parameter(nx=181,nz=121)
      parameter(n=nx,l=nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
compute available storage in //
      parameter (nl=n*l)
      dimension px(n,l),pz(n,l)

      dxil=.5*dxi
      dzil=.5*dzi

compute pressure derivatives everywhere
      do 18 k=1,l
      do 1 i=2,n-1
    1 px(i,k)=     dxil*(p(i+1,k)-p(i-1,k))
      px(1,k)=dxil*(p(2,k)-p(n-1,k))
      px(n,k)=dxil*(p(2,k)-p(n-1,k))
   18 continue
      do 38 i=1,n
      do 3 k=2,l-1
    3 pz(i,k)=dzil*(p(i,k+1)-p(i,k-1))
      pz(i,1)= dzi*(p(i,2)-p(i,1))
   38 pz(i,l)= dzi*(p(i,l)-p(i,l-1))

compute interior pressure forces
      do 21 i=1,n
      do 10 k=2,l-1
      u(i,k)=-px(i,k)
   10 w(i,k)=-pz(i,k)
      w(i,1)=0.
      w(i,l)=0.
      u(i,1)=-px(i,1)
   21 u(i,l)=-px(i,l)

      do 99 k=1,l
      coef=rho0(k)*gac(k)
      do 99 i=1,n
      u(i,k)=coef*u(i,k)
   99 w(i,k)=coef*w(i,k)

compute laplacian
      do 911 k=1,l
      do 91  i=2,n-1
   91 r(i,k)=dxil*(u(i+1,k)-u(i-1,k))
      r(1,k)=dxil*(u(2,k)-u(n-1,k))
      r(n,k)=dxil*(u(2,k)-u(n-1,k))
  911 continue
      do 931 i=1,n
      do 93  k=2,l-1
   93 r(i,k)=r(i,k)+dzil*(w(i,k+1)-w(i,k-1))
      r(i,1)=r(i,1)+dzi *(w(i,2)-w(i,1)) 
      r(i,l)=r(i,l)+dzi *(w(i,l)-w(i,l-1))
  931 continue
      do 94 i=1,n
      do 94 k=1,l
   94 r(i,k)=-r(i,k)/(rho0(k)*gac(k))

      return
      end

      subroutine prforc_1(p,pfx,pfz,u,w,n1,n3)
      dimension p(n1,n3),u(n1,n3),w(n1,n3),pfx(n1,n3),pfz(n1,n3)
      parameter(nx=181,nz=121)
      parameter(n=nx,l=nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      parameter (nl=n*l)
      dimension px(n,l),pz(n,l)

      dxil=.5*dxi
      dzil=.5*dzi

compute pressure derivatives everywhere
      do 18 k=1,l
      do 1 i=2,n-1
    1 px(i,k)=     dxil*(p(i+1,k)-p(i-1,k))
      px(1,k)=dxil*(p(2,k)-p(n-1,k))
      px(n,k)=dxil*(p(2,k)-p(n-1,k))
   18 continue
      do 38 i=1,n
      do 3 k=2,l-1
    3 pz(i,k)=dzil*(p(i,k+1)-p(i,k-1))
      pz(i,1)= dzi*(p(i,2)-p(i,1))
   38 pz(i,l)= dzi*(p(i,l)-p(i,l-1))

compute interior pressure forces
      do 21 i=1,n
      do 10 k=2,l-1
      pfx(i,k)=u(i,k)-px(i,k)
   10 pfz(i,k)=w(i,k)-pz(i,k)
      pfz(i,1)=0.
      pfz(i,l)=0.
      pfx(i,1)=u(i,1)-px(i,1)
   21 pfx(i,l)=u(i,l)-px(i,l)

      return
      end

      subroutine rhsdiv_1(u,w,d,r,n1,l1,iflg)
      dimension u(n1,l1),w(n1,l1),d(n1,l1),r(n1,l1)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

      n=n1
      l=l1
      nl=n*l

      do 200 i=1,nl
  200 r(i,1)=0.

      dxil=.5*dxi
      dzil=.5*dzi

      do 1 i=2,n-1
      do 1 j=1,l
    1 r(i,j)=dxil*(u(i+1,j)*d(i+1,j)-u(i-1,j)*d(i-1,j))
      do 11 j=1,l
      r(1,j)=dxil*(u(2,j)*d(2,j)-u(n-1,j)*d(n-1,j))
      r(n,j)=dxil*(u(2,j)*d(2,j)-u(n-1,j)*d(n-1,j))
   11 continue
      do 3 k=2,l-1
      do 3 i=1,n
    3 r(i,k)=r(i,k)
     3        +dzil*(w(i,k+1)*d(i,k+1)-w(i,k-1)*d(i,k-1))
      do 13 i=1,n
      r(i,1)=r(i,1)+dzi*(w(i,2)*d(i,2)-w(i,1)*d(i,1))
   13 r(i,l)=r(i,l)+dzi*(w(i,l)*d(i,l)-w(i,l-1)*d(i,l-1))

      if(iflg.ne.0) then
      do 4 i=1,nl
    4 r(i,1)=iflg*r(i,1)/d(i,1)
      endif

      return
      end

       subroutine minmax(a,n,an,ax)
       dimension a(n)
       an= 1.e15
       ax=-1.e15
       do i=1,n
       an=amin1(a(i),an)
       ax=amax1(a(i),ax)
       enddo
       return
       end

      subroutine diagno_1(ux,uy,uz,th,nx,nz,scr1,scr2,rho)
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),scr1(nx,nz),scr2(nx,nz)
      dimension uy(nx,nz)
      dimension rho(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
cc pressure solver diagnostics
      common/itero/ niter,nitsm,icount,eer,eem

      nxz=nx*nz
      do i=1,nx
      do k=1,nz
      scr2(i,k)=rho(i,k)
      enddo
      enddo

      print 200, time
 200  format(1x,' ****** analysis for time (min): ',f8.2)

      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max ux: ',2e12.4)

      call minmax(uy,nxz,amn,amx)
      print 301,amn,amx
 301  format(1x,' --> min, max uy: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max uz: ',2e12.4)

      cour=0.
      do i=1,nxz
      cour=amax1(cour,abs(ux(i,1))*dt/dx+abs(uz(i,1))*dt/dz)
      enddo
      print 302,cour
 302  format(1x,' --> max courant: ',e12.4)

      call minmax(th,nxz,amn,amx)
      print 203,amn,amx
 203  format(1x,' --> min, max th: ',2e12.4)
      call rhsdiv_1(ux,uz,scr2,scr1,nx,nz,1)

      call minmax(scr1,nxz,amn,amx)
      print 204,amn,amx
 204  format(1x,' --> min, max div: ',2e12.4)

      nitav=nitsm/max0(icount,1)
      print 205, eer,eem,niter,nitav
  205 format(1x,'            eer,eem:',2e11.4/
     .       1x,'       niter, nitav:',2i4)

       if(cour.gt.1.) then
       call clsgks
       stop 'courant'
       endif

       return
       end

      subroutine diagno_2(ux,uz,th,nx,nz,scr1,scr2,rho)
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),scr1(nx,nz),scr2(nx,nz)
      dimension rho(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
cc pressure solver diagnostics
      common/itero/ niter,nitsm,icount,eer,eem

      nxz=nx*nz

      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max qv: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max qc: ',2e12.4)

      call minmax(th,nxz,amn,amx)
      print 203,amn,amx
 203  format(1x,' --> min, max qr: ',2e12.4)

       return
       end

      subroutine plot_1(ux,uy,uz,th)
      parameter(nx=181,nz=121)
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),xx(nx),zz(nz)
      dimension uy(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)
      dimension pl(nx,nz),w1(nx,nz),w2(nx,nz)

      character*50 lhead

      xl=float(nx-1)*dx/1.e3
      zl=float(nz-1)*dz/1.e3

      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo
      rat=zl/xl

cc theta
      do i=1,nx
      do k=1,nz
      pl(i,k)=th(i,k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=288.
      amx=352.
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,5,2,5)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,290.01,301.01,0.5,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,951) time
  951 format('   theta (K)      time (min): ',f6.1)
c      CALL plchmq(.525,0.93, lhead(1:50), 0.016,0.,0)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame

cc ux
      do i=1,nx
      do k=1,nz
      pl(i,k)=ux(i,k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=-5.
      amx=5.
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,5,2,5)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,-38.,38.,2.,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,952) time
  952 format('     ux (m/s)      time (min): ',f6.1)
c      CALL plchmq(.525,0.93, lhead(1:50), 0.016,0.,0)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame
ccc uy
c      do i=1,nx
c      do k=1,nz
c      pl(i,k)=uy(i,k)
c      enddo
c      enddo
c      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
c      amn=-2.
c      amx=2.
c      sn=1.
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
c      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
c      call periml(2,5,2,5)
c      call cpsetc('ILT',' ')
c      call cpseti('LLP',0)
c      call cpcnrc(pl,nx,nx,nz,-38.,38.,2.,-1,-1,-682)
c      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      write(lhead,752) time
c  752 format('     uy (m/s)      time (min): ',f6.1)
c      CALL plchmq(.525,0.93, lhead(1:50), 0.016,0.,0)
c      CALL plchhq(.525,0.08, 'X (km)', 0.016,0.,0)
c      CALL plchhq(.02,0.525, 'Z (km)', 0.016,90.,0)
c      call frame
cc uz
      do i=1,nx
      do k=1,nz
      pl(i,k)=uz(i,k)*gac(k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=-1.
      amx=1.
      sn=1.e4
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call prof(uz,height,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,5,2,5)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,-28.,28.,0.5,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,953) time
  953 format('     uz (m/s)      time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame

      return
      end

      subroutine plot_2(th,qv,qc,qr)
      parameter(nx=181,nz=121)
      dimension qv(nx,nz),qc(nx,nz),qr(nx,nz),xx(nx),zz(nz)
      dimension sc(nx,nz),th(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      character*50 lhead
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      dimension pl(nx,nz),w1(nx,nz),w2(nx,nz)


      xl=float(nx-1)*dx/1.e3
      zl=float(nz-1)*dz/1.e3

      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo
      rat=zl/xl

cc qv
      do i=1,nx
      do k=1,nz
      pl(i,k)=qv(i,k)*1.e3
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=10.
      amx=18.
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,5,2,5)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,1.,30.,1.,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,951) time
  951 format('   qv (g/kg)      time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame
cc qc
      do i=1,nx
      do k=1,nz
      pl(i,k)=qc(i,k)*1.e3
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=0.
      amx=.1
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,5,2,5)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,.1,30.,.1,-1,-1,-682)
      call cpcnrc(pl,nx,nx,nz,.01,.011,.01,-1,-1,682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,952) time
  952 format('     qc (g/kg)      time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame
ccc qr
c      do i=1,nx
c      do k=1,nz
c      pl(i,k)=qr(i,k)*1.e3
c      enddo
c      enddo
c      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
c      amn=0.
c      amx=.1
c      sn=1.
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
c      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
c      call periml(2,5,2,5)
c      call cpsetc('ILT',' ')
c      call cpseti('LLP',0)
c      call cpcnrc(pl,nx,nx,nz,.1,30.,.1,-1,-1,-682)
c      call cpcnrc(pl,nx,nx,nz,.01,.011,.01,-1,-1,682)
c      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      write(lhead,953) time
c  953 format('     qr (g/kg)      time (min): ',f6.1)
c      CALL plchmq(.525,0.93, lhead(1:50), 0.016,0.,0)
c      CALL plchhq(.525,0.02, 'x (km)', 0.016,0.,0)
c      CALL plchhq(.02,0.525, 'Z (km)', 0.016,90.,0)
c      call frame

cc rh:
      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

      do k=1,nz
      thetme=th_e(k)/tm_e(k)
      pre=1.e5*thetme**e
      do i=1,nx
      tt=th(i,k)/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      qvsw=a * esw /(pre-esw)
      sc(i,k)=qv(i,k)/qvsw * 100.
      enddo
      enddo
      do i=1,nx
      do k=1,nz
      pl(i,k)=sc(i,k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=0.0
      amx=100.0
      sn=1
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,5,2,5)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,.0,120.,10.,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,993) time
  993 format('      rh (%)        time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame

      return
      end

      subroutine prof(aa,zz,nx,nz,sn,an,ax,zl)
      dimension zz(nz),aa(nx,nz)
      dimension z(400),prf(400)

      if(nz.gt.400) stop 'dim prof'

      do k=1,nz
      z(k)=zz(k)/1000.
      prf(k)=0.
      do i=1,nx-1
      prf(k)=prf(k)+aa(i,k)*sn/float(nx-1)
      enddo
      enddo
      call set(.30,.70,.15,.9, an,ax,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(1,10,4,5)
      call curved(prf,z,nz)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.5,0.93, 'profile of ____', 0.016,0.,0)
      call frame

      return
      end

      subroutine moist_init

create parametrs for the model:

cc mass, terminal velocity, diameter
      data ar,br,cr,dr /5.2e2,3.,130.,0.5/
      data as,bs,cs,ds /2.5e-2,2.,4.,0.25/
cc collection ef., alpha, beta
      data er,alphr,betr /0.8, 1., 2./
      data es,alphs,bets /0.2, .3, 3./
cc No
      data anor,anos /2*1.e7/
cc latent heats:
      data hlatv,hlats /2.53e6,2.84e6/
cc cloud droplet concentration (per cc)
      data dconc /100./ ! must be between 50 and 2000
c      data dconc /1000./ ! must be between 50 and 2000
cc limiting temperatures
c      data tup,tdn /268.,253./   
      data tup,tdn /168.,153./   
cc gammas:
      data gamb1r,gambd1r /6.0,11.7/
      data gamb1s,gambd1s /2.0,2.56/
cc reference temperature and saturated vapor pressure:
      data tt0,ee0 /273.16,611./

      common/rain_p0/ ar,br,cr,dr,er,alphr,betr,gamb1r,gambd1r,anor
      common/rain_p1/ dconc
      common/snow_p0/ as,bs,cs,ds,es,alphs,bets,gamb1s,gambd1s,anos
      common/temp_p/ tup,tdn
      common/latent/hlatv,hlats
      common/reference/ tt0,ee0

      common /const/ gg,cp,rg,rv
      data gg,cp,rg,rv   /9.81,1005.,287.,461./

      print 2075,anor,anos,dconc
2075  format(1x,' N0 in raindrop distr.: ',e15.3/
     1 1x,' N0 in snowflake distr.: ',e15.3/
     1 1x,' Berry parameters of cloud droplet spectrum:'/
     1 1x,' droplet conc.: ',f12.4)

       return
       end

      subroutine rain_fall(qr,tm_e,rho,uza)
cc modify vertical advective velocity for rain fallout
      parameter(nx=181,nz=121)
      dimension qr(nx,nz),uza(nx,nz+1),tm_e(nz),rho(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common/rain_p0/ ar,br,cr,dr,er,alphr,betr,gamb1r,gambd1r,anor
      common/snow_p0/ as,bs,cs,ds,es,alphs,bets,gamb1s,gambd1s,anos
      common/temp_p/ tup,tdn

      real lambdr,lambds,massr,masss

cc statement functions:
      alim01(fi)=amax1(0.,amin1(1.,fi))
      comb(tm,td,tu,ad,au)=
     1  alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad

      gc3=dt*dzi

      do k=2,nz
        do i=1,nx
               dens= 0.5*(rho(k)+rho(k-1))
               qrv = 0.5*(qr(i,k)+ qr(i,k-1))
                coe_l=comb(tm_e(k),tdn,tup,0.,1.)   ! liquid part
                 qpr=qrv*coe_l         ! divide between rain and snow
                 qps=qrv-qpr           ! divide between rain and snow
         lambdr=(ar*anor*gamb1r/dens/(qpr+1.e-6))**(1./(1.+br)) 
         lambds=(as*anos*gamb1s/dens/(qps+1.e-6))**(1./(1.+bs)) 

        vtr=cr*gambd1r/gamb1r / lambdr**dr  ! terminal velocity
        vts=cs*gambd1s/gamb1s / lambds**ds  ! terminal velocity

               vtf=coe_l*vtr+(1.-coe_l)*vts   ! TERMINAL VELOCITY

c                vtf = amin1(7.,vtf)

               uza(i,k)=uza(i,k)-vtf*dens*gc3
        end do
      end do
ccc
CC LOWER AND UPPER BOUNDARIES:
cc lower:
        do i=1,nx
               dens=1.5*rho(1)-.5*rho(2)
               qrv =amax1(0.,1.5*qr(i,1)-.5*qr(i,2))
                coe_l=comb(tm_e(1),tdn,tup,0.,1.)   ! liquid part
                 qpr=qrv*coe_l         ! divide between rain and snow
                 qps=qrv-qpr           ! divide between rain and snow
         lambdr=(ar*anor*gamb1r/dens/(qpr+1.e-6))**(1./(1.+br))
         lambds=(as*anos*gamb1s/dens/(qps+1.e-6))**(1./(1.+bs))

        vtr=cr*gambd1r/gamb1r / lambdr**dr  ! terminal velocity
        vts=cs*gambd1s/gamb1s / lambds**ds  ! terminal velocity

               vtf=coe_l*vtr+(1.-coe_l)*vts   ! TERMINAL VELOCITY

                vtf = amin1(5.,vtf)

               uza(i,1)=uza(i,1)-vtf*dens*gc3
        end do
cc upper:
        do i=1,nx
               dens=1.5*rho(nz)-.5*rho(nz-1)
               qrv =amax1(0.,1.5*qr(nz,1)-.5*qr(nz-1,2))
                coe_l=comb(tm_e(nz),tdn,tup,0.,1.)   ! liquid part
                 qpr=qrv*coe_l         ! divide between rain and snow
                 qps=qrv-qpr           ! divide between rain and snow
         lambdr=(ar*anor*gamb1r/dens/(qpr+1.e-6))**(1./(1.+br))
         lambds=(as*anos*gamb1s/dens/(qps+1.e-6))**(1./(1.+bs))

        vtr=cr*gambd1r/gamb1r / lambdr**dr  ! terminal velocity
        vts=cs*gambd1s/gamb1s / lambds**ds  ! terminal velocity

               vtf=coe_l*vtr+(1.-coe_l)*vts   ! TERMINAL VELOCITY

                vtf = amin1(5.,vtf)

               uza(i,nz+1)=uza(i,nz+1)-vtf*dens*gc3
        end do

        return
        end

      subroutine thermo(th,qv,qc,qr,fth,fqv,fqc,fqr)
      parameter(nx=181,nz=121)
      parameter(n=nx,l=nz)
      dimension th(n,l),qv(n,l),qc(n,l),qr(n,l),fth(n,l),
     1         fqv(n,l),fqc(n,l),fqr(n,l)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /const/ gg,cp,rg,rv

      common/rain_p0/ ar,br,cr,dr,er,alphr,betr,gamb1r,gambd1r,anor
      common/rain_p1/ dconc
      common/snow_p0/ as,bs,cs,ds,es,alphs,bets,gamb1s,gambd1s,anos
      common/temp_p/ tup,tdn
      common/latent/hlatv,hlats
      common/reference/ tt0,ee0
      real lambdr,lambds,massr,masss

cc statement functions:
      alim01(fi)=amax1(0.,amin1(1.,fi))
      comb(tm,td,tu,ad,au)=
     1  alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad

condensation/evaporation

      pi=3.1416
      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

      do 100 k=1,l
      thetme=th_e(k)/tm_e(k)
      do 100 i=1,n
      coe_l=comb(tm_e(k),tdn,tup,0.,1.)   ! liquid contribution
      pre=1.e5*thetme**e
      tt=th(i,k)/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      esi=ee0*exp(b * delt)
      qvsw=a * esw /(pre-esw)
      qvsi=a * esi /(pre-esi)
      qvs=coe_l*qvsw + (1.-coe_l)*qvsi
ccc linearized condensation rate is next:
      cf1=thetme/th(i,k)
      cf1=cf1*cf1
      cf1=c*cf1*pre/(pre-esw)*d
      delta=(qv(i,k)-qvs)/(1.+qvs*cf1)
c--->
ccc one Newton-Raphson iteration is next:
      thn=th(i,k)+c*thetme*delta
      tt=thn/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      esi=ee0*exp(b * delt)
      qvsw=a * esw /(pre-esw)
      qvsi=a * esi /(pre-esi)
      qvs=coe_l*qvsw + (1.-coe_l)*qvsi
      fff=qv(i,k)-delta-qvs
      cf1=thetme/thn
      cf1=cf1*cf1
      cf1=c*cf1*pre/(pre-esw)*d
      fffp=-1.-qvs*cf1
      delta=delta-fff/fffp
ccc end of the iteration; if required, it can be repeated
c--->
      delta=amin1( qv(i,k), amax1(-qc(i,k),delta) )
      qv(i,k)=qv(i,k)-delta
      qc(i,k)=qc(i,k)+delta
      th(i,k)=th(i,k)+c*thetme*delta

      delta=amin1( qv(i,k), amax1(-qc(i,k),delta) )
      fqv(i,k)=-delta*2.*dti
      fth(i,k)=-c*thetme*fqv(i,k)
      fqc(i,k)=-fqv(i,k)

  100 continue

ccc remove trace of water variables:
      nl=n*l
      do i=1,nl
      qc(i,1)=cvmgm(0.,qc(i,1),qc(i,1)-1.e-9)
      qr(i,1)=cvmgm(0.,qr(i,1),qr(i,1)-1.e-10)
      enddo


compute moist forces update
      do 300 k=1,l
      thetme=th_e(k)/tm_e(k)
      do 300 i=1,n
       tt=th(i,k)/thetme
       pre=1.e5*thetme**e
       coe_l=comb(tm_e(k),tdn,tup,0.,1.)   ! liquid contribution
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      esi=ee0*exp(b * delt)
      qvsw=a * esw /(pre-esw)
      qvsi=a * esi /(pre-esi)

       ssw=qv(i,k) / qvsw      ! saturation ratio
       ssi=qv(i,k) / qvsi      ! saturation ratio

      qpr=qr(i,k)*coe_l                ! divide between rain and snow
      qps=qr(i,k)-qpr                  ! divide between rain and snow
      qcc=qc(i,k)*coe_l                ! divide between ice and water
      qci=qc(i,k)-qcc                  ! divide between ice and water

      lambdr=(ar*anor*gamb1r/rho0(k)/(qpr+1.e-6))**(1./(1.+br)) ! lambda
      lambds=(as*anos*gamb1s/rho0(k)/(qps+1.e-6))**(1./(1.+bs)) ! lambda

CC AUTOCONVERSION:
cc Kogan 2013:
       autc=7.98e10 * qcc**4.22 * dconc**(-3.01)
cc snow G98:
       tc=tt-273.16
       times=amin1(1.e3,(3.56*tc+106.7)*tc+1.e3) ! time scale for
      auti=qci/times
      AUT = autc + auti

cccccccccccccccccccccccccccccccccccccc
ccc NO RAIN:
           AUT=0.
cccccccccccccccccccccccccccccccccccccc

CC GROWTH:
      conr=anor/lambdr ! concentration
      cons=anos/lambds ! concentration

      massr=rho0(k)*(qpr+1.e-7) / conr  ! mass
      masss=rho0(k)*(qps+1.e-7) / cons  ! mass

      diamr=(massr/ar)**(1./br) ! diameter
      diams=(masss/as)**(1./bs) ! diameter

      rer=cr*diamr**(dr+1.)/2.e-5  ! Reynolds number
      res=cs*diams**(ds+1.)/2.e-5  ! Reynolds number

      ventr=amax1(1.,.78+.27*sqrt(rer))  ! ventilation factor
      vents=amax1(1.,.65+.39*sqrt(res))  ! ventilation factor

      thfun=1.e-7/(2.2*tm_e(k)/esw+2.2e-2/tm_e(k))  ! thermodynamic fun.

      g_acc_r=pi/4.*cr*diamr**(2.+dr)*er*alphr*rho0(k)*qc(i,k) ! growth
      g_acc_s=pi/4.*cs*diams**(2.+ds)*es*alphs*rho0(k)*qc(i,k) ! growth

      g_dep_r=4.*pi*diamr/betr*(ssw-1.)*ventr*thfun   ! growth/evap
      g_dep_s=4.*pi*diams/bets*(ssi-1.)*vents*thfun   ! growth/evap

      acc_r=conr * g_acc_r * qpr / (qpr + 1.e-9)
      acc_s=cons * g_acc_s * qps / (qps + 1.e-9)

       ACC= acc_r + acc_s  ! growth by accretion

      dep_r=conr * g_dep_r * qpr / (qpr + 1.e-9)
      dep_s=cons * g_dep_s * qps / (qps + 1.e-9)

       DEP= dep_r + dep_s  ! growth by deposition

      dcol=2.*(AUT + ACC)
      dcol=amin1(dcol,  2.*dti*qc(i,k)+fqc(i,k))
      devp=2.*DEP
      devp=amax1(devp, -2.*dti*qr(i,k)-dcol)
cc
      fqr(i,k)=devp+dcol
      fqc(i,k)=fqc(i,k)-dcol
      fqv(i,k)=fqv(i,k)-devp
      fth(i,k)=fth(i,k)+c*devp*thetme

  300 continue

      return
      end

      subroutine absor(ux,uz,th,qv,qc,qr,fth,fqv,fqc,fqr)
      parameter(nx=181,nz=121)
      parameter(n=nx,l=nz)
      dimension ux(n,l),uz(n,l)
      dimension th(n,l),qv(n,l),qc(n,l),qr(n,l),fth(n,l),
     1         fqv(n,l),fqc(n,l),fqr(n,l)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common /prof_a/ tau(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

      do k=1,nz
      coe1=.5*dt*tau(k)
      coe2=1.+coe1
      do i=1,nx
      ux(i,k)=(ux(i,k)+coe1*ux_e(k))/coe2
      uz(i,k)=(uz(i,k)             )/coe2
      th(i,k)=(th(i,k)+coe1*th_e(k))/coe2
      qv(i,k)=(qv(i,k)+coe1*qv_e(k))/coe2
      qc(i,k)=(qc(i,k)             )/coe2
      qr(i,k)=(qr(i,k)             )/coe2

      fth(i,k)=fth(i,k) - tau(k)*(th(i,k)-th_e(k))
      fqv(i,k)=fqv(i,k) - tau(k)*(qv(i,k)-qv_e(k))
      fqc(i,k)=fqc(i,k) - tau(k)*(qc(i,k)-     0.)
      fqr(i,k)=fqr(i,k) - tau(k)*(qr(i,k)-     0.)
      enddo
      enddo
      return
      end

      subroutine integz(a,b,n1,n2)
      dimension a(n1,n2),b(n1,n2)

cc z smoothing:
      do k=2,n2-1
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do i=1,n1
      b(i,1 )=0.5*(a(i, 2)+a(i,   1))
      b(i,n2)=0.5*(a(i,n2)+a(i,n2-1))
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
      end

      subroutine integxz(a,b,n1,n2)
      dimension a(n1,n2),b(n1,n2)

cc z smoothing:
      do k=2,n2-1
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do i=1,n1
      b(i,1 )=0.5*(a(i, 2)+a(i,   1))
      b(i,n2)=0.5*(a(i,n2)+a(i,n2-1))
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

cc x smoothing:
      do i=1,n1
      im=i-1
      if(im.eq.0) im=n1-1
      ip=i+1
      if(ip.eq.n1+1) ip=2
      do k=1,n2
      b(i,k)=0.25*(a(im,k)+2.*a(i,k)+a(ip,k))
      enddo
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
      end

      subroutine prof_init_1layer  
      parameter(nx=181,nz=121)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)

      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)

      common/temp_p/ tup,tdn
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

c stability for the virtual potential temperature stab and relative humidity relhum
      parameter(stab=1.3e-5,relhum=0.2)

      dimension pres(nz)

      a=rg/rv
      c=hlatv/cp
      d=hlatv/rv
      e=-cp/rg

      cap=rg/cp
      capi=1./cap

cc surface data:
      zz=height(1)
      tm_e(1)=283.
      pres(1)=850.e2
      tt=tm_e(1)
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      qvs=a * esw /(pres(1)-esw)
      qv_e(1)=relhum*qvs

cc surface values for dry profiles:
      rho0(1)=pres(1)/rg/tm_e(1)
      th0(1)=tm_e(1)*(1.e5/pres(1))**(cap)

      th_e(1)=th0(1)*(1.+a*qv_e(1))

         k=1
       print*,'z,the,tme,qve,p: ',zz,th_e(k),tm_e(k),qv_e(k),pres(k)

cc move upward:
       do k=2,nz
       zz=height(k)

       th_e(k)=th_e(1)*exp(stab*zz)

c predictor:
       rhob=pres(k-1)/rg/(tm_e(k-1)*(1.+a*qv_e(k-1)))
       pres(k)=pres(k-1) - gg*rhob*dz
cc iteration for T and qv:
       qv_e(k)=qv_e(k-1)
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               do iter=1,4
          tt=tm_e(k)
          delt=(tt-tt0)/(tt*tt0)
          esw=ee0*exp(d * delt)
          qvs=a * esw /(pres(k)-esw)
          qv_e(k)=relhum*qvs
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               enddo
         
c corrector:
       rhon=pres(k)/rg/(tm_e(k)*(1.+a*qv_e(k)))
       pres(k)=pres(k-1) - gg*.5*(rhob+rhon)*dz
cc iteration for T and qv:
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               do iter=1,4
          tt=tm_e(k)
          delt=(tt-tt0)/(tt*tt0)
          esw=ee0*exp(d * delt)
          qvs=a * esw /(pres(k)-esw)
          qv_e(k)=relhum*qvs
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               enddo
         
       print*,'z,the,tme,qve,p: ',zz,th_e(k),tm_e(k),qv_e(k),pres(k)

              enddo

cc set constant-stability dry profiles:
        sum=0.
        do k=2,nz-1
          sum = sum + (th_e(k+1)-th_e(k-1))/th_e(k)
        enddo
      st=sum/(float(nz-2)*2.*dz)
       print*,'checking stability: ',stab,st
cc compute reference state vertical profiles
      cap=rg/cp
      capi=1./cap
      cs=gg/(cp*tm_e(1)*st)
      zz=0.
      print*,'z,th0,rho0: ',zz,th0(1),rho0(1)
      do k=2,nz
      zz=height(k)
      exs=exp(-st*zz)
      th0(k)=th0(1)/exs
      rho0(k)=rho0(1)*exs*(1.-cs*(1.-exs))**(capi-1.)
      print*,'z,th0,rho0: ',zz,th0(k),rho0(k)
      enddo

      return
      end

      subroutine prof_init_dry
      parameter(nx=181,nz=121)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)

      dimension xx(nx),zz(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)

      common/temp_p/ tup,tdn
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

      parameter(npin=2)
      dimension press(npin),temp(npin),zin(npin),relhum(npin),
     1          vap(npin),uu(npin),vv(npin),zinkm(npin)

       DATA ZIN   /  0. ,2.4e3 /
       DATA TEMP  /300., 300. /
       DATA UU /NPIN*0./
       DATA VV /NPIN*0./

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cc surface data:
      iisn=1
      zz(1)=0.
      tm_e(1)=temp(iisn)
      th_e(1)=tm_e(1) 
      qv_e(1)=0.
      ux_e(1)=uu(1)
      uy_e(1)=vv(1)
c      print*,'DETERMINED SURFACE DATA'
cc higher levels - interpolate:
c     print*,'INTERPOLATION TO HIGHER LEVELS'
      l=nz
      do k=2,l
      zz(k)=height(k)
c       print*,'k,z= ',k,zz(k)
        do kk=2,npin
          iisn=kk-1
          if(zin(kk).ge.zz(k)) go to 665
        enddo
c       print*,'INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.'
        stop 'SOUNDING'
 665    continue
c       print*,'iisn=',iisn
        coe2=(zz(k)-zin(iisn))/(zin(iisn+1)-zin(iisn))
        tm_e(k)=coe2*temp(iisn+1) + (1.-coe2)*temp(iisn)
        th_e(k)=tm_e(k)
        qv_e(k)=0.
        ux_e(k)=coe2*uu(iisn+1) + (1.-coe2)*uu(iisn)
        uy_e(k)=coe2*vv(iisn+1) + (1.-coe2)*vv(iisn)
      end do

compute th00,tt00,pr00,rh00 and average stability for base state profiles
      th00=th_e(1)
      tt00=tm_e(1)
      rh00=1.e5/(rg*tt00)
      print*,'th00,tt00,pr00,rh00,st: ',th00,tt00,pr00,rh00,st

compute reference state vertical profiles 
      do k=1,l
      th0(k)=th00
      rho0(k)=rh00
      enddo

      print*,'PROFILES'
      do k=1,l
       print 200,zz(k)/1.e3,th0(k),rho0(k),th_e(k),
     .             tm_e(k),qv_e(k)*1.e3,ux_e(k)
 200    format(1x,'z,th0,rho0,the,tme,qve,ue:',
     .          2f7.1,f5.2,2f7.1,e10.3,f6.1)
      enddo

      return
      end


      subroutine surfflux(theta,qv,fth,fqv,ux,uy,fx)
cc
cc this a routine that calculates surface flux and distributes
cc the flux in the vertical
cc    spatial fluctuations are allowed to exist, see the logic below...
cc
      parameter(nx=181,nz=121)
      dimension theta(nx,nz),qv(nx,nz),ux(nx,nz),uy(nx,nz)      
      dimension fth(nx,nz),fqv(nx,nz),fx(nx,nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)

      dimension ss1(nx,nz),ss2(nx,nz),scr(nx,nz)
      dimension pthsfl(nx),pqvsfl(nx),pusfl(nx),pvsfl(nx)
      common /bl_param/ blh

      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

      dimension ventfc(2)
      common /vent_hist/ ventfc
      dimension uua(nz),vva(nz),tha(nz),qva(nz)
      dimension thsflx(nx),qvsflx(nx)

c      double precision rand

cc
      if(nz.gt.300) stop 'surfflux'

cc orig (Klaasen/Clark):
c            thsfl = 150. /(rho0(1)*cp) 
c            qvsfl = 0. / (rho0(1)*hlatv)
cc mod1:
c            thsfl = 50. /(rho0(1)*cp) 
c            qvsfl = 100. / (rho0(1)*hlatv)
cc mod2:
            thsfl =350. /(rho0(1)*cp) 
            qvsfl =500. / (rho0(1)*hlatv)

      h_dis=300.
      xct=.5*float(nx-1)*dx
      sigma=700.

      do k=1,nz
cc note shift of flux positions:
      zz=(k-1)*dz+.5*dz
      expfun=exp(-zz/h_dis)
      fun=expfun

      do i=1,nx
      xx=float(i-1)*dx
      fun1=exp(- (xx-xct)**2/sigma**2 )
      ss1(i,k)=thsfl * fun * fun1 
      ss2(i,k)=qvsfl * fun * fun1 
      thsflx(i)=thsfl * fun1
      qvsflx(i)=qvsfl * fun1
      enddo
      enddo

      do i=1,nx-1
cc first level above the ground:
      fth(i,1)=fth(i,1)-2.*(ss1(i,1)-thsflx(i))/(0.5*dz)
      fqv(i,1)=fqv(i,1)-2.*(ss2(i,1)-qvsflx(i))/(0.5*dz)
cc higher levels:
      do k=2,nz-1
      fth(i,k)=fth(i,k)-2.*(ss1(i,k)-ss1(i,k-1))/dz
      fqv(i,k)=fqv(i,k)-2.*(ss2(i,k)-ss2(i,k-1))/dz
      enddo
      enddo

cc cyclicity:
       do k=1,nz
       fth(nx,k)=fth(1,k)
       fqv(nx,k)=fqv(1,k)
       enddo


      return
      end

      subroutine zstrtch(zz,nz1,dz)
cc define vertical grid for stretched coordinate, derive Jacobian
      dimension zz(nz1)
      parameter(nx=181,nz=121)
      common /strtch/ height(nz),gac(nz)
      dimension zb(400)
      if(nz.gt.400) stop 'grid formulation'

      top=float(nz-1)*dz

cc exponent of stretching function:
c      ex1=4./3.
c      ex1=4.1/2.9
      ex1=1.
cc
       if(ex1.ne.1.) stop 'fix surface flux for strtching'
cc

      aa=top / top**ex1

cc vertical coordinate:
      do k=1,nz
      zb(k)=aa*zz(k)**ex1
      height(k)=zb(k)
      enddo

cc jacobian:
      do k=1,nz
      if(k.eq.1 ) gac(k)=(zb(2)-zb(1))/dz
      if(k.eq.nz) gac(k)=(zb(nz)-zb(nz-1))/dz
      if(k.ne.1 .and. k.ne.nz) gac(k)=(zb(k+1)-zb(k-1))/(2.*dz)
      print*,k,zb(k),gac(k)
      enddo

cc check consistency (int gac ddzeta = H)
      sum1=.5*(gac(1)*dz + gac(nz)*dz)
      do i=2,nz-1
      sum1=sum1+gac(i)*dz
      enddo
      print*,'int Jacobian before adj: ',sum1

cc adjust numerical jacobian:
      coe=float(nz-1)*dz/sum1
      do i=1,nz
      gac(i)=gac(i)*coe
      enddo

cc check:
      sum1=.5*(gac(1)*dz + gac(nz)*dz)
      do i=2,nz-1
      sum1=sum1+gac(i)*dz
      enddo
      print*,'int Jacobian after adj: ',sum1

      return
      end


      function cvmgm(a,ab,ac)
       if(ac.lt.0) then
        cvmgm=a
       else
        cvmgm=ab
       endif
       return
       end

      subroutine gridint
     1 (ain,nx,nz,xx,zz,nx1,nz1,xx1,zz1,x0,z0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc this subroutine performs interpolation of the data given in the
cc input array ain(nx,nz) on the grid defined by xx(nx) and zz(nz)
cc into the grid given by xx1(nx1) and zz1(nz1). NIETHER OF THE
cc GRIDS HAS TO BE REGULAR. Data is returned in the ain(nx1,nz1)
cc part of the input array.
cc 
cc    levels in the input array are given in zz(nz), 
cc    levels in the output array are given in zz1(nz1)
cc      x-coordinate in the input array are in xx(nx)
cc      x-coordinate in the output array are in xx1(nx1)
cc        x0(nx,nz) and z0(nx,nz) are working arrays
cc
cc NOTE that nx1 (or nz1) must be smaller than nx (or nz) and xx1 (zz1)
cc  must be a subdomain of xx (zz)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension ain(nx,nz),zz(nz),xx(nx)
      dimension zz1(nz1),xx1(nx1) 
      dimension x0(nx,nz),z0(nx,nz)

c      check consistency of the input data
      ier=0
cc nx1,nz1 not larger than nx,nz
      if(nx1.gt.nx) ier=1
      if(nz1.gt.nz) ier=2
cc limits (zz1(nz1).le.zz(nz) ?) 
      if(zz1(1).lt.zz(1)) ier=3
      if(zz1(nz1).gt.zz(nz)) ier=4
cc limits (xx1(nx1).le.xx(nx) ?) 
      if(xx1(1).lt.xx(1)) ier=5
      if(xx1(nx1).gt.xx(nx)) ier=6
      if(ier.ne.0) then
      print 999,ier
 999  format(2x,' ** problems with input data. will stop.'/
     1 ' ier = ',i3,'. see code why stoped.')
      stop
      endif
cc
      nxz=nx*nz
      do 99 i=1,nxz
      z0(i,1)=1.
  99  x0(i,1)=1.
cc  map vertical grid positions:
      do 1 k1=1,nz1
      zzh=zz1(k1)
      do 2 k=1,nz
      kk=k
      if(zz(k).ge.zzh) go to 6
  2   continue
  6   kkm=max0(1,kk-1)
      z0(1,k1)=float(kkm)+(zzh-zz(kkm))/(zz(kk)-zz(kkm)+1.e-6)
  1   continue
      do 3 i1=2,nx1
      do 3 k1=1,nz1
  3   z0(i1,k1)=z0(1,k1)
c
cc  map horizontal grid positions:
      do 11 i1=1,nx1
      xxh=xx1(i1)
      do 12 i=1,nx
      ii=i
      if(xx(i).ge.xxh) go to 16
 12   continue
 16   iim=max0(1,ii-1)
      x0(i1,1)=float(iim)+(xxh-xx(iim))/(xx(ii)-xx(iim)+1.e-6)
 11   continue
      do 13 i1=1,nx1
      do 13 k1=2,nz1
 13   x0(i1,k1)=x0(i1,1)
cc
cc  call Piotr's interpolation routine
      call inter2(ain,x0,z0,nx,nz)
cc
      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INTER2(XF,XD1,XD2,NX,NZ)
C IOR=ORDER OF ACCURACY/2; ONLY EVEN ORDER TRMBACK SCHEMES ARE CONSIDERED

      PARAMETER(IOR=2)
      PARAMETER(LINER=0)
      DIMENSION XF(*),XD1(*),XD2(*)
CC  N1 - HORIZONTAL INDEX, N2 - VERTICAL INDEX
      PARAMETER (N1=181, N2=121,NN=N1*N2)
      DIMENSION Z(NN,-IOR:IOR)
      DATA  EP/ 1.E-10/
      PARAMETER(NONOS=1)
      REAL      MX,MN
      PARAMETER(IBC=0)
      COMMON // IG0(NN),JG0(NN),X(-IOR+1:N1+IOR,-IOR+1:N2+IOR)
C  next is for shavano: 
C      DONOR(Y1,Y2,A)=CVMGM(Y2,Y1,A)*A
C  next is for workstation:
      DONOR(Y1,Y2,A)=AMAX1(0.,A)*Y1 + AMIN1(0.,A)*Y2
      TR2(Y1,Y2,A)=A*.5*(Y1+Y2)-A**2*.5*(Y2-Y1)
      TR4(YM1,Y0,YP1,YP2,A)=A/12.*(7.*(YP1+Y0)-(YP2+YM1))
     1 -A**2/24.*(15.*(YP1-Y0)-(YP2-YM1))-A**3/12.*((YP1+Y0)
     2 -(YP2+YM1))+A**4/24.*(3.*(YP1-Y0)-(YP2-YM1))
      TR6(YM2,YM1,Y0,YP1,YP2,YP3,A)=-A/60.*(-YM2+8.*YM1-37.*Y0
     1                                     -37.*YP1+8.*YP2-YP3)
     2-A**2/360.*(-2.*YM2+25.*YM1-245.*Y0+245.*YP1-25.*YP2+2.*YP3)
     3-A**3/48.*(YM2-7.*YM1+6.*Y0+6.*YP1-7.*YP2+YP3)
     4-A**4/144.*(YM2-11.*YM1+28.*Y0-28.*YP1+11.*YP2-YP3)
     5-A**5/240.*(-YM2+3.*YM1-2.*Y0-2.*YP1+3.*YP2-YP3)
     6-A**6/720.*(-YM2+5.*YM1-10.*Y0+10.*YP1-5.*YP2+YP3)
      PP(XI)=AMAX1(0.,XI)
      PN(XI)=AMIN1(0.,XI)
C
CC CHECK COSISTENCY OF THE DATA:
      IF(NX.NE.N1.OR.NZ.NE.N2) THEN
      PRINT 777
 777  FORMAT(2X,'!!! CALLS TO INTER2 WITH NON-MATCHING DIMENSIONS.'
     1  ,' STOP.')
      STOP
      ENDIF
CC
      DO 1 K=1,NN
      IG0(K)=NINT(XD1(K))
    1 JG0(K)=NINT(XD2(K))

C  GRID EXTENSION FOR BC REMOVAL 
      DO 508 I=1,N1      
      DO 509 J=1,N2
      II=(J-1)*N1+I
  509 X(I,J)=XF(II) 
      DO 5091 IS=1-IOR,0
C     II=(1-1)*N1+I
 5091 X(I,IS)=XF(I)
      DO 5092 IS=1,IOR
      II=(N2-1)*N1+I
 5092 X(I,N2+IS)=XF(II)
  508 CONTINUE
      DO 507 J=-IOR+1,N2+IOR
      DO 5071 IS=-IOR+1,1
 5071 X(IS,J)=X(1,J)*(1-IBC)+IBC*X(N1+IS-1,J)
      DO 5072 IS=0,IOR
 5072 X(N1+IS,J)=X(N1,J)*(1-IBC)+IBC*X(1+IS,J)
  507 CONTINUE
C  END OF GRID EXTENSION
C                     
C
C  HERE STARTS REZIDUAL ADVECTION
C
                     DO 50 J=-IOR,IOR
C
      IF(LINER.EQ.1) THEN
      DO 211 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
  211 Z(II,J)=Y0-(FL1-FL0) 
      GO TO 50
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 311 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  311 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 321 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
  321 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 312 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  312 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 322 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
  322 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 313 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  313 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 323 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
  323 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
C
   50 CONTINUE
C  
      IF(LINER.EQ.1) THEN
      DO 212 II=1,NN
      U=JG0(II)-XD2(II)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
  212 XF(II)=Z(II,0)-(FL1-FL0) 
      RETURN
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 411 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  411 CONTINUE
        ELSE
      DO 421 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
  421 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 412 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  412 CONTINUE
        ELSE
      DO 422 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
  422 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 413 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  413 CONTINUE
        ELSE
      DO 423 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
  423 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF
      RETURN
      END   

