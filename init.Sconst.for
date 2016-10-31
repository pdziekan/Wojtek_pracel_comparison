c      program initialize
      subroutine prof_init_stab_th_v
c      parameter(nx=181,nz=121)
c      common/grid/ time,dt,dx,dz,dti,dxi,dzi
c      common /strtch/ height(nz),gac(nz)
c
c      dimension xx(nx),zz(nz)
c      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
c      common /prof_m/ qv_e(nz),tm_e(nz)
c
c      common/temp_p/ tup,tdn
c      common /const/ gg,cp,rg,rv
c      common/reference/ tt0,ee0
c      common/latent/hlatv,hlats

c calculation of atmospheric profiles for model initialization assuming constant
c stability for the virtual potential temperature and constant relative humidity

c grid in the vertical, here the domain is 2.4 km deep
      parameter(nz=121,dz=20.)
c      parameter(nz=2401,dz=1.)

c stability for the virtual potential temperature stab and relative humidity relhum
      parameter(stab=1.3e-5,relhum=0.2)

cc base state (rho0 and th0; dry) and environmnetal (th_e, tm_e, qv_e; moist) profiles
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      dimension pres(nz)

      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

cc constants - as in the model:
      data gg,cp,rg,rv   /9.81,1005.,287.,461./
cc reference temperature and saturated vapor pressure:
      data tt0,ee0 /273.16,611./
cc latent heats:
      data hlatv,hlats /2.53e6,2.84e6/

      a=rg/rv
      c=hlatv/cp
      d=hlatv/rv
      e=-cp/rg

      cap=rg/cp
      capi=1./cap

cc surface data:
      zz=0.
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
         print*,'z,the,tme,qve: ',zz,th_e(k),tm_e(k),qv_e(k)

cc move upward:
       do k=2,nz
       zz=(k-1)*dz

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
         
         print*,'z,the,tme,qve: ',zz,th_e(k),tm_e(k),qv_e(k)

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
      zz=(k-1)*dz
      exs=exp(-st*zz)
      th0(k)=th0(1)/exs
      rho0(k)=rho0(1)*exs*(1.-cs*(1.-exs))**(capi-1.)
      print*,'z,th0,rho0: ',zz,th0(k),rho0(k)
      enddo

      rho0(1)=pres(1)/rg/tm_e(1)
      th0(1)=tm_e(1)*(1.e5/pres(1))**(cap)
              stop
              end

