      program parcel
cc ! this program calculates thermodynamic parameters inside
cc ! a simple adiabatic parcel rising in the atmosphere
cc ! only condensation is considered
cc ! SI units are used throughout

cc BOMEX sounding...

C nz, dz - numer of levels in the vertical, vertical spacing
C note that vertical extent is (nz-1)*dz
      parameter(nz=241,dz=10.)

C arrays: environmental profiles: z is height, pre is pressure,
C         te is temperature, the is potential temperature,
C         qve is environmental water vapor mixing ratio
C         parcel parameters: th is potential temperature,
C         qv is water vapor mixing ratio, qc is cloud water
C         mixing ratio, t is temperature
      dimension z(nz),pre(nz),te(nz),the(nz),qve(nz)
      dimension t(nz),th(nz),qv(nz),qc(nz),tbuo(nz)
      dimension cape(nz)

cc for ploting (and output)
      dimension pl1(nz),pl2(nz),x(2),y(2),rho(nz)

C model constants: 
      data hlat,gg,rd,rv,cp /2.5e6,9.81,287.,461.,1005./

C base values for saturated vapor pressure:
c      data e00,t00 /2337.,293.16/
      data e00,t00 /611.,273.16/
c      data e00,t00 /1227.,283.16/

cc data to get Grabow/Clark JAS 1999 profiles
      parameter(npin=2)
      dimension press(npin),temp(npin),theta(npin),zin(npin),
     1          vap(npin),zinkm(npin)

       data zinkm  /
     1    0.0,2.4/
       data press  /
     1    850.,631.37/
       data theta  /
     1    296.78, 306.19/
       data vap  /
     1    1.8262, .8471 /

C saturated vapor pressure statement function:
      esat(t11)=e00*exp(hlat/rv*(1./t00 - 1./t11))

        do k=1,npin
        zin(k)=zinkm(k)*1.e3
        enddo

ccc pressure levels:
c        press(1)=1015.
c        rc=rd/cp
c        rci=1./rc
c        do k=2,npin
c          km=k-1
c          tempk =theta(k )*(1. + .61*vap(k )*1.e-3)
c          tempkm=theta(km)*(1. + .61*vap(km)*1.e-3)
c          thetavm=.5*(tempk+tempkm)
c          press(k)=( press(km)**rc -
c     .             gg*(1.e3)**rc*(zin(k)-zin(km))/(cp*thetavm) )**rci
c         enddo
c
ccc temperature; convert qv to SI units:
          do k=1,npin
          temp(k)=theta(k)*(press(k)/1.e3)**(rd/cp)
          vap(k)=vap(k)*1.e-3
          enddo
         
          do k=1,npin
          print*,'z,p,th,t,qv: ',zin(k),press(k),theta(k),temp(k),vap(k)
          enddo

C interpolate profiles
      z(1)=0.
      pre(1)=press(1)*1.e2
      te(1)=temp(1)
      the(1)=theta(1)
      qve(1)=vap(1)
          do k=2,nz
        z(k)=z(k-1)+dz
           do kk=2,npin
          iisn=kk-1
          if(zin(kk).ge.z(k)) go to 665
        enddo
        print*,'INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.'
        stop 'SOUNDING'
 665    continue
        coe2=(z(k)-zin(iisn))/(zin(iisn+1)-zin(iisn))
        te(k)=coe2*temp(iisn+1) + (1.-coe2)*temp(iisn)
        the(k)=coe2*theta(iisn+1) + (1.-coe2)*theta(iisn)
        pre(k)=1.e5*(te(k)/the(k))**(cp/rd)
        qve(k)=coe2*vap(iisn+1) + (1.-coe2)*vap(iisn)
        enddo

        do k=1,nz
        tv=te(k)*(1.+.61*qve(k))
        rho(k)=pre(k)/(rd*tv)
        enddo

C test printout:
         print*,' ***** profiles: *****'
         do k=1,nz
         print 102, z(k),te(k),pre(k),the(k),qve(k)
102      format(1x,5e14.6,'    z,t,p,th,qv')
         pl1(k)=qve(k)*1.e3
         enddo

cc plot theta and qv profiles:

c      call opngks
c      call gsclip(0)
c
c      call setusv('LW',2000)
c      call set(.3,.51,.55,.95,290.,310.,0.,2400.,1)
c      call labmod('(f5.0)','(f5.0)',5,5,2,2,28,20,0)
c      call perim (3,2,3,4)
c      call dashdc('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$', 3,12)
c      call curved(the,z,nz)
c
c      call set(.52,.73,.55,.95,0.,2.,0.,2400.,1)
c      call labmod('(f5.1)','(f5.0)',5,5,2,2,28,20,0)
c      call perim (2,5,3,4)
c      call dashdc('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$', 3,12)
c      call curved(pl1,z,nz)
c
c      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      CALL plchhq(.22,0.75,'height (km)',0.012,90.,0)
c      CALL plchhq(.26,0.55,'0.0',0.012,0.,0)
c      CALL plchhq(.26,0.95,'2.4',0.012,0.,0)
c
c      CALL plchhq(.40,0.92,
c     1   ':PGU:H:PRL: :PRU:(K)', 0.014,0.,0)
c      CALL plchhq(.31,0.53,'290',0.012,0.,0)
c      CALL plchhq(.50,0.53,'310',0.012,0.,0)
c
c      CALL plchhq(.65,0.92,'q:B:v:N: (g/kg)', 0.014,0.,0)
c      CALL plchhq(.53,0.53,'0',0.012,0.,0)
c      CALL plchhq(.72,0.53,'2',0.012,0.,0)

C first level parameters
         k0=81
         th(k0)=the(k0)
         t(k0)=te(k0)
         qv(k0)=qve(k0) / 0.2   ! RH=20%
         qc(k0)=0.

C main loop 
              do k=k0+1,nz

C first step: move one level up:
              th(k)=th(k-1)
              qv(k)=qv(k-1)
              qc(k)=qc(k-1)

          t(k)=th(k)*te(k)/the(k)  !<-- convert pot. temp into temp

C check if adjustement is required
         es=esat(t(k))
         qvs=0.621*es/(pre(k)-es)

C  no adjustement required if conditions below saturation and there is
C  no cloud water to evaporate
       if(qv(k).le.qvs.and.qc(k).lt.1.e-15) go to 100  

C  adjustement procedure:
C    first guess
       tempk=t(k)
c if the next line is open, the first guess is delta=0 each time
C however, it makes more sense to leave delta from previous level
c to be the first guess, so the line is commented out
       delta=0.
C    iterate:
                   do it=1,10
Ccccccccccccccccccccccccccccccc
Ccc below is printout to see convergence:
C          if(k.gt.9.and.k.lt.21) then
C            print*,'** k,it,delta: ',k,it,delta
C          endif
Cccccccccccccccccccccccccccccccccc
         es=esat(tempk)
         qvs=0.621*es/(pre(k)-es) 
          top=qvs-qv(k)+delta
          bottom=.621*pre(k)/(pre(k)-es)**2  * hlat*es/(rv*tempk**2) 
     1     * hlat/cp   +  1. 

         delta = delta - top/bottom
         tempk=t(k)+hlat/cp*delta
                 enddo

C final adjustement:
          delta = amax1(delta,-qc(k)) ! <-- limit for evaporation
          th(k)=th(k) +  hlat/cp * the(k)/te(k) *delta
          qv(k)=qv(k)-delta
          qc(k)=qc(k)+delta
         
          t(k)=th(k)*te(k)/the(k)  !<-- convert pot. temp into temp

 100   continue
              enddo   ! <-- end of main loop 

C test printout:
c         print*,' ***** parcel parameters: *****'
c         do k=1,nz
c         print 103, k,z(k),th(k),qv(k),qc(k)
c103      format(1x,'k,z,th,qv,qc: ',i4,4e12.4)
c         enddo

         do k=1,nz
         print 104, z(k),qc(k),rho(k)
104      format(1x,3e12.4,'    z,qc,rho: ')
         enddo

        do k=1,nz
        pl1(k)=qc(k)*1.e3
        enddo

c      call setusv('LW',2000)
c      call set(.1,.37,.10,.50,0.,5.,0.,2400.,1)
c      call labmod('(f5.0)','(f5.0)',5,5,2,2,28,20,0)
c      call perim (5,1,3,4)
c      call dashdc('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$', 3,12)
c      call curved(pl1(k0),z(k0),nz-k0)

cc buo:
      eps=rv/rd-1
      do k=k0,nz
      thve=the(k)*(1.+eps*qve(k))
      thv =th(k)*(1.+eps*qv(k)-qc(k))
      pl1(k)=thv-thve
      enddo
c      call setusv('LW',2000)
c      call set(.38,.65,.10,.50,-2.,2.,0.,2400.,1)
c      call labmod('(f5.0)','(f5.0)',5,5,2,2,28,20,0)
c      call perim (4,1,3,4)
c      call dashdc('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$', 3,12)
c      call curved(pl1(k0),z(k0),nz-k0)
c      call dashdc('$''''$$''''$$''''$$''''$$''''$$''''$', 6,12)
c      call setusv('LW',1000)
      x(1)=0.
      x(2)=0.
      y(1)=0.
      y(2)=2400.
c      call curved(x,y,2)
c      call setusv('LW',2000)
c      call dashdc('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$', 3,12)

cc cumulative cape:
      eps=rv/rd-1
      cape(1)=0.
      do k=k0+1,nz
      thve=the(k)*(1.+eps*qve(k))
      thv =th(k)*(1.+eps*qv(k)-qc(k))
      pl1(k)=max(0., (thv-thve)/thve )
      cape(k)=cape(k-1) + gg*pl1(k)*dz
      enddo
c      call setusv('LW',2000)
c      call set(.66,.93,.10,.50,0.,100.,0.,2400.,1)
c      call labmod('(f5.0)','(f5.0)',5,5,2,2,28,20,0)
c      call perim (2,5,3,4)
c      call dashdc('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$', 3,12)
c      call curved(cape(k0),z(k0),nz-k0)
c
c      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c      CALL plchhq(.02,0.30,'height (km)',0.012,90.,0)
c      CALL plchhq(.06,0.10,'0.0',0.012,0.,0)
c      CALL plchhq(.06,0.50,'2.4',0.012,0.,0)
c
c      CALL plchhq(.23,0.06,'q:B:c:N: (g/kg)', 0.014,0.,0)
c      CALL plchhq(.11,0.08,'0',0.012,0.,0)
c      CALL plchhq(.36,0.08,'5',0.012,0.,0)
c
c      CALL plchhq(.51,0.06,
cc     1 ':PGU:D:PGU:H:PRL::B:d :N: :PRU:(K)', 0.014,0.,0)
c     1 ':PGU:D:PGU:H:PRL::B:d:S::PRU::S:(a) :N: :PRU:(K)', 0.014,0.,0)
c      CALL plchhq(.39,0.08,'-2',0.012,0.,0)
c      CALL plchhq(.64,0.08,' 2',0.012,0.,0)
c
c      CALL plchhq(.79,0.06,'cCAPE (J/kg)', 0.014,0.,0)
c      CALL plchhq(.67,0.08,' 0',0.012,0.,0)
c      CALL plchhq(.92,0.08,'100',0.012,0.,0)

c        call frame



c         call clsgks
          stop
          end


 
