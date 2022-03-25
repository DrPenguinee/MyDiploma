
      PROGRAM eval_KATRIN_KNM1


c KATRIN KNM1 neutrino mass experiment result evaluation.
c Data point and intensity for KNM1 - from September 2019 final file 
c Parameter shift - SpRes shift
c NB! instead of trapping - background slope				 
c Probabilities are dependant on dE (cos(theta))
c Step value vs error bars check by MC-simulation
c step + source tickness + scattered electrons+resolution function modification
c dead time variation
c HEAVY NEUTRINO

      implicit double precision (a-h,o-z)
      external chi2
      common /const/em,alpha,pi

      em=510999.06
      alpha=1./137.0359895
c      pi=3.141592653589793238
	pi = 6.d0*dasin(0.5d0)


      print*,'Read of data from file input_KNM1.dat'
      open (4,file='input_KNM1.dat',status='old')
      open (7,file='fit_kat.dat',status='old')
      open (1,file='eval_kat.log',status='old')
      call mintio(4,1,7)
      call minuit(chi2,0)
      close (4)
      close (7)
      close (1)
      stop
      end



      subroutine chi2(npar,grad,fval,xval,iflag)
c beta-spectrum calculation
      implicit double precision (a-h,o-z)
      dimension grad(25),xval(25)
	dimension epoint(301),Edata(301),error(301),espnm(301)
      dimension espec(301),bacval(301),endef(301),eml(301)
     *,trspec(301),hnu(301),espnm2(301),fakesp(301)
	*,rancor(301),tspec(301)
	data ifl /0/					 
      common /const/em,alpha,pi
      common /source/ thick
      common /resfun/ shift,shift2
      common /lattice/ epoint,npoint
      save Edata,error
      e0=xval(1)+18570.d0
c                   e0 - spectrum endpoint
      s=xval(2)
c                   s - spectrum square
      snm=xval(3)
c                   snm - squared neutrino mass
      bac=xval(4)
c                   bac - background value
      bacpar=xval(5)
c                   bacpar - background parametr
      step=xval(6)
c                   step - value of endpoint effect
      eend=xval(7)+18000.
c                   eend - end point of step effect
      shift=xval(8)
c             shift - probability shift of excitation-ionization losess
      shift2=xval(9)
c             shift2 - main state probability shift
      thick=xval(10)
c                   thick - source thickness factor
      emin=xval(11)
c                   emin - lower limit of analysis interval
      emax=xval(12)
c                   emax - upper limit of analysis interval

      dtime=xval(13)
c             dtime - counting system deadtime
      dtime0=6.e-6
      ntype=int(xval(14)+0.1)
c             ntype=1 - calculation by splines
c             ntype=2 - direct calculation
      snm2=xval(15)**2
c                   snm2 - squared heavy neutrino mass
      hnupr=xval(16)
c                   hnupr - heavy neutrino probability
      unvis=xval(17)
c	              unvis - unvisible correction shift
	 cut=xval(18)
c				  cut - point "glitches" cut
	pos_ml=e0-xval(19)
c				  pos_ml - missing level energy
	prob_ml=xval(20)
c				  prob_ml - missing level probability
      acc=.0002

      if (ifl.eq.0) then
      ifl=ifl+1

	Num_MC=0	!MC simulations counter

c data file read------------------------------------
      
101   format (f8.1,f9.6,f11.8)
1011  format (i3,f10.1,f11.6,f13.8)

	npoint=27

	open (33,file='KNM1.dat',status='old') 
	do i=1,npoint
	 read (33,101) epoint(i),Edata(i),error(i)
	 epoint(i)= 1000.d0*epoint(i)			!keV -> eV transform
	 write (6,1011) i,epoint(i),Edata(i),error(i)
      end do
	close (33)
c	stop

c end of data point set-----------------------------------------

c Arrays calculation
      do i=1,100
      entr=dble(i)
      call transmission (entr,tran)
      call fstail(e0-entr-165.,e0,tail)
      call specint(e0-entr,vesp,e0,snm,vespnm,thick)
c      call trapsp(e0-entr,trap,e0)
      write (6,220) entr,tran,tail
 220  format (4f10.3)
       end do
c Arrays calculated

      goto  50     !Real spectrum analysis

c Faked spectrum formation___________________________________
c      w=1.38*26.       !26 - Katrin to Troitsk Run34 intensity ratio
c	 w=168.0   		!3.3 Hz at 18555V normalization
c	w=295.9   		!5.83 Hz at 18555V KATRIN 10m, 3y normalization	 
	w= 96.4   		!1.9  Hz at 18555V KATRIN - KNM1	  
      open (2,file='run_kat_MC2.dat',status='old')
	write (2,102) 'Simulation with Elow=',Emin
 102	format (a21,f7.1)
	write (2,103) 'HV','Freq','Err','Time'
 103	format (5x,a2,8x,a4,10x,a3,9x,a4) 
      do  i1=1,npoint
      e=epoint(i1)
      call specint(e,espec(i1),e0,snm,espnm(i1),thick)
      call specint(e,hnu(i1),e0,snm2,espnm2(i1),thick)
	call specmlint(e,eml(i1),pos_ml,e0)
c      call trapbacground(e,trspec(i1),e0)
      call endeffect (e,eend,step,endef(i1))
      call background (e,e0,bac,bacpar,bacval(i1))
c Background slope 3mHz per 40eV -switched OFF
	bacval(i1)=bacval(i1)
c	*-.001d0*(e0-e)/40.d0
      fakesp(i1)=w*(espec(i1)+hnu(i1)*hnupr
c	*+trspec(i1)*bacpar
     *+eml(i1)*prob_ml)+bacval(i1)+endef(i1)
	Edata(i1)= fakesp(i1)
c	error(i1)=dsqrt(fakesp(i1)-bacval(i1) + bacval(i1)*1.21)
c	* /tspec(i1)												   !Correction for non-Poissonic baground
c     error(i1)=fakesp(i1)/dsqrt(fakesp(i1)*tspec(i1))
 106	format (f9.1,2f12.4,e12.3)
      write (2,106) e,Edata(i1),error(i1)
c	,w*espec(i1),w*hnu(i1)*hnupr,
c	*w*trspec(i1)*bacpar,w*eml(i1)*prob_ml,bacval(i1),endef(i1)
     
	end do
	close (2)
      end if

c Faked spectrum formed___________________________________


      if (iflag.eq.5) then
      Num_MC=Num_MC+1
      do i1=1,npoint
      Edata(i1)=fakesp(i1)+error(i1)*gauss(-1.d0)
      rancor(i1)=(Edata(i1)-fakesp(i1))/error(i1)
 120  format (1x,8f8.2)
      end do
      write (6,120) (rancor(j),j=1,npoint)
      write (6,*) Num_MC
      end if

c Current random faked spectrum formed____________________

 50   continue

      if (iflag.eq.3) then
       open (8,file='graph_kat.dat',status='old')
       open (9,file='graph_katd.dat',status='old')
      end if

      fval=0.
      sespec=0.
	sml=0.
      strap=0.
      sbacval=0.
      sendef=0.
      do 2 i1=1,npoint
      e=epoint(i1)
      if (iflag.eq.3) goto 12
      if (e.lt.emin .or. e.gt.emax) goto 2
 12   continue
      if (ntype.eq.1) then
      call specint(e,espec(i1),e0,snm,espnm(i1),thick)
      call specint(e,hnu(i1),e0,snm2,espnm2(i1),thick)
      end if
      if (ntype.eq.2) then
      call expspectrum(e,espec(i1),e0,snm,ACC,6)
      call expspectrum(e,hnu(i1),e0,snm2,ACC,6)
      end if
	call specmlint(e,eml(i1),pos_ml,e0)
c      call trapbacground(e,trspec(i1),e0)
      call background (e,e0,bac,bacpar,bacval(i1))
      call endeffect (e,eend,step,endef(i1))
      if (e.lt.(emin-.01) .or. e.gt.emax) goto 2
      sespec=sespec+espec(i1)+hnu(i1)*hnupr
	sml=sml+eml(i1)*prob_ml
c      strap=strap+trspec(i1)*bacpar
      sbacval=sbacval+bacval(i1)
      sendef=sendef+endef(i1)
  2   continue
      w=(s-sbacval-sendef)/(sespec+strap+sml)
      do 4 i1=1,npoint
      e=epoint(i1)
      theo=w*(espec(i1)+hnu(i1)*hnupr
c	*+trspec(i1)*bacpar
     *+eml(i1)*prob_ml)+bacval(i1)+endef(i1)
      exper=Edata(i1)*(1-Edata(i1)*dtime0)/(1-Edata(i1)*dtime)
c      if (e.lt.18540.) exper=exper-unvis
      if (iflag.eq.3)	then
      write (8,104) i1,e-18000.
     *,(exper-w*(espec(i1)-espnm(i1)+hnupr*(hnu(i1)-espnm2(i1)))
     *-bacval(i1))/error(i1),1.
     *,endef(i1)/error(i1)
c     *,w*trspec(i1)*bacpar/error(i1)
     *,w*espnm(i1)/error(i1)
     *,w*espnm2(i1)*hnupr/error(i1)
      write (9,104) i1,e-18000.
     *,exper,w*(espec(i1)-espnm(i1)+hnupr*(hnu(i1)-espnm2(i1)))
     * -bacval(i1)
     * ,endef(i1)
c     * ,w*trspec(i1)*bacpar
     * ,w*espnm(i1)
     * ,w*espnm2(i1)*hnupr
     * ,w*eml(i1)*prob_ml
      end if

      if (e.lt.(emin-.01).or. e.gt.emax) goto 4
      dif=(exper-theo)/error(i1)
      if (dabs(dif) .gt. cut) goto 4
      fval=fval+dif**2
  4   continue
      if (iflag.eq.3) then
      close(8)
      close(9)
      end if

 104  format (i4,8e12.4)
 105  format (f8.2,7e12.4)
      return
      end

      function Gauss(x)
C Random Gauss value generator
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      vls=0.       
       do i1=1,20
       vl=g05caf(-1.d0)
       vls=vls+vl/20.
       end do
      vls=(vls-.5)/.068/0.949
      gauss=vls
      return      
      end

      DOUBLE PRECISION FUNCTION G05CAF(X)
C     * GENERATE PSEUDO-RANDOM NUMBER UNIFORMLY DISTRIBUTED
C                    BETWEEN 0 AND 1 *
C     ALGORITHM AS.183 APLL. STATISTICS (1982) V.31, NO.2, P.188
C     BY B.A.WICHMAN & I.D.HILL (NPL, TEDDINGTON, ENGLAND)
C     IMPLEMENTED ON HP-1000 BY E.L.KOSAREV  23.08.83
C     ( SEE FILE RANDOM FOR SEPARATE VERSION OF PROGRAM )
C     MODIFIED FOR NAF LIBRARY BY KOLPAKOW I.V. 20.04.84
C     X IS DUMMY ARGUMENT TO MATCH FORTRAN SYNTAX
C     ..
C     ..SCALAR ARGUMENTS..
      DOUBLE PRECISION X
C     ..
C     ..SCALARS IN COMMON..
      INTEGER IX,IY,IZ
      COMMON /CAG05A/ IX,IY,IZ
C
      IX = 171*MOD(IX,177) - 2*(IX/177)
      IY = 172*MOD(IY,176) - 2*(IY/176)
      IZ = 170*MOD(IZ,178) - 2*(IZ/178)
C
      IF (IX.LT.0) IX = IX+30269
      IF (IY.LT.0) IY = IY+30307
      IF (IZ.LT.0) IZ = IZ+30323
      G05CAF = DMOD(DFLOAT(IX)/30269.0D0+DFLOAT(IY)/30307.0D0+
     +         DFLOAT(IZ)/30323.0D0, 1.0D0)
      RETURN
      END
C
      BLOCK DATA
C     * INITIALIZATION OF COMMON USED IN G05 SECTION *
C     SET INITIAL VALUE TO CAG05A - G05CAF START STATE
C     SET INITIAL VALUE TO CAG05B - G05DDF AND G05DGF START STATE
C     KOLPAKOW I.V. 20.04.84
C     ..
C     ..SCALARS IN COMMON..
      INTEGER IX,IY,IZ
      DOUBLE PRECISION NORMAL,GAMMA
C     ..
      COMMON /CAG05A/ IX,IY,IZ
      COMMON /CAG05B/ NORMAL,GAMMA
C     ..
      DATA IX/1/,IY/255/,IZ/25555/
      DATA NORMAL/1.0D0/,GAMMA/-1.0D0/
      END


      SUBROUTINE specmlint(e,eml,e0,emax)
c Missimg level spectrum calculation (m=0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W1(301),W2(301),W3(301)
      DIMENSION DE(301),epoint(301)
      DIMENSION SLEV1(301),CM1(301)
      COMMON /CONST/EM,ALPHA,PI
c      common /lattice/ epoint,npoint
      common /source/ thick
      SAVE IFL,DE,SLEV1,CM1
      data IFL /0/

      if (ifl.eq.0) then	  !Missing level spline calculation
	ifl=ifl+1
      print*,'Missing level spectrum calculation'
C E-POINTS NODES LIST FORMATION------------------------------------
	do i=1,301
	de(i)=dble(i)*.2-2.d0
	enddo


c END OF E-POINT FORMATION___________________

c Spline calculation
       ac=.001

      DO IP=1,301
       Elow=Emax-de(IP)
       IF (Elow.GE.Emax) THEN
        SLEV1(IP)=0
       ELSE
	  snm=0.		!Regular spectrum for normalization
        call specint (Elow,ESP,Emax,SNM,espnm,thick)
	  accur=esp*ac	!absolute required accuracy calculated
c	write (6,*) ip,de (ip),esp,accur
        CALL expmlspectrum(Elow,SLEV1(IP),Emax,accur)
       END IF
c	write (6,*) ip,de (ip),slev1(ip),accur
	end do

      CALL SPLINE (301,DE,SLEV1,CM1,W1,W2,W3)

      print*,'Spline is calculated for missing level'
      open (unit=22,file='spn_MCml.dat',status='old')
      do j=1,301
      write (22,*) de(j),slev1(j),cm1(j)
      end do
      close (22)
C END OF SPLINES CALCULATION _____________________________

      end if	 !End of initialization

	if (e0-e.lt.0.) then
	 eml=0.
	else
       CALL SPLINT (301,E0-E,eml,DE,SLEV1,CM1)
	end if

	return
	end

      SUBROUTINE 	expmlspectrum(E,Specml,E0,ac)
c experimental value of missing level spectrum 
c     subroutine expspectrum(e,espec,e0,snm,ac,ntyp)
      implicit double precision (a-h,o-z)
      dimension w(2000),iw(260)
      external convol_ml
      common /energy_ml/ pe,pe0
      pe=e
      pe0=e0
      ifail=0
      ERabs=1.D-4
      call d01ajf(convol_ml,e,e0,ERabs,ac,specml,er,w,
     *2000,iw,260,ifail)
      return
      end


      function convol_ml(x)
      implicit double precision (a-h,o-z)
      common /energy_ml/ pe,pe0
      call truspectrum_ml(x,tspec,pe0)
      call transmission (x-pe,tran)
      convol_ml=tspec*tran
      return
      end

      subroutine truspectrum_ml(e,tspec,e0)
      implicit double precision (a-h,o-z)
      common /const/em,alpha,pi

      pe=sqrt(e*(e+2.*em))

c --- snm=0 formulae
      de=e0-e
      de2=de**2
      if (de.le.0.) then
	 tspec=0.
	else
       tspec=de2
      end if
c -------------------------------------------

      tspec=fermi(e)*(e+em)/em*pe*tspec/4.e9
      return
      end





      SUBROUTINE specint (E,ESP,E0,SNM,espnm,thick)
C LATTICE APPROXIMATION CALCULATION FOR THE SPECTRA
C 193 or 155 FINAL STATES
C NEUTRINO MASS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ELEV1(193),PROLEV1(193)
      DIMENSION W1(151),W2(151),W3(151)
      DIMENSION DE(151),Enode(151),epoint(301)
      DIMENSION SLEV1(151),CM1(151)
      DIMENSION VSNM(110),VNEUT1(151,110)
     *,cm201(151,110),cm221(151,110),cm021(151,110)
      character name(24)*12
      COMMON /LEVDAT/ ELEV(193),PROLEV(193),ist
      COMMON /CONST/EM,ALPHA,PI
      common /resfun/ shift,shift2																				  77yg
      common /lattice/ epoint,npoint
      SAVE IFL,DE,thickold,shiftold,shift2old,e0old
     *,VSNM,VNEUT1,cm201,cm221,cm021
     *,elev1,prolev1,enode,SLEV1,CM1

      DATA IFL /0/
      data e0old /0./
      data thickold /0./
      data shiftold /100./
      data shift2old /100./
      data name /'spn_KNM1.dat',
     *'spn_sX00.dat','spn_sX01.dat',
     *'spn_sX02.dat','spn_sX03.dat',
	*'spn_sX05.dat','spn_sX10.dat','spn_sX15.dat',
c     *'spn_dx00.dat','spn_dx01.dat',
c     *'spn_dx02.dat','spn_dx03.dat',
c     *'spn_5m05.dat','spn_5p02.dat','spn_5p10.dat',
     *'spn_5m02.dat','spn_5p05.dat','spn_t115.dat',
     *'spn_t120.dat','spn_t125.dat','spn_t130.dat',
     *'spn_2310.dat','spn_2311.dat','spn_2312.dat',
     *'spn_2313.dat','spn_2314.dat','spn_2315.dat',
     *'spn_2316.dat','spnov017.dat','spnov018.dat',
     *'spnov019.dat'/

      if (ifl.eq.0) then
C Read exitation level data---------------------
      OPEN (8,FILE='EXCITAT2.DAT',STATUS='OLD')
      i1=1
 200   READ (8,*,end=201) ELEV1(I1),PROLEV1(I1)
      PROLEV1(I1)=PROLEV1(I1)/100.
      PROLEV(I1)=PROLEV1(I1)
      eLEV(I1)=eLEV1(I1)
      i1=i1+1
      goto 200
c      END DO
 193  format (i6,f12.3,2f12.5)
 201  ist=i1-1
      CLOSE (8)
C-----------------------------------------------
      end if

      IF (dabs(shift2-shift2old).ge..000001) then
          print*,'Final states spectrum with shift=',shift2
      slev=0.
      sl1=.57408
      sl2=.42077
      open (22,file='FS_out.dat',status='old')
        DO I1=1,ist
      if (elev1(i1).le.5.)
     *    prolev(i1)=prolev1(i1)*(sl1-shift2/100.)/sl1
      if (elev1(i1).gt.5.)
     *    prolev(i1)=prolev1(i1)*(sl2+shift2/100.)/sl2
      slev=slev+prolev(i1)
          WRITE (22,193) i1,ELEV(I1),PROLEV(I1),slev
        END DO
	close (22)
      end if


      IFL=IFL+1
      IF (dabs(thick-thickold).ge..000001.or.
     *    dabs(shift2-shift2old).ge..000001.or.
     *    dabs(e0-e0old).gt.15.35.or.
     *    dabs(shift-shiftold).ge..000001) then
      e0old=e0
      thickold=thick
      shiftold=shift
      shift2old=shift2
c calculation of new splines------------------------
      print*,'Needs a spline'
      do ithick=1,1
c      ithick=1
      open (unit=2,file=name(ithick),status='old')
      read (2,*) e0file,thickfile,shiftfile,shift2file
      IF (dabs(thick-thickfile).le..0001.and.
	*    dabs(e0-e0file).le.15.25.and.
     *    dabs(shift2-shift2file).le..0001.and.
     *    dabs(shift-shiftfile).le..0001) then
      print*,'Reading file=',name(ithick),'with parameters:'
      print*,'Spectrum endpoit=',e0file
      print*,'Thickness factor=',thickfile
      print*,'Ex/ion shift, %=',shiftfile
      print*,'First state shift, %=',shift2file
      do j=1,151
      read (2,*) de(j),slev1(j),cm1(j)
      end do
      read (2,220) (vsnm(j),j=1,110)
      read (2,220) ((vneut1(k,j),k=1,151),j=1,110)
 220  format (3e24.16)
      close (2)
      CALL SPLINE (151,DE,SLEV1,CM1,W1,W2,W3)
      CALL spline2(151,110,VNEUT1,DE,VSNM,cm201,cm221,cm021)
      goto 33
      else
      close (2)
      end if
      end do
      print*,'New spline should be calculated'

      print*,e0
c	DO i=1,27
c	 WRITE (25,*) i,epoint(i)
c	end do
c	stop
C E-POINTS NODES LIST FORMATION------------------------------------
      DO i=1,9
       enode(i)=19621.d0-100.d0*dble(i)
	 WRITE (25,*) i,enode(i)
	END  DO
      DO i=11,37
       enode(i)= epoint(38-i)+0.5d0
	 WRITE (25,*) i,enode(i)
	END  DO
      DO i=38,48
       enode(i)=epoint(1) - 2.d0*(i-37)
	 WRITE (25,*) i,enode(i)
	END  DO
	DO i=49,150
	 enode(i)= enode(48)-100.d0*dble(i-48) 
	 WRITE (25,*) i,enode(i) 
	END DO
 
c      stop

c END OF E-POINT FORMATION___________________

c SNM POINT FORMATION________________________
      VSNM(36)=0.
      DO I=1,20
       VSNM(36+I)=.05*DBLE(I)
       VSNM(36-I)=-VSNM(36+I)
      END DO
      DO I=21,35
       VSNM(36+I)=.2*DBLE(I-20)+1.
       VSNM(36-I)=-VSNM(36+I)
      END DO
      do i=72,110
      vsnm(i)=4.+4.*dble(i-71)**2
      end do
c      do i=1,110
c	write (22,*) i,vsnm(i)
c	end do
c	stop
C END OF SNM POINT FORMATION

       ac=.0001

       DO IP=1,151
        de(ip)=e0-enode(IP)
        Elow=enode(IP)
        VNEUT1(IP,36)=0.
        IF (Elow.GE.E0-ELEV1(1)) THEN
        SLEV1(IP)=0
        ELSE
        snmt=0.
        CALL expspectrum(Elow,SLEV1(IP),E0,SNMt,ac,3)
        END IF
        DO IM=1,35
         snmt=VSNM(36+IM)
         IF (Elow.GE.E0-ELEV1(1)) THEN
         snmspec=0.
         ELSE
         CALL expspectrum(Elow,snmspec,E0,snmt,3.*ac*SLEV1(IP),7)
c         CALL expspectrum(Elow,snmspec,E0,snmt,ac*3.,7)
         END IF
         VNEUT1(IP,36+IM)=snmspec
         VNEUT1(IP,36-IM)=-snmspec
         write (6,*) ip, im,snmspec
        END DO
        DO IM=72,110
         snmt=VSNM(IM)
         IF (Elow.GE.E0-ELEV1(1)) THEN
         snmspec=0.
         ELSE
         CALL expspectrum(Elow,snmspec,E0,snmt,3.*ac*SLEV1(IP),7)
c         CALL expspectrum(Elow,snmspec,E0,snmt,ac*5.,7)
         END IF
        VNEUT1(IP,IM)=snmspec
        write (6,*) ip, im,snmspec
        END DO
       END DO

       CALL SPLINE (151,DE,SLEV1,CM1,W1,W2,W3)
       CALL spline2(151,110,VNEUT1,DE,VSNM,cm201,cm221,cm021)

      print*,'New splines are calculated with:'
      print*,'Endpoint energy=',e0
      print*,'Thickness factor=',thick
      print*,'Ex/ion shift, %=',shift
      print*,'First state shift, %=',shift2
      open (unit=22,file='spn_last.dat',status='old')
      write (22,*) e0,thick,shift,shift2
      do j=1,151
      write (22,*) de(j),slev1(j),cm1(j)
      end do
      write (22,220) (vsnm(j),j=1,110)
      write (22,220) ((vneut1(k,j),k=1,151),j=1,110)
      close (22)

C END OF SPLINES CALCULATION _____________________________
 33   continue
      end if

      CALL splint2 (151,110,VNEUT1,DE,VSNM,cm201,
     *cm221,cm021,E0-E,SNM,ESPnm,dx,dy,0,ifail)

      CALL SPLINT (151,E0-E,esp,DE,SLEV1,CM1)

      ESP=ESP+ESPnm

 100  format (2x,f10.2,4f12.6)

      RETURN
      END


      subroutine endeffect (e,eend,step,endef)
      implicit double precision (a-h,o-z)
      energy=eend-e
      call transmission (energy,tran)
      endef=tran*step
c      endef=0
c      if (energy.ge.0) endef=step
      return
      end

      subroutine background (e,e0,bac,bacpar,bacval)
      implicit double precision (a-h,o-z)
      bacval=bac+bacpar*(18575.d0-e)/40.d0 
      return
      end

      subroutine expspectrum(e,espec,e0,snm,ac,ntyp)
      implicit double precision (a-h,o-z)
      dimension w(2000),iw(260)
      external convol
      common /energy/ pe,pe0,psnm,ntype
      pe=e
      pe0=e0
      psnm=snm
      ntype=ntyp
      ifail=0
      if (ntyp.lt.7) then
      ERABS=1.D-10
      call d01ajf(convol,e,e0,ERABS,ac,espec,er_t,w,
     *2000,iw,260,ifail)
      else
      ERrel=1.D-10
      call d01ajf(convol,e,e0,ac,ERrel,espec,er_t,w,
     *2000,iw,260,ifail)
      end if
      return
      end


      function convol(x)
      implicit double precision (a-h,o-z)
      common /energy/ pe,pe0,psnm,ntype
      call truspectrum(x,tspec,pe0,psnm,ntype)
      call transmission (x-pe,tran)
      convol=tspec*tran
      return
      end


      subroutine trapbacground(e,trbac,e0)
      implicit double precision (a-h,o-z)
      dimension w(2000),iw(260)
      DIMENSION W1(98),W2(98),W3(98)
      DIMENSION DE(98),enode(98),epoint(301)
      DIMENSION SLEV1(98),CM1(98)
      external conlin
      common /conl/ pe,pe0
      common /lattice/ epoint,npoint
      SAVE IFL,DE,SLEV1,CM1
      DATA IFL /0/

      IFL=IFL+1
      IF (ifl.eq.1) then
c calculation of trapped electrons splines------------------------
      print*,'Start of trapped electrons background calculation'

C E-POINTS NODES LIST FORMATION------------------------------------

      DO i=1,npoint
      enode(i)=epoint(npoint+1-i)+0.25
	WRITE (6,*) i,enode(i) 
      END  DO
      do i=npoint+1,98
      enode(i)=epoint(1)-5.*dble(i-npoint)+0.25
	WRITE (6,*) i,enode(i) 
      end do

c END OF E-POINT FORMATION___________________

      ac=.0001

       DO IP=1,98
        DE(IP)=E0-enode(IP)
        Elow=enode(IP)
        IF (Elow.GE.E0) THEN
        SLEV1(IP)=0
        ELSE
      ac=.0001
      pe=elow
      pe0=e0
      ifail=0
      ERABS=1.D-10
      call d01ajf(conlin,elow,e0,ERABS,ac,SLEV1(IP),erest,w,
     *2000,iw,260,ifail)
      end if
       END DO

       CALL SPLINE (98,DE,SLEV1,CM1,W1,W2,W3)
c test:
      open (unit=22,file='trbac.dat',status='old')
      do i=7,98
      el=enode(i)
       CALL SPLINT (98,E0-El,trbac,DE,SLEV1,CM1)
        write (22,*) el,trbac
c        write (6,*) el,trbac
      end do
      close (22)
c     stop
c end of test

C END OF SPLINES CALCULATION _____________________________
      print*,'Trapped electrons background calculated'
      end if

      CALL SPLINT (98,E0-E,trbac,DE,SLEV1,CM1)
      return
      end

      function conlin(x)
      implicit double precision (a-h,o-z)
      common /conl/ pe,pe0
      call trapsp(x,trap,pe0)
      call transmission (x-pe,tran)
      conlin=trap*tran
      return
      end

      subroutine trapsp(e,trap,e0)
      implicit double precision (a-h,o-z)
      dimension w(2000),iw(260)
      DIMENSION W1(166),W2(166),W3(166)
      DIMENSION et(166),tt(166),gt(166)
      external ftrap
      common /trapco/ e1,e2
      common /temp/ tspec,scprob
      save ifl,et,tt,gt
      DATA IFL /0/

c trapspectrum spline calculation
      IFL=IFL+1
      IF (IFL.EQ.1) THEN
      print*,'Differential trapping spectrum calculation'
      do i=1,166
      e1=18581.-1.5*dble(i-1)
      et(i)=e0-e1
      if (e1.ge.e0-10.) then
       tt(i)=0.
      else
       ac=0.0001
       e2=e0
       ifail=0
       ERABS=1.D-10
       call d01ajf(ftrap,e1,e2,ERABS,ac,vtr,er,w,
     * 2000,iw,260,ifail)
       tt(i)=vtr
      end if
c      write (6,*) e1,e2,tt(i)
      end do


      CALL SPLINE (166,et,tt,gt,W1,W2,W3)
      print*,'Differential trapping spectrum calculated'
c test
      open (unit=22,file='trutrap.dat',status='old')
c      print*,'Trapped spectra'
      do i=1,166
      write (22,*) et(i),tt(i),gt(i)
c      write (6,*) et(i),tt(i),gt(i)
 100  format (f12.2,f12.6)
      end do
      close (22)
c end of spline calculation______________________________________
      end if

      if (e.ge.e0) then
       trap=0.
      else
       call splint (166,e0-e,trap,Et,Tt,GT)
      end if
      return
      end



      function ftrap(x)
      implicit double precision (a-h,o-z)
      common /trapco/ e1,e2
      common /temp/ tspec,scprob
      snmt=0.
      nt=6
      call truspectrum(x,tspec,e2,snmt,nt)
      call scatprob(e1,x,e2,scprob)
      ftrap=tspec*scprob
c      tspec=0.
c      scprob=0.
c      ftrap=1.
      return
      end

      subroutine scatprob(e,ep,e0,prob)			!Trapped electrons spectrum
      implicit double precision (a-h,o-z)
      DIMENSION En(221),S1(221),G1(221)			!Scattering spectrum spline
      DIMENSION W1(221),W2(221),W3(221)			!Spline working arrays
	Save en,s1,g1
      DATA IFL /0/
c Read of precalculated scattering probability spline 
      IFL=IFL+1
      IF (IFL.EQ.1) THEN
	 open (22,file='Nscat.dat',status='old')
	 do i=1,221
	 read(22,*) en(i),s1(i)
	 end do
	 call spline (221,en,s1,g1,w1,w2,w3)
c Spectrum restored
      end if
	if ((ep-e).lt.0.) then
	 prob=0.
	else
	 call splint(221,ep-e,prob,en,s1,g1)       !data interpolation by spline
      end if
	return
      end


      subroutine truspectrum(e,tspec,e0,snm,ntype)

      implicit double precision (a-h,o-z)
      common /const/em,alpha,pi
      COMMON /LEVDAT/ ELEV(193),PROLEV(193),ist

      pe=sqrt(e*(e+2.*em))
      tspec=0

      do 2 i1=1,ist
         if (ntype.eq.1) then
c ----mainz approximation
      de=e0-elev(i1)-e
      de2=de**2
      if (snm.ge.-1.d-8) then
      if (de.lt.0.0 .or. de2.lt.snm) goto 2
      fuck=de
      end if
      if (snm.lt.-1.d-8) then
      vmu=.76*dsqrt(-snm)
      if (de.lt.(-vmu)) goto 2
      fuck=de+vmu*dexp(-de/vmu-1.)
      end if
      tspec=tspec+prolev(i1)*fuck*dsqrt(de2-snm)
      end if
c -------------------------------------------
      if (ntype.eq.2) then
c --- stupid formulae
      de=e0-elev(i1)-e
      de2=de**2
      if (de.le.0.) goto 2
       tspec=tspec+prolev(i1)*(de2-snm/2.)
      end if
c -------------------------------------------
      if (ntype.eq.3) then
c --- snm=0 formulae
      de=e0-elev(i1)-e
      de2=de**2
      if (de.le.0.) goto 2
       tspec=tspec+prolev(i1)*de2
      end if
c -------------------------------------------


      if (ntype.eq.6) then
c --- basic symmetric formulae
      de=e0-elev(i1)-e
      de2=de**2
      if (de.le.0.) goto 2
      if (de2.gt.dabs(snm)) then
       dtspec=de*dsqrt(de2-dabs(snm))-de2
       else
       dtspec=-de2
      end if
      if (snm.ge..0) then
       tspec=tspec+prolev(i1)*(de2+dtspec)
       else
       tspec=tspec+prolev(i1)*(de2-dtspec)
       end if
      end if
c -------------------------------------------
      if (ntype.eq.7) then
c neutrino mass term from Troitsk approximation
      de=e0-elev(i1)-e
      de2=de**2
      if (de.le.0.) goto 2
      if (de2.gt.dabs(snm)) then
       dtspec=de*dsqrt(de2-dabs(snm))-de2
       else
       dtspec=-de2
      end if
      if (snm.ge.0.) then
       tspec=tspec+prolev(i1)*(+dtspec)
       else
       tspec=tspec+prolev(i1)*(-dtspec)
       end if
      end if
C -------------------------------------------
  2   continue
      if (ist.ne.193.or.ntype.eq.7) goto 3
      if ((e0-e).gt.165.) then
      call fstail(e,e0,tail)
      tspec=tspec+tail
      end if
  3    tspec=fermi(e)*(e+em)/em*pe*tspec/4.e9

      return
      end

      subroutine fstail(e,e0,tail)
C final states tail effect
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(2000),IW(260)
      DIMENSION Ep(200),S(200),GM1(200)
      DIMENSION w1(200),w2(200),w3(200)
      EXTERNAL ftail
      common /tail/Emax
      save ifl,Ep,S,GM1
      data ifl /0/

      if (ifl.eq.0) then

      ERABS=1.D-10
      ERREL=1.D-4
      ep(1)=165.
      s(1)=0.
      open (unit=22,file='fstail.dat',status='old')
      DO I=2,200
      Ep(I)=165.+5*dble(i-1)
      IFAIL=0
      EMIN=ep(1)
      EMAX=ep(i)
      CALL D01AJF(ftail,EMIN,EMAX,ERABS,ERREL,St,ER,W,
     *    2000,IW,260,IFAIL)
      s(i)=st
      write (22,*) ep(i),s(i)
      end do
      close (22)
      CALL SPLINE (200,Ep,S,GM1,w1,w2,w3)
      print*,'Final states tail spline calculated'
      end if
      ifl=ifl+1
      call splint (200,e0-e,tail,ep,s,gm1)

      return
      END

      FUNCTION ftail(E)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c Final states tail probability distribution from S.Jonsell,H.J.Monkhorst
c Phys.Rev.Let. 76,(1996),4476
      common /tail/Emax
      Ry=13.606
      pi=3.1415926
      ep=e-45.
      ftail=.511/0.137*64./(ep/ry+1.)**4
     *    /(1.-dexp(-4.*pi*dsqrt(ry/ep)))
     *    *dexp(-8.*dsqrt(ry/ep)*datan(dsqrt(ep/ry)))/ry
c	 end of tail formula, next is electron spectrum
     *    *(emax-e)**2
      RETURN
      END



      subroutine transmission (energy,tran)
      implicit double precision (a-h,o-z)
      dimension etrans(221),trans(221),gtrans(221)
      common /source/ thick
      common /resfun/ shift,shift2
      save thickold,shiftold,shift2old,etrans,trans,gtrans
      data thickold /0./
      data shiftold /100./
      if (dabs(thick-thickold).ge..000001.or.
     *dabs(shift-shiftold).ge..000001) then
      thickold=thick
      shiftold=shift
c calculation of transmission data------------------------
      print*,'New resolution function should be calculated'
      call RFCAL_KATRIN(thick,etrans,trans,gtrans)
c-----------------------------------------------
      end if
      if (energy.le.0.) then
      tran=0.
      goto 1
      end if
      if (energy.ge.300.) then
      tran=1.
      goto 1
      end if
      call splint (221,energy,tran,etrans,trans,gtrans)
  1   continue
      return
      end


      SUBROUTINE RFCAL_KATRIN(FAC,E,RF,GMR)
C Response function calculation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(4000),IW(1040)
c      DIMENSION S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION E(221),RF(221),GMR(221)
      DIMENSION ER(221)
      DIMENSION WR1(221),WR2(221),WR3(221)
      DIMENSION sProb(0:5)
      COMMON/DESPL/E0,sxt
	common /theta_max/ theta_max

      EXTERNAL clr
c Scattering probabilities calculation from X-value
      sx=1.d0			!Input value 3.64e(-18)cm2 times 5.0e(+17) T2/cm2
					!Information about thickness is shifted to the input file 
	Bs=2.52d0 
	Bpinch=4.24d0		!Magnetic fields for KNM1
      SXT=SX*FAC		!Variable correction

	cos_theta_max =	dsqrt(1.d0 - Bs/Bpinch)
	write (6,*)'cos_theta_max =	',cos_theta_max
c	stop

 170  format (7f8.4)
 	open (22,file='temp.dat',status='old')
	do i=1,100
	theta_tmp=1.d0-.37d0*G05CAF(1.d0)
	call ProbCalc_int(sxt,sProb,theta_tmp)
	write (22,170) theta_tmp,sProb
	end do
	close (22)
C scattering probalities calculated and checked________________________
	
      call LOSPAT(er,0.d0,s1,s2,s3,s4,s5)			  !Los patterns calculation

      DO I=1,221
      e(i)=er(i)
      END DO


      ERABS=1.d-5
      RF(1)=0.
      DO I=2,221
      IFAIL=0
      E1=0.
      E0=ER(I)
      ERREL=1.E-3
      CALL D01AJF(clr,E1,E0,ERABS,ERREL,VAL,ERpr,W,
     *    4000,IW,1040,IFAIL)
c     WRITE (6,*) I,Er(I),VAL
      Spres_point=SPRES(Er(I))
	call ProbCalc_int(sxt,sProb,theta_max)
        rf(i)=SPRES_point*SProb(0)+VAL
c      WRITE (66,120) I,Er(I),SPRES_point*SProb(0),VAL
120 	format (i6,f8.2,2e12.4)
      END DO

c      DO I=1,221
c      RF(I)=RF(I)/RF(221)
c      END DO
      CALL SPLINE (221,ER,RF,GMR,WR1,WR2,WR3)
c test_______________________________
      open (unit=9,file='rf.dat',status='old')
      do i=1,221
c      WRITE (6,*) er(i),rf(i),gmr(i)
      WRITE (9,*) er(i),rf(i),gmr(i)
      end do
      close (9)
c______________________________________
      RETURN
      END

      SUBROUTINE ProbCalc_int(x,sp,theta_max)
C scattering probabilities for source thickness x and theta_max
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION sp(0:5),P0(40),P1(40),P2(40),P3(40),P4(40),P5(40)
	DIMENSION theta(40),G0(40),G1(40),G2(40),G3(40),G4(40),G5(40)
	DIMENSION W1(40),W2(40),W3(40)
	save P0,P1,P2,P3,P4,P5,theta,G0,G1,G2,G3,G4,G5
      data ifl /0/

	if (ifl.eq.0) then		
	ifl=ifl+1

	P0(40)=0.d0					!Counters initialisation
	P1(40)=0.d0
	P2(40)=0.d0
	P3(40)=0.d0
	P4(40)=0.d0
	P5(40)=0.d0
	theta(40)=1.d0

	do i=39,1,-1					!Probabilities precalculation 
	 cos_th=1.-dble(40-i)/100.d0+.005d0
	 theta(i)=1.-dble(40-i)/100.d0
	 sx=x/cos_th			        !spirality correction
       call ProbCalc(sx,0,spr)		!Probility summation
	 P0(i)=P0(i+1)+spr		 
       call ProbCalc(sx,1,spr)
	 P1(i)=P1(i+1)+spr		 
       call ProbCalc(sx,2,spr)
	 P2(i)=P2(i+1)+spr		 
       call ProbCalc(sx,3,spr)
	 P3(i)=P3(i+1)+spr		 
       call ProbCalc(sx,4,spr)
	 P4(i)=P4(i+1)+spr		 
       call ProbCalc(sx,5,spr)
	 P5(i)=P5(i+1)+spr		 
	end do

	theta(40)=1.d0				
 	sx=x/theta(40)			        
      call ProbCalc(sx,0,spr)		!COS(theta)=1 data point
	  P0(40)=spr		 
      call ProbCalc(sx,1,spr)
	  P1(40)=spr		 
      call ProbCalc(sx,2,spr)
	  P2(40)=spr		 
      call ProbCalc(sx,3,spr)
	  P3(40)=spr		 
      call ProbCalc(sx,4,spr)
	  P4(40)=spr		 
      call ProbCalc(sx,5,spr)
	  P5(40)=spr		 

	do i=1,39
	 P0(i)=P0(i)/dble(40-i)		!Probability normalization
	 P1(i)=P1(i)/dble(40-i)
 	 P2(i)=P2(i)/dble(40-i)
	 P3(i)=P3(i)/dble(40-i)
	 P4(i)=P4(i)/dble(40-i)
	 P5(i)=P5(i)/dble(40-i)
	end do

	call spline(40,theta,P0,G0,w1,w2,w3)
	call spline(40,theta,P1,G1,w1,w2,w3)
	call spline(40,theta,P2,G2,w1,w2,w3)
	call spline(40,theta,P3,G3,w1,w2,w3)
	call spline(40,theta,P4,G4,w1,w2,w3)
	call spline(40,theta,P5,G5,w1,w2,w3)

	open (2,file='scat_prob.dat',status='old')
	write (2,*) 'KATRIN scattering probabilities'
	write (2,*) 'Cos_th_max','   P0  ','   P1  ','   P2  ','   P3  ',
	*'   P4  ','   P5  '

	do i=1,40
	write (2,100) theta(i),P0(i),P1(i),P2(i),P3(i),P4(i),P5(i)
 100	format (f10.3,6f7.4)
      end do
	close (2)
	print*,'Scattering probabilities are calculated'
c Probabilities splines are calculated________________________
	end if

	call splint(40,theta_max,spr,theta,P0,G0)
	sp(0)=spr
	call splint(40,theta_max,spr,theta,P1,G1)
	sp(1)=spr
	call splint(40,theta_max,spr,theta,P2,G2)
	sp(2)=spr
	call splint(40,theta_max,spr,theta,P3,G3)
	sp(3)=spr
	call splint(40,theta_max,spr,theta,P4,G4)
	sp(4)=spr
	call splint(40,theta_max,spr,theta,P5,G5)
	sp(5)=spr

	return
	end


      SUBROUTINE ProbCalc(x,n,sp)
C n-scattering probability for source thickness x calculation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION spr(0:5),pr(0:5)
	pr(0)=dexp(-x)          !Poisson probability of 0 scatterings
	spr(0)=(1.-pr(0))/x		!Mean over x probability of 0 scatterings
	do i=1,n
	pr(i)=pr(i-1)*x/dble(i)			!Poisson probability of i scatterings
	spr(i)=spr(i-1)-pr(i-1)/dble(i)	!Mean over x probability of i scatterings
	end do
	sp=spr(n)
	return
	end

      FUNCTION clr(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Ep(221),sProb(0:5)
	common /theta_max/ theta_max
      COMMON/DESPL/E0,sxt

      Spres_point=SPRES(E0-X)
	call ProbCalc_int(sxt,sProb,theta_max)

      CALL LOSPAT(ep,X,s1,s2,s3,s4,s5)
	POG=sProb(1)*s1+sProb(2)*s2+sProb(3)*s3+sProb(4)*s4+sProb(5)*s5

      clr=POG*SPRES(E0-X)
      RETURN
      END

      FUNCTION SPRES(E)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	common /theta_max/ theta_max
      common /resfun/ shift,shift2
      dE=0.93d0		!Spectrometer resolution
	Bs=6.d0			!Source magmetic field
	Bm=10.d0		!Pinch magnetic field
	dE=dE*(1.d0+shift/100.d0)
	IF (E.LE.0.) then 
	SPRES=0.d0
	theta_max=1.d0
	end if
      IF (E.GT.0..AND.E.LE.dE) then
	SPRES=(1.d0-dsqrt(1.d0-E/dE*Bs/Bm))/(1.d0-dsqrt(1.d0-Bs/Bm))
	theta_max=dsqrt(1.d0-E/dE*Bs/Bm)
	end if
      IF (E.GT.dE) then 
	SPRES=1.d0
	theta_max=dsqrt(1.d0-Bs/Bm)
	end if
      RETURN
      END


      subroutine LOSPAT(ep,e_point,s1p,s2p,s3p,s4p,s5p)
C CALCULATION OF 2-,3-times SCATTERED ELECTRON ENERGY SPECTRA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(4000),IW(1040)
c     DIMENSION W(2000),IW(260)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      DIMENSION w1(221),w2(221),w3(221)
      DIMENSION ep(221)
c	*,sp1(221),sp2(221),sp3(221),sp4(221),sp5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
c      common /resfun/ shift,shift2
      EXTERNAL C11,C12,C13,C14,SPf2,SPf3,SPf4,SPf5
	shift=0.d0		!shift parameter is used at SPRES 
	shift2=0.d0		!shift2 is not in use in this routine
      data ifl /0/
	if (ifl.eq.0) then		!Scattered electron spectra calculated
	ifl=ifl+1

      DO I=1,21
       E(I)=DBLE(I-1)/20.	 !0-1.0V
      END DO

      DO I=22,191
        E(I)=.2d0*DBLE(I-21)+1.d0 !2.0-35.0V
      END DO

      DO I=192,212
       E(I)=3.*DBLE(I-191)+35.d0  !35.-98.V
      END DO
        
	do i=213,221
	 e(i)=20.*dble(i-212)+100.0d0 !100-300V
	end do


c____________________________________________

C Single inelastic scattering spectrum
C Final version from V.N.Aseev e.a.,Eur.Phys,J.D10,39-52(2000)

      xc1=12.6
      xc2=14.30
      w_1=1.85
      w_2=12.5
      a1=.204
      a2=.0556


      do i=1,221
	 de=e(i)
       if (de.le.14.3) then
        s1(i)=a1*dexp(-2.*(de-xc1)**2/w_1**2)
       else
        s1(i)=a2*w_2**2/(w_2**2+4.*(de-xc2)**2)
       end if
c      if (e(i).le.14.3) s1(i)=s1(i)*(.51+shift/100.)/.51
c      if (e(i).gt.14.3) s1(i)=s1(i)*(.49-shift/100.)/.49
c          write (6,*) e(i),s1(i)
      end do

c  TEST
c	do i=1,200
c	de=dble(i)/10.
c	sp1=a1*dexp(-2.*(de-xc1)**2/w1**2) 
c	sp2=a2*w2**2/(w2**2+4.*(de-xc2)**2)
c	write (22,100) de,sp1,sp2
c	end do
c	stop
c100	format (3e12.4)
      call spline(221,e,s1,gm1,w1,w2,w3)
     

      OPEN (UNIT=22,FILE='losin_mo.DAT',STATUS='old')
      DO I=1,221
       WRITE (22,*) E(I),S1(I),GM1(I)
      END DO
      CLOSE (22)
 
      print*,'Single scattering spectrum by Ostroumov (Aseev)'
      print*,'is taken and NOT modified with ex/ion shift=',shift
	 
      IFAIL=0
c 2-scatterings											
      S2(1)=0.
      DO I=2,221
      Emax=E(I)
      ERABS=1.E-10
      ERREL=1.E-4
      EMIN=0.
      CALL D01AJF(C11,EMIN,Emax,ERABS,ERREL,S2(I),ER,W,
     *    4000,IW,1040,IFAIL)
      END DO

      CALL SPLINE (221,E,S2,GM2,w1,w2,w3)

      EMIN=0.
      EMAX=300.
      CALL D01AJF(SPf2,EMIN,EMAX,ERABS,ERREL,S2N,ER,W,
     *    4000,IW,1040,IFAIL)

      DO I=1,221
c      S2(I)=S2(I)/S2N
      END DO

      CALL SPLINE (221,E,S2,GM2,w1,w2,w3)

c 3-scatterings
      S3(1)=0.
      DO I=2,221
      Emax=E(I)
      CALL D01AJF(C12,EMIN,Emax,ERABS,ERREL,S3(I),ER,W,
     *    4000,IW,1040,IFAIL)
      END DO

      CALL SPLINE (221,E,S3,GM3,w1,w2,w3)

      CALL D01AJF(SPf3,EMIN,EMAX,ERABS,ERREL,S3N,ER,W,
     *    4000,IW,1040,IFAIL)

      DO I=1,221
c      S3(I)=S3(I)/S3N
      END DO

      CALL SPLINE (221,E,S3,GM3,w1,w2,w3)


c 4-scatterings
      S4(1)=0.
      DO I=2,221
      Emax=E(I)
      CALL D01AJF(C13,EMIN,Emax,ERABS,ERREL,S4(I),ER,W,
     *    4000,IW,1040,IFAIL)
      END DO

      CALL SPLINE (221,E,S4,GM4,w1,w2,w3)

      CALL D01AJF(SPf4,EMIN,EMAX,ERABS,ERREL,S4N,ER,W,
     *    4000,IW,1040,IFAIL)

      DO I=1,221
c      S4(I)=S4(I)/S4N
      END DO

      CALL SPLINE (221,E,S4,GM4,w1,w2,w3)


c 5-scatterings
      S5(1)=0.
      DO I=2,221
      Emax=E(I)
      CALL D01AJF(C14,EMIN,Emax,ERABS,ERREL,S5(I),ER,W,
     *    4000,IW,1040,IFAIL)
      END DO

      CALL SPLINE (221,E,S5,GM5,w1,w2,w3)

      CALL D01AJF(SPf5,EMIN,EMAX,ERABS,ERREL,S5N,ER,W,
     *    4000,IW,1040,IFAIL)

      DO I=1,221
c      S5(I)=S5(I)/S5N
      END DO

      CALL SPLINE (221,E,S5,GM5,w1,w2,w3)


      DO I=1,221
      ep(i)=e(i)
c      write (66,200) e(i),s1(i),s2(i),s3(i),s4(i),s5(i)
200   format (f10.3,5f8.5)
      END DO
c	stop
      print*,'Loss patterns for 2- 5-scatterings are calculated'
	end if	  !spectra splines are precalculated
c scattering patterns interpolation at e_point __________________
	call splint (221,e_point,s1p,E,S1,GM1)
	call splint (221,e_point,s2p,E,S2,GM2)
 	call splint (221,e_point,s3p,E,S3,GM3)
	call splint (221,e_point,s4p,E,S4,GM4)
	call splint (221,e_point,s5p,E,S5,GM5)

      return
      END


      FUNCTION SPf2(DE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,DE,SPf2,E,S2,GM2)
      RETURN
      END

      FUNCTION SPf3(DE)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,DE,SPf3,E,S3,GM3)
      RETURN
      END


      FUNCTION SPf4(DE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,DE,SPf4,E,S4,GM4)
      RETURN
      END


      FUNCTION SPf5(DE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,DE,SPf5,E,S5,GM5)
      RETURN
      END




      FUNCTION C11(DE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,Emax-DE,SPLU,E,S1,GM1)
      CALL SPLINT (221,DE,SPMI,E,S1,GM1)
      C11=SPLU*SPMI
      RETURN
      END


      FUNCTION C12(DE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,Emax-DE,SPLU,E,S2,GM2)
      CALL SPLINT (221,DE,SPMI,E,S1,GM1)
      C12=SPLU*SPMI
      RETURN
      END


      FUNCTION C13(DE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,Emax-DE,SPLU,E,S3,GM3)
      CALL SPLINT (221,DE,SPMI,E,S1,GM1)
      C13=SPLU*SPMI
      RETURN
      END


      FUNCTION C14(DE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(221),S1(221),S2(221),S3(221),S4(221),S5(221)
      DIMENSION GM1(221),GM2(221),GM3(221),GM4(221),GM5(221)
      COMMON/scatt/Emax,E,S1,S2,S3,S4,S5,GM1,GM2,GM3,GM4,GM5
      CALL SPLINT (221,Emax-DE,SPLU,E,S4,GM4)
      CALL SPLINT (221,DE,SPMI,E,S1,GM1)
      C14=SPLU*SPMI
      RETURN
      END



      function fermi(e)
      implicit double precision (a-h,o-z)
      common /const/em,alpha,pi
      beta=sqrt(1.-(em/(em+e))**2)
      y=4.*pi*alpha/beta
      fermi=y/abs(1-dexp(-y))*(1.002037-0.001427*beta)
      return
      end



      subroutine spline2(n,m,f,x,y,cm20,cm22,cm02)

*  Routine to make two-dimensional spline interpolation.
*  Origin: "Methods of Spline Functions", Yu.S.Zavialov, B.I.Kvasov and
*           V.L.Miroshnichenko, Moscow, "Nauka", 1980, page 135.

      implicit double precision (a-h, o-z)
      dimension f(n,m), x(n), y(m), cm20(n,m), cm22(n,m), cm02(n,m)
      dimension fy(200), g(400), r(400), s(400), cmy(400)

*   Calculation of m20

      do j = 1,m
      call spline(n,x,f(1,j),cm20(1,j),g,r,s)
      enddo

*   Calculation of m22

      do i = 1,n
      do j = 1,m
       fy(j) = cm20(i,j)
      enddo
      call spline(m,y,fy,cmy,g,r,s)
      do j = 1,m
       cm20(i,j) = cm20(i,j)/6.
       cm22(i,j) = cmy(j)/36.
      enddo
      enddo

*   Calculation of m02

      do i = 1,n
      do j = 1,m
       fy(j) = f(i,j)
      enddo
      call spline(m,y,fy,cmy,g,r,s)
      do j = 1,m
       cm02(i,j) = cmy(j)/6.
      enddo
      enddo
      end


      subroutine splint2(n,m,f0,xi,yj,cm20,cm22,cm02,x,y,
     1            sxy,dx,dy,idiff,ifail)

*     if (idiff.eq.1) program calculates the first derivatives
*                 dx = D(s(x,y))/Dx and dy = D(s(x,y))/Dy
*     pt(i) = psi(t), ptx(i) = D(psi)/Dt * Dt/Dx
*     pu(i) = psi(u), puy(i) = D(psi)/Du * Du/Dy

      implicit double precision (a-h, o-z)
      dimension f0(n,m), xi(n), yj(m), cm20(n,m), cm22(n,m), cm02(n,m)
      dimension pt(4), pu(4), f(4,4), ft(4), ptx(4), puy(4)

      ifail = 0
      i1 = 0
      if (x.gt.xi(n)) goto 88888
      if (x.le.xi(1)) goto 88888

10    i1 = i1+1
      if (x.gt.xi(i1)) goto 10
      i = i1-1
      h = xi(i1)-xi(i)
      t = (x-xi(i))/h

      j1 = 0
      if (y.gt.yj(m)) goto 88888
      if (y.le.yj(1)) goto 88888

20    j1 = j1+1
      if (y.gt.yj(j1)) goto 20
      j = j1-1
      al = yj(j1)-yj(j)
      u = (y-yj(j))/al

      hq = h*h
      t1 = t-1.
      pt(1)= -t1
      pt(2) = t
      pt(3) = t*t1*(2.-t)*hq
      pt(4) = t*(t*t-1.)*hq

      alq = al*al
      u1 = u-1.
      pu(1) = -u1
      pu(2) = u
      pu(3) = u*u1*(2.-u)*alq
      pu(4) = u*(u*u-1.)*alq

      if (idiff.eq.1) then
      ptx(1) = -1./h
      ptx(2) = -ptx(1)
      ptx(3) = (3.*t*(2.-t)-2.)*h
      ptx(4) = (3.*t*t-1.)*h
      puy(1) = -1./al
      puy(2) = -puy(1)
      puy(3) = (3.*u*(2.-u)-2.)*al
      puy(4) = (3.*u*u-1.)*al
      endif

      f(1,1) = f0(i,j)
      f(2,1) = f0(i1,j)
      f(3,1) = cm20(i,j)
      f(4,1) = cm20(i1,j)

      f(1,2) = f0(i,j1)
      f(2,2) = f0(i1,j1)
      f(3,2) = cm20(i,j1)
      f(4,2) = cm20(i1,j1)

      f(1,3) = cm02(i,j)
      f(2,3) = cm02(i1,j)
      f(3,3) = cm22(i,j)
      f(4,3) = cm22(i1,j)

      f(1,4) = cm02(i,j1)
      f(2,4) = cm02(i1,j1)
      f(3,4) = cm22(i,j1)
      f(4,4) = cm22(i1,j1)

      do i = 1,4
      ft(i) = 0.
      do j = 1,4
       ft(i) = ft(i)+f(i,j)*pu(j)
      enddo
      enddo

      sxy = 0.
      if (idiff.eq.1) dx = 0.
      do i = 1,4
      sxy = sxy+pt(i)*ft(i)
      if (idiff.eq.1) then
       dx = dx+ptx(i)*ft(i)
      endif
      enddo

      if (idiff.eq.1) then
      do i = 1,4
      ft(i) = 0.
      do j = 1,4
       ft(i) = ft(i)+f(i,j)*puy(j)
      enddo
      enddo
      dy = 0.
      do i = 1,4
      dy = dy+pt(i)*ft(i)
      enddo
      endif
      goto 99999

88888     ifail = 1
99999     continue
      end


      SUBROUTINE SPLINE (N,X,Y,CM,G,R,S)

C "SPLINE COEFFICIENTS CALCULATION"

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(N), Y(N), CM(N), G(N), R(N), S(N)

      D = X(2)-X(1)
      E = (Y(2)-Y(1))/D

      DO K = 2,N-1
      H = D
      D = X(K+1)-X(K)
      F = E
      E = (Y(K+1)-Y(K))/D
      G(K) = D/(D+H)
      R(K) = 1.-G(K)
      S(K) = 6.*(E-F)/(D+H)
      ENDDO

      DO K = 2,N-1
      P = 1./(R(K)*G(K-1)+2.)
      G(K) = -G(K)*P
      S(K) = (S(K)-R(K)*S(K-1))*P
      ENDDO

      CM(N) = 0.
      G(N-1) = S(N-1)
      CM(N-1) = G(N-1)

      DO K = N-2,1,-1
      G(K) = G(K)*G(K+1)+S(K)
      CM(K) = G(K)
      ENDDO
      END

      SUBROUTINE SPLINT (N,AR,FN,X,Y,GM)
C "FUNCTION CALCULATION VIA SPLINE APPROXIMATION"
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(N), Y(N), GM(N)

      I1 = 0
      IF (AR.GT.X(N)) GOTO 8
      IF (AR.LE.X(1)) GOTO 9
 3    I1 = I1+1
      IF (AR.GT.X(I1)) GOTO 3
      I = I1-1
      D = X(I1)-X(I)
      H = AR-X(I)
      R = X(I1)-AR
      P = D*D/6.
      FN = (GM(I)*R**3+GM(I1)*H**3)/6./D
      FN = FN+((Y(I)-GM(I)*P)*R+(Y(I1)-GM(I1)*P)*H)/D
      GOTO 10
 8    D = X(N)-X(N-1)
      FN = D*GM(N-1)/6.+(Y(N)-Y(N-1))/D
      FN = FN*(AR-X(N))+Y(N)
      GOTO 10
 9    D = X(2)-X(1)
      FN = -D*GM(2)/6.+(Y(2)-Y(1))/D
      FN = FN*(AR-X(1))+Y(1)
10    CONTINUE
      END
