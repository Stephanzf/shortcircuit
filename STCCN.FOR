	subroutine stccncal
	PARAMETER (NBU=70,MAXEVENT=5, NTR0=30  ,
     *           INTWKI=15, INTWKO=16, IZGNDO=17, IZLEQO=18)
	COMMON/S15/Y(NBU,NBU),BUSTYPE(NBU),   PG(NBU),
     *                QG(NBU), PL(NBU),QL(NBU),VSPEC(NBU),
     *                QMIN(NBU),QMAX(NBU) ,VMAG(NBU), VANG(NBU)
	COMMON/S151/  DF(NBU)  ,DX(NBU)   ,A(NBU,NBU)
	COMMON/S16 /B1P(NBU,NBU) , B2P(NBU,NBU) , 
     *                YIJ0(2*NBU), IJ0(2*NBU, 2)
	COMMON/S161/B0P(NBU,NBU), ZEROBUS(2,NBU)
	COMMON/S17/FLTINFM( MAXEVENT,5 ), FLTINDI( 5 ) 
	
	COMMON/S20/ FLTZIJ( MAXEVENT, 4)
	common/s21/ B1P0(NBU,NBU) , B2P0(NBU,NBU) ,B0P0(NBU,NBU) 
	complex     B1P0    , B2P0   , B0P0 

	INTEGER     BUSTYPE ,  NBUS     , NFAULT, NEVENT  , NBRANCH
	INTEGER     ZEROBUS ,  NBUSINZERO       , IJ0 
	COMPLEX     Y   , B1P   , B2P   , B0P   , A    , DX  , DF
	COMPLEX     YIJ , Zgnd  , CMPLX , Zflt  , YIJ0 , B0PIJ0
	COMPLEX     Zsu , Zkk0  , Zkk1  , ZKK2  , XIJ0 , XIJ1
	
	COMPLEX     FLTZIJ

	CHARACTER*20    FLTINDI  
	CHARACTER*75    STR

	OPEN(INTWKI,FILE='NTWK.INP' , STATUS = 'OLD' )
	OPEN(INTWKO,FILE='NTWK.OUT')
	OPEN(IZGNDO,FILE='ZGND.OUT')
	OPEN(IZLEQO,FILE='ZLEQ.OUT')
C
	WRITE(  *   , '(/A38/)' ) '   Reading datafile NTWK.INP ...'
	
	CALL STRING( INTWKI, INTWKO )
100     READ ( INTWKI, *, END = 110 ) NBUS, NFAULT, NEVENT
	WRITE( INTWKO, 99 )           NBUS, NFAULT, NEVENT
99      FORMAT(T10, I4, T17, I3 ,T26, I3)
110     DO  I = 1, NBUS
	DF (I) = ( 0.0, 0.0 ) 
	ZEROBUS(1, I ) = 0
	ZEROBUS(2, I ) = 0
	DO  J=1,NBUS
	B1P(I,J)  = CMPLX( 0.0, 0.0 )
	B2P(I,J)  = CMPLX( 0.0, 0.0 )
	B0P(I,J)  = CMPLX( 0.0, 0.0 )
	Y(I,J)    = CMPLX( 0.0, 0.0 )
	END DO
	END DO
	NZsu = 0
C   read fault information to FLTINFM( ):
	
	REA = 0.0
	CALL STRING( INTWKI, INTWKO )
	DO II = 1, NEVENT
119     READ(INTWKI,*, END=130) I, J, J1,
     *   ( FLTINFM(II,K),K=4,5 ) ,Zgnd,  Zflt
	FLTINFM(II,1) = I
	FLTINFM(II,2) = J
	FLTINFM(II,3) = J1
	FLTZIJ(II,1)    = Zgnd
	FLTZIJ(II,2)    = Zflt
130     WRITE(INTWKO ,132)I,J,J1, ( FLTINFM(II,K), K=4,5),
     +                        (FLTZIJ(II,K),K=1,2)
132     FORMAT(T5,I4,T11,I4,T18,I2, T26, F6.3,T34,  F6.3,
     +                              T41,2F7.2,T57, 2F7.2 )
	END DO
C   read fault type indicator
	II = 1
135     CALL STRING( INTWKI, INTWKO )
	READ(INTWKI, '(A22)', ERR = 140 ) STR
138     FORMAT( A22 )
	IF (STR.EQ.'3-phase-to-ground') THEN
	    FLTINDI(II)='3-phase-to-ground'
	    II = II + 1
	END IF
	IF (STR.EQ.'1-phase-to-ground')  THEN
	    FLTINDI(II)='1-phase-to-ground'
	    II = II + 1
	END IF
	IF (STR.EQ.'2-phase-to-ground') THEN
	    FLTINDI(II)='2-phase-to-ground'
	    II = II + 1
	END IF
       IF (STR.EQ.'phase-to-phase')   THEN
	   FLTINDI(II)='phase-to-phase'
	    II = II + 1
	END IF
	IF (STR.EQ.'1-phase-broken')   THEN
	    FLTINDI(II)='1-phase-broken'
	    II = II + 1
	END IF
	IF (STR.EQ.'2-phase-broken')   THEN
	    FLTINDI(II)='2-phase-broken'
	    II = II + 1
	END IF
	IF ( II.GT.NFAULT ) GOTO 140
	GOTO 135
140     CONTINUE
	WRITE(INTWKO , '(1X, A25)' ) ( FLTINDI(J),J=1, NFAULT )
C
C    positive (1) et negative (2) admittance matrices:
C   read each line data & assembly  Y, B1P, B2P matrix
C

162     CALL STRING( INTWKI, INTWKO )
165     READ(INTWKI, *, END=180 ) I, J, R1, R2, X1, X2, BK
180     IF (I.EQ.0) GOTO 5
	WRITE(INTWKO, 188 ) I, J, R1, R2, X1, X2, BK
188     FORMAT(4X, 2I4, 1X, 5F9.5)
C   line branch assembly:
	IF (J.GT.0) THEN
	  IF ( (I.EQ.ABS(FLTINFM(1,1))).AND.
     +         ( J.EQ.ABS(FLTINFM(1,2)) ) ) BI0 = BK
	  YIJ      = ( 1.0, 0.0 )/CMPLX( R1, X1 )
	  NUL      = 0.0
	  Y(I,I)   = Y(I,I) + YIJ + CMPLX( NUL, BK/2.0 )
	  Y(J,J)   = Y(J,J) + YIJ + CMPLX( NUL, BK/2.0 )
	  Y(I,J)   = Y(I,J) - YIJ
	  Y(J,I)   = Y(J,I) - YIJ
	  YIJ      = ( 1.0, 0.0 )/CMPLX( R2, X2 ) 
	  B2P(I,J) = B2P(I,J) - YIJ
	  B2P(J,I) = B2P(I,J)
	  B2P(I,I) = B2P(I,I) + YIJ + CMPLX( NUL, BK/2.0 ) 
	  B2P(J,J) = B2P(J,J) + YIJ + CMPLX( NUL, BK/2.0 ) 
	END IF
C   transformer line assembly:
	IF (J.LT.0) THEN
	  J=-J
	  YIJ      = ( 1.0, 0.0 )/CMPLX( R1,X1 )
	  Y(I,I)   = Y(I,I) + YIJ
	  Y(J,J)   = Y(J,J) + YIJ/BK/BK
	  Y(I,J)   = Y(I,J) - YIJ/BK
	  Y(J,I)   = Y(J,I) - YIJ/BK
	  YIJ      = ( 1.0, 0.0 )/CMPLX( R2, X2 ) 
	  B2P(I,J) = B2P(I,J) - YIJ/BK
	  B2P(J,I) = B2P(I,J)
	  B2P(I,I) = B2P(I,I) + YIJ 
	  B2P(J,J) = B2P(J,J) + YIJ/BK/BK  
	END IF
C   grounded line assembly:
	IF (J.EQ.0) THEN
	  IF (BK.EQ.1.0) THEN
	    Y(I,I)   = Y(I,I)   + ( 1.0, 0.0 )/CMPLX( R1, X1 ) 
	    B2P(I,I) = B2P(I,I) + ( 1.0, 0.0 )/CMPLX( R2, X2 ) 
	  END IF
	  IF (BK.EQ.-1.0) THEN
	    Y(I,I)   =  Y(I,I)   + CMPLX( R1, X1 )   
	    B2P(I,I) =  B2P(I,I) + CMPLX( R2, X2 )   
	  END IF
	END IF
C    go back to read next line data
	GOTO 165
5       CONTINUE
	DO 10  I = 1, NBUS
	DO 10  J = 1, NBUS
10      B1P(I,J) = Y(I,J)
C
C     zero sequence (0) admittance matrix:
C   read each line data & assembly B0P matrix
C
	CALL STRING( INTWKI, INTWKO )
	L = 1
2       READ(INTWKI, *, ERR =200, END =200 ) I, J, R0, X0, BK
200     IF (I.EQ.0) GOTO 15
	WRITE(INTWKO, 205 ) I, J, R0, X0, BK
205     FORMAT( 4X, 2I4,1X,3F9.5 )
C
C   line branch assembly:
	IF (J.NE.0) THEN
C   store non grounded line data for B0P modification
	  IF ( (I.EQ.ABS(FLTINFM(1,1))).AND.
     +         ( J.EQ.ABS(FLTINFM(1,2)) ) ) B0PIJ0 =CMPLX(R0,X0)
	  IF ( (J.EQ.ABS(FLTINFM(1,1))).AND.
     +         ( I.EQ.ABS(FLTINFM(1,2)) ) ) B0PIJ0 =CMPLX(R0,X0)
       
	   IJ0(L,1) = I
	   IJ0(L,2) = J
	  NUL       = 0.0
	  YIJ       = ( 1.0, 0.0 )/CMPLX( R0, X0 ) 
	  YIJ0(L)   = YIJ
	  L = L + 1
	  B0P(I,J)  = B0P(I,J) - YIJ
	  B0P(J,I)  = B0P(I,J)
	  B0P(I,I)  = B0P(I,I) + YIJ + CMPLX( NUL, BK/2.0 ) 
	  B0P(J,J)  = B0P(J,J) + YIJ + CMPLX( NUL, BK/2.0 ) 
	END IF
C   grounded line assembly:
	IF (J.EQ.0) THEN
	  IF (BK.EQ.1.0) THEN
	    B0P(I,I) = B0P(I,I) + ( 1.0, 0.0 )/CMPLX( R0, X0 ) 
	  END IF
	  IF (BK.EQ.-1.0) THEN
	    B0P(I,I) =  B0P(I,I) + CMPLX( R0, X0 )   
	  END IF
	END IF
C    go back to read next line data
	GOTO 2
15      CONTINUE
	NBRANCH = L - 1
C
C   equivalent load admittance added to seq networks
C
156     READ ( IZLEQO, 163, END =175 ) I, YIJ
163     FORMAT( I5, 1X, 2F15.7  )
	IF (I.EQ.0) GOTO 175
	Y (I,I) = Y (I,I ) + YIJ
	B0P(I,I ) = B0P( I,I ) + YIJ
	B1P(I,I ) = B1P( I,I ) + YIJ
	B2P(I,I ) = B2P( I,I ) + YIJ
	GOTO 156
175     CONTINUE
C   end of (+), (-) et (0) sequence Y matrices assembly.
C
C   read in bus numbers which are either 
C     in or not in zero seq network
C
	CALL STRING( INTWKI, INTWKO )
230     READ ( INTWKI, *, END = 242  ) NBUSINZERO
242     READ ( INTWKI, *, END = 233 )  ( ZEROBUS(1,I), I=1, NBUS )
233     READ ( INTWKI, *, END = 243  ) NBUSINZERO
243     READ ( INTWKI, *, END = 244 )  ( ZEROBUS(2,I), I=1, NBUS )
C        
C   given those buses which are in (1) the zero seq network 
C      find those which are not in (-1)
C
244     IF ( NBUSINZERO.EQ.1 ) THEN
	DO KK = 1, 2
	WRITE(INTWKO, 252 )  ( ZEROBUS(KK, I), I=1, NBUS )
	END DO

	DO KK = 1, 2
	  DO K1=1, NBUS
	  DX(K1)=0.0
	  END DO

	  II = 1
	  DO I = 1, NBUS
	     J = 0   
	     L = 1
245          IF ( ZEROBUS(KK,L).EQ.0 ) GOTO 250
	     IF ( I.EQ.ZEROBUS(KK, L) ) THEN
		J = 1
	     END IF
	     L = L + 1
	     GOTO 245                    
250          IF ( J.EQ.0 )  THEN
		DX( II ) = I
		II = II + 1
	     END IF
	  END DO
	  DO L = 1, NBUS
	  ZEROBUS(KK,L) = DX(L)
	  END DO
	END DO ! KK LOOP
	END IF
	WRITE( INTWKO,   '( 3X, I5)' )  - IABS( NBUSINZERO )
	DO KK = 1, 2
	WRITE(INTWKO, 252 )  ( ZEROBUS(KK, I), I=1, NBUS )
	END DO
252     FORMAT(3X, 10I6  )


	CALL STRING( INTWKI, INTWKO )

	
	WRITE(  *   , '( A35/A35)' ) 'Computing the Z eq impedance',
     +                               '    at the fault point ... '
c --- store the original (+),(-) and (0) seq matrices -----
	do i=1, nbus
	do j=1, nbus
	B1P0(I,J) = B1P(I,J)
	B2P0(I,J) = B2P(I,J)
	B0P0(I,J) = B0P(I,J)
	end do
	end do
c ------- Zeq for ievent=grouned and =open ciruit -----
	
	L0 = 1
	Zsu = (0.0,0.0)

	DO I0 = 1,3   ! I0 LOOP
	do i=1, nbus
	do j=1, nbus
	B1P(I,J) = B1P0(I,J)
	B2P(I,J) = B2P0(I,J)
	B0P(I,J) = B0P0(I,J)
	end do
	end do

	   IF ( FLTINFM(I0,3).NE.0.0) THEN
C        
C   modify B0P according to ZEROBUS and YIJ0 and IJ0
	DO J = 1, 2
	   DO I = 1, NBRANCH
	      DO K = 1, NBUS
	      IF ( ZEROBUS(L0,K).EQ.IJ0(I,J) ) THEN
		 I1 = IJ0(I,1)
		 J1 = IJ0(I,2)
		 B0P( I1, J1 ) = B0P( I1, J1 ) + YIJ0(I)
		 B0P( J1, I1 ) = B0P( J1, I1 ) + YIJ0(I)
		 B0P( I1, I1 ) = B0P( I1, I1 ) - YIJ0(I)
		 B0P( J1, J1 ) = B0P( J1, J1 ) - YIJ0(I)
	      END IF
	      END DO
	   END DO
	END DO

C       WRITE(INTWKO,5233) ((B0P(I,J),J=1,NBUS),I=1,NBUS)
5233    FORMAT(1X,'B0P '/,5(E12.4,F12.4,'J',2X))



	IF ( FLTINDI(1).EQ.'3-phase-to-ground') THEN
		 IF ( L0.EQ.1 ) THEN
		    Zsu = FLTZIJ(I0,2)               
		    FLTZIJ(I0,3) = Zsu
		    GOTO   38
		 END IF
	END IF
		 IF( ( L0.EQ.1 ).AND. 
     +               ( FLTINDI(1).EQ.'1-phase-to-ground')) THEN
		     IF (FLTINFM(I0,1).LT.0.0) II0= FLTINFM(I0,1)
		     IF (FLTINFM(I0,2).LT.0.0) II0= FLTINFM(I0,2)
		     DF(-II0) = ( 1.0, 0.0 )
		     idfp = -II0
		     idfn = 0
		 END IF
		 IF(( ( L0.EQ.2 ).AND. 
     +                ( FLTINDI(1).EQ.'1-phase-to-ground')).OR. 
     +              ( FLTINDI(1).EQ.'1-phase-broken'      )    )THEN
C       phase a of line i-j is switched off
		      I1 = ABS(FLTINFM(I0,1))
		      J1 = ABS(FLTINFM(I0,2))
c  --------- if I1 or J1 not in (0) seq network, Z(0)=infinite -----
	 DO KL = 1, NBUS
	 IF ((I1.EQ.ZEROBUS(L0,KL)).OR.( J1.EQ.ZEROBUS(L0,KL))) NZsu=1
	 END DO
C 
	XIJ1 = FLTINFM(I0,3)
	B0P( I1, I1 ) = B0P( I1, I1 )+B0P( I1, J1 )/XIJ1-BI0 
	B0P( J1, J1 ) = B0P( J1, J1 )+B0P( I1, J1 )/XIJ1-BI0
	B0P( I1, J1 ) = B0P( I1, J1 )*(1-1.0/XIJ1) 
	B0P( J1, I1 ) = B0P( I1, J1 )*(1-1.0/XIJ1) 
	B1P( I1, I1 ) = B1P( I1, I1 )+B1P( I1, J1 )/XIJ1-BI0
	B1P( J1, J1 ) = B1P( J1, J1 )+B1P( I1, J1 )/XIJ1-BI0
	B1P( I1, J1 ) = B1P( I1, J1 )*(1-1.0/XIJ1)  
	B1P( J1, I1 ) = B1P( I1, J1 )*(1-1.0/XIJ1) 
	B2P( I1, I1 ) = B2P( I1, I1 )+B2P( I1, J1 )/XIJ1-BI0 
	B2P( J1, J1 ) = B2P( J1, J1 )+B2P( I1, J1 )/XIJ1-BI0
	B2P( I1, J1 ) = B2P( I1, J1 )*(1-1.0/XIJ1)  
	B2P( J1, I1 ) = B2P( I1, J1 )*(1-1.0/XIJ1) 
	IF (FLTINFM(I0,1).LT.0.0) THEN
	  DF(I1)=(1.0,0.0)
	  DF(J1)=(-1.0,0.0)
	  idfp =I1
	  idfn =J1
	END IF
	IF (FLTINFM(I0,2).LT.0.0) THEN
	  DF(J1)=(1.0,0.0)
	  DF(I1)=(-1.0,0.0)
	  idfp =J1
	  idfn =I1
	END IF
		 END IF  ! L0=2 ..... LOOP
	L = 1
255     IF ( ZEROBUS(L0, L).EQ.0 ) GOTO 360
	DO I = ZEROBUS(L0, L) - L + 1, NBUS
	  DO J = 1, NBUS
	  B0P(I,J) = B0P(I+1,J)
	  END DO
	END DO  
	DO I = ZEROBUS(L0,L) - L + 1, NBUS
	  DO K = 1, NBUS
	  B0P(K,I) = B0P(K,I+1)
	  END DO
	DF(I) = DF(I+1)
	END DO
	  L = L + 1
	  GOTO  255
360     NBUSINZERO = NBUS - L + 1
C -------- zero SEQ NET -------------------------        
3333    CONTINUE

       WRITE(INTWKO,523) ((B0P(I,J),J=1,NBUSINZERO),I=1,NBUSINZERO)
523     FORMAT(1X,'B0P '/,21(E8.4,E8.4,'J',1X))
	IF (NZsu.EQ.1) THEN
	ZKK0 = (1.0E6,1.0E6)
	GOTO 444
	END IF

      write(intwko,*) '0,idfp,idfn=', idfp,idfn
	CALL FACTERIZE ( B0P, NBUSINZERO )
	
	JJ1 = 0
	DO J0 = 1, NBUSINZERO
	IF ( DF(J0).EQ.(1.0,0.0) )  II1=J0
	IF ( DF(J0).EQ.(-1.0,0.0) ) JJ1=J0
	END DO
	write(intwko,'(i3,2e12.4)') (i, df(i),i=1,nbusinzero)

	CALL AXEQB ( B0P, DX, DF, NBUSINZERO )
C        write(intwko,'(i3,2e12.4)') (i, dx(i),i=1,nbusinzero)


	IF ( JJ1.EQ.0 ) THEN     
	  YIJ = (0.0, 0.0)
	ELSE
	  YIJ =DX(JJ1)
	END IF
	    
	    Zkk0 = DX(II1)- YIJ
444         WRITE(*, 528)  Zkk0
	    WRITE(INTWKO,528) Zkk0

528     FORMAT(/1X,'Z (0) = ', E12.5, 1X, E12.5,'J' )

C   ------- positive, negative -------------------        
	DO I=1, NBUS
	DF(I) = (0.0,0.0)
	END DO

c print Y matrix:
	do im=1,nbus
	write(intwko, '(i4, 2e12.5)') im, b1p(im,im)
	end do

	CALL FACTERIZE ( B1P, NBUS )

	DF( idfp ) = ( 1.0, 0.0 )
	IF ( idfn.NE.0 ) DF( idfn ) = (-1.0, 0.0 )
	write(intwko,*) '+,idfp,idfn=', idfp,idfn
	write(intwko,'(i3,2e12.4)') (i, df(i),i=1,nbus)


	CALL AXEQB ( B1P, DX, DF, NBUS  )

	IF ( idfn.EQ.0 ) THEN     
	  YIJ = (0.0, 0.0)
	ELSE
	  YIJ =DX(idfn)
	END IF

	Zkk1 = DX(idfp) - YIJ

	WRITE(*,524)  Zkk1
	WRITE(INTWKO,524)  Zkk1 
524     FORMAT(1X,'Z (+) = ', E12.5, 1X, E12.5,'J'  )
	
	Zkk2=Zkk1

	STR =  FLTINDI(1) 
	IF  (STR.EQ.'1-phase-to-ground') THEN 
	    IF (L0.EQ.1) Zsu = Zkk2 + Zkk0 + 3.0 * FLTZIJ(I0,2)
	END IF
	IF(( ( L0.EQ.2 ).AND. 
     +       ( STR.EQ.'1-phase-to-ground')).OR. 
     +     ( STR.EQ.'1-phase-broken'      )    )THEN
		      I1 = ABS(FLTINFM(I0,1))
		      J1 = ABS(FLTINFM(I0,2))
	    XIJ1 =-FLTINFM(I0,3)/B1P0(I1,J1)
	    YIJ  =-FLTINFM(I0,3)/B1P0(I1,J1)

	    XIJ0 =-FLTINFM(I0,3)/B0PIJ0
	     Zsu = Zkk0 + XIJ0
	     Zsu = (Zkk2 + YIJ ) * Zsu/((Zkk2+YIJ)+Zsu)+XIJ1
	END IF
	
	WRITE( *, 260 ) Zsu
	WRITE(INTWKO, 260 ) Zsu
260     FORMAT(/2X, 'Z eq = ', E12.5, 1X, E12.5, 'J' /)
	FLTZIJ(I0,3) = Zsu

77      L0 = L0 + 1
	IF (L0.GT.2) GOTO 89
	    END IF  ! f3 = 2 OR 1
	
	END DO  !  I0 LOOP

89      CONTINUE

C        
C   output to file ZGND.OUT for future use
C
38      DO I = 1, NEVENT
	WRITE(IZGNDO,280)( FLTINFM(I,K),K=1,5), FLTZIJ(I,3),BI0
	END DO
	WRITE(IZGNDO , '( A20)' ) ( FLTINDI(J),J=1, NFAULT )

280     FORMAT( 5F10.4, 3E12.6 )
	return
	END

