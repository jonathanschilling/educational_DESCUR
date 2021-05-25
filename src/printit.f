      SUBROUTINE printit(rin,zin,rbc,zbs,rbs,zbc,rmnaxis,zmnaxis)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(*) :: rin, zin
      REAL(rprec), DIMENSION(0:mpol-1,-nphi2:nphi2) ::
     1   rbc, zbs, rbs, zbc
      REAL(rprec), DIMENSION(0:nphi2) :: rmnaxis, zmnaxis
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p25 = 0.25_dp
      INTEGER :: i, m, n, n11
      REAL(rprec) :: tol, rhobcc, rhobss, rhobcs, rhobsc
      REAL(rprec), DIMENSION(0:mpol-1, -nphi2:nphi2) :: 
     1             rhoc, rhos
      CHARACTER*250 :: form_string
!-----------------------------------------------

      OPEN(unit=10, file='plotout', status='unknown')
      WRITE (10, 1990) mpol, ntheta, nphi, mpol-1, nphi2, nfp, mpnt
 1990 FORMAT(7i10)
      DO i = 1, nphi*ntheta
         WRITE (10, 1995) rin(i), zin(i)
      END DO
 1995 FORMAT(1p,2e12.4)

      PRINT *, ' OUTPUTTING FOURIER COEFFICIENTS TO OUTCURVE FILE'
!**********************************************************************
!       This SUBROUTINE prints out data into the file "outcurve"
!       The FORMAT of the modes is compatible with input to VMEC
!**********************************************************************
      WRITE (3, 10)
   10 FORMAT(/'   MB  NB      RBC         RBS         ',
     1   'ZBC         ZBS        RAXIS       ZAXIS')
      tol = 1.E-6_dp*ABS(rbc(1,0))
      DO m = 0, mpol-1
         DO n = -nphi2, nphi2
            WRITE (10, 2000) rbc(m,n), zbs(m,n), rbs(m,n), zbc(m,n)
            IF (.not.(ABS(rbc(m,n)).lt.tol .and. ABS(zbs(m,n)).lt.tol
     1          .and. ABS(rbs(m,n)).lt.tol .and. ABS(zbc(m,n)).lt.tol))
     2      THEN
               IF (m.eq.0 .and. n.ge.0) THEN
                  WRITE (3, 30) m, n, rbc(m,n), rbs(m,n), zbc(m,n),
     1               zbs(m,n), rmnaxis(n), zmnaxis(n)
               ELSE
                  WRITE (3, 40) m, n, rbc(m,n), rbs(m,n), zbc(m,n),
     1               zbs(m,n)
               ENDIF
            ENDIF
         END DO
      END DO
   30 FORMAT(i5,i4,1p,6e12.4)
   40 FORMAT(i5,i4,1p,4e12.4)
 2000 FORMAT(1p,4e12.4)

!**********************************************************************
!     WRITE OUT IN FORMAT THAT CAN BE IMPORTED INTO VMEC
!**********************************************************************
      WRITE(3, *)
      DO m = 0, mpol-1
         DO n = -nphi2, nphi2
            IF (.not.(ABS(rbc(m,n)).lt.tol .and. ABS(zbs(m,n)).lt.tol
     1          .and. ABS(rbs(m,n)).lt.tol .and. ABS(zbc(m,n)).lt.tol))
     2      THEN
               n11 = nphi2/10 + 1
               IF (n .lt. 0) n11 = n11+1
               WRITE(form_string,'(a,4(a,i1,a,i1,a))')"(2x,",
     1           "'RBC(',i",n11,",',',i",m/10+1,",') = ',1p,e14.6,3x,",
     2           "'RBS(',i",n11,",',',i",m/10+1,",') = ',e14.6,3x,",
     3           "'ZBC(',i",n11,",',',i",m/10+1,",') = ',e14.6,3x,",
     4           "'ZBS(',i",n11,",',',i",m/10+1,",') = ',e14.6)"
               WRITE(3, form_string) n,m,rbc(m,n), n,m,rbs(m,n), 
     1               n,m,zbc(m,n), n,m,zbs(m,n)
            END IF
         END DO
      END DO

!**********************************************************************
!     WRITE OUT R0C, Z0S, RHOCC, RHOSS, RHOCS, RHOSS THAT CAN BE IMPORTED INTO V3FIT/VMEC
!**********************************************************************
      IF (HB_Parameter .NE. one) RETURN

      WRITE (3, *)
      WRITE(3, *) '    MB       T1         T2'
      DO m = 0, mrho-1
         WRITE(3, 100) m,t1m(m),t2m(m)
      DO n=-nphi2, nphi2
         rhobcc = p25*(rbc(m+1,n)+zbs(m+1,n))/t1m(m+1)
         IF (m .LE. 1) THEN
            rhobss = rhobcc
         ELSE
            rhobss = p25*(rbc(m-1,n)-zbs(m-1,n))/t2m(m-1)
         END IF
         rhoc(m,n) = rhobcc + rhobss
         rhobcs = p25*(rbs(m+1,n)-zbc(m+1,n))/t1m(m+1)
         IF (m .LE. 1) THEN
            rhobsc = rhobcs
         ELSE 
            rhobsc = p25*(rbs(m-1,n)+zbc(m-1,n))/t2m(m-1)
         END IF
         rhos(m,n) = rhobcs+rhobsc
      END DO
      END DO
      WRITE(3, 100) mrho,t1m(mrho),t2m(mrho)

 100  FORMAT(4x,i3,2(1p,e12.3))

      WRITE(3, *)
      DO m = 0, mrho-1
         DO n = 0, nphi2

            IF ((ABS(rhoc(m,n)).lt.tol) .AND. 
     1          (ABS(rhos(m,n)).lt.tol)) CYCLE

            n11 = nphi2/10 + 1
            IF (m .eq. 0) THEN
            END IF
            WRITE(form_string,'(a,2(a,i1,a,i1,a))')"(2x,",
     1          "'RHOC(',i",n11,",',',i",m/10+1,",') = ',1p,e14.6,3x,",
     2          "'RHOS(',i",n11,",',',i",m/10+1,",') = ',e14.6)"
            WRITE(3, form_string) n,m,rhoc(m,n), n,m,rhos(m,n) 
         END DO
      END DO

      END SUBROUTINE printit
