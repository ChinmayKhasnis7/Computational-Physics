!---------------------------------------------------------------
! Chinmay Khasanis
! Radial distribution function
!---------------------------------------------------------------
      
      MODULE parameters
      
      IMPLICIT NONE
      SAVE
      
      INTEGER, PARAMETER :: npart = 500
      INTEGER, PARAMETER :: nframe = 5000
      INTEGER, PARAMETER :: maxbin = 400     
     
     !  Calculate the box length
      double precision, parameter :: rho = 0.80d0
      double precision, parameter :: boxl = (npart/rho)**(1.0/3.0)
    
      DOUBLE PRECISION,PARAMETER :: pi = 3.141592653589d0
      DOUBLE PRECISION,PARAMETER :: delr = 0.02 !/AA
      
      
      END MODULE parameters
      
!     -----------------------------------------------------------------------
      
      PROGRAM gr
      
      USE parameters
      IMPLICIT NONE
      
      INTEGER :: imol, jmol, iframe
      INTEGER :: bin
      
      DOUBLE PRECISION,DIMENSION (npart) :: rx, ry, rz
      DOUBLE PRECISION, DIMENSION (maxbin) :: Hist, g_r
      DOUBLE PRECISION :: rxij, ryij, rzij, rij     
      DOUBLE PRECISION :: rlower, rupper, r, const, nideal

      CHARACTER (LEN = 6) :: junk
      
      open (unit = 1, file = "traj.xyz", action = "read")
      open (unit = 2, file = "gr.dat")
      
      Hist = 0.0
 
      DO iframe = 1, nframe
         
         READ(1,*)
         READ(1,*)
         DO imol = 1, npart
            READ(1,*) junk, rx(imol), ry(imol), rz(imol)
         END DO
        
         DO imol = 1, npart-1
            DO jmol = imol+1, npart
               
!     calculate the com separation vector
               rxij = rx(jmol) - rx(imol)
               ryij = ry(jmol) - ry(imol)
               rzij = rz(jmol) - rz(imol)
               
!     apply the periodic boundary conditions
               rxij = rxij - dnint(rxij/boxl)*boxl
               ryij = ryij - dnint(ryij/boxl)*boxl
               rzij = rzij - dnint(rzij/boxl)*boxl
               
               rij = dsqrt(rxij**2 + ryij**2 + rzij**2)
               bin = int(rij / delr) + 1
               IF (bin .LE. maxbin) then
                  Hist(bin) =  Hist(bin) + 2
               END IF
               
            END DO
         END DO
         
      END DO
      CLOSE(1)
     
      DO bin  = 1, maxbin
         rlower  = (bin - 1) * delr
         rupper  = rlower + delr
         nideal  = (4.0*pi*rho/3.0) * (rupper**3 - rlower**3)
         g_r(bin) = Hist(bin) / (nframe * npart * nideal)
         WRITE(2,*)rlower, g_r(bin)
      END DO
      
      END 
      
