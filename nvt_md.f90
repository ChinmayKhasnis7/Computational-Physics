!---------------------------------------------------------------
! Chinmay Khasanis
! MD simulation in canonical ensemble
! 500 Lennard-Jones particles
!---------------------------------------------------------------

module parameters
  implicit none
  save
  
  integer, parameter :: npart = 500             ! number of partciles
  integer, parameter :: nstep = 500000          ! total no. of timee steps..

  double precision, parameter :: Tb = 1.0d0  ! initial temperature of the system..
  
  double precision, parameter :: dens = 0.80d0     ! number density
  double precision, parameter :: dt = 0.0020d0     ! time step
  double precision, parameter :: dt2 = dt*dt       ! square of the time step..

  integer, parameter :: thermo_dump = 10 ! dump thermodynamic variables
  integer, parameter :: traj_dump = 100  ! dump trajectory 

  !  Calculate the box length 
  double precision, parameter :: boxl = (npart/dens)**(1.0/3.0)

  ! LJ potential cut-off
  double precision, parameter :: rc = 2.50d0, rc2 = rc*rc

  ! Berendsen thermostat time scale...
  double precision, parameter :: tau = 50*dt !( or 100*dt)

end module parameters

! ---------------------------------------------------
! Main program..
program md_nvt_basics
  
  use parameters
  implicit none

  integer :: i, istep, countn

  double precision :: rx(npart), ry(npart), rz(npart) 
  double precision :: vx(npart), vy(npart), vz(npart)
  double precision :: fx(npart), fy(npart), fz(npart)
  double precision :: pe, ke, tot_e, Tinst
  double precision :: t
  double precision :: lambda

  common / coords / rx, ry, rz
  common / vels /  vx, vy, vz
  common / forces / fx, fy, fz

  ! Open files to write output...
  open (unit = 101, file = "thermo.out", action="write")
  open (unit = 102, file = "temp.out", action="write")
  open (unit = 103, file = "traj.xyz", action="write")

  ! Set up the simulation ..
  call setup

  ! initialize time
  t = 0.0d0

    ! Main loop
  do istep = 1, nstep  ! loop over the time...

     print*, istep

     ! compute time
     t = t + (istep - 1)*dt

     ! Use velocity verlet algorithm to propage the trajectory...
     ! Implimentation..
     do i = 1, npart  ! Loop over the no. of particles...

       ! calculate velocity at half time step
        vx(i) =  vx(i) + 0.5*dt*fx(i)
        vy(i) =  vy(i) + 0.5*dt*fy(i)
        vz(i) =  vz(i) + 0.5*dt*fz(i)

        ! Calculate positions at the next time (Delta_T) step..
        rx(i) = rx(i) + vx(i)*dt 
        ry(i) = ry(i) + vy(i)*dt 
        rz(i) = rz(i) + vz(i)*dt

        ! Ensure that particles are inside the simulation box..
        if (rx(i) .lt. 0.0) rx(i) = rx(i) + boxl
        if (ry(i) .lt. 0.0) ry(i) = ry(i) + boxl
        if (rz(i) .lt. 0.0) rz(i) = rz(i) + boxl
        if (rx(i) .gt. boxl) rx(i) = rx(i) - boxl
        if (ry(i) .gt. boxl) ry(i) = ry(i) - boxl
        if (rz(i) .gt. boxl) rz(i) = rz(i) - boxl
     end do

     ! Update force after getting the new position...
     call comp_force(rx,ry,rz,fx,fy,fz,pe)

     ! Update the velocity for second half time step
     ke = 0.0
     do i = 1, npart
        vx(i) = vx(i) + 0.5*dt*fx(i)
        vy(i) = vy(i) + 0.5*dt*fy(i)
        vz(i) = vz(i) + 0.5*dt*fz(i)
        ke = ke + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
     end do
     ke = 0.5*ke

    ! -- This part is for fixing T

    ! Berendsen thermostat
    ! Q. Change the tau parameter and plot the temperature with time...
    tinst = 2.0*ke/(3.0*npart)
    lambda = dsqrt(1+(dt/tau)*((tb/tinst)-1.0))

    ke = 0.0d0
    do i = 1, npart
       vx(i) = vx(i)*lambda
       vy(i) = vy(i)*lambda
       vz(i) = vz(i)*lambda
       ke = ke + vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
    end do
    ke = 0.5*ke

    ! --

     tot_e = pe + ke

    if (mod(istep,thermo_dump) .eq. 0) then
       ! Compute instantaneous temperature
       Tinst = 2.0*ke/(3.0*npart)
       write(101,5) istep, ke, pe, tot_e
       write(102,*) istep, Tinst
    end if

    ! Print trajectories here.... 
    if ((istep .gt. 10000) .and. (mod(istep,traj_dump) .eq. 0)) then
       write(103,*) npart
       write(103,*) 
       do i = 1, npart
          write(103,*) "O", rx(i), ry(i), rz(i)
       end do
    end if

  end do
  close(101)
  close(102)
  close(103)

5 format(i8,2x,3(f15.6,2x))

end program md_nvt_basics

!  =========================================================
!     setup the positions and other simulation parameters

subroutine setup
  
  use parameters
  implicit none
  
  integer :: ipart
  
  double precision, dimension (npart) :: rx, ry, rz, vx, vy, vz
  double precision, dimension (npart) :: fx, fy, fz
  double precision :: gauss, dummy
  double precision :: ke, pe, Tinst
  
  common / coords / rx, ry, rz
  common / vels /  vx, vy, vz
  common / forces / fx, fy, fz

  ! Open file to write the initial configuration in xyz format...
  open (unit = 11, file = "init_config.xyz", action="write")

  !  Generate inital configuration box length between 0 and 1
  call fcc
  
  ! scale coordinates with the box length
  do ipart = 1, npart
     rx(ipart) = rx(ipart) * boxl
     ry(ipart) = ry(ipart) * boxl
     rz(ipart) = rz(ipart) * boxl
  end do

  ! Write Coordinate in the XYZ format..
  write(11,*) npart
  write(11,*) "inital configuration"
  do ipart = 1, npart
     write(11,*) "O", rx(ipart), ry(ipart), rz(ipart)
  end do
  close(11)

  ! generate distribution of velocity components
  ! Uniformaly distributed with zero average and velocity distribution consistent with the temp..
  call gen_vel
  
  ! Calculate total kinetic energy
  ! mass is 1 in reduced unit.
  ke = 0.0
  do ipart = 1, npart
     ke = ke + vx(ipart)*vx(ipart)+vy(ipart)*vy(ipart)+vz(ipart)*vz(ipart)
  end do
  ke = 0.5*ke
 
  ! Calculate potential energy and force ... 
  call comp_force(rx, ry, rz, fx, fy, fz, pe)
 
  ! Calculate instantaneous temperature..
  ! Remember equipartition theorem..

  Tinst = 2.0*ke/(3.0*npart)

  write(*,*) "intial kinetic energy = ", ke
  write(*,*) "intial potential energy = ", pe
  write(*,*) "Total energy = ", pe + ke
  write(*,*) "intial system temperature = ", Tinst
  
  write(*,*) " setup complete ......."
 
  ! Check the following...
  ! Q. Instantaneous temp is (approx.) equal to Tb
  ! Q. Average velocity is zero..

  return
end subroutine setup
 
!     ******************************************************************
!     Subroutine to compute force
!     Q. Derive the equations coded here for the Lennard-Jones interaction 

subroutine comp_force(rx, ry, rz, fx, fy, fz, pe)
  
  use parameters
  implicit none
  
  integer :: ipart, jpart
  
  double precision, dimension (npart) :: rx, ry, rz
  double precision, dimension (npart) :: fx, fy, fz
  double precision :: dx, dy, dz, r2, r6i, f, pe

  pe    = 0.0d0
  fx(:) = 0.0d0
  fy(:) = 0.0d0
  fz(:) = 0.0d0

  do ipart = 1, npart-1
     do jpart = ipart + 1, npart
        dx = rx(ipart) - rx(jpart)
        dy = ry(ipart) - ry(jpart)
        dz = rz(ipart) - rz(jpart)
        ! Periodic boundary condition...
        dx = dx - dnint(dx/boxl)*boxl
        dy = dy - dnint(dy/boxl)*boxl
        dz = dz - dnint(dz/boxl)*boxl
        r2 = dx*dx + dy*dy + dz*dz
        ! Use the potential cut-off..
        if (r2 .lt. rc2) then   
           r6i = 1.0 / (r2*r2*r2)
           ! pe computation 
           pe  = pe + 4*(r6i*r6i - r6i)  ! Lennard-Jones Potential
           ! force computation
           f   = 48*(r6i*r6i-0.5*r6i)
           fx(ipart) = fx(ipart) + dx*f/r2  ! LJ force - X 
           fx(jpart) = fx(jpart) - dx*f/r2
           fy(ipart) = fy(ipart) + dy*f/r2
           fy(jpart) = fy(jpart) - dy*f/r2 
           fz(ipart) = fz(ipart) + dz*f/r2
           fz(jpart) = fz(jpart) - dz*f/r2
        end if
     end do
  end do
  
 return
end subroutine comp_force

! ==================================================================
! generate velocity and take away CoM velocities.. - uniformaly distributed..

subroutine gen_vel

  use parameters
  implicit none
  
  integer :: i
  
  double precision :: dummy, ranf
  double precision :: sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2
  double precision :: sumv2, fs

  external ranf

  double precision, dimension (npart) :: vx, vy, vz

  common / vels /  vx, vy, vz
 
  ! Iniitialize the variables... VERY IMPORTANT!  
  sumvx = 0.0d0
  sumvy = 0.0d0
  sumvz = 0.0d0
  sumv2 = 0.0d0
  do i = 1, npart
     vx(i) = (ranf(dummy) - 0.5)
     vy(i) = (ranf(dummy) - 0.5)
     vz(i) = (ranf(dummy) - 0.5)
     sumvx = sumvx + vx(i)
     sumvy = sumvy + vy(i)
     sumvz = sumvz + vz(i)
     sumv2 = sumv2 + (vx(i)**2 + vy(i)**2 + vz(i)**2)
  end do
  ! velocity of the center of mass (Not mass is taken as unity, 1)
   sumvx = sumvx / npart
   sumvy = sumvy / npart
   sumvz = sumvz / npart
   ! mean-squared velocity..
   sumv2 = sumv2 / npart
   ! Scale factor of the velocity  (in reduced units, K_B = 1)
   fs = sqrt(3*Tb/sumv2)

  ! Set desired kinetic energy and set velocity center of mass to zero..
  do i = 1, npart
     vx(i) = (vx(i) - sumvx)*fs
     vy(i) = (vy(i) - sumvy)*fs
     vz(i) = (vz(i) - sumvz)*fs
  end do
 
  ! Q. Check that the instantaneous temp is (approx.) equal to Tb 
  !  and average velocity is zero.. 

  return
end subroutine gen_vel

!     ********************************************************
! setup the initial coordinate..

subroutine fcc
  
  use parameters
  implicit none
  
  integer :: i, ix, iy, iz, nuc, m, iref
  
  double precision, dimension (npart) :: rx, ry, rz
  double precision :: ucl, uclh
  
  common / coords / rx, ry, rz
  
  !     number of unit cells
  i = 1
  do while (4*i**3 .lt. npart)
     i = i + 1
  end do
  nuc = i
  
  !     unit cell length
  ucl    = 1.0d0 / nuc
  uclh   = 0.5 * ucl
  
  !     build the unit cell
  rx(1) =  0.d0
  ry(1) =  0.d0
  rz(1) =  0.d0
  
  rx(2) =  uclh
  ry(2) =  uclh
  rz(2) =  0.d0

  rx(3) =  0.d0
  ry(3) =  uclh
  rz(3) =  uclh

  rx(4) =  uclh
  ry(4) =  0.d0
  rz(4) =  uclh

  m = 0
  do iz = 1, nuc
     do iy = 1, nuc
        do ix = 1, nuc
           do iref = 1, 4
              rx(iref+m) = rx(iref) + ucl * (ix - 1)
              ry(iref+m) = ry(iref) + ucl * (iy - 1)
              rz(iref+m) = rz(iref) + ucl * (iz - 1)
           end do
           m = m + 4
        end do
     end do
  end do
  
  return
end subroutine fcc

! ==========================================================
! function to generate random numbers between 0 and 1
! Q. Read about random number generators.. 

double precision function ranf(dummy)

  integer   l, c, m
  parameter (l = 1029, c = 221591, m = 1048576 )
  integer     seed
  double precision dummy
  save        seed
  data        seed / 0 /

  seed = mod(seed * l + c, m)
  ranf = dfloat(seed) / dfloat(m)
  
  return
end function ranf

!     ---------------------  end  --------------------
