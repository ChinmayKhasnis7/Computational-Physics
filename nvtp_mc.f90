!---------------------------------------------------------------
! Chinmay Khasanis
! Monte-Carlo samples in canonical ensemble to calculate pressure
! 500 Lennard-Jones particles
!---------------------------------------------------------------

module parameters      
  implicit none      
  save 
  
  integer, parameter :: npart = 500           ! number of partciles
  integer, parameter :: nstep = npart*50000   ! total no. of MC steps = nstep / npart

  double precision, parameter :: temp = 1.0d0  ! bath temperature.
  double precision, parameter :: invt = 1.0d0/temp   ! Boltzmann factor in reduced unit
  double precision, parameter :: dens = 0.80d0    ! number density

  integer, parameter :: thermo_dump = 10*npart ! dump thermodynamic variables
  integer, parameter :: traj_dump = 100*npart  ! dump trajectory

  !  calculate the box length
  double precision, parameter :: vol  = npart*1.0/dens
  double precision, parameter :: boxl = vol**(1.0d0/3)

  ! lj potential cut-off
  double precision, parameter :: rc = 2.50d0, rc2 = rc*rc

  ! MC 
  double precision, parameter :: disp_p = 1.50d-01     ! maximum disp of a particle
  ! Adjust this maximum disp to get the acceptance between 0.3 to 0.4 (30 to 40 %). 

end module parameters

!     ------------------------------------------------------------------

program mc_basics 

  use parameters
  implicit none
  
  integer :: istep,ipart
  integer :: moves_accepted, move_attempted
  integer :: mc_step
  double precision, dimension (npart) :: fx, fy, fz
  double precision :: rx(npart), ry(npart), rz(npart)
  double precision :: pe , virial , pressure
  
  logical  accept
  
  common / coords  / rx, ry, rz
  common / cnfg_dat2  / pe
  common / forces / fx, fy, fz

 
   ! Open files to write output...
 ! open (unit = 101, file = "thermo.out", action="write")
  open (unit = 102, file = "pe.out", action="write")
  open (unit = 103, file = "traj.xyz", action="write") 
  open (unit = 104, file = "pressure.dat", action="write")

!   pe    = 0.0d0
!   fx(:) = 0.0d0
!   fy(:) = 0.0d0
!   fz(:) = 0.0d0
!   virial = 0.0d0    
  call setup
  
  !     equilibrium part
  
  write (*,*) "# equilibriation start #"
  
  move_attempted = 0
  moves_accepted = 0

  ! Propagate the initial configuratio.. in the configuration space.. 
   do istep = 1, nstep
      mc_step = istep/npart
      call displacement(accept)

        ! Calculate potential energy and force ... 
      call comp_force(rx, ry, rz, fx, fy, fz, pe, virial)
      pressure = (dens*temp) + (virial)/vol

      if (mod(istep, thermo_dump) .eq. 0) then
         write(104,*) pressure
         flush(104)
         print*, pressure , virial/vol
      end if

      move_attempted = move_attempted + 1
      if (accept) moves_accepted = moves_accepted + 1
      ! write mc-step, volume, energy, density, probabilty of accepted moves
      if (mod(istep, thermo_dump) .eq. 0) write (*,5) mc_step, pe/npart, &
         npart/vol, moves_accepted*1.0/move_attempted
      ! write mc-step, and energy/particle 
      if (mod(istep, thermo_dump) .eq. 0) write(102,*) mc_step, pe/npart

      
      ! Print trajectories here....
      if ((mc_step .gt. 10000) .and. (mod(istep,traj_dump) .eq. 0)) then
         write(103,*) npart
         write(103,*)
         do ipart = 1, npart
          write(103,*) "O", rx(ipart), ry(ipart), rz(ipart)
       end do
    end if
   
   end do
   close(103)
   close(104)

5 format(i8,2x,3(f15.6,2x))
      
end program mc_basics
  
!     ****************************************************************
!     Setup the positions and other simulation parameters

subroutine setup

  use parameters
  implicit none
  
  integer :: ipart, jpart
  
  double precision :: lj
  double precision  :: rx(npart), ry(npart), rz(npart)
  double precision  :: pe, virial
  double precision, dimension (npart) :: fx, fy, fz
  external  lj
  
  common / coords / rx, ry, rz
  common / cnfg_dat2  / pe
  common / forces / fx, fy, fz
  
  !  generate inital configuration box length between 0 and 1

  call fcc

  ! scale coordinates with the box length
  do ipart = 1, npart
     rx(ipart) = rx(ipart) * boxl
     ry(ipart) = ry(ipart) * boxl
     rz(ipart) = rz(ipart) * boxl
  end do

  ! compute initial potential energy..
  pe = 0.0d0
  do ipart = 1, npart-1
     do jpart = ipart+1, npart
        pe = pe + lj(rx(ipart),ry(ipart),rz(ipart),rx(jpart),ry(jpart),rz(jpart))
     end do
  end do
    ! Calculate potential energy and force ... 
  call comp_force(rx, ry, rz, fx, fy, fz, pe, virial)
  
  write (*,*) "initial potential energy = ", pe 
  
  return
end subroutine setup

subroutine comp_force(rx, ry, rz, fx, fy, fz, pe, virial)
  
   use parameters
   implicit none
   
   integer :: ipart, jpart
   
   double precision, dimension (npart) :: rx, ry, rz
   double precision, dimension (npart) :: fx, fy, fz
   double precision :: dx, dy, dz, r2, r6i, f, pe, virial
 
   pe    = 0.0d0
   fx(:) = 0.0d0
   fy(:) = 0.0d0
   fz(:) = 0.0d0
   virial = 0.0d0
 
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
            !adding the virial computation 
            virial = virial + 48*(r6i*r6i - 0.5*r6i)
         end if
      end do
      
   end do
   virial = virial * 1.0d0/3.0d0
   ! print*, virial
   
  return
  
 end subroutine comp_force

!     ********************************************************************
!     Subroutine to randomly displace a particle and calcualte if the move is accepted

subroutine displacement(accept)

  use parameters
  implicit none
  
  integer :: p, i
  
  double precision :: rx_o, ry_o, rz_o, rx_n, ry_n, rz_n
  double precision :: pe_o, pe_n, de, lj, prob
  double precision :: ranf, dummy
  double precision :: dx, dy, dz
  double precision :: rx  double precision :: rx(npart), ry(npart), rz(npart)
  (npart), ry(npart), rz(npart)  double precision :: rx  double precision :: rx(npart), ry(npart), rz(npart)

  double precision, dimension (npart) :: fx, fy, fz

  double precision :: pe

  logical :: accept
  
  external ranf, lj
  
  common / coords / rx, ry, rz
  common / cnfg_dat2 / pe
  common / forces / fx, fy, fz
  
  !  select randomly a particle
  p = int(ranf(dummy)*npart) + 1
  
  ! save the old position.. 
  rx_o = rx(p)
  ry_o = ry(p)
  rz_o = rz(p)
 
  ! get new position  
  ! generate a new trial position..
  dx = disp_p * (2*ranf(dummy)-1.0d0)
  dy = disp_p * (2*ranf(dummy)-1.0d0)
  dz = disp_p * (2*ranf(dummy)-1.0d0)
  rx_n = rx_o + dx
  ry_n = ry_o + dy
  rz_n = rz_o + dz
  ! ensure that particles are inside the simulation box.. )
  if (rx_n .lt. 0.0) rx_n = rx_n + boxl
  if (ry_n .lt. 0.0) ry_n = ry_n + boxl
  if (rz_n .lt. 0.0) rz_n = rz_n + boxl
  if (rx_n .gt. boxl) rx_n = rx_n - boxl
  if (ry_n .gt. boxl) ry_n = ry_n - boxl
  if (rz_n .gt. boxl) rz_n = rz_n - boxl
  
  !     calculate the old and the new "particle" energies
  pe_o = 0.0d0
  pe_n = 0.0d0
  do i = 1, npart
     if (i .ne. p) then
        pe_o = pe_o + lj(rx_o,ry_o,rz_o,rx(i),ry(i),rz(i))
        pe_n = pe_n + lj(rx_n,ry_n,rz_n,rx(i),ry(i),rz(i))
     end if
  end do
  
  de      = pe_n - pe_o
  prob    = min(1.0d0, exp(-invt*de))

  accept  = .true.
  
  if (ranf(dummy) .gt. prob) accept = .false.
  
  if (.not. accept) then    ! revert to the old configuration
     return
  end if
  
  !     if accepted then update the configuration data
  pe  = pe + de
  rx(p) = rx_n
  ry(p) = ry_n
  rz(p) = rz_n
  

  return
end subroutine displacement
      
!     ******************************************************************
!     Function to calculate lj potential energy 

function lj(rx1, ry1, rz1, rx2, ry2, rz2)
  
  use parameters
  implicit none
  
  double precision :: rx1, ry1, rz1, rx2, ry2, rz2, lj
  double precision :: dx, dy, dz, r2, r6i
  
  dx = rx1 - rx2
  dy = ry1 - ry2
  dz = rz1 - rz2
  ! periodic boundary condition.. 
  dx = dx - dnint(dx/boxl)*boxl
  dy = dy - dnint(dy/boxl)*boxl
  dz = dz - dnint(dz/boxl)*boxl
  r2 = (dx*dx) + (dy*dy) + (dz*dz)
  ! use the potential cut-off..
  if (r2 .lt. rc2) then
     r6i = 1.0 / (r2*r2*r2)
     ! pe computation
     lj  = 4*(r6i*r6i - r6i)  ! lennard-jones potential
  else 
     lj = 0.0d0     
  end if

  return
end function lj

!     ********************************************************
! Subroutine to generate the intial configuration as FCC solid

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

!     ********************************************************
! function to generate random numbers between 0 and 1

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
