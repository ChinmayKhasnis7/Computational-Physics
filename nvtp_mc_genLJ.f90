!---------------------------------------------------------------
! Chinmay Khasanis
! Monte-Carlo (MC) for NVT ensemble to calculate Pressure
! 500 particles in Genetralised Lennard-Jones potential
!---------------------------------------------------------------

module parameters_q2b      
  implicit none      
  save 
  integer, parameter :: npart = 500       
  integer, parameter :: nstep = npart*5000   
  double precision, parameter :: temp = 1.0d0  
  double precision, parameter :: invt = 1.0d0/temp   !Boltzmann factor in reduced unit
  double precision, parameter :: dens = 0.80d0    !number density
  integer, parameter :: thermo_dump = 10*npart 
  integer, parameter :: traj_dump = 10*npart 
  double precision, parameter :: vol  = npart*1.0/dens
  double precision, parameter :: boxl = vol**(1.0d0/3)
  double precision, parameter :: rc = 1.50d0, rc2 = rc*rc !npot potential cut-off
  double precision, parameter :: alpha = 50.0d0
  double precision, parameter :: disp_p = 1.50d-01  !maximum displacement of a particle
end module parameters_q2b

!     ------------------------------------------------------------------

program mc_nvtp 
  use parameters_q2b
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

  open (unit = 102, file = "pe.out", action="write")
  open (unit = 103, file = "traj.xyz", action="write") 
  open (unit = 104, file = "pressure.dat", action="write")
  call setup
  
  move_attempted = 0
  moves_accepted = 0

  ! Propagate the initial configuration in the configuration space.. 
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
      if (mod(istep, thermo_dump) .eq. 0) write (*,5) mc_step, pe/npart, &
         npart/vol, moves_accepted*1.0/move_attempted
      if (mod(istep, thermo_dump) .eq. 0) then
         write(102,*) mc_step, pe
         flush(102)
      end if

      if ((mc_step .gt. 500) .and. (mod(istep,traj_dump) .eq. 0)) then
         do ipart = 1, npart
            write(103,*) rx(ipart), ry(ipart), rz(ipart)
            flush(103)
         end do
      end if
   
   end do
   close(103)
   close(104)

5 format(i8,2x,3(f15.6,2x))
      
end program mc_nvtp
  
subroutine setup
  use parameters_q2b
  implicit none
  integer :: ipart, jpart
  double precision :: npot
  double precision  :: rx(npart), ry(npart), rz(npart)
  double precision  :: pe, virial
  double precision, dimension (npart) :: fx, fy, fz
  external  npot
  
  common / coords / rx, ry, rz
  common / cnfg_dat2  / pe
  common / forces / fx, fy, fz
  
  !  generate inital configuration box length between 0 and 1
  call fcc
  do ipart = 1, npart
     rx(ipart) = rx(ipart) * boxl
     ry(ipart) = ry(ipart) * boxl
     rz(ipart) = rz(ipart) * boxl
  end do

  ! compute initial potential energy
  pe = 0.0d0
  do ipart = 1, npart-1
     do jpart = ipart+1, npart
        pe = pe + npot(rx(ipart),ry(ipart),rz(ipart),rx(jpart),ry(jpart),rz(jpart))
     end do
  end do
    ! Calculate potential energy and force 
  call comp_force(rx, ry, rz, fx, fy, fz, pe, virial)
  
  write (*,*) "initial potential energy = ", pe 
  return
end subroutine setup

subroutine comp_force(rx, ry, rz, fx, fy, fz, pe, virial)
   use parameters_q2b
   implicit none
   integer :: ipart, jpart
   double precision, dimension (npart) :: rx, ry, rz
   double precision, dimension (npart) :: fx, fy, fz
   double precision :: dx, dy, dz, r2, r2p , r6i, f, pe, virial
 
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
         ! Periodic boundary condition
         dx = dx - dnint(dx/boxl)*boxl
         dy = dy - dnint(dy/boxl)*boxl
         dz = dz - dnint(dz/boxl)*boxl
         r2 = dx*dx + dy*dy + dz*dz
         ! potential cut-off
         if (r2 .lt. rc2) then   
            r2p = 1.0 / (r2-1)
            r6i = r2p**3
            pe  = pe + (4.0/(alpha**2))*(r6i*r6i - alpha*r6i)  ! given potential
            ! force
            f   = (48*sqrt(r2)/(alpha**2))*((r2p**7)-(0.5*alpha*r2p**4))
            fx(ipart) = fx(ipart) + dx*f/r2 
            fx(jpart) = fx(jpart) - dx*f/r2
            fy(ipart) = fy(ipart) + dy*f/r2
            fy(jpart) = fy(jpart) - dy*f/r2 
            fz(ipart) = fz(ipart) + dz*f/r2
            fz(jpart) = fz(jpart) - dz*f/r2
            !------------virial -------------------------
            virial = virial +  (48*sqrt(r2)/(alpha**2))*((r2p**7)-(0.5*alpha*r2p**4))
         end if
      end do
   end do
   virial = virial * 1.0d0/3.0d0   
  return
  
 end subroutine comp_force

!----Subroutine to randomly displace a particle and calcualte if the move is accepted

subroutine displacement(accept)
  use parameters_q2b
  implicit none
  integer :: p, i
  double precision :: rx_o, ry_o, rz_o, rx_n, ry_n, rz_n
  double precision :: pe_o, pe_n, de, npot, prob
  double precision :: ranf, dummy
  double precision :: dx, dy, dz
  double precision :: rx(npart), ry(npart), rz(npart)
  double precision, dimension (npart) :: fx, fy, fz
  double precision :: pe
  logical :: accept
  external ranf, npot
  
  common / coords / rx, ry, rz
  common / cnfg_dat2 / pe
  common / forces / fx, fy, fz
  
  !  select randomly a particle
  p = int(ranf(dummy)*npart) + 1
  
  ! save the old position
  rx_o = rx(p)
  ry_o = ry(p)
  rz_o = rz(p)
 
  ! generate a new trial position..
  dx = disp_p * (2*ranf(dummy)-1.0d0)
  dy = disp_p * (2*ranf(dummy)-1.0d0)
  dz = disp_p * (2*ranf(dummy)-1.0d0)
  rx_n = rx_o + dx
  ry_n = ry_o + dy
  rz_n = rz_o + dz

  if (rx_n .lt. 0.0) rx_n = rx_n + boxl
  if (ry_n .lt. 0.0) ry_n = ry_n + boxl
  if (rz_n .lt. 0.0) rz_n = rz_n + boxl
  if (rx_n .gt. boxl) rx_n = rx_n - boxl
  if (ry_n .gt. boxl) ry_n = ry_n - boxl
  if (rz_n .gt. boxl) rz_n = rz_n - boxl
  
  pe_o = 0.0d0
  pe_n = 0.0d0
  do i = 1, npart
     if (i .ne. p) then
        pe_o = pe_o + npot(rx_o,ry_o,rz_o,rx(i),ry(i),rz(i))
        pe_n = pe_n + npot(rx_n,ry_n,rz_n,rx(i),ry(i),rz(i))
     end if
  end do
  
  de      = pe_n - pe_o
  prob    = min(1.0d0, exp(-invt*de))
  accept  = .true.
  if (ranf(dummy) .gt. prob) accept = .false.
  if (.not. accept) then    ! revert to the old configuration
     return
  end if
  !if accepted then update the configuration data
  pe  = pe + de
  rx(p) = rx_n
  ry(p) = ry_n
  rz(p) = rz_n
  return
end subroutine displacement
      
!-------------Function to calculate given interaction potential energy--------------------

function npot(rx1, ry1, rz1, rx2, ry2, rz2)
  use parameters_q2b
  implicit none
  double precision :: rx1, ry1, rz1, rx2, ry2, rz2, npot
  double precision :: dx, dy, dz, r2, r6i
  
  dx = rx1 - rx2
  dy = ry1 - ry2
  dz = rz1 - rz2
  ! periodic boundary condition 
  dx = dx - dnint(dx/boxl)*boxl
  dy = dy - dnint(dy/boxl)*boxl
  dz = dz - dnint(dz/boxl)*boxl
  r2 = (dx*dx) + (dy*dy) + (dz*dz)
  !potential cutoff
  if (r2 .lt. rc2) then
     r6i = 1.0 / (r2-1)**3
     npot  = (4.0/(alpha*alpha))*(r6i*r6i - alpha*r6i) 
  else 
     npot = 0.0d0     
  end if

  return
end function npot

!  FCC solid init config

subroutine fcc
  use parameters_q2b
  implicit none
  integer :: i, ix, iy, iz, nuc, m, iref
  double precision, dimension (npart) :: rx, ry, rz
  double precision :: ucl, uclh

  common / coords / rx, ry, rz
  i = 1
  do while (4*i**3 .lt. npart)
     i = i + 1
  end do
  nuc = i

  !unit cell length
  ucl    = 1.0d0 / nuc
  uclh   = 0.5 * ucl

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

! function to generate random numbers between 0 and 1

double precision function ranf(dummy)
  integer   l, c, m
  parameter (l = 1029, c = 221591, m = 1048576 )
  integer     seed
  double precision dummy
  save        seed
  data        seed / 0 /
  dummy = 0
  seed = mod(seed * l + c, m)
  ranf = dfloat(seed) / dfloat(m)
  return
end function ranf
