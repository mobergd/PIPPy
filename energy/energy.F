      program energy
c     DiNT:  Adiabatix and Nonadiabatic Trajectories
c     Main program.  This routine calls HEADER and READIN once, then
c     calls PREMOL once per AG to precalculate some constants.
c     DRIVER is called NTRAJ times, once for each trajectory.
c     FINALSTATE is called after each trajectory to write some data.

      implicit none
      include 'param.f'
      include 'c_sys.f'
!      include 'c_ran.f'
      include 'c_traj.f'
!#include <sprng_f.h>

      integer im,initrng,minitrng,i,ia
      double precision t1,t2,t3,cputime,t3old,r

      call header

      write(6,'("...entering READIN...",/)')
      call readin
      write(6,'("...exiting READIN...",/)')

      call timing(t1)
      write(6,'("...entering PREMOL...",/)')
c     prepare each molecule
c     compute quantities that are the same for all trajectories
c     initialize electronic coordinates
      do im = 1, nmol
        call premol(im)
      enddo
      write(6,'("...exiting PREMOL...",/)')
      call timing(t2)
      write(6,'(" CPU time in PREMOL is ",f10.5," s",/)')(t2-t1)

      t3 = 0.d0
      cputime = 0.d0

      write(6,'("...entering INITMOL...",/)')
c     generate specific initial coordinates for each molecule
      do im = 1, nmol
        call initmol(im)
      enddo
      write(6,'("...exiting INITMOL...",/)')

      write(6,'("...entering DRIVER...",/)')
      call driver
      write(6,'("...exiting DRIVER...",/)')
      t3old = t3
      call timing(t3)
      cputime = cputime + (t3-t3old)
      write(6,'(" CPU time for this configuration is ",f10.5," s",/)')
     & (t3-t3old)

!      write(6,'(" CPU time for all configurations is ",f10.5," s or ",
!     & f10.5," s per trajectory.",/)')
!     & cputime,cputime/dble(nconf)

 999  write(6,'("Well done!")')

      write(6,'("PE of each AG and total")')
      write(6,*)(peagi(i)*autoev,i=1,nmol)
      write(6,*) pei*autoev

      end
