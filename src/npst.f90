!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine: NPST
! Language:   Fortran 90/95
! This code can be compiled to a shared library (DLL on Windows) and
! then be interfaced via R (www.r-project.org) to derive the empirical
! distribution of a generalization of Hewitt's and Rogerson's
! nonparametric seasonality tests.
!
! This code has been successfully compiled  on Windows XP, GNU/Linux (Ubuntu) and FreeBSD
! using 'gfortran' from the GNU Compiler Collection (gcc)
!
! Copyright (C) 2007-2010  Roland Rau
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine npst(lng, pk, rpts, mxrksum, resvec, basedata, basedata2)
  implicit none
  integer               :: lng, pk, rpts, mxrksum
  integer               :: resvec(mxrksum)
  real                  :: x
  integer               :: basedata(lng)
  integer               :: basedata2(lng*2)
  integer               :: tempvalue, posis, maxconsecsum, temp

  integer               :: i,j

  integer               :: ii, n, clock

  integer, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  
  CALL SYSTEM_CLOCK(COUNT=clock)
  
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
  
  DEALLOCATE(seed)

  i = 0
 
  ! the main loop
  do j=1,rpts,1
     call random_number(x)
  
     ! generating a random permutation of 'basedata' and 'basedata2'
     ! 
     do i = lng, 1, -1
        posis = floor(1 + x*i)
        tempvalue = basedata(posis)
        basedata(posis)=basedata(i)
        basedata(i) = tempvalue
        basedata2(i) = basedata(i)
        basedata2(i+lng) = basedata(i)
     end do
     
    
     maxconsecsum = sum(basedata2(1:(1+pk-1)))
     temp =  sum(basedata2(1:(1+pk-1)))
      do i =2, lng, 1
         temp = temp - basedata2(i-1) + basedata2(i+pk-1)
         if (temp > maxconsecsum) maxconsecsum = temp
      end do
     
      resvec(maxconsecsum) = resvec(maxconsecsum) + 1

   end do
   !  deallocate(iseed)
end subroutine npst
