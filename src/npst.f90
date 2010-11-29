!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine: NPST
! Language:   Fortran 90/95
!
! This subroutine estimates via Monte-Carlo simulation the empirical
! cumulative distribution function for a generalized version of
! Hewitt's et al. (1971) and Rogerson's (1996) nonparametric
! seasonality tests. Please see the associated documentation in
! package 'npst' for further details.
!
! Copyright (C) 2007-2010  Roland Rau (roland.rau@gmail.com)
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

subroutine npst(lng, pk, rpts, mxrksum, resvec, basedata, basedata2, seed)
  implicit none
  integer, intent(in)   :: lng, pk, rpts, mxrksum
  integer, intent(out)  :: resvec(mxrksum)
  real                  :: x
  integer, intent(in)   :: seed(8)
  integer               :: basedata(lng)
  integer               :: basedata2(lng*2)
  integer               :: tempvalue, posis, maxconsecsum, temp
  integer               :: i,j

  call random_seed(PUT=seed)
 
  ! the main loop
  do j=1,rpts,1
     call random_number(x)
     
     ! generating a random permutation of 'basedata' and 'basedata2'
     ! (a "shuffle")
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
     ! little trick (?) coming up here: instead of summing, for
     ! instance from 2 to 7, 3 to 8, ... --- we take the previous sum
     ! (e.g. from 1 to 6), subtract element 1 and add element 7. I
     ! played around with both approaches and this one gave a
     ! significant speed boost (at least 10% faster now)
     do i = 2, lng, 1
        temp = temp - basedata2(i-1) + basedata2(i+pk-1)
        if (temp > maxconsecsum) maxconsecsum = temp
     end do
     
     resvec(maxconsecsum) = resvec(maxconsecsum) + 1
     
  end do

end subroutine npst
