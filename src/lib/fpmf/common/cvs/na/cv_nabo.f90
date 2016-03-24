!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2016 Ivo Durnik, 
!    Copyright (C) 2016 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, 
!    Boston, MA  02110-1301  USA
!===============================================================================

module cv_nabo

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeNABO
    contains
        procedure :: load_cv        => load_nabo
        procedure :: calculate_cv   => calculate_nabo
end type CVTypeNABO

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_nabo
!===============================================================================

subroutine load_nabo(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common

    implicit none
    class(CVTypeNABO)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'NABO'
    cv_item%unit          = AngleUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 6
    call cv_common_init_groups(cv_item,prm_fin)

    ! read group a,b ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)
    call cv_common_read_group(cv_item,prm_fin,2)

    ! read group c,d ----------------------------------
    write(PMF_OUT,60)
    call cv_common_read_group(cv_item,prm_fin,3)
    call cv_common_read_group(cv_item,prm_fin,4)

    ! read group e,f ----------------------------------
    write(PMF_OUT,70)
    call cv_common_read_group(cv_item,prm_fin,5)
    call cv_common_read_group(cv_item,prm_fin,6)

    return

50 format('   == Helical axis ===============================')
60 format('   == Base #1 direction ==========================')
70 format('   == Base #2 direction ==========================')

end subroutine load_nabo

!===============================================================================
! Subroutine:  calculate_nabo
!===============================================================================

subroutine calculate_nabo(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNABO)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: ai,m
    real(PMFDP)         :: d1(3),d2(3),d3(3),d4(3),d5(3),d6(3),h(3)
    real(PMFDP)         :: v1(3),v2(3),v1p(3),v2p(3),hp(3)
    real(PMFDP)         :: totmass1,totmass2,totmass3,totmass4,totmass5,totmass6,amass
    real(PMFDP)         :: h2,t1,t2,arg,v1p2,v2p2,v1pv,v2pv,scal
    ! --------------------------------------------------------------------------

    ! calculate actual value
    totmass1 = 0.0d0
    d1(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d1(:) = d1(:) + x(:,ai)*amass
        totmass1 = totmass1 + amass
    end do
    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_nabo!')
    end if
    d1(:) = d1(:) / totmass1
 
    totmass2 = 0.0d0
    d2(:) = 0.0d0
    do  m = cv_item%grps(1) + 1 , cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d2(:) = d2(:) + x(:,ai)*amass
        totmass2 = totmass2 + amass
    end do
    if( totmass2 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass2 is zero in calculate_nabo!')
    end if
    d2(:) = d2(:) / totmass2

    totmass3 = 0.0d0
    d3(:) = 0.0d0
    do  m = cv_item%grps(2) + 1 , cv_item%grps(3)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d3(:) = d3(:) + x(:,ai)*amass
        totmass3 = totmass3 + amass
    end do
    if( totmass3 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass3 is zero in calculate_nabo!')
    end if
    d3(:) = d3(:) / totmass3

    totmass4 = 0.0d0
    d4(:) = 0.0d0
    do  m = cv_item%grps(3) + 1 , cv_item%grps(4)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d4(:) = d4(:) + x(:,ai)*amass
        totmass4 = totmass4 + amass
    end do
    if( totmass4 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass4 is zero in calculate_nabo!')
    end if
    d4(:) = d4(:) / totmass4

    totmass5 = 0.0d0
    d5(:) = 0.0d0
    do  m = cv_item%grps(4) + 1 , cv_item%grps(5)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d5(:) = d5(:) + x(:,ai)*amass
        totmass5 = totmass5 + amass
    end do
    if( totmass5 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass5 is zero in calculate_nabo!')
    end if
    d5(:) = d5(:) / totmass5

    totmass6 = 0.0d0
    d6(:) = 0.0d0
    do  m = cv_item%grps(5) + 1 , cv_item%grps(6)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d6(:) = d6(:) + x(:,ai)*amass
        totmass6 = totmass6 + amass
    end do
    if( totmass6 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass6 is zero in calculate_nabo!')
    end if
    d6(:) = d6(:) / totmass6

    h(:) = d1(:) - d2(:)

    v1(:) = d3(:) - d4(:)
    v2(:) = d5(:) - d6(:)

    h2 = h(1)**2 + h(2)**2 + h(3)**2

    t1 = ( v1(1)*h(1) + v1(2)*h(2) + v1(3)*h(3) ) / h2
    v1p(:) = v1(:) - t1*h(:)

    t2 = ( v2(1)*h(1) + v2(2)*h(2) + v2(3)*h(3) ) / h2
    v2p(:) = v2(:) - t2*h(:)

    v1p2 = v1p(1)**2 + v1p(2)**2 + v1p(3)**2
    v2p2 = v2p(1)**2 + v2p(2)**2 + v2p(3)**2

    v1pv = sqrt(v1p2);
    v2pv = sqrt(v2p2);

    arg = (v1p(1)*v2p(1) + v1p(2)*v2p(2) + v1p(3)*v2p(3)) / (v1pv * v2pv)

    hp(1) = v1p(2)*v2p(3) - v2p(2)*v1p(3)
    hp(2) = v1p(3)*v2p(1) - v2p(3)*v1p(1)
    hp(3) = v1p(1)*v2p(2) - v2p(1)*v1p(2)

    scal = h(1)*hp(1) + h(2)*hp(2) + h(3)*hp(3)

    ! coordinate ctx%CVsValues(cv_item%idx) ------------------------------------------------------------

    if ( arg .gt.  1.0 ) then
        arg =  1.0
    else if ( arg .lt. -1.0 ) then
        arg = -1.0
    end if

    ctx%CVsValues(cv_item%idx) = sign(1.0d0,scal)*acos( arg )
 
!     ctx%CVsValues(cv_item%idx) = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
! 
!     ! ------------------------------------------------
! 
!     if( ctx%CVsValues(cv_item%idx) .gt. 1e-7 ) then
!         sc = 1.0d0 / ctx%CVsValues(cv_item%idx)
!     else
!         sc = 0.0d0
!     end if
! 
!     ! warning - groups can overlap - it is therefore important to add gradients
! 
!     do  m = 1, cv_item%grps(1)
!         ai = cv_item%lindexes(m)
!         amass = mass(ai)
!         ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + sc*dx(:)*amass/totmass1
!     end do
! 
!     do  m = cv_item%grps(1) + 1 , cv_item%grps(2)
!         ai = cv_item%lindexes(m)
!         amass = mass(ai)
!         ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - sc*dx(:)*amass/totmass2
!     end do

 return

end subroutine calculate_nabo

!===============================================================================

end module cv_nabo
