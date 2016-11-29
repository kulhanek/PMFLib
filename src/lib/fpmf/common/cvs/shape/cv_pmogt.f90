!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
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

! Principal Moments Of Gyration Tensor - PMOGT
! gyration tensor is mass weighted

module cv_pmogt

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePMOGT

    integer :: axis

    contains
        procedure :: load_cv        => load_pmogt
        procedure :: calculate_cv   => calculate_pmogt
end type CVTypePMOGT

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_pmogt
!===============================================================================

subroutine load_pmogt(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePMOGT)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(1)                        :: cdir
    ! -----------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'PMOGT'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1
    call cv_common_read_groups(cv_item,prm_fin)

    ! rest of controls ------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'axis',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'axis is not specified!')
    end if

    write(PMF_OUT,60) cdir

    select case(cdir)
        case('x')
            cv_item%axis = 1
        case('y')
            cv_item%axis = 2
        case('z')
            cv_item%axis = 3
        case default
            call pmf_utils_exit(PMF_OUT,1,'axis has to be x, y, or z!')
    end select

    return

60 format('   ** axis               : ',A)

end subroutine load_pmogt

!===============================================================================
! Subroutine:  calculate_pmogt
!===============================================================================

subroutine calculate_pmogt(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypePMOGT)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: ai,m,info
    real(PMFDP)    :: com(3),amass,totmass,itotmass,fac,dv
    real(PMFDP)    :: t11,t22,t33,t12,t13,t23,t(3,3),work(26*3),eigenvalues(3)
    ! --------------------------------------------------------------------------

    ! calculate centre of mases --------------------
    totmass = 0.0d0
    com(:) = 0.0
    do  m = 1, cv_item%natoms
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        com(:) = com(:) + x(:,ai)*amass
        totmass = totmass + amass
    end do
    if( totmass .gt. 0 ) then
        itotmass = 1.0d0 / totmass
        com(:) = com(:) * itotmass
    else
        itotmass = 0.0d0
    end if

    ! calculate tensor of gyration
    t11 = 0.0d0
    t22 = 0.0d0
    t33 = 0.0d0
    t12 = 0.0d0
    t13 = 0.0d0
    t23 = 0.0d0

    ! construct matrix
    do m = 1, cv_item%natoms
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        t11 = t11 + amass*(x(1,ai) - com(1))*(x(1,ai) - com(1))
        t22 = t22 + amass*(x(2,ai) - com(2))*(x(2,ai) - com(2))
        t33 = t33 + amass*(x(3,ai) - com(3))*(x(3,ai) - com(3))
        t12 = t12 + amass*(x(1,ai) - com(1))*(x(2,ai) - com(2))
        t13 = t13 + amass*(x(1,ai) - com(1))*(x(3,ai) - com(3))
        t23 = t23 + amass*(x(2,ai) - com(2))*(x(3,ai) - com(3))
    end do

    t(1,1) = t11 * itotmass
    t(2,2) = t22 * itotmass
    t(3,3) = t33 * itotmass

    t(1,2) = t12 * itotmass
    t(1,3) = t13 * itotmass
    t(2,3) = t23 * itotmass

    t(2,1) = t12 * itotmass
    t(3,1) = t13 * itotmass
    t(3,2) = t23 * itotmass

    ! calculate eignevalues and eigenvectors of matrix
    eigenvalues(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 3, t, 3, eigenvalues, work, 26*3, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize T matrix in calculate_pmogt!')
    end if

    ctx%CVsValues(cv_item%idx) = sqrt(eigenvalues(cv_item%axis))

    dv = 1.0 / ctx%CVsValues(cv_item%idx)

    ! calculate derivatives -------------------------
    do m = 1, cv_item%natoms
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        fac = dv * amass * itotmass
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                                + fac*((x(1,ai) - com(1))*t(1,cv_item%axis)*t(1,cv_item%axis) &
                                +  (x(2,ai) - com(2))*t(1,cv_item%axis)*t(2,cv_item%axis) &
                                +  (x(3,ai) - com(3))*t(1,cv_item%axis)*t(3,cv_item%axis) )

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                                + fac*((x(1,ai) - com(1))*t(2,cv_item%axis)*t(1,cv_item%axis) &
                                +  (x(2,ai) - com(2))*t(2,cv_item%axis)*t(2,cv_item%axis) &
                                +  (x(3,ai) - com(3))*t(2,cv_item%axis)*t(3,cv_item%axis) )

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                                + fac*((x(1,ai) - com(1))*t(3,cv_item%axis)*t(1,cv_item%axis) &
                                +  (x(2,ai) - com(2))*t(3,cv_item%axis)*t(2,cv_item%axis) &
                                +  (x(3,ai) - com(3))*t(3,cv_item%axis)*t(3,cv_item%axis) )
    end do

    return

end subroutine calculate_pmogt

!===============================================================================

end module cv_pmogt

