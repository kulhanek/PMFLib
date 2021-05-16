!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cst_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_core_main_lf
!===============================================================================

subroutine cst_core_main_lf

    use cst_constraints
    use cst_lambdas
    use cst_output
    use cst_restart

    implicit none
    ! --------------------------------------------------------------------------

    call cst_constraints_increment
    call cst_lambdas_calculate
    call cst_core_analyze
    call cst_core_enthalpy_entropy_lf
    call cst_output_write
    call cst_restart_update
    call cst_restart_trajectory_write_snapshot

end subroutine cst_core_main_lf

!===============================================================================
! Subroutine:  cst_core_main_vv_shake
!===============================================================================

subroutine cst_core_main_vv_shake

    use cst_lambdas
    use cst_output
    use cst_velocities

    implicit none
    ! --------------------------------------------------------------------------

    call cst_lambdas_calculate
    call cst_velocities_correct_a
    call cst_core_analyze
    call cst_output_write

end subroutine cst_core_main_vv_shake

!===============================================================================
! Subroutine:  cst_core_main_vv_rattle
!===============================================================================

subroutine cst_core_main_vv_rattle()

    use cst_constraints
    use cst_velocities

    implicit none
    ! --------------------------------------------------------------------------

    call cst_constraints_increment
    call cst_velocities_correct_b

end subroutine cst_core_main_vv_rattle

!===============================================================================
! Subroutine:  cst_core_analyze
!===============================================================================

subroutine cst_core_analyze

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,ci,j,cj,k,info
    real(PMFDP)            :: fzv,isrz,lam,mu,dval,dval2,invn
    ! --------------------------------------------------------------------------

    ! reset accumulators ---------------------------------------------------
    if ( faccurst .eq. 0 ) then
        faccumulation   = 0
        flfsteps        = 0
        misrz           = 0
        m2isrz          = 0
        faccurst        = -1

        mlambda(:)      = 0.0d0
        m2lambda(:)     = 0.0d0

        if( has_lambdav ) then
            mlambdav(:) = 0.0d0
            m2lambdav(:) = 0.0d0
        end if

        fentaccu    = 0.0d0

        if( fenthalpy .or. fentropy ) then
            mlambdae(:) = 0.0d0
            lambda0(:)  = 0.0d0
            lambda1(:)  = 0.0d0
            lambda2(:)  = 0.0d0
            cds_hp(:)   = 0.0d0
            cds_hk(:)   = 0.0d0
        end if

        metot       = 0.0d0
        m2etot      = 0.0d0
        mepot       = 0.0d0
        m2epot      = 0.0d0
        mekin       = 0.0d0
        m2ekin      = 0.0d0

        epothist0   = 0.0d0
        epothist1   = 0.0d0

        CONList(:)%sdevtot = 0.0d0

        write(CST_OUT,'(A)') '#-------------------------------------------------------------------------------'
        write(CST_OUT,'(A)') '# INFO: ALL ACCUMULATORS WERE RESETED                                           '
        write(CST_OUT,'(A)') '#       PRODUCTION STAGE OF ACCUMULATION IS STARTED                             '
        write(CST_OUT,'(A)') '#-------------------------------------------------------------------------------'
    end if

    ! accumulate results -----------------------------------------------------
    if( faccurst .gt. 0 ) then
        faccurst = faccurst - 1
    end if

    ! it is not production part
    if( fstep .le. 0 ) return

    faccumulation = faccumulation + 1
    flfsteps = flfsteps + 1

    ! this will occur for velocity verlet algorithm
    if( faccumulation .le. 0 ) return

    ! calculate final constraint deviations
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        CONList(i)%deviation = get_deviation(CONList(i)%cv,CVContextP%CVsValues(ci),CONList(i)%value)
        CONList(i)%sdevtot = CONList(i)%sdevtot + CONList(i)%deviation**2
    end do

    ! calculate Z matrix --------------------------------------------------------
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfCONs
            cj = CONList(j)%cvindx
            fzv = 0.0
            do k=1,NumOfLAtoms
                fzv = fzv + MassInv(k)*dot_product(CVContextP%CVsDrvs(:,k,ci),CVContextP%CVsDrvs(:,k,cj))
            end do
            fz(i,j) = fzv
        end do
    end do

    ! calculate Z determinant ------------------------------------
    if( NumOfCONs .gt. 1 ) then
        ! LU decomposition
        call dgetrf(NumOfCONs,NumOfCONs,fz,NumOfCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_main!')
        end if
        fzdet = 1.0d0
        ! and finaly determinant
        do i=1,NumOfCONs
            if( indx(i) .ne. i ) then
                fzdet = - fzdet * fz(i,i)
            else
                fzdet = fzdet * fz(i,i)
            end if
        end do
    else
        fzdet = fz(1,1)
    end if

    invn = 1.0d0/real(faccumulation,PMFDP)

    ! calculate metric tensor correction --------------------------------------------------------
    isrz    = 1.0d0/sqrt(fzdet)
    dval    = isrz - misrz
    misrz   = misrz  + dval * invn
    dval2   = isrz - misrz
    m2isrz  = m2isrz + dval*dval2

    do i=1,NumOfCONs
        ! lambda ----------------------------------
        lam             = lambda(i)
        dval            = lam - mlambda(i)
        mlambda(i)      = mlambda(i) + dval * invn
        dval2           = lam - mlambda(i)
        m2lambda(i)     = m2lambda(i) + dval * dval2

        if( has_lambdav ) then
            mu              = lambdav(i)
            dval            = mu - mlambdav(i)
            mlambdav(i)     = mlambdav(i) + dval * invn
            dval2           = mu - mlambdav(i)
            m2lambdav(i)    = m2lambdav(i) + dval * dval2
        end if
    end do

end subroutine cst_core_analyze

!===============================================================================
! Subroutine:  cst_core_enthalpy_entropy
!===============================================================================

subroutine cst_core_enthalpy_entropy_lf

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i
    real(PMFDP)            :: invn, epot, ekin, etot
    real(PMFDP)            :: detot1, detot2
    real(PMFDP)            :: depot1, depot2
    real(PMFDP)            :: dekin1, dekin2
    real(PMFDP)            :: dlam1
    ! --------------------------------------------------------------------------

    if( .not. (fenthalpy .or. fentropy) ) return

    ! record history of lambda
    lambda0(:) = lambda1(:)   ! t-dt
    lambda1(:) = lambda2(:)   ! t
    lambda2(:) = lambda(:)    ! t+dt

    ! record history of Epot
    epothist0  = epothist1             ! t-dt
    epothist1  = PotEne + fepotoffset  ! t

    ! do we have enough samples?
    if( flfsteps .le. 2 ) return

    fentaccu = fentaccu + 1
    invn = 1.0d0/real(fentaccu,PMFDP)

    epot = epothist0
    ekin = KinEne + fekinoffset         ! t-dt

    etot = epot + ekin

    ! total energy
    detot1 = etot - metot
    metot  = metot  + detot1 * invn
    detot2 = etot - metot
    m2etot = m2etot + detot1 * detot2

    ! potential energy
    depot1 = epot - mepot
    mepot  = mepot  + depot1 * invn
    depot2 = epot - mepot
    m2epot = m2epot + depot1 * depot2

    ! potential energy
    dekin1 = ekin - mekin
    mekin  = mekin  + dekin1 * invn
    dekin2 = ekin - mekin
    m2ekin = m2ekin + dekin1 * dekin2

    ! entropy
    do i=1,NumOfCONs
        dlam1 = lambda0(i) - mlambdae(i)
        mlambdae(i) = mlambdae(i)  + dlam1 * invn
        cds_hp(i)   = cds_hp(i) + dlam1 * depot2
        cds_hk(i)   = cds_hk(i) + dlam1 * dekin2
    end do


end subroutine cst_core_enthalpy_entropy_lf

!===============================================================================

end module cst_core

