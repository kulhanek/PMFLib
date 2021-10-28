!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more detajls.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module abf_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_core_main
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine abf_core_main

    use abf_trajectory
    use abf_restart
    use abf_output
    use abf_client
    use abf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------

    select case(fmode)
        case(1)
            ! simplified ABF algorithm
            call abf_core_force_2p
        case(2)
            ! original ABF algorithm
            call abf_core_force_4p
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    call abf_client_exchange_data(.false.)

end subroutine abf_core_main

!===============================================================================
! Subroutine:  abf_core_force_2p
! this is leap-frog ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_force_2p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,e
    ! --------------------------------------------------------------------------

! shift accuvalue history
    cvhist(:,1) = cvhist(:,2)

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,2) = CVContext%CVsValues(ci)
    end do

! shift epot ene
    epothist(1) = epothist(2)
    if( fenthalpy ) then
        epothist(2) = PotEne - fepotaverage
    else
        epothist(2) = 0.0d0
    end if

    ersthist(1) = ersthist(2)
    if( fenthalpy ) then
        ersthist(2) = PMFEne
    else
        ersthist(2) = 0.0d0
    end if

! shift etot ene
    if( fentropy ) then
        ekinhist(1) = KinEne - fekinaverage ! shifted by -dt
    else
        ekinhist(1) = 0.0d0
    end if

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ki = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zd1(m,j,i) = v
            end do
        end do
    end do

    do i=1,NumOfABFCVs
        v = 0.0d0
        e = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! zd0 in t-dt
                ! Vel in t-1/2dt
                ! v0 (OldVel) in t-3/2dt
                ! a <- Vel(m,j)-v0(m,j))/fdtx in t-dt, do not use forces (Frc) because of SHAKE
                v = v + zd0(m,j,i)*(Vel(m,j)-v0(m,j))
                ! zd1 in t
                ! zd0 in t-dt
                ! vel in t-1/2dt
                e = e + (zd1(m,j,i)-zd0(m,j,i))* Vel(m,j)
            end do
        end do
        pxi0(i) = v / fdtx ! in t-dt
        pxip(i) = e / fdtx ! in t-1/2dt
    end do

    if( fstep .ge. 4 ) then

        ! complete ICF in t-dt
        ! pxi0 in t-dt
        ! pxi1 - old ABF forces in t-dt
        ! pxip in t-1/2dt
        ! pxim in t-3/2dt
        pxi0(:) = pxi0(:) - pxi1(:)
        pxim(:) = 0.5d0*(pxim(:)+pxip(:))

        ! total ABF force
        pxi0(:) = pxi0(:) + pxim(:)

        ! add data to accumulator
        if( fsmooth_enable ) then
            call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                call abf_accu_add_data_ksmooth(cvhist(:,1),pxi0,epothist(1))
            call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
        else
            call abf_accu_add_data_online(cvhist(:,1),pxi0,epothist(1))
        end if

        call abf_accu_add_data_entropy(cvhist(:,1),fzinv0,pxi0,pxi1,epothist(1),ersthist(1),ekinhist(1))
    end if

    ! backup to the next step
    zd0     = zd1
    pxim    = pxip
    v0      = Vel
    fzinv0  = fzinv

    ! apply ABF bias
    la(:) = 0.0d0

    ! apply force filters
    if( fapply_abf ) then
        ! calculate abf force to be applied
        select case(feimode)
            case(0)
                call abf_accu_get_data(cvhist(:,2),la)
            case(1)
                call abf_accu_get_data_lramp(cvhist(:,2),la)
            case(2)
                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                    call abf_accu_get_data_ksmooth(cvhist(:,2),la)
                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if

    ! keep ABF forces to subtract them in the next step
    pxi1 = la

    return

end subroutine abf_core_force_2p

!===============================================================================
! Subroutine:  abf_core_force_4p
! this is leap-frog ABF version, original implementation
!===============================================================================

subroutine abf_core_force_4p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer     :: i,j,k,m
    integer     :: ci,ck
    real(PMFDP) :: v
    ! --------------------------------------------------------------------------

! calculate acceleration in time t for all pmf atoms
    do i=1,NumOfLAtoms
        a1(:,i) = MassInv(i)*Frc(:,i)
    end do

! shift accuvalue history
    cvhist(:,1) = cvhist(:,2)
    cvhist(:,2) = cvhist(:,3)
    cvhist(:,3) = cvhist(:,4)

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,4) = CVContext%CVsValues(ci)
    end do

! apply force filters
    if( fapply_abf ) then
    ! calculate abf force to be applied -------------
        select case(feimode)
                case(0)
                    call abf_accu_get_data(cvhist(:,4),la)
                case(1)
                    call abf_accu_get_data_lramp(cvhist(:,4),la)
                case(2)
                    call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                        call abf_accu_get_data_ksmooth(cvhist(:,4),la)
                    call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select
    ! project abf force along coordinate ------------
        do i=1,NumOfABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                a1(:,j) = a1(:,j) + la(i) * MassInv(j) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    else
        ! we do not want to apply ABF force - set the forces to zero
        la(:) = 0.0d0
   end if

! rest of ABF stuff -----------------------------

    ! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

    ! pxip = zd0(t-dt)*[v(t-dt/2)/2 - dt*a1(t)/12]
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) - fdtx*a1(m,j)/12.0)
            end do
        end do
        pxip(i) = v
    end do

    ! ZD0(3, NumOfLAtoms, NumOfABFCVs)   <-- 1/dt * m_ksi grad ksi(r0)
    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ck = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k) * CVContext%CVsDrvs(m,j,ck)
                end do
                zd0(m,j,i) = v / fdtx
            end do
        end do
    end do

    ! pxim = zd0(t)*[v(t-dt/2)/2 + dt*a0(t-dt)/12]
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) + fdtx*a0(m,j)/12.0)
            end do
        end do
        pxim(i) = v
    end do

    !a0 <-- a1
    a0(:,:) = a1(:,:)

    ! pxi0 <--- pxi0 + pxip
    pxi0(:) = pxi0(:) + pxip(:)

    ! update accumulator - we need at least four samples
    if( fstep .ge. 4 ) then
        ! calculate coordinate values at time t-3/2dt
        do i=1,NumOfABFCVs
            cvave(i) = ABFCVList(i)%cv%get_average_value(cvhist(i,2),cvhist(i,3))
        end do

        ! add data to accumulator
        if( fsmooth_enable ) then
            call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                call abf_accu_add_data_ksmooth(cvhist(:,1),pxi0,0.0d0)
            call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
        else
            call abf_accu_add_data_online(cvhist(:,1),pxi0,0.0d0)
        end if
    end if

    ! pxi0 <--- -pxip + pxim + pxi1 - la/2
    pxi0(:) = -pxip(:) + pxim(:) + pxi1(:) - 0.5d0*la(:)

    ! pxi1 <--- -pxim - la/2
    pxi1(:) = -pxim(:) - 0.5d0*la(:)

    ! project updated acceleration back
    do i=1,NumOfLAtoms
        Frc(:,i) = a1(:,i)*Mass(i)
    end do

    return

end subroutine abf_core_force_4p

!===============================================================================
! subroutine:  abf_core_calc_Zmat
!===============================================================================

subroutine abf_core_calc_Zmat(ctx)

    use pmf_utils
    use abf_dat

    implicit none
    type(CVContextType) :: ctx
    integer             :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfABFCVs
            cj = ABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(ctx%CVsDrvs(:,k,ci),ctx%CVsDrvs(:,k,cj))
            end do
            fzinv(i,j) = fz(i,j)            ! we need this for LAPACK
        end do
    end do

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfABFCVs .gt. 1) then
        call dgetrf(NumOfABFCVs,NumOfABFCVs,fzinv,NumOfABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat!')
        end if

        call dgetri(NumOfABFCVs,fzinv,NumOfABFCVs,indx,vv,NumOfABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_core_calc_Zmat!')
        end if
    else
        fzinv(1,1)=1.0d0/fz(1,1)
    end if

    return

end subroutine abf_core_calc_Zmat

!===============================================================================

end module abf_core
