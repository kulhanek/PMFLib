!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_init_method
!===============================================================================

subroutine abf_init_method

    use abf_output
    use abf_restart
    use abf_trajectory
    use abf_client
    use abf_cvs

    implicit none
    ! --------------------------------------------------------------------------

    call abf_init_print_header
    call abf_init_arrays
    call abf_init_print_summary
    call abf_output_open
    call abf_restart_read
    call abf_trajectory_open
    call abf_client_register
    call abf_client_get_initial_data
    call abf_output_write_header
    call abf_cvs_init_values

end subroutine abf_init_method

!===============================================================================
! Subroutine:  abf_init_dat
!===============================================================================

subroutine abf_init_dat

    use abf_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode           = 0
    fsample         = 5000
    frestart        = .false.
    frstupdate      = 5000
    ftrjsample      = 0

    fapply_mask     = .false.
    fapply_abf      = .true.
    fupdate_abf     = .true.
    ficfsample      = 1

    fenthalpy       = .false.
    fenthalpy_der   = 0
    fentropy        = .false.
    fentdecomp      = .false.
    fenesample      = 1

    finclude_pv     = 0

    ftds_ekin_src   = 1
    ftds_add_bias   = .false.

    fepotaverage    = 0.0d0
    fekinaverage    = 0.0d0

    fepotsmooth     = 0
    ferstsmooth     = 0
    fekinsmooth     = 0

    feimode         = 1
    fhramp_min      = 20000
    fhramp_max      = 30000

    fusmode         = .false.
    falignbias      = .false.

    NumOfABFCVs         = 0
    NumOfABFCVs   = 0

    fserver_enabled = .false.
    fserverkey      = ''
    fserver         = ''
    fserverupdate   = 20000
    fconrepeats     = 0
    fabortonmwaerr  = .true.
    fmwamode        = 1

    client_id       = -1
    failure_counter = 0

    insidesamples   = 0
    outsidesamples  = 0

    fsmooth_kernel  = 0
    fswitch2zero    = .false.

    fene_step       = 0

    abf_p2_vx = 7
    abf_p2_px = 7

end subroutine abf_init_dat

!===============================================================================
! Subroutine:  abf_print_header
!===============================================================================

subroutine abf_init_print_header

    use pmf_constants
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' *********************** ADAPTIVE BIASING FORCE METHOD ************************ '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)

120 format(A)

end subroutine abf_init_print_header

!===============================================================================
! Subroutine:  abf_init_print_summary
!===============================================================================

subroutine abf_init_print_summary

    use prmfile
    use pmf_constants
    use pmf_dat
    use abf_dat
    use abf_cvs
    use pmf_utils

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' ABF mode (fmode)                        : ', fmode
    select case(fmode)
    case(1)
    write(PMF_OUT,120)  '      |-> ABF algorithm (2pX)'
    write(PMF_OUT,130)  '          Velocity order (abf_p2_vx)     : ', abf_p2_vx
    write(PMF_OUT,130)  '          Momenta order (abf_p2_px)      : ', abf_p2_px
    case(2)
    write(PMF_OUT,120)  '      |-> ABF algorithm (2pV)'
    write(PMF_OUT,130)  '          Velocity order (abf_p2_vx)     : ', abf_p2_vx
    write(PMF_OUT,130)  '          Momenta order (abf_p2_px)      : ', abf_p2_px

    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown fmode in abf_init_print_summary!')
    end select
    write(PMF_OUT,125)  ' Coordinate definition file (fabfdef)    : ', trim(fabfdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfABFCVs
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Control'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Apply ABF force (fapply_abf)            : ', prmfile_onoff(fapply_abf)
    write(PMF_OUT,125)  ' Update ABF force (fupdate_abf)          : ', prmfile_onoff(fupdate_abf)
    write(PMF_OUT,125)  ' ABF mask mode (fapply_mask)             : ', prmfile_onoff(fapply_mask)
    write(PMF_OUT,125)  ' ABF mask file (fabfmask)                : ', trim(fabfmask)
    write(PMF_OUT,130)  ' Sampling for ICF (ficfsample)           : ', ficfsample

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' US-ABF Control'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' US-ABF enable (fusmode)                 : ', prmfile_onoff(fusmode)
    write(PMF_OUT,125)  ' Align bias by a bin (falignbias)        : ', prmfile_onoff(falignbias)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Interpolation/Extrapolation '
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Switch ICF to zero (fswitch2zero)       : ', prmfile_onoff(fswitch2zero)
    write(PMF_OUT,130)  ' Extra/interpolation mode (feimode)      : ', feimode
    select case(feimode)
    case(0)
    write(PMF_OUT,120)  '      |-> Disabled'
    case(1)
    write(PMF_OUT,120)  '      |-> Min/Max linear ramp'
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case(2)
    write(PMF_OUT,120)  '      |-> Kernel smoother'
    write(PMF_OUT,130)  '          Kernel type (fsmooth_kernel)   : ', fsmooth_kernel
    select case(fsmooth_kernel)
    case(0)
    write(PMF_OUT,120)  '          |-> Epanechnikov (parabolic)'
    case(1)
    write(PMF_OUT,120)  '          |-> Triweight'
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown kernel in abf_init_print_summary!')
    end select
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case(3)
    write(PMF_OUT,120)  '      |-> Linear interpolation'
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown extrapolation/interpolation mode in abf_init_print_summary!')
    end select

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Enthalpy/entropy options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Accumulate enthalpy (fenthalpy)         : ', prmfile_onoff(fenthalpy)
    write(PMF_OUT,125)  ' Accumulate entropy (fentropy)           : ', prmfile_onoff(fentropy)
    write(PMF_OUT,125)  ' Decompose entropy (fentdecomp)          : ', prmfile_onoff(fentdecomp)
    write(PMF_OUT,125)  ' Use ABF bias for -TdS (ftds_add_bias)   : ', prmfile_onoff(ftds_add_bias)

    write(PMF_OUT,130)  ' Kinetic energy source (ftds_ekin_src)   : ', ftds_ekin_src
    select case(ftds_ekin_src)
    case(1)
    write(PMF_OUT,120)  '      |-> VV (Velocity Verlet KE)'
    case(2)
    write(PMF_OUT,120)  '      |-> LF (Leap-Frog KE)'
    case(3)
    write(PMF_OUT,120)  '      |-> HA (Harmonic Approximation Verlet KE)'
    case(4)
    write(PMF_OUT,120)  '      |-> LF (Leap-Frog KE), shifted by 0.5'
   ! case default
   ! call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown kinetic energy source in abf_init_print_summary!')
    end select
    write(PMF_OUT,150)  ' Potential energy offset (fepotaverage)  : ', pmf_unit_get_rvalue(EnergyUnit,fepotaverage),  &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,150)  ' Kinetic energy offset (fekinaverage)    : ', pmf_unit_get_rvalue(EnergyUnit,fekinaverage), &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,130)  ' Pot energy smoothing mode (fepotsmooth) : ', fepotsmooth
    write(PMF_OUT,130)  ' Rst energy smoothing mode (ferstsmooth) : ', ferstsmooth
    write(PMF_OUT,130)  ' Kin energy smoothing mode (fekinsmooth) : ', fekinsmooth
    write(PMF_OUT,130)  ' Sampling for -TdS and ENT (fenesample)  : ', fenesample

    write(PMF_OUT,130)  ' Include pV term (finclude_pv)           : ', finclude_pv
    select case(finclude_pv)
    case(0)
    write(PMF_OUT,120)  '      |-> none'
    case(1)
    write(PMF_OUT,120)  '      |-> p0V'
    case(2)
    write(PMF_OUT,120)  '      |-> pV'
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown pV mode in abf_init_print_summary!')
    end select

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fabfrst)                  : ', trim(fabfrst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Restart file update (frstupdate)        : ', frstupdate

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fabfout)                   : ', trim(fabfout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample
    write(PMF_OUT,125)  ' Trajectory file (fabftrj)               : ', trim(fabftrj)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' MWA server options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Server communication is                      : ', prmfile_onoff(fserver_enabled)
    if( fserver_enabled ) then
    write(PMF_OUT,125)  ' Server key file name (fserverkey)            : ', trim(fserverkey)
    else
    write(PMF_OUT,125)  ' Server key file name (fserverkey)            : ', 'none'
    end if
    write(PMF_OUT,130)  ' Server update interval (fserverupdate)       : ', fserverupdate
    write(PMF_OUT,130)  ' Number of connection repeats (fconrepeats)   : ', fconrepeats
    write(PMF_OUT,125)  ' Abort on MWA failure (fabortonmwaerr)        : ', prmfile_onoff(fabortonmwaerr)
    write(PMF_OUT,130)  ' MWA mode (fmwamode)                          : ', fmwamode
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of ABF collective variables'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfABFCVs
    write(PMF_OUT,140) i
    call abf_cvs_cv_info(ABFCVList(i))
    write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

 return

120 format(A)
125 format(A,A)
130 format(A,I6)
150 format(A,F10.1,1X,A)

140 format(' == Collective variable #',I4.4)

end subroutine abf_init_print_summary

!===============================================================================
! Subroutine:  abf_init_arrays
!===============================================================================

subroutine abf_init_arrays

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accu

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

! init accumulator
    call abf_accu_init

! general arrays --------------------------------
    allocate(                                   &
            la(NumOfABFCVs),                    &
            vint(3,NumOfLAtoms),                &
            pxia(NumOfABFCVs),                  &
            pxif(NumOfABFCVs),                  &
            picf(NumOfABFCVs),                  &
            sfac(NumOfABFCVs),                  &
            fz(NumOfABFCVs,NumOfABFCVs),        &
            fzinv(NumOfABFCVs,NumOfABFCVs),     &
            indx(NumOfABFCVs),                  &
            vv(NumOfABFCVs),                    &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for arrays used in ABF calculation!')
    end if

    la(:)       = 0.0d0
    vint(:,:)   = 0.0d0
    pxia(:)     = 0.0d0
    pxif(:)     = 0.0d0
    picf(:)     = 0.0d0
    fz(:,:)     = 0.0d0
    fzinv(:,:)  = 0.0d0
    sfac(:)     = 1.0d0

! history buffers ------------------------------------------
    select case(fmode)
        case(1)
            hist_len = 15
        case(2)
            hist_len = 8
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_init_arrays!')
    end select

    allocate(                                               &
            cvhist(NumOfABFCVs,hist_len),                   &
            micfhist(NumOfABFCVs,hist_len),                 &
            icfhist(NumOfABFCVs,hist_len),                  &
            icfphist(NumOfABFCVs,hist_len),                 &
            xphist(NumOfABFCVs,hist_len),                   &
            vhist(3,NumOfLAtoms,hist_len),                  &
            fhist(3,NumOfLAtoms,hist_len),                  &
            fzinvhist(NumOfABFCVs,NumOfABFCVs,hist_len),    &
            cvderhist(3,NumOfLAtoms,NumOfABFCVs,hist_len),  &
            zdhist(3,NumOfLAtoms,NumOfABFCVs,hist_len),     &
            epothist(hist_len),                             &
            ersthist(hist_len),                             &
            ekinhist(hist_len),                             &
            ekinlfhist(hist_len),                           &
            epvhist(hist_len),                              &
            volhist(hist_len),                              &
            enevalidhist(hist_len),                         &
            fdetzhist(hist_len),                            &
            fziihist(NumOfABFCVs,hist_len),                 &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for buffers used in ABF calculation!')
    end if

    cvhist(:,:)         = 0.0d0
    micfhist(:,:)       = 0.0d0
    icfhist(:,:)        = 0.0d0
    icfphist(:,:)       = 0.0d0
    vhist(:,:,:)        = 0.0d0
    fhist(:,:,:)        = 0.0d0
    epothist(:)         = 0.0d0
    ersthist(:)         = 0.0d0
    ekinhist(:)         = 0.0d0
    epvhist(:)          = 0.0d0
    ekinlfhist(:)       = 0.0d0
    xphist(:,:)         = 0.0d0
    fzinvhist(:,:,:)    = 0.0d0
    cvderhist(:,:,:,:)  = 0.0d0
    zdhist(:,:,:,:)     = 0.0d0
    volhist(:)          = 0.0d0
    fdetzhist(:)        = 0.0d0
    fziihist(:,:)       = 0.0d0
    enevalidhist(:)     = .false.

! other setup ----------------------------------------------

! sanity checks
    if( feimode .eq. 2 ) then
        call abf_init_snb_list
    end if

    if( (feimode .eq. 3) .and. (NumOfABFCVs .gt. 1) ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] feimode == 3 can be used only with one CV!')
    end if

end subroutine abf_init_arrays

!===============================================================================
! Subroutine:  abf_init_arrays
!===============================================================================

subroutine abf_init_snb_list

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accu

    implicit none
    integer         :: i,j,k,idx,alloc_failed
    real(PMFDP)     :: dx,u2,fac
    ! --------------------------------------------------------------------------

! use bigger distance buffer, fac is square of this buffer
    fac = 2**2 ! 2^2 = 4

! calculate the number of required pairs
    max_snb_size = 0
    do i=1,abfaccu%PMFAccuType%tot_nbins
        do j=1,abfaccu%PMFAccuType%tot_nbins
            u2 = 0.0d0
            do k=1,abfaccu%PMFAccuType%tot_cvs
                dx = abfaccu%PMFAccuType%sizes(k)%cv%get_deviation(abfaccu%binpos(k,i),abfaccu%binpos(k,j)) &
                   / (ABFCVList(k)%wfac * abfaccu%PMFAccuType%sizes(k)%bin_width)
                u2 = u2 + dx**2
            end do
            if( u2 .le. fac ) then
                max_snb_size = max_snb_size + 1
            end if
        end do
    end do

    max_snb_size = max_snb_size + 1 ! terminating null

    allocate(                                                       &
            snb_list(max_snb_size,abfaccu%PMFAccuType%tot_nbins),   &
            sweights(max_snb_size),                                 &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for arrays for kernel smoothing in abf_init_snb_list!')
    end if

    sweights(:)   = 0.0d0
    snb_list(:,:) = 0

    do i=1,abfaccu%PMFAccuType%tot_nbins
        idx = 0
        do j=1,abfaccu%PMFAccuType%tot_nbins
            u2 = 0.0d0
            do k=1,abfaccu%PMFAccuType%tot_cvs
                dx = abfaccu%PMFAccuType%sizes(k)%cv%get_deviation(abfaccu%binpos(k,i),abfaccu%binpos(k,j)) &
                   / (ABFCVList(k)%wfac * abfaccu%PMFAccuType%sizes(k)%bin_width)
                u2 = u2 + dx**2
            end do
            if( u2 .le. fac ) then
                idx = idx + 1
                if( idx .gt. max_snb_size ) then
                    call pmf_utils_exit(PMF_OUT,1, &
                                '[ABF] Max index into snb_list overflow in abf_init_snb_list!')
                end if
                snb_list(idx,i) = j
            end if
        end do
    end do

end subroutine abf_init_snb_list

!===============================================================================

end module abf_init
