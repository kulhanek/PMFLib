!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module usabf_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  usabf_accu_init
!===============================================================================

subroutine usabf_accu_init()

    use usabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    usabfaccu%tot_cvs = NumOfUSABFCVs
    allocate(usabfaccu%sizes(usabfaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator!')
    endif

    tot_nbins = 1
    do i=1,usabfaccu%tot_cvs
        usabfaccu%sizes(i)%min_value  = USABFCVList(i)%min_value
        usabfaccu%sizes(i)%max_value  = USABFCVList(i)%max_value
        usabfaccu%sizes(i)%nbins      = USABFCVList(i)%nbins
        usabfaccu%sizes(i)%width      = abs(usabfaccu%sizes(i)%max_value - usabfaccu%sizes(i)%min_value)
        usabfaccu%sizes(i)%bin_width  = usabfaccu%sizes(i)%width / usabfaccu%sizes(i)%nbins
        usabfaccu%sizes(i)%cv         => USABFCVList(i)%cv
        tot_nbins = tot_nbins * usabfaccu%sizes(i)%nbins
    end do

    usabfaccu%tot_nbins = tot_nbins

    ! ABF force arrays
    allocate(  &
            usabfaccu%nsamples(usabfaccu%tot_nbins), &
            usabfaccu%micf(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
            usabfaccu%m2icf(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
            stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (icf)!')
    endif

! enthalpy ---------------------------------------------------------------------
    if( fenthalpy .or. fentropy ) then
        allocate(  &
                usabfaccu%metot(usabfaccu%tot_nbins), &
                usabfaccu%m2etot(usabfaccu%tot_nbins), &

                usabfaccu%mepot(usabfaccu%tot_nbins), &
                usabfaccu%m2epot(usabfaccu%tot_nbins), &

                usabfaccu%merst(usabfaccu%tot_nbins), &
                usabfaccu%m2erst(usabfaccu%tot_nbins), &
                stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (emthalpy)!')
        endif
    end if

! entropy ----------------------------------------------------------------------
    if( fentropy ) then
        if( fblock_size .eq. 0 ) then
            allocate(  &
                    usabfaccu%c11hh(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
                    stat = alloc_failed)

            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (entropy)!')
            endif
        else
            allocate(  &
                    usabfaccu%mcovhh(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
                    usabfaccu%m2covhh(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
                    stat = alloc_failed)

            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (entropy)!')
            endif
        end if
    end if

! pre-blocking ------------------------------------------------------------------
    if( fblock_size .gt. 0 ) then

        ! ABF force arrays
        allocate(  &
                usabfaccu%block_nsamples(usabfaccu%tot_nbins), &
                usabfaccu%block_micf(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
                stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (block_icf)!')
        endif

        if( fenthalpy .or. fentropy ) then
            allocate(  &
                    usabfaccu%block_metot(usabfaccu%tot_nbins), &
                    usabfaccu%block_mepot(usabfaccu%tot_nbins), &
                    usabfaccu%block_merst(usabfaccu%tot_nbins), &
                    stat = alloc_failed)

            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (block_emthalpy)!')
            endif
        end if

        if( fentropy ) then
            allocate(  &
                    usabfaccu%block_c11hh(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
                    stat = alloc_failed)

            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (block_entropy)!')
            endif
        end if
    end if

    call usabf_accu_clear()

    return

end subroutine usabf_accu_init

!===============================================================================
! Subroutine:  usabf_accu_clear
!===============================================================================

subroutine usabf_accu_clear()

    use usabf_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    usabfaccu%nsamples(:)       = 0

    usabfaccu%micf(:,:)         = 0.0d0
    usabfaccu%m2icf(:,:)        = 0.0d0

    if( fenthalpy .or. fentropy ) then
        usabfaccu%metot(:)      = 0.0d0
        usabfaccu%m2etot(:)     = 0.0d0

        usabfaccu%mepot(:)      = 0.0d0
        usabfaccu%m2epot(:)     = 0.0d0

        usabfaccu%merst(:)      = 0.0d0
        usabfaccu%m2erst(:)     = 0.0d0
    end if

    if( fentropy ) then
        if( fblock_size .eq. 0 ) then
            usabfaccu%c11hh(:,:)    = 0.0d0
        else
            usabfaccu%mcovhh(:,:)   = 0.0d0
            usabfaccu%m2covhh(:,:)  = 0.0d0
        end if
    end if

    if( fblock_size .gt. 0 ) then
        usabfaccu%block_nsamples(:) = 0
        usabfaccu%block_micf(:,:)   = 0.0d0

        if( fenthalpy .or. fentropy ) then
            usabfaccu%block_metot(:)    = 0.0d0
            usabfaccu%block_mepot(:)    = 0.0d0
            usabfaccu%block_merst(:)    = 0.0d0
        end if

        if( fentropy ) then
            usabfaccu%block_c11hh(:,:)   = 0.0d0
        end if
    end if

end subroutine usabf_accu_clear

!===============================================================================
! Subroutine:  usabf_accu_write
!===============================================================================

subroutine usabf_accu_write(iounit)

    use usabf_dat

    implicit none
    integer                     :: iounit
    !---------------------------------------------------------------------------

    usabfaccu%method = 'US-ABF'
    call pmf_accu_write_header(usabfaccu%PMFAccuType,iounit)
    call pmf_accu_write_ibuf_B(usabfaccu%PMFAccuType,iounit,'NSAMPLES',     'AD',usabfaccu%nsamples)
    call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'MICF',         'WA',usabfaccu%micf)
    call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'M2ICF',        'M2',usabfaccu%m2icf,'MICF')

    if( fenthalpy .or. fentropy ) then
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'METOT',    'WA',usabfaccu%metot)
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'M2ETOT',   'M2',usabfaccu%m2etot,'METOT')

        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'MEPOT',    'WA',usabfaccu%mepot)
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'M2EPOT',   'M2',usabfaccu%m2epot,'MEPOT')

        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'MERST',    'WA',usabfaccu%merst)
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'M2ERST',   'M2',usabfaccu%m2erst,'MERST')
    end if

    if( fentropy ) then
        if( fblock_size .eq. 0 ) then
            call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'C11HH',    'CO',usabfaccu%c11hh,'MICF','METOT')
        else
            call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'MCOVHH',   'WA',usabfaccu%mcovhh)
            call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'M2COVHH',  'M2',usabfaccu%m2covhh, 'MCOVHH')
        end if
    end if

end subroutine usabf_accu_write

!===============================================================================
! Subroutine:  usabf_accu_add_data_online
!===============================================================================

subroutine usabf_accu_add_data_online(cvs,gfx,epot,ekin,erst)

    use usabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: ekin
    real(PMFDP)    :: erst
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, etot, icf
    real(PMFDP)    :: detot1, detot2
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: derst1, derst2
    real(PMFDP)    :: dicf1, dicf2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(usabfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    usabfaccu%nsamples(gi0) = usabfaccu%nsamples(gi0) + 1
    invn = 1.0d0 / real(usabfaccu%nsamples(gi0),PMFDP)

    if( fenthalpy .or. fentropy ) then
        ! total energy
        etot = epot + ekin + erst
        detot1 = etot - usabfaccu%metot(gi0)
        usabfaccu%metot(gi0)  = usabfaccu%metot(gi0)  + detot1 * invn
        detot2 = etot - usabfaccu%metot(gi0)
        usabfaccu%m2etot(gi0) = usabfaccu%m2etot(gi0) + detot1 * detot2

        ! potential energy
        depot1 = epot - usabfaccu%mepot(gi0)
        usabfaccu%mepot(gi0)  = usabfaccu%mepot(gi0)  + depot1 * invn
        depot2 = epot - usabfaccu%mepot(gi0)
        usabfaccu%m2epot(gi0) = usabfaccu%m2epot(gi0) + depot1 * depot2

        ! restraint energy
        derst1 = erst - usabfaccu%merst(gi0)
        usabfaccu%merst(gi0)  = usabfaccu%merst(gi0)  + derst1 * invn
        derst2 = erst - usabfaccu%merst(gi0)
        usabfaccu%m2erst(gi0) = usabfaccu%m2erst(gi0) + derst1 * derst2

    end if

    do i=1,NumOfUSABFCVs
        icf = gfx(i)

        dicf1 = - icf - usabfaccu%micf(i,gi0)
        usabfaccu%micf(i,gi0)  = usabfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = - icf -  usabfaccu%micf(i,gi0)
        usabfaccu%m2icf(i,gi0) = usabfaccu%m2icf(i,gi0) + dicf1 * dicf2

        if( fentropy ) then
            usabfaccu%c11hh(i,gi0)  = usabfaccu%c11hh(i,gi0) + dicf1     * detot2
        end if
    end do

end subroutine usabf_accu_add_data_online

!===============================================================================
! Subroutine:  usabf_accu_add_data_block
!===============================================================================

subroutine usabf_accu_add_data_block(cvs,gfx,epot,ekin,erst)

    use usabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: ekin
    real(PMFDP)    :: erst
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, etot, icf, c11hh, block_invn
    real(PMFDP)    :: detot1, detot2
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: dicf1, dicf2
    real(PMFDP)    :: derst1, derst2
    real(PMFDP)    :: dc11hh1, dc11hh2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(usabfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    usabfaccu%block_nsamples(gi0) = usabfaccu%block_nsamples(gi0) + 1
    invn = 1.0d0 / real(usabfaccu%block_nsamples(gi0),PMFDP)

    if( fenthalpy .or. fentropy ) then
        ! total energy
        etot = epot + ekin + erst
        detot1 = etot - usabfaccu%block_metot(gi0)
        usabfaccu%block_metot(gi0)  = usabfaccu%block_metot(gi0)  + detot1 * invn
        detot2 = etot - usabfaccu%block_metot(gi0)

        ! potential energy
        depot1 = epot - usabfaccu%block_mepot(gi0)
        usabfaccu%block_mepot(gi0)  = usabfaccu%block_mepot(gi0)  + depot1 * invn
        depot2 = epot - usabfaccu%block_mepot(gi0)

        ! restraint energy
        derst1 = erst - usabfaccu%block_merst(gi0)
        usabfaccu%block_merst(gi0)  = usabfaccu%block_merst(gi0)  + derst1 * invn
        derst2 = erst - usabfaccu%block_merst(gi0)
    end if

    do i=1,NumOfUSABFCVs
        icf = gfx(i)

        dicf1 = - icf - usabfaccu%block_micf(i,gi0)
        usabfaccu%block_micf(i,gi0)  = usabfaccu%block_micf(i,gi0)  + dicf1 * invn
        dicf2 = - icf -  usabfaccu%block_micf(i,gi0)

        if( fentropy ) then
            usabfaccu%block_c11hh(i,gi0)  = usabfaccu%block_c11hh(i,gi0) + dicf1     * detot2
        end if
    end do

    ! do we have enough samples?
    if( usabfaccu%block_nsamples(gi0) .lt. fblock_size ) return

    ! process them
    block_invn = invn

    ! increase number of samples
    usabfaccu%nsamples(gi0) = usabfaccu%nsamples(gi0) + 1
    invn = 1.0d0 / real(usabfaccu%nsamples(gi0),PMFDP)

    if( fenthalpy .or. fentropy ) then
        ! total energy
        etot = usabfaccu%block_metot(gi0)
        detot1 = etot - usabfaccu%metot(gi0)
        usabfaccu%metot(gi0)  = usabfaccu%metot(gi0)  + detot1 * invn
        detot2 = etot - usabfaccu%metot(gi0)
        usabfaccu%m2etot(gi0) = usabfaccu%m2etot(gi0) + detot1 * detot2

        ! potential energy
        epot = usabfaccu%block_mepot(gi0)
        depot1 = epot - usabfaccu%mepot(gi0)
        usabfaccu%mepot(gi0)  = usabfaccu%mepot(gi0)  + depot1 * invn
        depot2 = epot - usabfaccu%mepot(gi0)
        usabfaccu%m2epot(gi0) = usabfaccu%m2epot(gi0) + depot1 * depot2

        ! restraint energy
        erst = usabfaccu%block_merst(gi0)
        derst1 = erst - usabfaccu%merst(gi0)
        usabfaccu%merst(gi0)  = usabfaccu%merst(gi0)  + derst1 * invn
        derst2 = erst - usabfaccu%merst(gi0)
        usabfaccu%m2erst(gi0) = usabfaccu%m2erst(gi0) + derst1 * derst2
    end if

    do i=1,NumOfUSABFCVs
        icf = usabfaccu%block_micf(i,gi0)

        dicf1 = icf - usabfaccu%micf(i,gi0)
        usabfaccu%micf(i,gi0)  = usabfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = icf -  usabfaccu%micf(i,gi0)
        usabfaccu%m2icf(i,gi0) = usabfaccu%m2icf(i,gi0) + dicf1 * dicf2

        if( fentropy ) then
            c11hh = usabfaccu%block_c11hh(i,gi0) * block_invn
            dc11hh1 = c11hh - usabfaccu%mcovhh(i,gi0)
            usabfaccu%mcovhh(i,gi0)  = usabfaccu%mcovhh(i,gi0)  + dc11hh1 * invn
            dc11hh2 = c11hh -  usabfaccu%mcovhh(i,gi0)
            usabfaccu%m2covhh(i,gi0) = usabfaccu%m2covhh(i,gi0) + dc11hh1 * dc11hh2
        end if
    end do

    usabfaccu%block_nsamples(gi0)       = 0
    usabfaccu%block_micf(:,gi0)         = 0.0d0

    if( fenthalpy .or. fentropy ) then
        usabfaccu%block_metot(gi0)      = 0.0d0
        usabfaccu%block_mepot(gi0)      = 0.0d0
        usabfaccu%block_merst(gi0)      = 0.0d0
    end if

    if( fentropy ) then
        usabfaccu%block_c11hh(:,gi0)    = 0.0d0
    end if

end subroutine usabf_accu_add_data_block

!===============================================================================

end module usabf_accu
