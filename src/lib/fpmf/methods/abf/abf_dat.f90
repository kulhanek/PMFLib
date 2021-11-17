!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module abf_dat

use pmf_sizes
use pmf_dat
use pmf_accu

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable ABF
                            ! 1 - enable ABF (simple ABF algorithm - 2p)
                            ! 2 - enable ABF (original ABF algorithm - 4p)
integer     :: fsample      ! output sample period in steps
logical     :: frestart     ! 1 - restart job with previous data, 0 - otherwise not
integer     :: frstupdate   ! how often is restart file written
integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"
logical     :: fapply_mask  ! off - disable ABF mask, on - enable ABF mask
logical     :: fapply_abf   ! on - apply ABF, off - do not apply ABF
integer     :: feimode      ! interpolation/extrapolation mode
                            ! 0 - disabled
                            ! 1 - linear ramp

! fmode == 3
integer     :: fsgframelen
integer     :: fsgorder
logical     :: fsgsmoothall

logical     :: fenthalpy    ! collect data for enthalpy calculation
logical     :: fentropy     ! collect data for entropy calculation
logical     :: frecord      ! record time progress

real(PMFDP) :: fepotaverage
real(PMFDP) :: fekinaverage

! linear ramp
integer     :: fhramp_min
integer     :: fhramp_max

! kernel smoothing
integer     :: fsmooth_kernel

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails
integer                 :: fmwamode             ! 0 - all is transferred
                                                ! 1 - only MICF is transferred

! abf server -----------------
integer                 :: client_id            ! abf walker client ID
integer                 :: failure_counter      ! current number of MWA failures

! item list --------------------------------------------------------------------
type CVTypeABF
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data
    integer                 :: set              ! coordinate set index
    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins
    real(PMFDP)             :: wfac             ! smoothing factor in number of bins
    logical                 :: nobias           ! do not use for bias and do not include in PMFAccumulator
end type CVTypeABF

! ----------------------
integer                     :: NumOfABFCVs          ! number of ALL CVs in a group
type(CVTypeABF),allocatable :: ABFCVList(:)         ! definition of CVs
integer                     :: NumOfBiasedABFCVs    ! number of biased CVS

! ----------------------

type,extends(PMFAccuType) :: ABFAccuType

    real(PMFDP),pointer    :: binpos(:,:)               ! position of grids
    real(PMFDP),pointer    :: weights(:)                ! mask weights
    real(PMFDP),pointer    :: nsamples(:)               ! number of hits into bins

! free energy
    real(PMFDP),pointer    :: micf(:,:)                 ! mean ICF
    real(PMFDP),pointer    :: m2icf(:,:)                ! M2 of ICF

! enthalpy
    real(PMFDP),pointer    :: mepot(:)                  ! mean of pot energy
    real(PMFDP),pointer    :: m2epot(:)                 ! M2 of pot energy
    real(PMFDP),pointer    :: merst(:)                  ! mean of rst energy
    real(PMFDP),pointer    :: m2erst(:)                 ! M2 of rst energy

! entropy
    real(PMFDP),pointer    :: metot(:)                  ! mean of tot energy
    real(PMFDP),pointer    :: m2etot(:)                 ! M2 of tot energy
    real(PMFDP),pointer    :: mpp(:,:)                  ! ICF+ETOT
    real(PMFDP),pointer    :: m2pp(:,:)                 !
    real(PMFDP),pointer    :: mpn(:,:)                  ! ICF-ETOT
    real(PMFDP),pointer    :: m2pn(:,:)                 !

! time recording for post-processing
    real(PMFDP),pointer    :: tcvs(:,:)
    real(PMFDP),pointer    :: tzinv(:,:,:)
    real(PMFDP),pointer    :: ticf(:,:)
    real(PMFDP),pointer    :: tmicf(:,:)
    real(PMFDP),pointer    :: tepot(:)
    real(PMFDP),pointer    :: terst(:)
    real(PMFDP),pointer    :: tekin(:)

! applied ICF - this is stored in accu but ignored
    real(PMFDP),pointer    :: bnsamples(:)              ! number of hits into bins
    real(PMFDP),pointer    :: bmicf(:,:)                ! applied MICF

! ABF force - incremental part for ABF-server
    real(PMFDP),pointer    :: inc_nsamples(:)           ! number of hits into bins
    real(PMFDP),pointer    :: inc_micf(:,:)             ! accumulated mean ICF
    real(PMFDP),pointer    :: inc_m2icf(:,:)            ! accumulated M2 of ICF
    real(PMFDP),pointer    :: inc_mepot(:)              ! accumulated mean of pot energy
    real(PMFDP),pointer    :: inc_m2epot(:)             ! accumulated M2 of pot energy
end type ABFAccuType

! ----------------------
type(ABFAccuType)           :: abfaccu                  ! accumulated forces
integer                     :: insidesamples
integer                     :: outsidesamples
! ----------------------

! global variables for force calculation ---------------------------------------
real(PMFDP)                 :: fdtx                 ! timestep

! global variables for abf - results -------------------------------------------
real(PMFDP),allocatable     :: fz(:,:)              ! Z matrix              in t
real(PMFDP),allocatable     :: fzinv(:,:)           ! inverse of Z matrix   in t
real(PMFDP),allocatable     :: vv(:)                ! for LU decomposition
integer,allocatable         :: indx(:)              ! for LU decomposition
real(PMFDP),allocatable     :: fzinv0(:,:)          ! inverse of Z matrix   in t-dt

! helper arrays -------
real(PMFDP),allocatable     :: a1(:,:)          ! acceleration in current step (t)
real(PMFDP),allocatable     :: a0(:,:)          ! acceleration from previous step (t-dt)
real(PMFDP),allocatable     :: v0(:,:)          ! velocity in previous step (t-dt)

real(PMFDP),allocatable     :: la(:)            ! ABF force in CV direction
real(PMFDP),allocatable     :: zd0(:,:,:)       ! ZD0
real(PMFDP),allocatable     :: zd1(:,:,:)       ! ZD1
real(PMFDP),allocatable     :: pxi0(:)          !
real(PMFDP),allocatable     :: pxi1(:)          !
real(PMFDP),allocatable     :: pxip(:)          !
real(PMFDP),allocatable     :: pxim(:)          !

real(PMFDP),allocatable     :: cvave(:)         ! CV values or their derivatives
real(PMFDP),allocatable     :: cvcur(:)
real(PMFDP),allocatable     :: cv1dr(:)
real(PMFDP),allocatable     :: cv2dr(:)

! ------------------------------------------------------------------------------

integer                     :: hist_len
real(PMFDP),allocatable     :: cvhist(:,:)          ! history of CV values (nCVS,hist_len)
real(PMFDP),allocatable     :: xihist(:,:)          ! history of CV values (hist_len,nCVs)
real(PMFDP),allocatable     :: micfhist(:,:)        ! history of ABF bias
real(PMFDP),allocatable     :: epothist(:)          ! history of Epot
real(PMFDP),allocatable     :: ersthist(:)          ! history of Erst
real(PMFDP),allocatable     :: ekinhist(:)          ! history of Ekin
real(PMFDP),allocatable     :: zinvhist(:,:,:)      ! history of zinv

! ------------------------------------------------------------------------------

real(PMFDP),allocatable     :: sg_c0(:)
real(PMFDP),allocatable     :: sg_c1(:)
real(PMFDP),allocatable     :: sg_c2(:)
real(PMFDP)                 :: invfdtx

! ------------------------------------------------------------------------------

! smoothing facility
integer                     :: max_snb_size
integer,allocatable         :: snb_list(:,:)
real(PMFDP),allocatable     :: sweights(:)

! FIXME - test
real(PMFDP)                 :: fzdet
real(PMFDP)                 :: fzdet0

real(PMFDP)                 :: fzdetA
real(PMFDP)                 :: fzdetA0

! ------------------------------------------------------------------------------

end module abf_dat

