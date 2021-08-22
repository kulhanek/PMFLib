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

module usabf_dat

use pmf_sizes
use pmf_dat
use pmf_accu

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable US-ABF
                            ! 1 - enable ABF (simple ABF algorithm - 2p)
                            ! 2 - enable ABF (original ABF algorithm - 4p)
integer     :: fsample      ! output sample period in steps
integer     :: frstupdate   ! how often is restart file written

integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"

logical     :: fenthalpy    ! collect data for enthalpy calculation
logical     :: fentropy     ! collect data for entropy calculation
real(PMFDP) :: fepotoffset
real(PMFDP) :: fekinoffset

integer     :: fblock_size  ! pre-blocking

! item list --------------------------------------------------------------------
type CVTypeUSABF
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data
    integer                 :: set              ! coordinate set index
    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins

    real(PMFDP)             :: target_value     ! required value of restraint
    real(PMFDP)             :: force_constant   ! sigma value

    real(PMFDP)             :: deviation        ! deviation between real and actual value
    real(PMFDP)             :: energy
end type CVTypeUSABF

! ----------------------
integer                         :: NumOfUSABFCVs        ! number of CVs in a group
type(CVTypeUSABF),allocatable   :: USABFCVList(:)       ! definition of CVs

real(PMFDP)                     :: TotalUSABFEnergy
! ----------------------

type,extends(PMFAccuType)   :: USABFAccuType

    integer,pointer         :: nsamples(:)              ! number of hits into bins

    ! MICF
    real(PMFDP),pointer     :: micf(:,:)                ! mean ICF - total
    real(PMFDP),pointer     :: m2icf(:,:)               ! M2 of ICF - total

    ! ENTHALPY
    real(PMFDP),pointer     :: metot(:)                 ! mean of total energy
    real(PMFDP),pointer     :: m2etot(:)                ! M2 of total energy
    real(PMFDP),pointer     :: mepot(:)                 ! mean of potential energy
    real(PMFDP),pointer     :: m2epot(:)                ! M2 of potential energy
    real(PMFDP),pointer     :: merst(:)                 ! mean of restraint energy
    real(PMFDP),pointer     :: m2erst(:)                ! M2 of restraint energy

    ! ENTROPY
    real(PMFDP),pointer     :: c11hh(:,:)               ! c11 - total/total

    real(PMFDP),pointer     :: mcovhh(:,:)
    real(PMFDP),pointer     :: m2covhh(:,:)

    ! pre-blocking
    integer,pointer         :: block_nsamples(:)        ! number of hits into bins
    real(PMFDP),pointer     :: block_micf(:,:)          ! mean ICF - total
    real(PMFDP),pointer     :: block_metot(:)           ! mean of total energy
    real(PMFDP),pointer     :: block_mepot(:)           ! mean of potential energy
    real(PMFDP),pointer     :: block_merst(:)           ! mean of restraint energy
    real(PMFDP),pointer     :: block_c11hh(:,:)         ! c11 - total/total

end type USABFAccuType

! ----------------------
type(USABFAccuType)         :: usabfaccu                ! accumulated forces
integer                     :: insidesamples
integer                     :: outsidesamples
! ----------------------

! global variables for force calculation ---------------------------------------
real(PMFDP)                 :: fdtx                 ! timestep

! global variables for abf - results -------------------------------------------
real(PMFDP),allocatable     :: fz(:,:)              ! Z matrix
real(PMFDP),allocatable     :: fzinv(:,:)           ! inverse of Z matrix
real(PMFDP),allocatable     :: vv(:)                ! for LU decomposition
integer,allocatable         :: indx(:)

! helper arrays -------
real(PMFDP),allocatable     :: a0(:,:)        ! acceleration from previous step (t-dt)
real(PMFDP),allocatable     :: a1(:,:)        ! acceleration in current step (t)
real(PMFDP),allocatable     :: v0(:,:)        ! velocity in previous step

real(PMFDP),allocatable     :: la(:)          ! ABF force in coordinate direction
real(PMFDP),allocatable     :: zd0(:,:,:)     ! ZD0
real(PMFDP),allocatable     :: zd1(:,:,:)     ! ZD1
real(PMFDP),allocatable     :: pxi0(:)        !
real(PMFDP),allocatable     :: pxi1(:)        !
real(PMFDP),allocatable     :: pxip(:)        !
real(PMFDP),allocatable     :: pxim(:)        !
real(PMFDP),allocatable     :: avg_values(:)  ! average values of coordinates at t - 3/2dt

real(PMFDP),allocatable     :: cvaluehist0(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist1(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist2(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist3(:)   ! history of coordinate values

real(PMFDP)                 :: epothist0   ! history of Epot
real(PMFDP)                 :: epothist1   ! history of Epot
real(PMFDP)                 :: epothist2   ! history of Epot
real(PMFDP)                 :: epothist3   ! history of Epot

real(PMFDP)                 :: ekinhist0   ! history of Ekin
real(PMFDP)                 :: ekinhist1   ! history of Ekin
real(PMFDP)                 :: ekinhist2   ! history of Ekin
real(PMFDP)                 :: ekinhist3   ! history of Ekin

real(PMFDP)                 :: ersthist0   ! history of Erst (PMFEne)
real(PMFDP)                 :: ersthist1   ! history of Erst
real(PMFDP)                 :: ersthist2   ! history of Erst
real(PMFDP)                 :: ersthist3   ! history of Erst

! ------------------------------------------------------------------------------

end module usabf_dat
