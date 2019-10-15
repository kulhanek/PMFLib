#ifndef ABFIntegratorGPRH
#define ABFIntegratorGPRH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// =============================================================================

#include <PMFMainHeader.hpp>
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <FortranMatrix.hpp>

//------------------------------------------------------------------------------

class CABFAccumulator;
class CEnergySurface;

// how to invert the matrix

enum EGPRINVMethod {
    EGPRINV_LU      = 1,    // DGETRI/DGETRF (via LU factorization)
    EGPRINV_SVD     = 2,    // via SVD factorization, divide and conquer driver
    EGPRINV_SVD2    = 3,    // via SVD factorization, simple driver
};

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing gaussian process
*/

class PMF_PACKAGE CABFIntegratorGPR {
public:
// constructor and destructor -------------------------------------------------
    CABFIntegratorGPR(void);
    virtual ~CABFIntegratorGPR(void);

// setup methods --------------------------------------------------------------
    /// set input ABF accumulator, only ABF forces are integrated
    void SetInputABFAccumulator(CABFAccumulator* p_accu);

    /// set output free energy surface
    void SetOutputFESurface(CEnergySurface* p_surf);

    /// multiply of bin sizes
    void SetWFac1(double wfac);

    /// multiply of bin sizes, if zero use wfac1
    void SetWFac2(double wfac);

    /// set sigmaf2
    void SetSigmaF2(double sigf2);

    /// set include error
    void SetIncludeError(bool set);

    /// skip energy calculation, it also disables errors
    void SetNoEnergy(bool set);

    /// set algorithm for LLS
    void SetINVMehod(EGPRINVMethod set);

    /// set rcond for SVD
    void SetRCond(double rcond);

    /// should we include glued area to energy calculation?
    void IncludeGluedAreas(bool set);

    /// set position of global minimum
    void SetGlobalMin(const CSmallString& spec);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout);

    /// get root mean square residuals
    double GetRMSR(int cv);

    /// get log of Marginal Likelihood
    double GetLogMarginalLikelihood(void);

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name);

// this destroys the state of the integrator
    /// remove mean force outliers from ABF data
    void FilterByMFZScore(double maxzscore,CVerboseStr& vout);

// section of private data ----------------------------------------------------
private:
    CABFAccumulator*    Accumulator;
    CEnergySurface*     FES;

    // GPR data
    int                     GPRSize;
    int                     NCVs;
    double                  WFac1;
    double                  WFac2;
    CSimpleVector<double>   CVLengths2;
    double                  SigmaF2;
    CFortranMatrix          K;          // kernels
    double                  logdetK;
    CSimpleVector<double>   Y;          // derivatives
    CSimpleVector<double>   GPRModel;   // weights
    bool                    NoEnergy;
    bool                    IncludeError;
    bool                    IncludeGluedBins;
    EGPRINVMethod           Method;
    bool                    GlobalMinSet;

    CSimpleVector<double>   ipos;
    CSimpleVector<double>   jpos;
    CSimpleVector<double>   gpos;   // global position
    CSimpleVector<double>   rk;
    CSimpleVector<double>   lk;
    CSimpleVector<double>   ik;

    // SVD setup
    double                  RCond;

    bool TrainGP(CVerboseStr& vout);
    void CalculateErrors(CSimpleVector<double>& gpos,CVerboseStr& vout); // gpos - position of global minimum

    double GetValue(const CSimpleVector<double>& position);
    double GetMeanForce(const CSimpleVector<double>& position,int icoord);
    double GetCov(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos);
    double GetVar(CSimpleVector<double>& lpos);
};

//------------------------------------------------------------------------------

#endif
