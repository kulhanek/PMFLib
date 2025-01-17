// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <ABFProxy_dG.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFProxy_dG::CABFProxy_dG(void)
{
    SetType(ABF_MICF);

    Requires.push_back("ABF");
}

//------------------------------------------------------------------------------

CABFProxy_dG::~CABFProxy_dG(void)
{
}

//------------------------------------------------------------------------------

bool CABFProxy_dG::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "ABF" ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

void CABFProxy_dG::SetType(EABFdGType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(ABF_MICF):
            Provide = "ABF dG(x)";
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CABFProxy_dG::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples = 0.0;
    double  micf     = 0.0;
    double  m2icf    = 0.0;
    double  ncorr    = Accu->GetNCorr();

    switch(Type){
    // -------------------
        case(ABF_MICF):
            nsamples = Accu->GetData("NSAMPLES",ibin);
            micf     = Accu->GetData("MICF",ibin,icv);
            m2icf    = Accu->GetData("M2ICF",ibin,icv);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
// mean force
        // -------------------
        case(E_PROXY_VALUE):
            return( micf );
        // -------------------
        case(E_PROXY_SIGMA):
            return( sqrt(m2icf / nsamples) );
        // -------------------
        case(E_PROXY_ERROR):
            return( sqrt(m2icf * ncorr) / nsamples );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



