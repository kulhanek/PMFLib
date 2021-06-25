// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "PMFAccuInfoOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAccuInfoOptions::CPMFAccuInfoOptions(void)
{
    SetShowMiniUsage(true);
    SetAllowProgArgs(true);
}

//------------------------------------------------------------------------------

int CPMFAccuInfoOptions::CheckOptions(void)
{
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CPMFAccuInfoOptions::FinalizeOptions(void)
{
    bool ret_opt = false;

    if(GetOptHelp() == true) {
        PrintUsage();
        ret_opt = true;
    }

    if(GetOptVersion() == true) {
        PrintVersion();
        ret_opt = true;
    }

    if(ret_opt == true) {
        printf("\n");
        return(SO_EXIT);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CPMFAccuInfoOptions::CheckArguments(void)
{
    if( (GetNumberOfProgArgs() <= 0) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: at least one argument is expected, but %d is provided\n",
                (const char*)GetProgramName(),GetNumberOfProgArgs());
        IsError = true;
    }
    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
