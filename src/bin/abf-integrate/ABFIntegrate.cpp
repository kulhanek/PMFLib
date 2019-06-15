// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Martin Petrek, petrek@chemi.muni.cz
//                       Petr Kulhanek, kulhanek@enzim.hu
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

#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <ABFIntegratorRFD.hpp>
#include <ABFIntegratorRBF.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include "ABFIntegrate.hpp"
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegrate::CABFIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFIntegrate::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-integrate (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABFAccuName() != "-") {
        vout << "# ABF accu file (in)    : " << Options.GetArgABFAccuName() << endl;
    } else {
        vout << "# ABF accu file (in)    : - (standard input)" << endl;
    }
    if(Options.GetArgFEOutputName() != "-") {
        vout << "# Free energy file (out): " << Options.GetArgFEOutputName() << endl;
    } else {
        vout << "# Free energy file (out): - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Limit                 : all bins will be taken into account" << endl;
    } else {
        vout << "# Limit                 : " << Options.GetOptLimit() << endl;
    }
    vout << "# Integration method    : " << Options.GetOptMethod() << endl;
    vout << "# FD order              : " << Options.GetOptOrder() << endl;
    vout << "# Integration offset    : " << Options.GetOptOffset() << endl;
    vout << "# Periodicity           : " << bool_to_str(Options.GetOptPeriodicity()) << endl;
    vout << "# Integrate errors      : " << bool_to_str(Options.GetOptErrors()) << endl;
    vout << "# ------------------------------------------------" << endl;
    vout << "# Output FES format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# No header to output   : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# X format              : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format              : " << Options.GetOptOEFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABFAccuName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgFEOutputName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::Run(void)
{
   CABFAccumulator accumulator;

// load accumulator
    vout << endl;
    vout << "1) Loading ABF accumulator: " << Options.GetArgABFAccuName() << endl;
    try {
        accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done" << endl;

    // print CVS info
    accumulator.PrintCVSInfo(vout);

    // print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        fprintf(OutputFile,"# data integrated by    : ");
        if(Options.GetOptMethod() == "rfd" ) {
            fprintf(OutputFile," RFD (reverse finite difference)\n");
        } else if( Options.GetOptMethod() == "rbf" ){
            fprintf(OutputFile," RBF (radial basis functions)\n");
        } else {
            INVALID_ARGUMENT("method - not implemented");
        }
        fprintf(OutputFile,"# Number of coordinates : %d\n",accumulator.GetNumberOfCoords());
        fprintf(OutputFile,"# Total number of bins  : %d\n",accumulator.GetNumberOfBins());
        fprintf(OutputFile,"# Sample limit          : %d\n",Options.GetOptLimit());
        fprintf(OutputFile,"# Periodicity           : %s\n",(const char*)bool_to_str(Options.GetOptPeriodicity()));
        if( Options.GetOptWithErrors() ) {
        fprintf(OutputFile,"# Integrate errors      : energy+errors\n");
        } else {
        fprintf(OutputFile,"# Integrate errors      : %s\n",(const char*)bool_to_str(Options.GetOptErrors()));
        }
        fprintf(OutputFile,"# RFD/RBF order         : %d\n",Options.GetOptOrder());
    }

    // simple mode
    if( Options.GetOptWithErrors() == false ){
        return(RunSimple(accumulator)); // simple integration - either energy or error
    } else {
        return(RunWithErrors(accumulator)); // energy + errors
    }
}

//------------------------------------------------------------------------------

bool CABFIntegrate::RunSimple(CABFAccumulator& accumulator)
{

// prepare accumulator --------------------------
    vout << endl;
    vout << "2) Preparing ABF accumulator for integration"<< endl;
    PrepareAccumulator(accumulator,Options.GetOptErrors());
    vout << "   Done" << endl;

// integrate data ------------------------------
    vout << endl;
    vout << "3) ABF accumulator integration"<< endl;

    CEnergySurface     fes;

    if(Options.GetOptMethod() == "rfd" ) {
        CABFIntegratorRFD   integrator;

        integrator.SetVerbosity(Options.GetOptVerbose());
        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDOrder(Options.GetOptOrder());

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate() == false) {
            ES_ERROR("unable to prepare ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptMethod() == "rbf" ){
        CABFIntegratorRBF   integrator;

        integrator.SetVerbosity(Options.GetOptVerbose());
        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetGaussianWidth(Options.GetOptOrder());

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to prepare ABF accumulator");
            return(false);
        }
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    vout << "   Done" << endl;

// apply offset
    fes.ApplyOffset(Options.GetOptOffset() - fes.GetGlobalMinimumValue());

// print result ---------------------------------
    vout << endl;
    vout << "4) Writing results to file: " << Options.GetArgFEOutputName() << endl;
    CESPrinter printer;

    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());
    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    } else if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    } else if(Options.GetOptOutputFormat() == "fes") {
        printer.SetOutputFormat(EESPF_PMF_FES);
    } else {
        INVALID_ARGUMENT("output format - not implemented");
    }

    if(Options.GetOptPrintWithLimit()) {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetPrintedES(&fes);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output free energy file");
        return(false);
    }
    vout << "   Done" << endl;

    return(true);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::RunWithErrors(CABFAccumulator& accumulator)
{

// prepare accumulator --------------------------
    vout << endl;
    vout << "2) Preparing ABF accumulator for integration (energy)"<< endl;
    PrepareAccumulator(accumulator,false);
    vout << "   Done" << endl;

// integrate data ------------------------------
    vout << endl;
    vout << "3) ABF accumulator integration (energy)"<< endl;

    CEnergySurface     fes;

    if(Options.GetOptMethod() == "rfd" ) {
        CABFIntegratorRFD   integrator;

        integrator.SetVerbosity(Options.GetOptVerbose());
        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDOrder(Options.GetOptOrder());

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(false) == false) {
            ES_ERROR("unable to prepare ABF accumulator");
            return(false);
        }
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    vout << "   Done" << endl;

 // apply offset
    fes.ApplyOffset(Options.GetOptOffset() - fes.GetGlobalMinimumValue());

// -----------------------------------------
    CABFAccumulator error_accu;

 // load accumulator
     vout << endl;
     vout << "4) Loading ABF accumulator: " << Options.GetArgABFAccuName() << endl;

     if( fseek(InputFile, 0L, SEEK_SET) != 0 ) {
         ES_ERROR("unable to rewind input ABF accumulator file - is it a regular file?");
         return(false);
     }

     try {
         error_accu.Load(InputFile);
     } catch(...) {
         ES_ERROR("unable to load the input ABF accumulator file");
         return(false);
     }
     vout << "   Done" << endl;

 // prepare accumulator --------------------------
     vout << endl;
     vout << "5) Preparing ABF accumulator for integration (errors)"<< endl;
     PrepareAccumulator(error_accu,true);
     vout << "   Done" << endl;

 // integrate data ------------------------------
     vout << endl;
     vout << "6) ABF accumulator integration (errors)"<< endl;

     if(Options.GetOptMethod() == "rfd" ) {
         CABFIntegratorRFD   integrator;

         integrator.SetVerbosity(Options.GetOptVerbose());
         integrator.SetPeriodicity(Options.GetOptPeriodicity());
         integrator.SetFDOrder(Options.GetOptOrder());

         integrator.SetInputABFAccumulator(&error_accu);
         integrator.SetOutputFESurface(&fes);

         if(integrator.Integrate(true) == false) {
             ES_ERROR("unable to prepare ABF accumulator");
             return(false);
         }
     } else {
         INVALID_ARGUMENT("method - not implemented");
     }

     vout << "   Done" << endl;

     // final tuning
     fes.AdaptErrorsToGlobalMinimum();

// print result ---------------------------------
    vout << endl;
    vout << "6) Writing results to file: " << Options.GetArgFEOutputName() << endl;
    CESPrinter printer;

    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());
    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    } else if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    } else if(Options.GetOptOutputFormat() == "fes") {
        printer.SetOutputFormat(EESPF_PMF_FES);
    } else {
        INVALID_ARGUMENT("output format - not implemented");
    }

    if(Options.GetOptPrintWithLimit()) {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetIncludeErrors(true);
    printer.SetPrintedES(&fes);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output free energy file");
        return(false);
    }
    vout << "   Done" << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero
//    b) accumulated sums/sum squares will be recalculated to averages

void CABFIntegrate::PrepareAccumulator(CABFAccumulator& accumulator,bool errors)
{

    for(int icoord=0; icoord < accumulator.GetNumberOfCoords(); icoord++) {
        for(int ibin=0; ibin < accumulator.GetNumberOfBins(); ibin++) {

            double value = 0.0;
            double sum = 0.0;
            double sum_square = 0.0;
            int    nsamples = 0;
            int    newsamples = 0;

            nsamples = accumulator.GetNumberOfABFSamples(ibin);
            // abf force
            sum = accumulator.GetABFForceSum(icoord,ibin);
            sum_square = accumulator.GetABFForceSquareSum(icoord,ibin);

            if((nsamples > 0) && (nsamples > Options.GetOptLimit())) {
                if( errors == false) {
                    // calculate average
                    value = sum / nsamples;
                } else {
                    // calculate error of average
                    double sq = nsamples*sum_square - sum*sum;
                    if(sq > 0) {
                        sq = sqrt(sq) / nsamples;
                    } else {
                        sq = 0.0;
                    }
                    value = sq / sqrt((double)nsamples);
                }
                newsamples = nsamples;
            }
            accumulator.SetABFForceSum(icoord,ibin,value);
            accumulator.SetNumberOfABFSamples(ibin,newsamples);
        }
    }

    // calculate sampled area
    double maxbins = accumulator.GetNumberOfBins();
    double sampled = accumulator.GetNumberOfBinsWithABFLimit(Options.GetOptLimit());
    if( maxbins > 0 ){
        vout << "   Sampled area: " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegrate::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-integrate terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

