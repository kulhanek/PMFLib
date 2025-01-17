#ifndef StringAdmOptionsH
#define StringAdmOptionsH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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
// ===============================================================================

#include <SimpleOptions.hpp>
#include <PMFMainHeader.hpp>

//------------------------------------------------------------------------------

class CStringAdmOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CStringAdmOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "stm-admin"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The utility that controls the behaviour and state of the server implementing the string method."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,Command)
    // options ------------------------------
    CSO_OPT(CSmallString,ServerKey)
    CSO_OPT(CSmallString,Password)
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Command,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "stm://server[:port]/command",                        /* parametr name */
                "The server can be specified as either the DNS name, the IP address of the server, or the keyword 'serverkey' or 'key'. "
                "In the latter case, the server information is read from the server key file. "
                "The listening port number can be optionally provided. Finally, the command specifies an administrative task, which can be one of the following: \n"
                "info                    = display information about registered clients\n"
                "flush                   = save the accumulated STM path on the server side\n"
                "get?file=output.path    = retrieve the accumulated STM path and save it locally to the file 'output.path'\n"
                "release?id=bead_id      = unregister a client associated with the specified bead ID\n"
                "terminate               = terminate the server at the end of the STM cycle\n"
                "shutdown                = stop the server execution\n"
                "errors                  = display errors from the server stack\n")   /* argument description */
                // description of options -----------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                ServerKey,                        /* option name */
                "server.key",                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "serverkey",                      /* long option name */
                "FILE",                           /* parametr name */
                "Name of file containing the server key. The server key contains the server name, port, and password. This option is mutually exclusive with 'password' option.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                Password,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                'p',                           /* short option name */
                "password",                        /* long option name */
                "FILE",                           /* parametr name */
                "Name of file containing the server magic word. If the pasword is not provided via this option or via the server key then it is read interactively from the keyboard. This option is mutually exclusive with 'serverkey' option.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Verbose,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'v',                           /* short option name */
                "verbose",                      /* long option name */
                NULL,                           /* parametr name */
                "Increase output verbosity.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Version,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "version",                      /* long option name */
                NULL,                           /* parametr name */
                "Output version information and exit.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Help,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'h',                           /* short option name */
                "help",                      /* long option name */
                NULL,                           /* parametr name */
                "Display this help and exit.")   /* option description */
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
    virtual int CheckArguments(void);
};

//------------------------------------------------------------------------------

#endif
