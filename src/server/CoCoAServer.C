//   Copyright (c)  2007-2010  John Abbott and Anna Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/library.H"
#include "CoCoA5io.H"  // #include "CoCoA4io.H"
#include "GlobalIO.H"
#include "ServerOp.H"
#include "RegisterServerOps.H"
#include "SocketStream.H"

using namespace CoCoA;

// #include <string>
using std::string;
// #include <iostream>
using std::endl;
using std::clog;
using std::hex;
using std::dec;
#include <memory>
using std::unique_ptr;
#include <cstdlib>
// using exit

// print CoCoA4 style info and result on GlobalOutput()
static bool CoCoA5Output;

// print CoCoA4 style result on clog
static bool ServerDebug = false;

//------------------------
void ChangeDefaultMempoolDebugLevel()
{
  MemPoolDebug::ourDefaultMarginSize=0; // 4
  MemPoolDebug::ourInitialDebugLevel=0; // 2
  MemPoolDebug::ourInitialVerbosityLevel = 0;
}


void ChangeDefaultStatLevel(int level)
{  // to make this mechanism work, use SetVerbosityLevel
//   if (level<0)  // default
//     GReductor::ourDefaultStatLevel = -1;     // <-- for the server (with timings)
//   if (level==0)
//     GReductor::ourDefaultStatLevel = 0;  // <-- for the benchmarks with nr of reductions
//   if (level>0)
//     GReductor::ourDefaultStatLevel = level;  // <-- debugging
}

/*
  GReductor::ourDefaultStatLevel
  -1 Gives nothing.  This is the default for the LIBRARY and SERVER
  0 adds some final stats.  This is the default for the BENCHMARKS.
  1 adds starting and final stats
  2 same as 1
  3 adds pair by pair and deg by deg minimal stats
  4 adds reduction by reduction stats
  5 adds the printing of input polys
  6 adds some poly by poly and deg by deg final data
*/

//----------------------------------------------------------------------
void PrintTime(double T)
{
  PrintTimeToLog(T);
  if (CoCoA5Output) PrintTimeToCoCoA5(T);
}

//----------------------------------------------------------------------
void program()
{
  CheckOperationRegistry();
  ChangeDefaultMempoolDebugLevel();
  GlobalManager CoCoAFoundations;
  //  debug_new::verbose_obj trace_freestore; // merely activates logging of new/delete

  string ExName = "NoName";
  string OperationString;
  string tag;
  double PkgVersion = 1.0;    // careful, PkgVersion is a double
  const string V103 = "    1.03: CoCoAServer waits for <end_of_session/>;";
  const string V102 = "    1.02: LT5 returns an ideal;";
  const string V101 = "    1.01: added <number_arguments>;";
  
  ChangeDefaultStatLevel(-1);

  // Exchange greetings
  SkipTag(GlobalInput(), "<Greetings>");
  SkipTag(GlobalInput(), "CoCoAServer?");
  if (CoCoA5Output)  PrintVersionToCoCoA5();
  GlobalInput() >> tag;
  while (tag!="</Greetings>")
  {
    if      (tag == "nolimits")     /* unused option */;
    else if (tag == "version")      GlobalInput() >> PkgVersion;
    else if (tag == "example_name") GlobalInput() >> ExName;
    //    else if (tag == "stat_level")   GlobalInput() >> StatLevel;
    else ThrowInputError(tag);
    GlobalInput() >> tag;
  }
//   if (PkgVersion < 1.02)
//   {
//     GlobalOutput() << 
//       "Block " << 
//       " Print   \"  Old version of cocoa5.cpkg: \", $cocoa5.Version();\n" <<
//       " PrintLn \"  (expected >= 1.02)\";\n"           << endl;
//   /*if (PkgVersion < 1.02)*/GlobalOutput() << "PrintLn \"" << V102 << "\";\n";
//     if (PkgVersion < 1.01)  GlobalOutput() << "Print \""   << V101 << "\";\n";
//     GlobalOutput() << "EndBlock;" << endl;
//   }
  OperationString = ReadOperationString(GlobalInput());
  ServerOp op = GetOperation(OperationString);
  if ( IsVoidOperation(op) )
    CoCoA_THROW_ERROR("No operation defined for: " + OperationString, "CoCoAServer");
  //2019  else    if (CoCoA5Output)      EndOfTransmissionToCoCoA5();
  GlobalLogput() << "-- INPUT: " << ExName << endl << "-- " << op << endl;
  SkipTag(GlobalInput(), "<number_arguments>");
  int NumArgs;
  GlobalInput() >> NumArgs;
  SkipTag(GlobalInput(), "</number_arguments>");
  op->myClear(); // in case some rubbish was left around after an exception was thrown
  op->myReadArgs(GlobalInput(), NumArgs);
  double t0 = CpuTime();
  op->myCompute();
  PrintTime(CpuTime() - t0);
  if (ServerDebug)  op->myWriteResult(clog);
  if (CoCoA5Output) op->myWriteResult(GlobalOutput());
  GlobalOutput() << "SERVER HAS FINISHED" << endl;
  op->myClear();
  if (PkgVersion > 1.02  &&  CoCoA5Output)
  {
    //2019    EndOfTransmissionToCoCoA4();
    // Wait for CoCoA4 to end the session before quitting
    // Must do this otherwise have a race condition: if we close the
    // socket too early then CoCoA4 cannot read all the data on some Linuxes.
    SkipTag(GlobalInput(), "<end_of_session/>");
  }
  GlobalLogput() << endl;
}


void SendErrorAndQuit(const string& s)
{
  GlobalOutput() <<
    "Catch  Error(\"CoCoAServer: " << s << "\") In " <<
    ServerOpBase::ourVarName4 << " EndCatch;" << endl;
}


//----------------------------------------------------------------------
int main(int argc, char**argv)
{
  bool UsingSocket = true;
  const unsigned short DefaultPort = 0xc0c0; // 49344
  unsigned short port = DefaultPort;

  unique_ptr<SocketStream> SockPtr;
  if (argc > 1)
  {
    if (argv[1]==string("-h") || argv[1] == string("--help"))
    {
      GlobalErrput() << "Options should be:" << endl;
      GlobalErrput() << "CoCoAServer -h --- print options" << endl;
      GlobalErrput() << "CoCoAServer    --- default port" << endl;
      GlobalErrput() << "CoCoAServer -d --- debugging (stdin/out)" << endl;
      GlobalErrput() << "CoCoAServer -t --- timing (stdin/out & reduced output)" << endl;
      GlobalErrput() << "CoCoAServer -p <port number>" << endl;
      exit(0);
    }
    if (argv[1]==string("-t"))
    {
      UsingSocket = false;
      CoCoA5Output = false;
    }
    else if (argv[1]==string("-d"))
    {
      UsingSocket = false;
      CoCoA5Output = true;
    }
    else if (argv[1]==string("-p"))
    {
      if (argc!=3)
      {
        GlobalErrput() << "Error: option should be -p <port number>" << endl;
        exit(1);
      }
      port = atoi(argv[2]);  // BUG: should check that has a sane value!!!
    }
    else
    {
      GlobalErrput() << "Unrecognized option: " << argv[1] << endl;
      GlobalErrput() << "Use option --help for a list of options" << endl;
      exit(1);
    }
  }

  if (UsingSocket)
  {
    CoCoA5Output = true;
    GlobalLogput() << "------[   Starting CoCoAServer on port " << port
                   << " (0x" << hex << port << dec << ")   ]------" << endl << endl;
    //    PrintOperations(GlobalLogput()); // print all operations (with library)
    PrintLibraries(GlobalLogput()); // print all (sub)libraries
    try { SockPtr.reset(new SocketStream(port)); } // this could throw (e.g. socket already in use)
    catch (const CoCoA::ErrorInfo& err)
    {
      GlobalErrput() << endl
                     << "***ERROR*** Failed to start CoCoAServer because " << message(err) << endl
                     << endl
                     << "Perhaps there is a CoCoAServer already running?" << endl;
      exit(1);
    }
    SetGlobalInput(*SockPtr);
    SetGlobalOutput(*SockPtr);
  }

  try
  {
    program();
    if (CoCoA5Output)
    {
      GlobalLogput() << "SERVER HAS FINISHED" << endl;
      //2019      EndOfTransmissionToCoCoA4();
    }
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    GlobalErrput() << "***ERROR***  UNCAUGHT CoCoAError" << endl;
    ANNOUNCE(GlobalErrput(), err);
    if (CoCoA5Output) SendErrorAndQuit(message(err));
  }
  catch (const std::exception& exc)
  {
    GlobalErrput() << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
    if (CoCoA5Output) SendErrorAndQuit(exc.what());
  }
  catch(...)
  {
    GlobalErrput() << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
    if (CoCoA5Output) SendErrorAndQuit("UNKNOWN EXCEPTION");
  }
  return 1;
}


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/CoCoAServer.C,v 1.47 2022/02/23 09:49:23 abbott Exp $
// $Log: CoCoAServer.C,v $
// Revision 1.47  2022/02/23 09:49:23  abbott
// Summary: Removed ResgisterServerOpsDo.H
//
// Revision 1.46  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.45  2021/06/12 13:13:24  abbott
// Summary: Improved handling when given unknown cmd line flags
//
// Revision 1.44  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.43  2019/03/20 16:20:55  bigatti
// -- updates for the new server
//
// Revision 1.42  2019/03/04 11:36:16  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.41  2017/04/18 16:23:58  bigatti
// -- commented out body of ChangeDefaultStatLevel
//
// Revision 1.40  2015/06/26 15:01:18  abbott
// Summary: Revised after changing interface to ErroInfo (now inherits from CoCoA::exception)
// Author: JAA
//
// Revision 1.39  2014/05/15 12:31:34  abbott
// Summary: Now using new files server/GlobalIO.HC (previously in CoCoA/io.H)
// Author: JAA
//
// Revision 1.38  2013/05/27 12:59:05  abbott
// Revised include directives after moving all
// server-related code into src/server/
//
// Revision 1.37  2011/10/04 13:03:03  bigatti
// -- new logo for gui
//
// Revision 1.36  2011/09/21 14:49:36  abbott
// Added missing "else" so that "-t" command line flag is recognized.
//
// Revision 1.35  2010/12/17 16:11:39  abbott
// Consequential change to call to ANNOUNCE.
//
// Revision 1.34  2010/05/28 15:51:27  bigatti
// -- updated copyright date
//
// Revision 1.33  2010/02/03 14:19:00  bigatti
// -- removed RegisterServerOpsFrobby.H (declaration moved into RegisterServerOps.H)
// -- added RegisterServerOpsUser declaration in RegisterServerOps.H
//
// Revision 1.32  2009/11/03 13:43:27  bigatti
// -- added SendErrorAndQuit(s)
// -- added const to strings describing versions
//
// Revision 1.31  2009/11/02 15:42:02  bigatti
// -- using the SkipTag function to check <end_of_session/> from CoCoA-4
// -- Moved ANNOUNCE(err) before writing error to GlobalOutput (in case
//    GlobalOutput, socket, get closed too early)
//
// Revision 1.30  2009/10/30 11:54:34  bigatti
// -- added check for cocoa5.cpkg version for latest change on communication
//
// Revision 1.29  2009/10/29 19:14:21  abbott
// Added EndOfTransmissionToCoCoA4 function.
// Used this new function in CoCoAServer to implement "over-and-out" handshake
// at the end of a CoCoA4-CoCoAServer session.
// New version number to mark this change -- not backward compatible!
//
// Revision 1.28  2009/10/27 15:15:23  bigatti
// -- added "sleep" at the end of "program" for problems on linux
//
// Revision 1.27  2009/10/27 13:48:12  bigatti
// -- added decimal number for the server port
//
// Revision 1.26  2009/05/22 15:37:30  bigatti
// -- changed:  Eof --> char(1)  for end-of-transmission
//
// Revision 1.25  2009/02/11 15:09:39  bigatti
// -- added: check the cocoa5.cpkg version number
//
// Revision 1.24  2009/01/30 17:12:31  bigatti
// -- frobby integration
//
// Revision 1.23  2009/01/26 15:59:54  bigatti
// -- added "-d" for debugging mode
//
// Revision 1.22  2008/11/24 17:10:21  abbott
// Server now explicitly clears an operation object before reading the args
// (this may be useful if an exception occurred).
//
// Revision 1.21  2008/09/22 16:07:02  bigatti
// -- tested (and fixed) communication with cocoa-4 (number of arguments
//    passed to myReadArgs   )
//
// Revision 1.20  2008/09/19 15:02:39  bigatti
// -- new mechanism for passing verbosity level
//
// Revision 1.19  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.18  2007/10/11 16:32:31  bigatti
// -- cleaned up output for benchmarks
//
// Revision 1.17  2007/09/27 10:31:42  bigatti
// -- added function PrintLibraries
//
// Revision 1.16  2007/09/25 16:32:30  abbott
// Several minor changes to silence gcc-4.3:
//    more #includes,
//    and fixed a template problemm in RegisterServerOps.C
//
// Revision 1.15  2007/09/25 16:28:31  bigatti
// -- ServerOp includes infos about the library it is defined in
// -- CoCoAServer no longer prints its own version
// -- CoCoAServer prints all offered operations
// -- CoCoAServer can change stat_level (Max's verbosity) at runtime
//
// Revision 1.14  2007/06/04 12:57:27  abbott
// Improved error message when server fails to create a listening socket.
//
// Revision 1.13  2007/05/18 11:17:25  bigatti
// -- fixed unpredictable behaviour when calling a non-existent function
//    (activated by cocoa5.cpkg version 1.1 officially from CoCoA-4.7.1)
//
// Revision 1.12  2007/04/27 14:54:22  bigatti
// -- content of CoCoAServer.C split into dedicated files
// -- new registration mechanism (through include "RegisterServerOps.H")
//
// Revision 1.11  2007/03/28 09:13:41  bigatti
// -- just a newline before returning error in AssertTag
//
// Revision 1.10  2007/03/27 21:02:23  bigatti
// -- added istream arg to Read* functions
//
// Revision 1.9  2007/03/27 17:02:12  bigatti
// -- removed 1234 from IsTree names
//
// Revision 1.8  2007/03/27 16:22:15  bigatti
// -- updated IsTree functions
// -- removed obsolete server code
//
// Revision 1.7  2007/03/27 14:21:31  bigatti
// -- removed "damaging" semicolons
//
// Revision 1.6  2007/03/27 12:39:57  bigatti
// -- new design applied to Numerical
//
// Revision 1.5  2007/03/21 14:57:32  bigatti
// -- probably the final version for CoCoA-4.7
//
// Revision 1.3  2007/03/20 13:58:42  bigatti
// -- new design for server using inheritanche and map
//
// Revision 1.2  2007/03/09 16:26:20  abbott
// Renamed some functions, and removed some fn arg names because it led
// to linker problems (?!?).
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.33  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.32  2007/03/08 15:28:50  bigatti
// -- more work on numerical
//
// Revision 1.31  2007/03/08 13:07:00  bigatti
// -- updated ParamNames
//
// Revision 1.30  2007/03/08 12:14:10  bigatti
// -- simplified WriteMatrix
// Revision 1.29  2007/03/08 11:45:12  bigatti
// -- cleaner WriteMatrix (was WriteQMatrix)
// Revision 1.28  2007/03/08 11:10:49  bigatti
// -- fixed bug in WriteQMatrix
//
// Revision 1.25  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
// Revision 1.23  2007/03/07 14:59:37  bigatti
// -- renamed complex -> FacetComplex
//
// Revision 1.21  2007/03/07 11:17:18  bigatti
// -- changes in numerical
// -- some CoCoA_ASSERT --> CoCoA_THROW_ERROR
//
// Revision 1.20  2007/03/06 13:38:32  bigatti
// -- cleaned ReadQMatrix (no longer transposed)
// -- added ReadZMatrix
//
// Revision 1.19  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.18  2007/03/05 16:17:11  bigatti
// -- clenup for numerical code (and 3 new error codes)
//
// Revision 1.17  2007/03/04 14:39:49  bigatti
// -- fixed function ReadQMatrix
// -- improved handling of MemPool debugging in ChangeDefaultStatLevel
//
// Revision 1.16  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
// Revision 1.15  2007/03/02 21:47:37  bigatti
// -- added new "foundations"
//
// Revision 1.14  2007/03/01 13:49:51  bigatti
// -- complex --> CoCoA::complex  (should be SimplicialComplex)
//
// Revision 1.13  2007/02/12 18:40:37  bigatti
// -- added numerical operations
//
// Revision 1.12  2007/02/08 22:34:22  cocoa
// Changed BuildInfo: only BuildInfo version is publicly visible.
// Only BuildInfo needs complicated compile-time flags, so several changes
// to Makefiles etc.  Added a new example: ex-BuildInfo.C
//
// Revision 1.11  2007/01/16 14:00:24  bigatti
// -- mainly a cleanup of TmpFGLM and TmpLESystemSolver (by S.Kaspar)
// -- added tests
//
// Revision 1.10  2007/01/08 16:56:17  cocoa
// -- added F5 code (by A.Arri)
//
// Revision 1.9  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.8  2006/11/27 13:10:58  cocoa
// -- fixed: CoCoAServer was a bit too silent at default level...
//
// Revision 1.7  2006/11/27 09:40:29  cocoa
// -- added Kaspar's FGLM and BBasis
//
// Revision 1.6  2006/11/20 14:42:30  cocoa
// -- removed STAT variable
//
// Revision 1.5  2006/10/27 19:02:11  cocoa
// Corrected WritePolyAsList when the poly is zero (added a missing return).
//
// Revision 1.4  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.3  2006/06/13 16:40:55  cocoa
// -- simpler code in CoCoAServer to print length of the result and time
//    (easier for benchmarcks)
//
// Revision 1.2  2006/06/13 10:59:11  cocoa
// -- changed: now the time is printed by the Server (was done by GReductor)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.26  2006/05/22 15:52:16  cocoa
// Added preprocess-disg algorithm to ApproxPts.
// Sundry minor improvements.
//
// Revision 1.25  2006/05/16 09:04:14  cocoa
// -- added IsTree operation on complexes
//
// Revision 1.24  2006/05/12 16:59:23  cocoa
// -- removed error about multiple parameters
//
// Revision 1.23  2006/05/12 13:16:30  cocoa
// Added functions for preprocessing approximate points.
//
// Revision 1.22  2006/05/12 10:49:41  cocoa
// -- changed: ElimVariable --> ElimIndet
// -- changed: using SparsePolyIter in WritePolyAsList
// -- changed: input errors throws instead of exiting
//
// Revision 1.19  2006/05/09 14:07:19  cocoa
// -- new communication style with CoCoA4: MEMORY.PKG.CoCoA5.Result
//    now contains a list of polys or vectors built using Append
//
// Revision 1.18  2006/05/05 16:58:49  cocoa
// -- Time sent to CoCoA4 via "$cocoa5.PrintTime"
//
// Revision 1.16  2006/04/26 09:50:40  cocoa
// -- added GReductor::ourDefaultStatLevel variable to allow CoCoAServer
//    to set statistics level
//
// Revision 1.15  2006/04/21 16:51:03  cocoa
// -- new approach for coefficients in FractionField/GCDDomain
// -- enum for ComputationType and OrderingType (cleaning process
//    started...)
// -- cleaner handling for communication of rings with parameters
// -- new TimeToCoCoA4 function
//
// Revision 1.14  2006/04/11 16:34:29  cocoa
// Completed overhaul of MemPool (incl documentation).
// Modified GMPAllocator so that you can specify the slice size
// in the ctor -- useful for experimentation.
//
// Revision 1.13  2006/03/17 18:17:16  cocoa
// -- changed: use of ReductionCog for reduction (major cleanup)
//
// Revision 1.12  2006/03/07 17:03:03  cocoa
// -- added call to GRI.mySetCoeffRingType(FrFldOfGCDDomain)
//
// Revision 1.10  2006/03/01 14:13:01  cocoa
// -- new management of "Catch" for CoCoA-4 (in "main")
// -- name of CoCoA-4 variable to contain CoCoAServer output is global (VarName4)
//
// Revision 1.9  2006/02/13 13:59:26  cocoa
// -- changes by Max (GRingInfo)
//
// Revision 1.8  2006/02/13 12:08:04  cocoa
// Fixed a problem with some missing assignment ops for certain PPMonoidElems.
// Fixed a bug in RingDistrMPoly::myIndetPower.
//
// Revision 1.7  2006/01/18 15:58:20  cocoa
// -- new changes by Max
// Revision 1.6  2006/01/17 15:44:56  cocoa
// -- changes by Max for operations with modules
//
// Revision 1.5  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.3  2005/10/26 14:38:33  cocoa
// No longer leaves GlobalLogput() in hexadecimal mode.
// All logging output now on GlobalLogput() instead of cerr.
//
// Revision 1.2  2005/10/18 12:06:36  cocoa
// Cleaned Makefiles, and fixed them so they should work wherever
// CoCoALib is unpacked.
//
// Replaced VERSION cpp macro with COCOA_VERSION.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.6  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.5  2005/10/12 13:17:14  cocoa
// -- set default: IsPosTo = false
//
// Revision 1.4  2005/10/07 17:22:39  cocoa
// -- set default "STAT" to 0
// -- changed the way the timing is printed by cocoa-4 (using $cocoa4.PrintTime)
//
// Revision 1.3  2005/09/28 11:50:34  cocoa
// -- new code for graded modules
//
// Revision 1.2  2005/05/05 16:31:47  cocoa
// -- added: "Greetings" session in communication
// -- changed: improved communication with CoCoA-4 for input errors
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.3  2005/05/02 13:46:04  cocoa
// -- new syntax: CoCoAServer [-h , -t , -p <port>]
// Revision 1.2  2005/04/29 14:07:15  cocoa
// -- "CoCoAServer 0" will run the server on port c0c0
//    "CoCoAServer n" will run the server on port n
//    "CoCoAServer" will run the server on stdin/stdout
//
// Revision 1.4  2005/03/10 16:43:22  cocoa
// -- updated: use of RefRingElem and ConstRefRingElem
// -- changed: iterator via myMoveLM and myAppendClear in WritePolyAsList
//             to be fixed when we have proper iterators
//
// Revision 1.3  2005/02/22 16:38:57  cocoa
// -- changed: macro (USE_SOCKET, defined) --> (DONT_USE_SOCKET, undefined)
//
// Revision 1.2  2005/02/09 11:16:59  cocoa
// -- easier switch between socket-file I/O
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.20  2004/11/29 16:53:50  cocoa
// -- change: PolyRing --> AsPolyRing
// -- removed definition of adjoint matrix (no longer needed)
//
// Revision 1.19  2004/11/11 15:33:37  cocoa
// -- changed: unsigned int --> size_t
// -- changed: LT --> LPP
// -- some function cleaning and tidying
//
// Revision 1.18  2004/11/09 16:05:18  cocoa
// -- changed: input functions take input from GlobalInput()
// -- changed: "out" argument is now first in all ouput functions
//
// Revision 1.17  2004/11/05 16:49:26  cocoa
// -- removed code for order matrices modulo 32749 (computed internally)
// -- cleaned up code for choosing socket or file I/O
//
// Revision 1.16  2004/11/03 18:23:31  cocoa
// -- started cleaning of ordering matrix code
//
// Revision 1.15  2004/11/02 15:15:24  cocoa
// -- changes for new PPOrdering syntax (matrix ordering not tested yet)
//
// Revision 1.13  2004/07/28 14:00:33  cocoa
// -- new syntax for parametric computations
//
// Revision 1.12  2004/07/27 16:15:38  cocoa
// -- added  SkipTag
//
// Revision 1.10  2004/07/27 14:08:33  cocoa
// -- changed "groebner_basis" into "groebner" (more OpenMath-like)
// -- changed line about socket port
//
// Revision 1.9  2004/07/23 13:37:43  cocoa
// -- tags are now more openMath-like
// -- checks on tags validity
//
// Revision 1.8  2004/07/23 09:09:00  cocoa
// -- for characteristic<60000  use  RingDistrMPolyInlFpPP
//    otherwise                      RingDistrMPolyInlPP
// Revision 1.7  2004/07/21 16:39:52  cocoa
// -- call for using RingDistrMPolyInlFpPP
//
// Revision 1.6  2004/07/16 10:52:20  cocoa
// -- more refined code for communication with cocoa-4:
//    .all output sent through GlobalLogput() or GlobalOutput()
//    .better error handling (assigned to result variable)
// -- all char* changed into string
//
// Revision 1.5  2004/06/30 16:47:11  cocoa
// Changed port number to 49344 (which is C0C0 in hexdecimal :-).
//
// Revision 1.4  2004/06/30 15:49:46  cocoa
// -- new communication protocol with cocoa-4
//
/*
  Changed MakeCoeffRing - now it takes both Characteristic and FloatPrecision.
  If FloatPrecision<>0, this is float. If not, is Z, Zp according to
  Characteristic. In the input file, only one of this two fields can be present
*/


// Revision 1.13  2004/06/09 12:51:43  cocoa
// -- change: ifstream/ofstream arguments into istream/ostream
// -- some cleaning
//
// Revision 1.12  2004/05/28 13:17:35  cocoa
// -- ReadPolyList: changed "unsigned int" into "size_t" and "int" to fix
//    input bug found by the new compiler on vector
// -- removed "ANNA_DEBUG" variable
//
// Revision 1.11  2004/05/27 17:00:17  cocoa
// -- slight change for testing RingDistrMPoly
//
// Revision 1.10  2004/05/25 10:33:39  cocoa
// -- corrected calls to NewPolyRing --> NewPolyRing_DMPI
//
// Revision 1.9  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.8  2004/04/02 16:24:36  cocoa
// -- input functions now have ifstream as argument instead of string
//    for file name  [ReadPolyList is much cleaner]
// -- all "FromFile" in function names have been removed
// -- same for output functions
// -- MakeCoeffRing: negative characteristic is interpreted as
//    precision for RingFloat, e.g. -128 will call NewRingFloat(128)
//
// Revision 1.7  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.6  2004/03/11 13:05:52  cocoa
// -- new: DenseMatrixOrderingMod46337 which uses the inverse of the
//    ordering matrix mod 46337 (returning exponents < 46337)
//
// Revision 1.5  2004/03/04 11:39:30  cocoa
// -- updated code for Borel reductors:
//    Borel reductors are used if the variable UBR is set to UseBorel
//
// Revision 1.4  2003/11/21 15:41:42  cocoa
// -- added function MakeCoeffRing
// -- added possible choice of DistrMPoly in MakePolyRing
//
// Revision 1.2  2003/09/25 13:49:38  cocoa
// - working syntax (at CoCoALib-0.0.0 time)
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.3  2003/05/29 17:11:18  bigatti
// - new code for orderings, modules, parameters
//
