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
