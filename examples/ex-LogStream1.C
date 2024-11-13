// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <fstream>
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a short example showing how to use LogStream, and how to  \n"
  "make it refer to a different output stream.                       \n";

const string LongDescription =
  "This is a short example showing how to use LogStream, and how to  \n"
  "make it refer to a different output stream.                       \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Uncomment next line to tell CoCoALib to log messages on std::clog:
//    LogStreamForThisBlock MainProgram(std::clog);
    LogStream() << "Output on default CoCoA::LogStream, namely std::cout" << endl;

    {
      // Use a different logging file in this block...
      ofstream BitBucket("/dev/null");
      LogStreamForThisBlock InnerBlock(BitBucket);
      LogStream() << "Message sent to BitBucket" << endl;
    }
    
    LogStream() << "Outside inner block; another message to CoCoA::LogStream" << endl;
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
