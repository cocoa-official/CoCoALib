// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example program showing how to see which external  \n"
  "libraries have been compiled into CoCoALib.                   \n";

const string LongDescription =
  "This is an example program showing how to see which external  \n"
  "libraries have been compiled into CoCoALib.                   \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    const int NumLibs = len(ExternalLibs());

    cout << "A short summary of external libs compiled into CoCoALib:" << endl;
    for (int i=0; i < NumLibs; ++i)
      cout << ExternalLibs()[i].myName << "(v." << ExternalLibs()[i].myVersion << ")" << endl;


    cout << endl;
    cout << "All information about external libs compiled into CoCoALib:" << endl;
    for (int i=0; i < NumLibs; ++i)
      cout << ExternalLibs()[i] << endl;
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
