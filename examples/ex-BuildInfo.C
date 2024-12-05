// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a very short example showing what you can do with BuildInfo. \n";

const string LongDescription =
  "BuildInfo::PrintAll gives important information in the (enormously \n"
  "unlikely :-) event that you need to report a bug in CoCoALib. \n";
//----------------------------------------------------------------------


namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // You can get the library version number as a C++ string.
    // The string has the format  "A.bcde"  where "A" is the major version,
    // "bc" the minor version, and "de" the patch level.
    const string& version = BuildInfo::version();
    cout << "The version of CoCoALib used to produce this executable was "
         << version << endl
         << endl;

    // You can print out a message containing all the build information like this:
    cout << "Here is a summary of all BuildInfo for the library:";
    BuildInfo::PrintAll(cout);
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
