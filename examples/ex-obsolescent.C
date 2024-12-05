// Copyright (c) 2016  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
// *** BELOW: EXTRA INCLUDE DIRECTIVE TO MAKE OBSOLESCENT FNS VISIBLE ***
#include "CoCoA/obsolescent.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows how to call an obsolescent CoCoALib function. \n";

const string LongDescription =
  "This example program shows how to call an obsolescent CoCoALib function. \n"
  "You must do two things: \n"
  "(1) include the header file CoCoA/obsolescent.H \n"
  "(2) give the option AllowObsolescentFns to GlobalManager";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    // *** SEE EXTRA INCLUDE DIRECTIVE AT TOP OF FILE ***
    GlobalManager CoCoAFoundations(AllowObsolescentFns); // *** GIVE OPTION TO GLOBALMANAGER ***
    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    PPMonoid PPM = NewPPMonoidOv(symbols("x,y"), lex);
    PPMonoidElem t = indet(PPM, 0);
    CoCoA_ASSERT_ALWAYS(IsRadical(t)); // the call to IsRadical will print a logging message.
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
