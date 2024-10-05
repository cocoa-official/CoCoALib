// Copyright (c) 2011  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program makes simple use of the \"floating point\" approximation \n"
  "functions: MantissaAndExponent10 and MantissaAndExponent2, FloatApprox.";

const string LongDescription =
  "Example of use of MantissaAndExponent10 and related MantExp10 structure. \n"
  "Example of use of MantissaAndExponent2 and related MantExp2 structure.   \n"
  "Example of use of FloatApprox (similatr to MantissaAndExponent2).        \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    const BigRat v = -(1-BigRat(5,100000)); // -0.99995
    cout << endl << "True rational value is " << v << endl;

    cout << "\nDecimal approx:" << endl;
    for (int sigfig = 2; sigfig <= 6; ++sigfig)
    {
      const MantExp10 ME = MantissaAndExponent10(v, sigfig);
      cout << "With " << sigfig << " significant figures: " << ME << endl;
    }

    cout << "\nBinary approx:" << endl;
    for (int sigbits = 2; sigbits <= 20; sigbits += 4)
    {
      const MantExp2 ME = MantissaAndExponent2(v, sigbits);
      cout << "With " << sigbits << " significant bits: " << ME << endl;
    }

    cout << "\nBinary approx using FloatApprox:" << endl;
    for (int sigbits = 2; sigbits <= 20; sigbits += 4)
    {
      cout << "With " << sigbits << " bits value is: " << FloatApprox(v, sigbits) << endl;
    }
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
