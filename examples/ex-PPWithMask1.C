// Copyright (c) 2010  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows the use of PPWithMasks for \n"
  "testing divisibility between PPs.                     \n"
  "Compare with ex-DivMask1.C.                           \n";

const string LongDescription =
  "We show how to use PPWithMasks for testing divisibility.      \n"
  "This program is very similar to ex-DivMask1.C.  The main      \n"
  "differences are that this program is shorter and clearer, and \n"
  "it does a proper job of testing divisibility (even when the   \n"
  " useless null rule is used).                                  \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void DoSomeOps(const DivMaskRule& DMR)
  {
    cout << "Doing some computations with DivMaskRule=" << DMR << endl;

    // Create a PPMonoid for 4 indets called x,y,z,t -- a bit laborious!
    //  const int NumIndets = 4;
    //  PPOrdering PPO = NewStdDegRevLexOrdering(NumIndets);
    PPMonoid PPM = NewPPMonoidEv(symbols("x,y,z,t"), StdDegRevLex);
    const PPMonoidElem x = indet(PPM, 0);
    const PPMonoidElem y = indet(PPM, 1);
    const PPMonoidElem z = indet(PPM, 2);
    const PPMonoidElem t = indet(PPM, 3);

    // Fill array pp with some PPs whose divisibility we'll test.
    vector<PPWithMask> ppwm;
    ppwm.push_back(PPWithMask(one(PPM), DMR));
    ppwm.push_back(PPWithMask(x*y, DMR));
    ppwm.push_back(PPWithMask(z*t, DMR));
    ppwm.push_back(PPWithMask(x*y*z*t, DMR));
    ppwm.push_back(PPWithMask(power(x,4)*power(y,3)*power(z,2)*t,DMR));
    ppwm.push_back(PPWithMask(x*power(y,2)*power(z,3)*power(t,4),DMR));
    ppwm.push_back(PPWithMask(power(x*y*z,5), DMR));
    ppwm.push_back(PPWithMask(power(x*y*z,10), DMR));
    const int NumPPs = ppwm.size();

    // Print out divisibility table:
    //   N -- not a factor.
    //   Y -- is a factor.
    //   = -- the two PPs are equal (just to help read the table)
    cout << "Divisibility table according to DivMask:" << endl;
    for (int i=0; i < NumPPs; ++i)
    {
      for (int j=0; j < NumPPs; ++j)
      {
        if (i == j) { cout << "= "; continue; }
        if (IsDivisibleFast(ppwm[j], ppwm[i]))
          cout << "Y "; 
        else
          cout << "N ";
      }
      cout << "     " << PP(ppwm[i]) << endl;
    }

    cout << "---------------------------------" << endl
         << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    DoSomeOps(NewDivMaskNull());
    DoSomeOps(NewDivMaskSingleBit());
    DoSomeOps(NewDivMaskSingleBitWrap());
    DoSomeOps(NewDivMaskEvenPowers());
    DoSomeOps(NewDivMaskHashing());
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PPWithMask1.C,v 1.7 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-PPWithMask1.C,v $
// Revision 1.7  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.5  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.4  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2010/03/18 18:02:54  bigatti
// -- improved syntax for PPOrdering argument
//
// Revision 1.1  2010/02/01 22:36:59  abbott
// Added new examples for DivMasks and PPWithMask.
//
