// Copyright (c) 2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows the use of DivMasks.  \n"
  "See also the example program ex-PPWithMask1.C    \n";

const string LongDescription =
  "We show how to use DivMasks directly -- it is a bit tedious!             \n"
  "See also the example ex-PPWithMask1.C, which has a friendlier interface. \n"
  "We show how the various different DivMaskRules can be selected.          \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void DoSomeOps(const DivMaskRule& DMR)
  {
    cout << "Doing some computations with DivMaskRule=" << DMR << endl;

    // Create a PPMonoid for indets called x,y,z,t.
    PPMonoid PPM = NewPPMonoid(symbols("x,y,z,t"), StdDegRevLex);
    const PPMonoidElem x = indet(PPM, 0);
    const PPMonoidElem y = indet(PPM, 1);
    const PPMonoidElem z = indet(PPM, 2);
    const PPMonoidElem t = indet(PPM, 3);

    // Fill array pp with some PPs whose divisibility we'll test.
    vector<PPMonoidElem> pp;
    pp.push_back(one(PPM));
    pp.push_back(x*y);
    pp.push_back(z*t);
    pp.push_back(x*y*z*t);
    pp.push_back(power(x,4)*power(y,3)*power(z,2)*t);
    pp.push_back(x*power(y,2)*power(z,3)*power(t,4));
    pp.push_back(power(x*y*z,5));
    pp.push_back(power(x*y*z,10));

    // Now put into the array dm the DivMask of each PP.
    const int NumInds = NumIndets(PPM);
    const int NumPPs = pp.size();
    vector<SmallExponent_t> exps(NumInds);
    vector<long> expl(NumInds);
    vector<DivMask> dm(NumPPs);
    for (int i=0; i < NumPPs; ++i)
    {
      exponents(expl, pp[i]);
      for (int k=0; k < NumInds; ++k) exps[k] = expl[k];
      DMR->myAssignFromExpv(dm[i], &exps[0], NumInds);
      cout << dm[i] << " is DivMask for " << pp[i] << endl;
    }

    // Print out divisibility table:
    //   N -- definitely not a factor.
    //   Y -- DivMask says might be a factor, & it really is.
    //   ! -- DivMask says might be a factor, BUT it is not.
    //   = -- the two PPs are equal (just to help read the table)
    cout << "Divisibility table according to DivMask:" << endl;
    for (int i=0; i < NumPPs; ++i)
    {
      for (int j=0; j < NumPPs; ++j)
      {
        if (i == j) { cout << "= "; continue; }
        if (!IsSubset(dm[i], dm[j])) { cout << "N "; continue; }
        if (IsDivisible(pp[j], pp[i])) { cout << "Y "; continue; }
        cout << "! ";
      }
      cout << "     " << pp[i] << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-DivMask1.C,v 1.6 2022/02/13 09:56:56 abbott Exp $
// $Log: ex-DivMask1.C,v $
// Revision 1.6  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2015/07/28 11:54:18  bigatti
// -- new syntax for indets in NewPolyRing
//
// Revision 1.4  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.3  2010/12/17 16:07:55  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2010/02/03 16:10:22  abbott
// Changed to more compact PPMonoid pseudo ctor
//
// Revision 1.1  2010/02/01 22:36:59  abbott
// Added new examples for DivMasks and PPWithMask.
//
