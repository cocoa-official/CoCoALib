// Copyright (c) 2010  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program performs some divisibility speed tests   \n"
  "using PPWithMasks; it shows the difference in speed which can \n"
  "be achieved using various DivMaskRules.                       \n";

const string LongDescription =
  "This example program test the speed of the divisibility test  \n"
  "on values of type PPwithMask.  The main aim is to show that   \n"
  "different DivMaskRules can produces differing behaviour, and  \n"
  "which rule is best depends on the problem (e.g. few indets    \n"
  "and high degrees, or many indets and low degrees).  Here we   \n"
  "see also that unused indets can affect the speed.             \n"
  "In case you're interested, the program ex-DivMask2 measures   \n"
  "the effectiveness of the various DivMaskRules.                \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void TimingTest(PPMonoidElem x1, PPMonoidElem x2, PPMonoidElem x3, PPMonoidElem x4, DivMaskRule DMR)
  {
    // Fill vector ppwm with some power products.
    const int MaxExp = 6;
    vector<PPWithMask> ppwm;
    for (int e1=0; e1 <= MaxExp; ++e1)
      for (int e2=0; e2 <= MaxExp; ++e2)
        for (int e3=0; e3 <= MaxExp; ++e3)
          for (int e4=0; e4 <= MaxExp; ++e4)
          {
            ppwm.push_back(PPWithMask(power(x1,e1)*power(x2,e2)*power(x3,e3)*power(x4,e4), DMR));
          }

    const long NumPPs = len(ppwm);

    // Consider all pairs of elements of ppwm...
    // count how many are divisible, and how many not.
    const double StartTime = CpuTime();
    long CountDiv = 0;
    long CountNotDiv = 0;
    for (long i=0; i < NumPPs; ++i)
      for (long j=0; j < NumPPs; ++j)
      {
        if (IsDivisibleFast(ppwm[j], ppwm[i]))
          ++CountDiv;
        else
          ++CountNotDiv;
      }

    cout << "[" << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "]"
         << "  Time = " << CpuTime()-StartTime
         << "   #divisible = " << CountDiv
         << "   #not-divisible = " << CountNotDiv << endl;

  }


  void DoTimingTests(const DivMaskRule& DMR)
  {
    cout << "Doing timing tests with DivMaskRule=" << DMR << endl;

    {
      PPMonoid PPM_4 = NewPPMonoidEv(SymbolRange("x",0,3), StdDegRevLex);

      // First timing test in PPMonoid with 4 indets...
      TimingTest(indet(PPM_4, 0),
                 indet(PPM_4, 1),
                 indet(PPM_4, 2),
                 indet(PPM_4, 3),
                 DMR);

      PPMonoid PPM_10 = NewPPMonoidEv(SymbolRange("y",0,9), StdDegRevLex);

      // Second timing test in PPMonoid with 10 indets, using the first 4 of them...
      TimingTest(indet(PPM_10, 0),
                 indet(PPM_10, 1),
                 indet(PPM_10, 2),
                 indet(PPM_10, 3),
                 DMR);

      // Third timing test in PPMonoid with 10 indets, using the last 4 of them...
      TimingTest(indet(PPM_10, 6),
                 indet(PPM_10, 7),
                 indet(PPM_10, 8),
                 indet(PPM_10, 9),
                 DMR);


      PPMonoid PPM_50 = NewPPMonoidEv(SymbolRange("z",0,49), StdDegRevLex);

      // Fourth timing test in PPMonoid with 50 indets, using the first 4 of them...
      TimingTest(indet(PPM_50, 0),
                 indet(PPM_50, 1),
                 indet(PPM_50, 2),
                 indet(PPM_50, 3),
                 DMR);

      // Fifth timing test in PPMonoid with 50 indets, using the last 4 of them...
      TimingTest(indet(PPM_50, 46),
                 indet(PPM_50, 47),
                 indet(PPM_50, 48),
                 indet(PPM_50, 49),
                 DMR);

    }


    cout << "---------------------------------" << endl
         << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    DoTimingTests(NewDivMaskNull());
    DoTimingTests(NewDivMaskSingleBit());
    DoTimingTests(NewDivMaskSingleBitWrap());
    DoTimingTests(NewDivMaskEvenPowers());
    DoTimingTests(NewDivMaskHashing());
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PPWithMask2.C,v 1.9 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-PPWithMask2.C,v $
// Revision 1.9  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.7  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.6  2014/07/09 14:21:29  abbott
// Summary: Reduced size of test to make it faster
// Author: JAA
//
// Revision 1.5  2012/06/29 15:13:11  abbott
// Improved alignment.
//
// Revision 1.4  2012/05/04 20:01:19  abbott
// Minor modifications to reduce execution time (usu. by reducing number of iterations).
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2010/03/18 18:03:32  bigatti
// -- improved syntax for PPOrdering argument
//
// Revision 1.1  2010/02/01 22:36:59  abbott
// Added new examples for DivMasks and PPWithMask.
//
