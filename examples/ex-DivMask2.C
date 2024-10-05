// Copyright (c) 2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program performs some divisibility tests on PPs     \n"
  "using DivMasks; it measures the effectiveness of the various     \n"
  "DivMaskRules.  See ex-PPWithMask2 for some similar timing tests. \n";

const string LongDescription =
  "This example program illustrates the effectiveness of the various    \n"
  "DivMaskRules when using DivMasks as a quick non-divisibility check.  \n"
  "Using DivMasks to check divisibility we can obtain two answers:      \n"
  "'Surely not' or 'Maybe'.  The intention here is to see how often     \n"
  "'Maybe' really means 'Yes, divisible'.  We see that the effectiveness\n"
  "depends also on the presence of unused indeterminate in the PPMonoid.\n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void EffectivenessTest(PPMonoidElem x1, PPMonoidElem x2, PPMonoidElem x3, PPMonoidElem x4, DivMaskRule DMR)
  {
    // Fill array with PPs we shall use for divisibility testing
    const int MaxExp = 6;
    vector<PPMonoidElem> pp;
    for (int e1=0; e1 <= MaxExp; ++e1)
      for (int e2=0; e2 <= MaxExp; ++e2)
        for (int e3=0; e3 <= MaxExp; ++e3)
          for (int e4=0; e4 <= MaxExp; ++e4)
          {
            pp.push_back(power(x1,e1)*power(x2,e2)*power(x3,e3)*power(x4,e4));
          }

    const int NumPPs = pp.size();

    // Now put into array dm the DivMask of each PP.
    // This part is ugly & laborious.
    const int Nvars = NumIndets(owner(x1));
    vector<SmallExponent_t> exps(Nvars);
    vector<long> expl(Nvars);
    vector<DivMask> dm(NumPPs);
    for (int i=0; i < NumPPs; ++i)
    {
      exponents(expl, pp[i]);
      for (int k=0; k < Nvars; ++k) exps[k] = expl[k];
      DMR->myAssignFromExpv(dm[i], &exps[0], Nvars);
    }

    // Try all divisibility tests, and count when DivMask misleads.
    int CountDiv = 0;
    int CountNotDiv = 0;
    long CountMistakes = 0;
    for (int i=0; i < NumPPs; ++i)
      for (int j=0; j < NumPPs; ++j)
      {
        if (!IsSubset(dm[i], dm[j]))
        {
          ++CountNotDiv;
          continue;
        }
        if (IsDivisible(pp[j], pp[i]))
          ++CountDiv;
        else
          ++CountMistakes;
      }

    cout << "[" << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "]"
//                  << "   #divisible = " << CountDiv
//                  << "   #not-divisible = " << CountNotDiv
//                  << "   #mistakes = " << CountMistakes
         << "   effectiveness=" << double(CountNotDiv)/(CountNotDiv+CountMistakes) << endl;

  }


  void DoEffectivenessTests(const DivMaskRule& DMR)
  {
    cout << "Doing effectiveness tests with DivMaskRule=" << DMR << endl;

    {
      //    PPOrdering PPO_4 = NewStdDegRevLexOrdering(4);
      //    PPMonoid PPM_4 = NewPPMonoidEv(SymbolRange("x",0,3), PPO_4);
      PPMonoid PPM_4 = NewPPMonoidEv(SymbolRange("x",0,3), StdDegRevLex);

      // First test in PPMonoid with 4 indets...
      EffectivenessTest(indet(PPM_4,0),
                        indet(PPM_4,1),
                        indet(PPM_4,2),
                        indet(PPM_4,3),
                        DMR);

      //-------------------------------------------------------
      //    PPOrdering PPO_10 = NewStdDegRevLexOrdering(10);
      PPMonoid PPM_10 = NewPPMonoidEv(SymbolRange("y",0,9), StdDegRevLex);

      // Second test in PPMonoid with 10 indets, using the first 4 of them...
      EffectivenessTest(indet(PPM_10,0),
                        indet(PPM_10,1),
                        indet(PPM_10,2),
                        indet(PPM_10,3),
                        DMR);

      // Third test in PPMonoid with 10 indets, using the last 4 of them...
      EffectivenessTest(indet(PPM_10,6),
                        indet(PPM_10,7),
                        indet(PPM_10,8),
                        indet(PPM_10,9),
                        DMR);

      //-------------------------------------------------------
      //    PPOrdering PPO_50 = NewStdDegRevLexOrdering(50);
      PPMonoid PPM_50 = NewPPMonoidEv(SymbolRange("z",0,49), StdDegRevLex);

      // Fourth test in PPMonoid with 50 indets, using the first 4 of them...
      EffectivenessTest(indet(PPM_50,0),
                        indet(PPM_50,1),
                        indet(PPM_50,2),
                        indet(PPM_50,3),
                        DMR);

      // Fifth test in PPMonoid with 50 indets, using the last 4 of them...
      EffectivenessTest(indet(PPM_50,46),
                        indet(PPM_50,47),
                        indet(PPM_50,48),
                        indet(PPM_50,49),
                        DMR);
    }


    cout << "---------------------------------" << endl
         << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    DoEffectivenessTests(NewDivMaskNull());
    DoEffectivenessTests(NewDivMaskSingleBit());
    DoEffectivenessTests(NewDivMaskSingleBitWrap());
    DoEffectivenessTests(NewDivMaskEvenPowers());
    DoEffectivenessTests(NewDivMaskHashing());
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-DivMask2.C,v 1.7 2022/02/13 09:56:56 abbott Exp $
// $Log: ex-DivMask2.C,v $
// Revision 1.7  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2015/11/04 10:14:20  abbott
// Summary: Made example faster
//
// Revision 1.5  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.4  2012/05/04 20:01:19  abbott
// Minor modifications to reduce execution time (usu. by reducing number of iterations).
//
// Revision 1.3  2010/12/17 16:07:55  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2010/03/18 18:06:11  bigatti
// -- improved syntax for PPOrdering argument
//
// Revision 1.1  2010/02/01 22:36:59  abbott
// Added new examples for DivMasks and PPWithMask.
//
