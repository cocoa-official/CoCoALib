// Copyright (c)  2009  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to interpret the result of a factorization. \n";

const string LongDescription =
  "This example shows how to interpret the result of a factorization. \n"
  "It creates a ring element (belonging to a polynomial ring), and    \n"
  "factorizes it.  The result is a \"factorization object\".  We show \n"
  "how to access/use the various fields in this object.               \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    SparsePolyRing P = NewPolyRing(RingQQ(), symbols("x,y")); // QQ[x,y];
    RingElem f = RingElem(P, "x^96 - y^96");

    const factorization<RingElem> FacInfo = factor(f);

    // These are convenient aliases for the 3 fields in the factorization:
    const RingElem&         content   = FacInfo.myRemainingFactor();
    const vector<RingElem>& IrredFacs = FacInfo.myFactors();
    const vector<long>&     mult      = FacInfo.myMultiplicities();
    
    // Print out the factorization in a "nice" way:
    cout << "The factors of " << f << " are:" << endl;
    if (!IsOne(content))
      cout << "content: " << content << endl;
    const int NumIrredFacs = len(IrredFacs);
    for (int i = 0; i != NumIrredFacs; ++i)
    {
      cout << IrredFacs[i] << "  with multiplicity  " << mult[i] << endl;
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-factor1.C,v 1.12 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-factor1.C,v $
// Revision 1.12  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.11  2021/10/17 18:34:19  abbott
// Summary: Corrected typo in comment
//
// Revision 1.10  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.9  2016/07/20 08:43:29  abbott
// Summary: Minor improvements to comments
//
// Revision 1.8  2015/07/27 11:50:49  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.7  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.6  2015/04/27 13:04:17  bigatti
// Summary: using ReadExpr for input poly
//
// Revision 1.5  2014/03/24 12:09:20  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.4  2012/10/05 09:29:43  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.3  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2009/09/23 14:09:36  abbott
// First example for (polynomial) factorization.
//
//
