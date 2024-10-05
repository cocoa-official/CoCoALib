// Copyright (c) 2005 John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing how to iterate through \"sparse\" polynomials.\n";

const string LongDescription =
  "Unlike in CoCoA-5, there no functions in CoCoALib to extract the      \n"
  "support, the coeff list or the monomials from a polynomial.  Instead  \n"
  "CoCoALib offers an \"iterator\" over the terms in a sparse polynomial:\n"
  "you can access the coeff and PP directly from the iterator.           \n";

// ----------------------------------------------------------------------

namespace CoCoA
{

  void SimpleDemo(const SparsePolyRing& P)
  {
    if (NumIndets(P) < 4) return;       // Do nothing if there are fewer than 4 indets...
    if (!IsInvertible(RingElem(P,3))) return; // ... or if 3 is not invertible.

    const vector<RingElem>& x = indets(P);
    RingElem f = 2*x[0] + 4*x[1]*x[1] - 2*x[0]*x[3]/3;
    cout << "Our poly f = " << f << "   element of " << owner(f) << endl;

    // ---- ITERATION over f: we just print out the coeffs and the PPs ----
    cout << "The terms in f are as follows\n";
    for (SparsePolyIter iter=BeginIter(f); !IsEnded(iter); ++iter)
    {
      cout << "coeff: " << coeff(iter) << "\t   element of " << owner(coeff(iter)) << endl
           << "PP: " << PP(iter) << "\t   element of " << owner(PP(iter)) << endl
           << endl;
    }
    cout << "Polynomial iterators are read-only, so f is unchanged: " << f << endl
         << "---------------------------------" << endl << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring Fp = NewZZmod(32003);          // Coefficient ring.
    // Now create 3 poly rings over Fp; they are all isomorphic to ZZ/(32003)[a,b,c],
    // but specify different internal implementations (and use different indet names).
    SparsePolyRing Fpx = NewPolyRing(Fp, SymbolRange("x",0,3));
    SparsePolyRing Fpy = NewPolyRing_DMPI(Fp, SymbolRange("y",0,3));
    SparsePolyRing Fpz = NewPolyRing_DMPII(Fp, SymbolRange("z",0,3));

    // Now run the demo in each of these rings:
    SimpleDemo(Fpx);
    SimpleDemo(Fpy);
    SimpleDemo(Fpz);
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PolyIterator1.C,v 1.10 2022/03/17 14:37:17 abbott Exp $
// $Log: ex-PolyIterator1.C,v $
// Revision 1.10  2022/03/17 14:37:17  abbott
// Summary: Added long descr; costmetic improvements
//
// Revision 1.9  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2016/11/18 18:03:28  abbott
// Summary: Changed name of loop variable (from i to iter) which is of type SparsePolyIter
//
// Revision 1.7  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.6  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.5  2012/02/08 17:41:29  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2008/10/07 12:06:12  abbott
// Removed useless #include.
//
// Revision 1.2  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.6  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/26 15:50:25  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 15:45:15  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.4  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.3  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1  2005/05/04 16:35:06  cocoa
// -- new examples
//
