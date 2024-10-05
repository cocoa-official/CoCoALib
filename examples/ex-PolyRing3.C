// Copyright (c) 2013  John Abbott, Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows how to extract the \"coefficients\" \n"
  "of a polynomial, and also how to compute some other operations \n"
  "on the coefficients                                            \n";

const string LongDescription =
  "This example program shows how to extract the \"coefficients\" \n"
  "of a polynomial, and also how to compute some other operations \n"
  "on the coefficients (e.g. content).  CoCoALib offers two notions\n"
  "of coefficient: one is the natural one dictated by the ring,   \n"
  "and the other comes from identifying e.g. k[x,y,z] and k[x,z][y]\n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    SparsePolyRing ZZxy = NewPolyRing(RingZZ(), symbols("x,y"));
    const RingElem x = indet(ZZxy,0);
    const RingElem y = indet(ZZxy,1);

    RingElem f = 2*x*x*y - 4*y*y + 6*x*x + 36;

    cout << "In the following we have f = " << f << endl << endl;

    cout << "Product of indets actually in f is " << IndetsProd(f) << endl << endl;

    // Accessing coeffs via SparsePolyIter:
    cout << "Using a SparsePolyIter we decompose f as follows:" << endl;
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      cout << "coeff = " << coeff(it) << "  in ring " << owner(coeff(it)) << endl;
      cout << "PP    = " << PP(it)    << "  in " << owner(PP(it)) << endl;
      cout << endl;
    }

    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << endl;

    // Regard f as a poly in just "x" or just "y" we obtain:
    cout << "Coefficients with respect to certain indeterminates:" << endl << endl;
    cout << "Case (1) considering f as a *univariate* polynomial in..." << endl;
    cout << "...the indeterminate x the coeffs are " << CoeffVecWRT(f, x) << endl;
    cout << "...the indeterminate y the coeffs are " << CoeffVecWRT(f, y) << endl;
    cout << endl;

    cout << "Case (2) considering f as a sparse multivariate polynomial in..." << endl;
    cout << "...the indet x its structure is " << CoefficientsWRT(f, x) << endl;
    cout << "...the indet y its structure is " << CoefficientsWRT(f, y) << endl;
    vector<long> XandY; XandY.push_back(0); XandY.push_back(1);
    cout << "...the indets x & y its structure is " << CoefficientsWRT(f, XandY) << endl;
    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "The content of a polynomial" << endl << endl;

    // Content of f
    RingElem ContF = content(f);
    cout << "The \"numerical\" content of f is " << ContF << "  -- an element of ring " << owner(ContF) << endl;
    cout << endl;
    RingElem ContWRTx = ContentWRT(f, x);
    cout << "Content WRT x is " << ContWRTx << "  -- element of " << owner(ContWRTx) << endl;

    RingElem ContWRTy = ContentWRT(f, y);
    cout << "Content WRT y is " << ContWRTy << "  -- element of " << owner(ContWRTy) << endl;

    RingElem ContWRTxy = ContentWRT(f, XandY);
    cout << "Content WRT x & y is " << ContWRTxy << "  -- element of " << owner(ContWRTxy) << endl;

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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PolyRing3.C,v 1.6 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-PolyRing3.C,v $
// Revision 1.6  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2020/02/13 14:16:36  abbott
// Summary: Added example of use of (new fn) IndetsProd
//
// Revision 1.4  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.3  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.2  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.1  2013/05/28 13:29:00  abbott
// New example for accessing coeffs of a poly.
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
