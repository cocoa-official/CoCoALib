// Copyright (c) 2005  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example showing operations on RingElem for a ring or a PolyRing.\n";

const string LongDescription =
  "This is a long list of function calls from different rings.  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Simple arithmetic operations, and equality test.
  void RingElemArith(RingElem a, const RingElem& b)
  {
    const ring& R = owner(a);
    if (owner(b) != R)
      CoCoA_THROW_ERROR("RingElems should belong to same ring", "TestRing");

    cout << endl << " --RingElemArith(a,b)-- ";
    cout << "a = " << a << ";  b = " << b << endl;

    cout << "a == b gives " << (a == b) << endl;
    cout << "a != b gives " << (a != b) << endl;
    cout << "a * b  gives " << a*b << endl;
    cout << "-a     gives " << -a << endl;
    cout << "a += b gives a = " << (a += b) << endl;
    // For division, see RingElemTests (below)
    cout << "power(a, 8)  gives "  << power(a, 8) << endl; // !!! CANNOT use a^8
    cout << endl;
  }


  // Numerical comparisons: greater than, less than etc.
  void RingElemComparisons(const RingElem& a, const RingElem& b)
  {
    const ring& R = owner(a);
    if (owner(b) != R)
      CoCoA_THROW_ERROR("RingElems should belong to same ring", "TestRing");

    if (!IsOrderedDomain(R))
    { cout << "Ring is not ordered, so comparisons not possible." << endl; return; }

    cout << " --RingElemComparisons(a,b)-- ";
    cout << "a = " << a << ";  b = " << b << endl;

    // Equality and not-equality are == and !=  (see RingElemArith above)
    cout << "a < b gives " << (a < b) << endl;
    cout << "a <= b gives " << (a <= b) << endl;
    cout << "a > b gives " << (a > b) << endl;
    cout << "a >= b gives " << (a >= b) << endl;
    cout << endl;
  }


  // Tests on RingElems, incl IsDivisible
  void RingElemTests(const RingElem& a, const RingElem& b)
  {
    const ring& R = owner(a);
    if (owner(b) != R) CoCoA_THROW_ERROR("RingElems should belong to same ring", "TestRing");
    const RingElem one(R, 1);

    cout << " --RingElemTests(a,b)-- ";
    cout << "a = " << a << ";  b = " << b << endl;

    // Equality:
    cout << "a == b gives " << (a == b) << endl;
    cout << "a != b gives " << (a != b) << endl;

    cout << "IsZero(one)     gives " << IsZero(one) << endl;
    cout << "IsOne(one)      gives " << IsOne(one) << endl;
    cout << "IsMinusOne(one) gives " << IsMinusOne(one) << endl;
    cout << endl;

    cout << "IsDivisible(a, b) gives " << IsDivisible_AllowFields(a, b) << endl;
    if (!IsDivisible_AllowFields(a, b)) 
      cout << "  so we CANNOT compute b / a" << endl;
    else 
      cout << "  b / a gives " << b/a << endl;
    cout << endl;
  }


  // Some "special" operations work only if the ring elems are
  // in the right sort of ring.  Here we test the ring type, and
  // then exhibit the special functions.
  void RingElemSpecial(RingElem a, const RingElem& b)
  {
    const ring& R = owner(a);
    if (owner(b) != R) CoCoA_THROW_ERROR("RingElems should belong to same ring", "TestRing");
    const RingElem one(R, 1);

    cout << " --RingElemSpecial(a,b)-- ";
    cout << "a = " << a << ";  b = " << b << endl;

    // FRACTION FIELDS
    cout << "IsFractionField(R) gives " << IsFractionField(R) << endl;
    if (!IsFractionField(R)) 
      cout << "  so we CANNOT compute num(a), den(a)" << endl;
    else 
    {
      cout << "  num(a)   gives " << num(a) << endl;
      cout << "  den(a)   gives " << den(a) << endl;
    }
    cout << endl;

    // (TRUE) GCD DOMAINS
    cout << "IsTrueGCDDomain(R) gives " << IsTrueGCDDomain(R) << endl;
    if (!IsTrueGCDDomain(R)) 
      cout << "  so we CANNOT compute gcd(a,b)" << endl;
    else 
      cout << "  gcd(a, b)   gives " << gcd(a, b) << endl;
    cout << endl;

    // POLYNOMIALS
    cout << "IsPolyRing(R) gives " << IsPolyRing(R) << endl;
    if (!IsPolyRing(R)) 
      cout << "  so we CANNOT compute deg(a) or StdDeg(a)" << endl;
    else 
    {
      cout << "  NB deg(a) and StdDeg(a) are synonyms;" << endl
           << "  StdDeg is a more precise name but is also more cumbersome" << endl;
      cout << "  deg(a)    gives " << deg(a) << endl;
      cout << "  StdDeg(a) gives " << StdDeg(a) << endl;
    }
    cout << endl;

    // (SPARSE) POLYNOMIALS
    cout << "IsSparsePolyRing(R) gives " << IsSparsePolyRing(R) << endl;
    if (!IsSparsePolyRing(R)) 
      cout << "  so we CANNOT compute LPP(a), wdeg(a)" << endl;
    else 
    {
      cout << "  LPP(a)  gives " << LPP(a) << endl;
      cout << "  wdeg(a) gives " << wdeg(a) << endl;
      cout << "  NB wdeg(a) agrees with deg(a) only if the ring is standard graded." << endl;
      if (GradingDim(R)>0)
      {
        cout << "  LF(a) gives " << LF(a) << endl;
        cout << "  CutLF(a) gives " << CutLF(a) << "  --> same as LF(a), EXCEPT that" << endl
             << "  a has now become a = " << a << " --> its leading form was `cut off'" << endl;
      }
    }
    cout << endl;
  }



  //-- main --------------------------------------------------------------
  // we run the above fns on elements of some rings

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools are printed as "true" and "false"

    //------------------------------------------------------------
    // ZZ
    ring ZZ = RingZZ();
    {
      RingElem a(ZZ), b(ZZ);
    
      a = 5;
      b = 7;
      cout << "--------------------" << endl;
      cout << "|||  ring is ZZ  |||" << endl;
      cout << "--------------------" << endl;
      RingElemArith(a, b);
      RingElemComparisons(a, b);
      RingElemTests(a, b);
      RingElemSpecial(a, b);
    }
  
    //------------------------------------------------------------
    // QQ
    ring QQ = RingQQ();  
    {
      RingElem a(QQ), b(QQ);
    
      a = 5; a /= 2;
      b = 7; b /= 3;
      cout << "--------------------" << endl;
      cout << "|||  ring is QQ  |||" << endl;
      cout << "--------------------" << endl;
      RingElemArith(a, b);
      RingElemComparisons(a, b);
      RingElemTests(a, b);
      RingElemSpecial(a, b);
    }

    //------------------------------------------------------------
    // QQ[x,y]

    PolyRing P = NewPolyRing(QQ, symbols("x,y"));
    RingElem f = RingElem(P, "15*y + y^3");
    RingElem g = RingElem(P, "4*y - 3*y^3 + x^2");
    cout << "-------------------------" << endl;
    cout << "|||  ring is QQ[x,y]  |||" << endl;
    cout << "-------------------------" << endl;
    RingElemArith(f, g);
    RingElemComparisons(f, g);
    RingElemTests(f, g);
    RingElemSpecial(f, g);
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingElem2.C,v 1.5 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingElem2.C,v $
// Revision 1.5  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2021/11/05 14:31:32  abbott
// Summary: Improved comments
//
// Revision 1.3  2021/11/05 14:30:28  bigatti
// Summary: mini-titles: more readable
//
// Revision 1.2  2021/11/03 16:45:30  abbott
// Summary: Restructured
//
// Revision 1.1  2021/09/13 13:54:10  abbott
// Summary: Renamed from ex-RingElem1
//
// Revision 1.14  2020/10/14 20:00:27  abbott
// Summary: Use IsDivisible_AllowFields instead of IsDivisible
//
// Revision 1.13  2017/07/03 19:51:54  abbott
// Summary: Added CutLF example
//
// Revision 1.12  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.11  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.10  2015/07/01 16:31:35  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.9  2015/06/29 15:47:58  bigatti
// -- code in namespace CoCoA
//
// Revision 1.8  2014/07/08 12:47:26  abbott
// Summary: Removed AsPolyRing, AsSparsePolyRing, AsQuotientRing
// Author: JAA
//
// Revision 1.7  2014/03/21 18:16:00  bigatti
// -- now using ReadExpr
//
// Revision 1.6  2012/10/05 10:21:39  bigatti
// -- added LF (leading form)
//
// Revision 1.5  2012/05/22 10:02:38  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.4  2012/05/20 09:43:54  abbott
// Corrected layout.
//
// Revision 1.3  2012/02/08 17:43:14  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.5  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
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
// Revision 1.3  2006/03/01 14:26:41  cocoa
// -- removed some useless "include"
//
// Revision 1.2  2006/02/13 12:08:04  cocoa
// Fixed a problem with some missing assignment ops for certain PPMonoidElems.
// Fixed a bug in RingDistrMPoly::myIndetPower.
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.3  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.2  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/12/09 15:08:42  cocoa
// -- added log info
//
