// Copyright (c) 2005 John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example program illustrating an implementation of Normal Remainder \n"
  "wrt a list of polynomials.                                         \n"
  "If the list is a Groebner Basis, NR returns the Normal Form.       \n";

const string LongDescription =
  "This is just an example!  If you want to compute Normal Forms \n"
  "you should use the library function \"NF\".                   \n";
//----------------------------------------------------------------------


namespace CoCoA
{

  long FindReducer(ConstRefRingElem f, const vector<RingElem>& v)
  {
    if (IsZero(f)) return -1;

    for (long j=0; j<len(v); ++j)  // 0 to len-1
      if (IsDivisible(LPP(f), LPP(v[j])))
        return j;
    return -1;
  }
  

  RingElem NRLM(ConstRefRingElem f, const vector<RingElem>& v)
  {
    const SparsePolyRing P = owner(f);
    RingElem m(P);
    RingElem r =f;
  
    long j = FindReducer(r, v);
    while (j != -1)
    {
      //   m = LM(r)/LM(v[i]);  
      P->myDivLM(raw(m), raw(r), raw(v[j])); // no checks
      r -= m * v[j];
      if (IsZero(r)) return zero(P);
      j = FindReducer(r, v);
    }
    return r;
  }


  RingElem NormalRemainder(ConstRefRingElem f, const vector<RingElem>& v)
  {
    if (IsZero(f)) return f;
    const SparsePolyRing P = owner(f);
    RingElem ansNR(P);
    RingElem tmpNR =f;
  
    tmpNR = NRLM(f, v);
    while (!IsZero(tmpNR))
    {
      P->myMoveLMToBack(raw(ansNR), raw(tmpNR));
      tmpNR = NRLM(tmpNR, v);
    }
    return ansNR;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring P = NewPolyRing(RingQQ(), symbols("x,y,z")); // QQ[x,y,z]

    vector<RingElem> g = RingElems(P, "x^3 -z^3,  x^2 -y");
    ideal I = ideal(g);
    cout << "gens(I)   = " << gens(I) << endl;
    vector<RingElem> GB = GBasis(I);  // it is the Groebner Basis of I
    cout << "GBasis(I) = " << GB << endl;
    cout << "When I is an ideal of polynomials \"TidyGens(I)\" returns its Groebner Basis." << endl << endl;

    RingElem f = RingElem(P, "x^12 + y^6 + z^6");
    cout << "-- f = " << f << endl;
    cout << "NormalRemainder(f, gens(I)) = " <<  NormalRemainder(f, gens(I)) << endl;
    cout << "NormalRemainder(f, GB)      = " <<  NormalRemainder(f, GB) << endl;
    cout << "NF(f, I)                    = " <<  NF(f, I) << endl;
    cout << endl;

    const vector<RingElem>& x = indets(P);
    RingElem h = g[0] + x[0]*g[1];
    cout << "-- h = " << h << endl;
    cout << "NormalRemainder(h, gens(I)) = " <<  NormalRemainder(h, gens(I)) << endl;
    cout << "NormalRemainder(h, GB)      = " <<  NormalRemainder(h, GB) << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-NF.C,v 1.12 2022/02/13 09:56:57 abbott Exp $
// $Log: ex-NF.C,v $
// Revision 1.12  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.11  2021/10/08 14:22:35  bigatti
// Summary: modernized
//
// Revision 1.10  2017/04/27 15:24:52  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.9  2015/06/29 15:27:25  bigatti
// -- minor fix
//
// Revision 1.8  2015/06/29 12:20:43  bigatti
// -- moved in cocoa namespace
//
// Revision 1.7  2015/04/27 12:50:13  bigatti
// Summary: updated (using ReadExpr and myMoveLMToBack)
//
// Revision 1.6  2014/07/08 12:47:26  abbott
// Summary: Removed AsPolyRing, AsSparsePolyRing, AsQuotientRing
// Author: JAA
//
// Revision 1.5  2012/02/08 17:40:23  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2010/12/17 16:07:55  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/07/14 13:00:28  bigatti
// -- cleaner code
//
// Revision 1.2  2009/08/10 14:46:51  abbott
// Cleaned up a call to ideal -- no longer need to specify the owning ring.
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
// Revision 1.4  2007/02/26 15:49:08  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 15:31:57  bigatti
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
// Revision 1.4  2006/02/14 16:22:20  cocoa
// -- defined "operator<<" for vector<RingElem>&  in ring.H/C
//
// Revision 1.3  2006/01/23 12:07:09  cocoa
// Further cleaning of Makefiles in response to a reported problem
// about spurious entries in src/Makefile_dependencies.
// Small change to doc for SmallFpImpl.
//
// Revision 1.2  2006/01/17 11:15:40  cocoa
// -- just a comment to say TidyGens is the Groebner Basis
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
// Revision 1.1  2005/05/04 16:35:06  cocoa
// -- new examples
//
