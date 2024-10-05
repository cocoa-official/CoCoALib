// Copyright (c) 2005,2008  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows how to create various types of ring, and \n"
  "several operations one can perform on rings (e.g. querying their    \n"
  "properties).                                                        \n";

const string LongDescription =
  "This example creates several different sorts of ring,            \n"
  "and then calls PrintRingInfo on each one.                        \n"
  "PrintRingInfo calls various functions to obtain information      \n"
  "about each ring passed to it.  Naturally, some query functions   \n"
  "make sense only for certain types of ring (e.g. NumIndets(R)).   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Function to print out some information about a ring.
  // Here you see how various query functions of a ring can be used.
  // If the ring is actually PolyRing, we create a "view" of it as
  // a PolyRing, and this permits us to make further queries which
  // are valid only for PolyRings.

  void PrintRingInfo(const ring& R)
  {
    cout << "--------========--------" << endl;
    cout << "R                   = " << R << endl; 
    cout << "its zero and one elements are " << zero(R) << " and " << one(R) << endl;

    // NOTE result is a BigInt not an int:
    cout << "char(R)             = " << characteristic(R) << endl; 

    // NB: "Is.." functions ask about implementation (not isomorphism classes)
    // "IsBla" should be read as "Is internally implemented as Bla"
    cout << "IsTrueGCDDomain(R)  = " << IsTrueGCDDomain(R) << endl;
    cout << "IsIntegralDomain(R) = " << IsIntegralDomain(R) << endl;
    cout << "IsField(R)          = " << IsField(R) << endl;

    // The next four queries tell how the ring is implemented:
    cout << "IsZZ(R)             = " << IsZZ(R) << endl;
    cout << "IsFractionField(R)  = " << IsFractionField(R) << endl;
    cout << "IsPolyRing(R)       = " << IsPolyRing(R) << endl;
    cout << "IsQuotientRing(R)   = " << IsQuotientRing(R) << endl;

    // If R is a PolyRing, we print out some more details...
    if (IsPolyRing(R))
    {
      cout << "........................." << endl;
      cout << "Additional information specific to a polynomial ring:" << endl;

      // P and R are two views of the same underlying ring; we can "see"
      // more detail by querying P since it "knows" that the ring is actually
      // a polynomial ring.
      cout << "NumIndets(R)        = " << NumIndets(R) << endl; // NB CANNOT do NumIndets(R)
      cout << "The indeterminates are: ";
      for (long i=0; i < NumIndets(R); ++i)
        cout << indet(R, i) << "  ";
      cout << endl;
      cout << "The information for the coefficient ring is:" << endl;
      PrintRingInfo(CoeffRing(R));
      cout << "-- end coefficient ring --" << endl;
    }

    cout << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    const int SmallPrime = 13;
    const BigInt LargePrime = NextProbPrime(power(10,20) + 1);
    const BigInt composite = power(7,510);  // a large composite

    ring ZZ = RingZZ();                  // the integers
    ring QQ = RingQQ();                  // the rationals
    ring Fsmallp = NewZZmod(SmallPrime); // small finite field
    ring Fbigp = NewZZmod(LargePrime);   // large finite field

    // We may also create quotient rings:
    ring ZZmodN = NewQuotientRing(ZZ, ideal(RingElem(ZZ, composite)));
    ring ZZmod0 = NewQuotientRing(ZZ, ideal(zero(ZZ))); // NB isomorphic to ZZ but not implemented as ZZ!

    // This is a ring of floating point values.
    ring R = NewRingTwinFloat(128);  // guarantees 128 bit of precision (otherwise throws).
  
    // Here we create two polynomial rings: see ex-PolyRing1 for more details.
    ring P = NewPolyRing(QQ, symbols("x,y,z"), lex);  // QQ[x,y,z]
    ring P30 = NewPolyRing(NewZZmod(30), SymbolRange("x",1,3));  // ZZ/(30)[x[1..3]]

    cout << "-- Info for ring ZZ" << endl; 
    PrintRingInfo(ZZ);

    cout << "-- Info for ring QQ" << endl; 
    PrintRingInfo(QQ);

    cout << "-- Info for ring F_13" << endl; 
    PrintRingInfo(Fsmallp);

    cout << "-- Info for ring F_p (for some large prime p)" << endl; 
    PrintRingInfo(Fbigp);

    cout << "-- Info for ring ZZ/(n) (for some large composite n)" <<endl;
    PrintRingInfo(ZZmodN);

    cout << "-- Info for ring ZZ/(0) (NB: different from ZZ)" << endl; 
    PrintRingInfo(ZZmod0);

    cout << "-- Info for ring R (twin floats)" << endl; 
    PrintRingInfo(R);

    cout << "-- Info for ring QQ[x,y,z]" << endl; 
    PrintRingInfo(P);

    cout << "-- Info for ring (ZZ/(30))[x[0],x[1],x[2]]" << endl; 
    PrintRingInfo(P30);
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-ring2.C,v 1.20 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-ring2.C,v $
// Revision 1.20  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.19  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.18  2015/07/27 11:50:49  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.17  2015/06/29 15:52:45  bigatti
// *** empty log message ***
//
// Revision 1.16  2015/06/29 13:17:29  bigatti
// -- code inside namespace CoCoA
//
// Revision 1.15  2014/07/07 11:59:56  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.14  2012/07/04 12:25:33  abbott
// Improved several printed messages.
//
// Revision 1.13  2012/05/22 10:02:38  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.12  2012/03/30 09:28:45  bigatti
// -- just aligned
//
// Revision 1.11  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.10  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.9  2011/05/24 14:58:25  abbott
// Consequential changes from removal of several ctors for principal ideals.
//
// Revision 1.8  2011/03/10 17:16:12  abbott
// Replaced size_t by long.
//
// Revision 1.7  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.6  2010/03/22 13:57:16  bigatti
// -- removed ring Z from ideal ctor (no longer necessary)
//
// Revision 1.5  2010/03/05 21:35:00  abbott
// Added some consts.  Changed an example to use PPOrderingCtor when creating a PolyRing.
//
// Revision 1.4  2009/06/05 12:14:55  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from BigInt.H/C
//   removed gcd from MachineInteger.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
// Revision 1.3  2008/07/21 07:58:51  abbott
// Several cosmetic changes to examples.
//
// Revision 1.2  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2007/03/08 17:28:47  bigatti
// -- improved NewPolyRing
//
// Revision 1.1  2007/03/07 14:30:38  bigatti
// -- was ex-ring1
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/26 17:40:34  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/22 17:26:23  bigatti
// -- added printing of separator after coeff ring
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
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.3  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
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
// Revision 1.4  2004/11/22 17:10:04  cocoa
// -- changed  "PolyRing" --> "AsPolyRing",  "FractionField" --> "AsFractionField"
//
