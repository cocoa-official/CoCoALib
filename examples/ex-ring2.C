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
