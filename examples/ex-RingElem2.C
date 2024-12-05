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
