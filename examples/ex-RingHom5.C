// Copyright (c) 2005  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "The example in this file shows how to create and use some \n"
  "homomorphisms between rings.                              \n";

const string LongDescription =
  "We compute these polynomials (with parameters) in some rings: \n"
  "f = (2*a/3-1)*x[0] + 1/a;  g = x[0]-a;                        \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // A simple function to test embedding ring homs in CoCoALib.
  // We know that Kx is of the form R(a)[x[0..N]]
  void SimpleTest(const PolyRing& Kx) 
  {
    cout << "Kx is " << Kx << endl;
    // Give a handy name to the indeterminate x[0]
    RingElem x0 = indet(Kx, 0);
    cout << "indet(Kx,0) is " << indet(Kx,0) << endl;

    // the following 6 lines show how to work with several rings and RingHoms

    ring K  = CoeffRing(Kx); // this is R(a)
    ring Ra = BaseRing(K);  // this is R[a]
    RingHom RaToK = CanonicalHom(Ra, K);     // K = FractionField(R[a])
    RingHom KToKx = CanonicalHom(K, Kx);
    RingElem a_R = indet(Ra, 0); // a as RingElem of R[a]
    RingElem a   = KToKx(RaToK(a_R));        // a as RingElem of Kx

    // for this particular example this would have been much simpler:
    //   RingElem a(Kx, symbol("a"));

    RingElem f = (2*a/3-1)*x0 + 1/a;         // This is a RingElem of Kx.
    RingElem g = x0-a;                       // This is another one.

    cout << "f = " << f << endl;
    cout << "f*g = " << f*g << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // this changes "cout" so that bools are printed true/false

    ring QQ = RingQQ();
    ring Fp = NewZZmod(32003);

    //---------------------------------------------------------------------------
    {
      cout << "  --- Coeffs in " << Fp << endl;
      PolyRing Fpa = NewPolyRing(Fp, symbols("a"));              // Fpa is Fp[a]  where Fp = ZZ/(32003)
      ring K = NewFractionField(Fpa);                            // K   is Fp(a)
      PolyRing Kxyzt = NewPolyRing(K, symbols("x,y,z,t")); // Kxyzt is Fp(a)[x,y,z,t]
      SimpleTest(Kxyzt);
      cout << endl;
    }

    //---------------------------------------------------------------------------
    {
      cout << "Coeffs in " << RingZZ() << endl;
      PolyRing ZZa = NewPolyRing(RingZZ(), symbols("a"));    // ZZ[a]
      PolyRing Qax4 = NewPolyRing(NewFractionField(ZZa), NewSymbols(4)); // QQ(a)[x[0..3]]
      SimpleTest(Qax4);
      cout << endl << endl;
    }

    //---------------------------------------------------------------------------
    {
      cout << "  --- Coeffs in " << QQ << endl;
      PolyRing Qa = NewPolyRing(QQ, symbols("a"));           // QQ[a]
      PolyRing Qax4 = NewPolyRing(NewFractionField(Qa), NewSymbols(4));  // QQ(a)[x[0..3]]
      SimpleTest(Qax4);
      cout << endl;
    }

    cout << "  --- Now we supply an unsuitable input" << endl;
    cout << "      (it is unsuitable because CoeffRing is not K(a))" << endl;
    try
    {
      SimpleTest(NewPolyRing(Fp, symbols("a")));
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      if (err != ERR::ReqPolyRing) throw;
      cout << "\nOK!  Our unsuitable input produced the expected exception." << endl;
    }
    cout << endl;
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
