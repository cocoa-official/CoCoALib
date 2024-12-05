// Copyright (c) 2005  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Operations between elements of different rings are not allowed \n"
  "but we can use homomorphisms to map the elements into the      \n"
  "same ring. \n";

const string LongDescription =
  "The example in this file shows how to create and use some      \n"
  "homomorphisms between rings.  In particular, it gives a simple \n"
  "example of mixed ring arithmetic: the user must map all values \n"
  "into a single ring before combining them arithmetically.       \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create some coefficient rings
    ring ZZ = RingZZ();
    ring QQ = RingQQ();
    ring Fp = NewZZmod(32003);
    ring P = NewPolyRing(ZZ, symbols("alpha,beta,gamma"));   // P is ZZ[alpha,beta,gamma]

    RingElem c(ZZ), f(P), q(QQ), a(Fp);
    c = 3;
    f = 100;
    q = 5;
    a = -1;
  
    cout << "c = " << c << "   \t in " << owner(c) << endl;
    cout << "f = " << f << "   \t in " << owner(f) <<  endl;
    cout << "q = " << q << "   \t in " << owner(q) <<  endl;
    cout << "a = " << a << "   \t in " << owner(a) <<  endl;
    cout << endl;

    // THE WRONG WAY TO DO IT...
    try
    {
      cout << "c * f = " << c * f;
      // Mixed ring arithmetic triggers an exception ==> WE NEVER REACH HERE.
    }
    catch (const CoCoA::ErrorInfo& err) 
    {
      if (err != ERR::MixedRings) throw; // rethrow any unexpected exceptions
      cout << " !!!  As expected we get this error:" << endl
           << "      " << message(err) << endl
           << "The correct way is to use ring homomorphisms as follows:" << endl;
    }

    // THE RIGHT WAY TO DO IT...
    RingHom emb = CanonicalHom(ZZ, P);
    // or, equivalently,  CoeffEmbeddingHom(AsPolyRing(P));
    cout << "emb(c) * f = " << emb(c) * f << endl << endl;

    RingHom frac = CanonicalHom(ZZ, QQ);
    // or, equivalently,  EmbeddingHom(AsFractionField(QQ));
    cout << "frac(c)/4 + q = " << frac(c)/4 + q << endl << endl;

    RingHom mod = CanonicalHom(ZZ, Fp);
    // or, equivalently,  QuotientingHom(AsQuotientRing(Fp));
    cout << "mod(c) - a = " << mod(c) - a << endl << endl;
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
