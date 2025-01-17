// Copyright (c) 2005,2007  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Some simple computations with rational numbers.                       \n"
  "This example illustrates how to create the field of rational numbers. \n"
  "It shows that we can compute 7/3 in QQ but not in ZZ.                 \n"
  "It shows how to map an integer into a rational number.                \n";

const string LongDescription =
  "Familiarize yourself with the example ex-RingZZ1.C before proceeding. \n"
  "As C++ does not natively have any rings, we must construct them from  \n"
  "scratch.                                                              \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Creation of the rings we shall be working in.
    ring ZZ = RingZZ();                // ZZ is the ring of integers.
    ring QQ = RingQQ();                // Create field of rationals directly.
    ring K  = NewFractionField(ZZ);    // Result is identical to QQ.
    cout << "RingQQ() == NewFractionField(ZZ)  is " << (QQ == K) << endl;
  
    // Now we have created the rings we can start using their elements:
    RingElem n(ZZ);       // n is a RingElem belonging to Z; initial value is 0.
    RingElem p(QQ);       // p is a RingElem belonging to Q; initial value is 0.

    //----------------------------------------------------------------------
    // Calculate and print 7/3 as an element of Q...
    p = 7;   // p is the rational number 7
    p = p/3; // now is the rational number 7/3
    cout << "7 divided by 3 is " << p 
         << " (an element of " << QQ << ")" << endl;
    // Why did I not simply write "p = 7/3;"?
    // What value would this assign to p?  And why?
    // I cannot use 7.0/3.0 either as C++ would compute only an approximation
    // to the true value.

    //----------------------------------------------------------------------
    // We cannot compute 7/3 as an element of ZZ as the division is not exact.
    // n = 7;
    // n /= 3; // <-- this would throw a CoCoA::ErrorInfo with code ERR::BadQuot
    cout << "But in " << ZZ << " we cannot compute 7 divided by 3." << endl;

    //----------------------------------------------------------------------
    // We can map n into QQ using a homomorphism (here a canonical ringhom).
    // See ex-RingHom* for more general examples of ringhoms.
    RingHom phi = CanonicalHom(ZZ, QQ);
    RingHom psi = EmbeddingHom(QQ); // same homomorphism
    cout << "n = " << n << " in " << owner(n) << endl;
    cout << "phi(n) = " << phi(n) << " in " << owner(phi(n)) << endl;
    cout << "psi(n) = " << psi(n) << " in " << owner(psi(n)) << endl;
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
