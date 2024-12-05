// Copyright (c) 2018  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a program showing how to use VerificationLevel.  \n";

const string LongDescription =
  "This is a program showing how to use VerificationLevel.  \n"
  "We construct a polynomial which confuses MinPolyQuot, then\n"
  "show that MinPolyQuot gives the wrong answer if verification\n"
  "is 0; but it get the right answer if verification is 2.\n"
  "The guaranteed correct answer can be obtained by passing \"guaranteed()\"";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    const ring P = NewPolyRing(RingQQ(), symbols("x"));
    const RingElem& x = indet(P,0);
    const long p0 = NextPrime(1024*1024*1024);
    const long p1 = NextPrime(p0);
    const long p2 = NextPrime(p1);

    const ideal I = ideal(x*x- p1*p2); // deliberately constructed poly
    const RingElem m_verif0 = MinPolyQuot(x,I,x,  VerificationLevel(0));
    const RingElem m_verif2 = MinPolyQuot(x,I,x,  VerificationLevel(2));

    cout << "Without verification, we get the wrong result: " << m_verif0 << endl;
    cout << "With verification, we get the right result: " << m_verif2 << endl;

    // Now compute the guaranteed correct answer.
    // Surely correct, but generally takes rather longer.
    const RingElem m = MinPolyQuot(x,I,x);
    const RingElem m2 = MinPolyQuot(x,I,x, guaranteed()); // same as line above
    CoCoA_ASSERT_ALWAYS(m == m2);
    CoCoA_ASSERT_ALWAYS(m == m_verif2);
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
  // catch (const CoCoA::InterruptReceived& intr)
  // {
  //   cerr << endl
  //        << "------------------------------" << endl
  //        << ">>>  CoCoALib interrupted  <<<" << endl
  //        << "------------------------------" << endl
  //        << "-->>  " << intr << "  <<--" << endl;
  //   return 2;
  // }
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
