// Copyright (c) 2015  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example illustrates basic use of SmallFpImpl for arithmetic in a\n"
  "small prime finite field.  ex-SmallFp3 shows how to compute more     \n"
  "efficiently, but it is not for the faint-hearted!                    \n";

const string LongDescription =
  "This example illustrates basic use of SmallFpImpl for arithmetic in a\n"
  "small prime finite field.  ex-SmallFp3 shows how to compute more     \n"
  "efficiently, but it is not for the faint-hearted!                    \n"
  "Here we see how to create a SmallFpImpl object, and its use for basic\n"
  "arithmetic (add, sub, mul, div, power) on finite field elements.     \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    const int p = NextPrime(4096);
    SmallFpImpl ModP(p);
    SmallFpImpl::value a,b,c; // initially 0
    a = ModP.myReduce(99);    // a = 99; mod p
    b = ModP.myReduce(28);    // b = 28; mod p
    c = ModP.myReduce(-3);    // c = -3; mod p
    a = ModP.myAdd(b,c);      // a = b+c;  also mySub, myMul, myDiv
    a = ModP.myNegate(a);     // a = -a;
    a = ModP.myRecip(a);      // a = 1/a;
    // Now verify Fermat's Little theorem
    b = ModP.myPower(a,p-1);  // b = a^(p-1);  where ^ means "power"
    if (!IsZero(a) && !IsOne(b)) cerr << "Fermat's Little Theorem failed!" << endl;

    a = ModP.myAddMul(a,b,c); // a = a+b*c;

    const long B = ModP.myExport(b); // deliver value of b as a long (see ex-SmallFp2)
    if (B != 1) cerr << "Exported value " << B << " should have been 1." << endl;
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
