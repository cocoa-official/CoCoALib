// Copyright (c) 2015  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates use of CRTMill to build a large integer \n"
  "from its residues modulo various different primes.               \n";

const string LongDescription =
  "This program illustrates use of CRTMill to build a large integer \n"
  "from its residues modulo various different primes.               \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // Daft example: we shall construct N from its modular images.
    // We pretend we know only the UPB, and have a way to compute
    // N modulo p for any prime p (in this case we just compute N%p).
    const BigInt N = power(10,100);
    const BigInt UPB = 2*N+1;

    CRTMill crt;
    int p = 101;
    while (true)
    {
      p = NextPrime(p);
      crt.myAddInfo(N%p, p); // tell crt the new residue-modulus pair
      if (CombinedModulus(crt) >= UPB) break;
    }

    // Since we already know the answer, we can check it is correct.
    if (CombinedResidue(crt) != N)
      CoCoA_THROW_ERROR("Wrong answer", "CoCoA::Program");
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-NumTheory2.C,v 1.8 2022/02/13 09:56:57 abbott Exp $
// $Log: ex-NumTheory2.C,v $
// Revision 1.8  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.7  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.6  2015/11/21 19:14:39  abbott
// Summary: Corrected "greater than" to "greater or equal"
//
// Revision 1.5  2015/06/29 15:27:25  bigatti
// -- minor fix
//
// Revision 1.4  2015/06/29 12:20:43  bigatti
// -- moved in cocoa namespace
//
// Revision 1.3  2015/06/25 16:08:05  abbott
// Summary: New example for CRTMill (moved old ex-NumTheory2 --> ex-NumTheory4)
// Author: JAA
//
//
