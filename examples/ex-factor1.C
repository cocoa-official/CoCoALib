// Copyright (c)  2009  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to interpret the result of a factorization. \n";

const string LongDescription =
  "This example shows how to interpret the result of a factorization. \n"
  "It creates a ring element (belonging to a polynomial ring), and    \n"
  "factorizes it.  The result is a \"factorization object\".  We show \n"
  "how to access/use the various fields in this object.               \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    SparsePolyRing P = NewPolyRing(RingQQ(), symbols("x,y")); // QQ[x,y];
    RingElem f = RingElem(P, "x^96 - y^96");

    const factorization<RingElem> FacInfo = factor(f);

    // These are convenient aliases for the 3 fields in the factorization:
    const RingElem&         content   = FacInfo.myRemainingFactor();
    const vector<RingElem>& IrredFacs = FacInfo.myFactors();
    const vector<long>&     mult      = FacInfo.myMultiplicities();
    
    // Print out the factorization in a "nice" way:
    cout << "The factors of " << f << " are:" << endl;
    if (!IsOne(content))
      cout << "content: " << content << endl;
    const int NumIrredFacs = len(IrredFacs);
    for (int i = 0; i != NumIrredFacs; ++i)
    {
      cout << IrredFacs[i] << "  with multiplicity  " << mult[i] << endl;
    }
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
