// Copyright (c) 2010   John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows what can be computed in CoCoALib using \n"
  "the C++ library Frobby\n";

const string LongDescription =
  "  \"Frobby - Computations With Monomial Ideals\"  \n"
  "Web page: http://www.broune.com/frobby/";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;
    cout << boolalpha; // prints true/false for bool

#ifndef CoCoA_WITH_FROBBY
    cout << "Frobby library is not linked." << endl << endl;
#else
    PolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    RingElem x(indet(P,0));
    RingElem y(indet(P,1));
    RingElem z(indet(P,2));
  
    //I := Ideal(x^2, x*y, y^2, z^2);
    ideal I = ideal(x*x, x*y, y*y, z*z);
    cout << "I = " << I << endl << endl;
  
    cout << "FrbAlexanderDual(I, LPP(x*x*y*y*z*z)) = "
         <<  FrbAlexanderDual(I, LPP(x*x*y*y*z*z)) << endl << endl;

    cout << "FrbMaximalStandardMonomials(I) = "
         <<  FrbMaximalStandardMonomials(I) << endl << endl;

    vector<ideal> ID;
    FrbIrreducibleDecomposition(ID, I);
    cout << "FrbIrreducibleDecomposition(ID, I):" << endl
         <<  "  ID[0] = " << ID[0] << endl
         <<  "  ID[1] = " << ID[1] << endl << endl;

    vector<ideal> PD;
    FrbPrimaryDecomposition(PD, I);
    cout << "FrbPrimaryDecomposition(PD, I):" << endl
         <<  "  PD[0] = " << PD[0] << endl << endl;

    vector<ideal> AP;
    FrbAssociatedPrimes(AP, I);
    cout << "FrbAssociatedPrimes(AP, I):" << endl
         <<  "  AP[0] = " << AP[0] << endl << endl;

    cout << "FrbDimension(ideal(x, y)) = "
         <<  FrbDimension(ideal(x, y)) << endl << endl;

    cout << "FrbMultigradedHilbertPoincareNumerator(ideal(x, y)) = "
         <<  FrbMultigradedHilbertPoincareNumerator(ideal(x, y)) << endl << endl;

    cout << "FrbTotalDegreeHilbertPoincareNumerator(ideal(x,y), x + 2 * y) = "
         <<  FrbTotalDegreeHilbertPoincareNumerator(ideal(x,y), x + 2 * y) << endl << endl;
  
#endif // CoCoA_WITH_FROBBY
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
