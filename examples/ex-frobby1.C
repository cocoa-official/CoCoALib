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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-frobby1.C,v 1.10 2024/04/15 17:05:16 abbott Exp $
// $Log: ex-frobby1.C,v $
// Revision 1.10  2024/04/15 17:05:16  abbott
// Summary: Corrected comment
//
// Revision 1.9  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.7  2016/07/20 08:44:11  abbott
// Summary: Changed ifdef into ifndef, and changed order of the two blocks
//
// Revision 1.6  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.5  2013/06/28 12:06:36  abbott
// Updated Frobby fn names.
//
// Revision 1.4  2012/10/12 13:23:24  bigatti
// -- added description
//
// Revision 1.3  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.2  2011/07/29 14:58:59  bigatti
// -- added (temporarily?) "Frobby" suffix to Frobby functions
//
// Revision 1.1  2011/07/20 09:13:51  bigatti
// -- first import
//
