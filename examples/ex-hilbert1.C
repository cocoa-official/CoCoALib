// Copyright (c) 2006   John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example is just for testing the Hilbert code.                   \n"
  "It might disappear as soon as HilbertSeries is included in CoCoALib. \n";

const string LongDescription =
  "This code also shows how to create the \"chess examples\".           \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Includes from the standard C++ library
#include <algorithm>  // using std::min;
  // #include <iostream>   // using std::endl;


  // ---- chess tests -----------------------------------------------------

  ConstRefRingElem CsbSquareIndet(SparsePolyRing P, long l, long sq1, long sq2)
  {
    CoCoA_ASSERT( l*l <= NumIndets(P) );
    CoCoA_ASSERT( sq1 <= l && sq2 <= l );
    return indet(P, (sq1-1)*l + (sq2-1));
  }


  ideal NewQueenMovesFrom(SparsePolyRing P, long Csb, long sq1, long sq2)
  {
    ConstRefRingElem x = CsbSquareIndet(P, Csb, sq1, sq2);
    vector<RingElem> g;
    for ( long i=sq2+1 ; i<=Csb ; ++i )
      g.push_back(x * CsbSquareIndet(P, Csb, sq1, i));
    for ( long i=sq1+1 ; i<=Csb ; ++i )
      g.push_back(x * CsbSquareIndet(P, Csb, i, sq2));
    for ( long i=min(Csb-sq1,Csb-sq2) ; i>0 ; --i )
      g.push_back(x * CsbSquareIndet(P, Csb, sq1+i, sq2+i));
    for ( long i=min(Csb-sq1, sq2-1) ; i>0 ; --i )
      g.push_back(x * CsbSquareIndet(P, Csb, sq1+i, sq2-i));
    return ideal(P, g);
  }


  ideal NewQueenIdeal(SparsePolyRing P, long Csb)
  {
    ideal I = ideal(zero(P));
    for ( long sq1=1 ; sq1<=Csb ; ++sq1 )
      for ( long sq2=1 ; sq2<=Csb ; ++sq2 )
        I += NewQueenMovesFrom(P, Csb, sq1, sq2);
    return I;
  }


  // ----------------------------------------------------------------------

  void program()
  {
    GlobalManager CoCoAFoundations;
    SignalWatcher MonitorInterrupt(SIGINT);

    cout << ShortDescription << endl;

    ring QQ = RingQQ();

    SparsePolyRing P = NewPolyRing(QQ, symbols("x,y,z,t"));
    const vector<RingElem>& x = indets(P);
    ideal I = ideal(x[1], x[2], x[3]);
    cout << "gens(I) = " << gens(I) << endl;
    cout << "TidyGens(I) = " << TidyGens(I) << endl;
    cout << "HilbertNumQuot(I) = "
         << HilbertNumQuot(I) << endl;

    SparsePolyRing CsbRing7 = NewPolyRing(QQ, SymbolRange("x", 1,49));
    ideal Queen7 = NewQueenIdeal(CsbRing7, 7);

    double T;
    T=CpuTime();
    TidyGens(Queen7);
    cout << "TidyGens time = " << CpuTime()-T << endl;
    T=CpuTime();  
    cout << "HilbertNumQuot(Queen7) = "
         << HilbertNumQuot(Queen7) << endl;
    cout << "Hilbert time = " << CpuTime()-T << endl;

    SparsePolyRing CsbRing = NewPolyRing(QQ, SymbolRange("x", 1,64));
    ideal Queen8 = NewQueenIdeal(CsbRing, 8);
    T=CpuTime();
    TidyGens(Queen8);
    cout << "TidyGens time = " << CpuTime()-T << endl;
    T=CpuTime();  
    RingElem HSNum = HilbertNumQuot(Queen8);
    cout << "Hilbert time = " << CpuTime()-T << endl;
    cout << "HilbertNumQuot(Queen8) = "
         << HSNum << endl;
    ring HRing = owner(HSNum);
    cout << "leading monomial of HSNum is "
         << monomial(HRing, LC(HSNum), LPP(HSNum)) << endl;

    FractionField HPSRing = NewFractionField(owner(HSNum));
    RingElem t = indet(owner(HSNum), 0);
    RingHom phi = EmbeddingHom(HPSRing);
    cout << "HilbertSeries (in FractionField) = "
         << phi(HSNum)/phi(power(1-t,64)) << endl;
    cout << "---------------------------------" << endl;
    cout << "HilbertSeries(NewQuotientRing(RingOf(I),I)) = "
         << HilbertSeries(NewQuotientRing(RingOf(I),I)) << endl;
    HPSeries HS = HilbertSeriesQuot(I);
    cout << "HS = HilbertSeriesQuot(I) = " << endl;
    cout << "HS = " << HS << endl;
    //    cout << "HSSimplified(HS) = "
    //         << HSSimplified(HS) << endl;
    cout << "num(HS) = " << num(HS) << endl;
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
  catch (const CoCoA::InterruptReceived& intr)
  {
    cerr << endl
         << "------------------------------" << endl
         << ">>>  CoCoALib interrupted  <<<" << endl
         << "------------------------------" << endl
         << "-->>  " << intr << "  <<--" << endl;
    return 2;
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
