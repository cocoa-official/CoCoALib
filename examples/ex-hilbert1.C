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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-hilbert1.C,v 1.16 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-hilbert1.C,v $
// Revision 1.16  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.15  2019/11/14 17:46:48  abbott
// Summary: Added interrupt handling
//
// Revision 1.14  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.13  2017/04/10 15:59:00  bigatti
// -- small updates
//
// Revision 1.12  2015/06/29 15:51:43  bigatti
// -- code in namespace CoCoA
//
// Revision 1.11  2014/07/08 13:13:51  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.10  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.9  2011/05/24 14:58:25  abbott
// Consequential changes from removal of several ctors for principal ideals.
//
// Revision 1.8  2011/04/27 09:40:40  bigatti
// -- syntax update
//
// Revision 1.7  2011/04/08 14:07:52  bigatti
// -- renamed HilbertNumeratorMod into HilbertNumQuot
//
// Revision 1.6  2011/03/11 12:01:13  bigatti
// -- changed size_t --> long
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2010/10/29 09:32:11  bigatti
// -- updated syntax for ideal ctor
//
// Revision 1.3  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/10/18 11:34:00  bigatti
// -- fixed: using DenseUPolyRing
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/08 17:27:20  bigatti
// -- improved NewPolyRing
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/02/28 15:14:14  bigatti
// -- minor: Q7 --> Queen7
//
// Revision 1.6  2007/02/26 17:40:34  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.5  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.4  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2006/11/17 18:16:02  cocoa
// -- added timings
//
// Revision 1.2  2006/11/16 17:39:51  cocoa
// -- fixed some calls to rum: lots of code could be deleted, but this is
//    all much more fragile than I would like
//
// Revision 1.1  2006/10/09 16:48:58  cocoa
// -- first import
//
