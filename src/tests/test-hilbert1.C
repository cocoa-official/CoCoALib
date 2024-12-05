//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-hilbert.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/time.H"


#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// First test for Hilbert.
//----------------------------------------------------------------------
namespace CoCoA
{

// ---- chess tests -----------------------------------------------------

  const RingElem& CsbSquareIndet(const SparsePolyRing& P, long l, long sq1, long sq2)
  {
    CoCoA_ASSERT_ALWAYS( l*l <= NumIndets(P) );
    CoCoA_ASSERT_ALWAYS( sq1 <= l && sq2 <= l );
    return indet(P, (sq1-1)*l + (sq2-1));
  }


  ideal NewQueenMovesFrom(const SparsePolyRing& P, long Csb, long sq1, long sq2)
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
    return ideal(P, g); // ideal(P, g) because g might be empty
  }


  ideal NewQueenIdeal(const SparsePolyRing& P, long Csb)
  {
    ideal I = ideal(zero(P));
    for ( long sq1=1 ; sq1<=Csb ; ++sq1 )
      for ( long sq2=1 ; sq2<=Csb ; ++sq2 )
        I += NewQueenMovesFrom(P, Csb, sq1, sq2);
    return I;
  }


  RingElem HilbertNumQuot_C_time(const ideal& I)
  {
    //    double t0 = CpuTime();
    RingElem HNum = HilbertNumQuot_C(I);
    //    std::cout << "HilbertNumQuot_C time: " << CpuTime() -t0 << endl;    
    return HNum;
  }
  
  
  RingElem HilbertNumQuot_time(const ideal& I)
  {
    //    double t0 = CpuTime();
    RingElem HNum = HilbertNumQuot(I);
    //    std::cout << "HilbertNumQuot   time: " << CpuTime() -t0 << endl;    
    return HNum;
  }
  
  
  RingElem MGHilbertNumQuot_time(const ideal& I)
  {
    //    double t0 = CpuTime();
    RingElem HNum = MGHilbertNumQuot(I);
    //    std::cout << "MGHilbertNumQuot time: " << CpuTime() -t0 << endl;    
    return HNum;
  }
  
  
// ----------------------------------------------------------------------

  void program()
  {
    GlobalManager CoCoAFoundations;

    const PolyRing QQt = RingQQt(1);     // owner(HNum)  is  QQt
    const SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",1,4));
    const vector<RingElem>& x = indets(P);

    CoCoA_ASSERT_ALWAYS(HilbertNumQuot_time(ideal(zero(P))) == 1);
  
    ideal I = ideal(x[1], x[2], x[3]);
    RingElem HNum = HilbertNumQuot_time(I);
    CoCoA_ASSERT_ALWAYS(HNum == RingElem(QQt, "-t^3 +3*t^2 -3*t +1"));
    CoCoA_ASSERT_ALWAYS(HilbertNumQuot_C_time(I) == HNum);
    CoCoA_ASSERT_ALWAYS(MGHilbertNumQuot_time(I) == HNum);
  
    //std::cout << "  --====Chessboard-examples====--" << std::endl;
    SparsePolyRing CsbRing = NewPolyRing(RingQQ(), SymbolRange("x",1,9*9));

    ideal Q3 = NewQueenIdeal(CsbRing, 3);
    RingElem HN3 = HilbertNumQuot_C_time(Q3);
    RingElem HN3_CPP = HilbertNumQuot_time(Q3);
    RingElem MGHN3_CPP = MGHilbertNumQuot_time(Q3);
    CoCoA_ASSERT_ALWAYS(HN3 == HN3_CPP);
    CoCoA_ASSERT_ALWAYS(HN3 == MGHN3_CPP);
  
    ideal Q4 = NewQueenIdeal(CsbRing, 4);
    RingElem HN4 = HilbertNumQuot_C_time(Q4);
    RingElem HN4_CPP = HilbertNumQuot_time(Q4);
    RingElem MGHN4_CPP = MGHilbertNumQuot_time(Q4);
    CoCoA_ASSERT_ALWAYS(HN4 == HN4_CPP);
    CoCoA_ASSERT_ALWAYS(HN4 == MGHN4_CPP);

    ideal Q6 = NewQueenIdeal(CsbRing, 6);
    RingElem HN = HilbertNumQuot_C_time(Q6);
    CoCoA_ASSERT_ALWAYS(LC(HN) == 19);
    CoCoA_ASSERT_ALWAYS(deg(HN) == 36);
    RingElem HN_CPP = HilbertNumQuot_time(Q6);
    RingElem MGHN_CPP = MGHilbertNumQuot_time(Q6);
    CoCoA_ASSERT_ALWAYS(HN == HN_CPP);
    CoCoA_ASSERT_ALWAYS(HN_CPP == MGHN_CPP);
  
    /*
      ideal Q8 = NewQueenIdeal(CsbRing, 8);
      ideal Q9 = NewQueenIdeal(CsbRing, 9);
      RingElem HN9 = HilbertNumQuot_C_time(Q9);
      RingElem HN9_CPP = HilbertNumQuot_time(Q9);
      CoCoA_ASSERT_ALWAYS(HN9 == HN9_CPP);
    */

    // this line is just to avoid mempool complaints with -DCoCoA_MEMPOOL_DEBUG
    EndPoincare_C();
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
