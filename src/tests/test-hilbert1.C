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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-hilbert1.C,v 1.33 2024/03/15 19:07:53 abbott Exp $
// $Log: test-hilbert1.C,v $
// Revision 1.33  2024/03/15 19:07:53  abbott
// Summary: Minor cleaning
//
// Revision 1.32  2024/03/01 17:33:55  bigatti
// Summary: cleaned-up timings code using dedicated functions (Blah_time(I))
//
// Revision 1.31  2024/02/23 11:08:29  abbott
// Summary: Commented out unusued variable t0
//
// Revision 1.30  2024/02/05 12:34:54  bigatti
// Summary: just some timings commented out
//
// Revision 1.29  2022/02/11 09:54:00  abbott
// Summary: Updated copyright (redmine 855)
//
// Revision 1.28  2018/09/28 15:54:06  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.27  2018/04/10 15:03:40  bigatti
// -- fixed includes
//
// Revision 1.26  2016/09/16 16:36:41  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.25  2015/05/07 15:00:40  abbott
// Summary: Moved test code into namespace CoCoA
// Author: JAA
//
// Revision 1.24  2014/09/05 16:13:58  abbott
// Summary: Commented out unused variable t0
// Author: JAA
//
// Revision 1.23  2014/07/07 13:40:28  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.22  2013/07/30 17:40:26  bigatti
// -- commented out temporary (printing) test
//
// Revision 1.21  2013/07/30 16:53:26  bigatti
// -- simplified input for hilbert
//
// Revision 1.20  2012/04/02 16:45:09  abbott
// Changed CoCoA_ASSERT into CoCoA_ASSERT_ALWAYS.
//
// Revision 1.19  2012/02/10 11:57:12  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.18  2012/02/08 17:36:47  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.17  2011/05/26 16:34:18  bigatti
// -- added test for ideal(0)
//
// Revision 1.16  2011/05/24 14:56:58  abbott
// Consequential changes from removal of several ctors for principal ideals.
//
// Revision 1.15  2011/04/27 09:48:37  bigatti
// -- updated syntax
//
// Revision 1.14  2011/04/26 10:32:21  bigatti
// -- added tests for multigraded case
//
// Revision 1.13  2011/04/08 14:07:37  bigatti
// -- renamed HilbertNumeratorMod into HilbertNumQuot
//
// Revision 1.12  2011/03/10 17:58:33  bigatti
// -- using long instead of size_t
//
// Revision 1.11  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.10  2010/10/29 09:43:09  bigatti
// -- manually freeing global memory for C implementation
//    to removed complaints when compiled with -DCoCoA_MEMPOOL_DEBUG
//
// Revision 1.9  2010/10/08 08:19:51  bigatti
// -- RingDistrMPoly.H --> RingDistrMPolyClean.H
//
// Revision 1.8  2009/07/30 15:41:25  bigatti
// -- now using the new nice constructors for ideals
//
// Revision 1.7  2007/10/30 17:14:05  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.6  2007/10/19 10:04:23  bigatti
// -- RingDenseUPolyClean now allow to specify the MinCapacity for all
//    coeff vectors (to avoid too many reallocations)
//
// Revision 1.5  2007/10/18 11:37:33  bigatti
// -- added Q9 for timings tests
//
// Revision 1.4  2007/10/15 12:59:58  bigatti
// -- added include
//
// Revision 1.3  2007/10/15 12:45:59  bigatti
// -- HP computed in Z[lambda] instead of Q[lambda]
//
// Revision 1.2  2007/10/10 14:40:48  bigatti
// -- added: test for comparing old C code with development version using CoCoAlib
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.9  2007/03/08 17:41:36  bigatti
// -- improved NewPolyRing
//
// Revision 1.8  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.7  2007/03/03 14:13:21  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.6  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.5  2007/02/26 17:11:58  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/10 18:44:02  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2007/01/17 17:38:11  bigatti
// -- moved all cocoa-4 code for hilbert into src/TmpHilbertDir
//
// Revision 1.2  2006/12/07 17:25:32  cocoa
// -- minimal set of #include's instead of library.H
//
// Revision 1.1  2006/11/16 18:16:10  cocoa
// -- added test for HilbertNumQuot
//
