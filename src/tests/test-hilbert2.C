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
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-hilbert.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/matrix.H"
#include "CoCoA/time.H"
#include "CoCoA/verbose.H"  // for SetVerbosityLevel


#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for Hilbert with weights
//----------------------------------------------------------------------
namespace CoCoA
{

  RingElem HilbertNumQuot_C_time(const ideal& I)
  {
    VerboseLog VERBOSE("HilbertNumQuot_C_time");
    //    double t0 = CpuTime();
    RingElem HNum = HilbertNumQuot_C(I);
    //    VERBOSE(20) << CpuTime() -t0 << endl;
    VERBOSE(30) << HNum << endl;    
    return HNum;
  }
  
  
  RingElem HilbertNumQuot_time(const ideal& I)
  {
    VerboseLog VERBOSE("HilbertNumQuot_time");
    //    double t0 = CpuTime();
    RingElem HNum = HilbertNumQuot(I);
    //    VERBOSE(20) << "  " << CpuTime() -t0 << endl;
    VERBOSE(30) << "  " << HNum << endl;    
    return HNum;
  }
  
  
  RingElem MGHilbertNumQuot_time(const ideal& I)
  {
    VerboseLog VERBOSE("MGHilbertNumQuot_time");
    //    double t0 = CpuTime();
    RingElem HNum = MGHilbertNumQuot(I);
    //    VERBOSE(20) << CpuTime() -t0 << endl;
    VERBOSE(30) << HNum << endl;    
    return HNum;
  }
  
  
// ----------------------------------------------------------------------

  void program()
  {
    GlobalManager CoCoAFoundations;
    VerboseLog VERBOSE("TestHilbert2");
    SetVerbosityLevel(0); // 30

    const PolyRing QQt = RingQQt(1);     // owner(HNum)  is  QQt
    RingElem t = indet(QQt, 0);
    
    const SparsePolyRing P = NewPolyRing(RingQQ(), symbols("x,y"));
    ideal I = ideal(indets(P));
    
    RingElem HNum = power((1-t), NumIndets(P));
    CoCoA_ASSERT_ALWAYS(num(HilbertSeriesQuot(I)) == HNum);
    CoCoA_ASSERT_ALWAYS(HilbertNumQuot_time(I) == HNum);
    CoCoA_ASSERT_ALWAYS(HilbertNumQuot_C_time(I) == HNum);
    CoCoA_ASSERT_ALWAYS(MGHilbertNumQuot_time(I) == HNum);

    vector<RingElem> W(2, RingElem(RingQQ(),3));
    const SparsePolyRing Pw = NewPolyRingWeights(RingQQ(), symbols("x,y"), RowMat(W));
    ideal I3 = ideal(indets(Pw));
    
    HNum = power((1-t*t*t), NumIndets(Pw));
    VERBOSE(30) << "---- IsStdGraded(P)?  " << IsStdGraded(P) << std::endl;
    VERBOSE(30) << "---- IsStdGraded(Pw)?  " << IsStdGraded(Pw) << std::endl;
    VERBOSE(30) << "-----------  " << HNum << std::endl;
    
    CoCoA_ASSERT_ALWAYS(num(HilbertSeriesQuot(I3)) == HNum);// OK
    CoCoA_ASSERT_ALWAYS(MGHilbertNumQuot_time(I3) == HNum); // OK  

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
