//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/ring.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/RingDenseUPolyClean.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;

namespace CoCoA
{

  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    GlobalManager CoCoAFoundations;

    ring R1 = RingZZ();
    ring R2 = RingQQ();
    SparsePolyRing SPR1 = NewPolyRing(R2, symbols("x"));
    SparsePolyRing SPR2 = NewPolyRing(R2, symbols("y"));
    DenseUPolyRing DUP1 = NewPolyRing_DUP(R2, symbol("x"));    
    DenseUPolyRing DUP2 = NewPolyRing_DUP(R2, symbol("y"));
    FractionField FrF1 = RingQQ();
    FractionField FrF2 = NewFractionField(SPR1);
    PolyRing PR1 = NewPolyRing(R2, symbols("z"));
    PolyRing PR2 = NewPolyRing(R2, symbols("w"));
    QuotientRing QR1 = NewZZmod(3);
    QuotientRing QR2 = NewZZmod(5);
    RingTwinFloat RTF1 = NewRingTwinFloat(16);
    RingTwinFloat RTF2 = NewRingTwinFloat(32);

    // Assign to ring
    ring R = R2;
    R = R1;
    R = SPR1;
    R = DUP1;
    R = FrF1;
    R = PR1;
    R = QR1;
    R = RTF2;

    // Assign to PolyRing
    PolyRing PR = PR1;
    PR = PR2;
    PR = SPR1;
    PR = DUP1;

    // Assign to FractionField
    FractionField FrF = FrF1;
    FrF = FrF2;

    // Assign to QuotientRing
    QuotientRing QR = QR1;
    QR = QR2;

    // Assign to SparsePolyRing
    SparsePolyRing SPR = SPR1;
    SPR = SPR2;

    // Assign to DenseUPolyRing
    DenseUPolyRing DUP = DUP1;
    DUP = DUP2;

    // Assign to RingTwinFloat
    RingTwinFloat RTF = RTF1;
    RTF = RTF2;
  }

} // end of namespace CoCoA


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
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
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
