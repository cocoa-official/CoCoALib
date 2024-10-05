#ifndef CoCoA_NumTheory_RatReconstruct_H
#define CoCoA_NumTheory_RatReconstruct_H

//   Copyright (c)  2009,2012  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/BigRat.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/NumTheory-CRT.H"

namespace CoCoA
{
  
///  class MachineInt; // fwd decl -- defined in MachineInt.H

  //////////////////////////////////////////////////////////////////
  // Heuristic fault-tolerant rational reconstruction

  class RatReconstructByContFrac
  {
  public:
    explicit RatReconstructByContFrac(MachineInt LogEps=20);
//    explicit RatReconstructByContFrac(const BigInt& threshold);
    void myAddInfo(const MachineInt& r, const MachineInt& m);
    void myAddInfo(const BigInt& R, const BigInt& M);
    friend const BigRat& ReconstructedRat(const RatReconstructByContFrac& reconstructor);
    friend bool IsConvincing(const RatReconstructByContFrac& reconstructor);
    friend const BigInt& BadMFactor(const RatReconstructByContFrac& reconstructor);
    friend std::ostream& operator<<(std::ostream& out, const RatReconstructByContFrac& reconstructor);
  private: // data members
    CRTMill myCRT;
//    const BigInt myThresholdValue;
    const long myLogEps;
    mutable bool myResultIsUpToDate;
    mutable bool myResultIsConvincing;
    mutable BigRat myResult;
    mutable BigInt myBadFactor;
  private: // impl details
    BigInt myThreshold(const BigInt& M) const;
    static BigInt myCheckThreshold(const BigInt& threshold);
    void myUpdateResult() const;
  };

  const BigRat& ReconstructedRat(const RatReconstructByContFrac& ModularImages);
  bool IsConvincing(const RatReconstructByContFrac& ModularImages);
  const BigInt& BadMFactor(const RatReconstructByContFrac& ModularImages);
  std::ostream& operator<<(std::ostream& out, const RatReconstructByContFrac& reconstructor);

  //////////////////////////////////////////////////////////////////
  // Heuristic fault-tolerant rational reconstruction

  class RatReconstructByLattice
  {
  public:
    explicit RatReconstructByLattice(const MachineInt& SafetyFactor); // SafetyFactor==0 --> use default value
    explicit RatReconstructByLattice(const BigInt& SafetyFactor);
    void myAddInfo(const MachineInt& r, const MachineInt& m);
    void myAddInfo(const BigInt& r, const BigInt& m);
    friend const BigRat& ReconstructedRat(const RatReconstructByLattice& ModularImages);
    friend bool IsConvincing(const RatReconstructByLattice& ModularImages);
    friend const BigInt& BadMFactor(const RatReconstructByLattice& ModularImages);
    friend std::ostream& operator<<(std::ostream& out, const RatReconstructByLattice& reconstructor);
  private: // data members
    CRTMill myCRT;
    const BigInt mySafetyFactor;
    mutable bool myResultIsUpToDate;
    mutable bool myResultIsConvincing;
    mutable BigRat myResult;
    mutable BigInt myBadFactor;
  private: // impl details
    static const long ourDefaultSafetyFactor;
    static BigInt myCheckSafetyFactor(const BigInt& SafetyFactor);
    void myUpdateResult() const;
  };

  const BigRat& ReconstructedRat(const RatReconstructByLattice& ModularImages);
  bool IsConvincing(const RatReconstructByLattice& ModularImages);
  const BigInt& BadMFactor(const RatReconstructByLattice& ModularImages);
  std::ostream& operator<<(std::ostream& out, const RatReconstructByLattice& reconstructor);

  BigRat RatReconstructWithBounds(long e, const BigInt& P, const BigInt& Q, const std::vector<long>& res, const std::vector<long>& mod);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-RatReconstruct.H,v 1.5 2022/02/18 14:11:46 abbott Exp $
// $Log: NumTheory-RatReconstruct.H,v $
// Revision 1.5  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.4  2021/01/07 14:17:52  abbott
// Summary: Corrected copyright
//
// Revision 1.3  2020/01/26 14:24:33  abbott
// Summary: Revised includes after splitting NumTheory
//
// Revision 1.2  2019/03/18 11:58:09  abbott
// Summary: Added missing includes
//
// Revision 1.1  2019/03/18 11:24:19  abbott
// Summary: Split NumTheory into several smaller files
//
//
//

#endif
