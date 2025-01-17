#ifndef CoCoA_NumTheory_CRT_H
#define CoCoA_NumTheory_CRT_H

//   Copyright (c)  2009,2012  John Abbott  &  Anna M Bigatti

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

#include "CoCoA/BigInt.H"


namespace CoCoA
{
  
  class MachineInt; // fwd decl -- defined in MachineInt.H


  class CRTMill
  {
  public:
    CRTMill();
    // default copy ctor and dtor are fine
    enum class CoprimeFlag { CheckCoprimality=0, CoprimeModulus=1};
    static constexpr CoprimeFlag CheckCoprimality = CoprimeFlag::CheckCoprimality;
    static constexpr CoprimeFlag CoprimeModulus = CoprimeFlag::CoprimeModulus;
    void myAddInfo(const MachineInt& r, const MachineInt& m, CoprimeFlag check=CheckCoprimality);
    void myAddInfo(const BigInt& r, const BigInt& m, CoprimeFlag check=CheckCoprimality);
//    const BigInt& myResidue() const { return myR; }
//    const BigInt& myModulus() const { return myM; }
    friend const BigInt& CombinedResidue(const CRTMill& CRT) noexcept { return CRT.myR; }
    friend const BigInt& CombinedModulus(const CRTMill& CRT) noexcept { return CRT.myM; }
    friend std::ostream& operator<<(std::ostream& out, const CRTMill& CRT);
  private: // data members
    BigInt myR;
    BigInt myM;
  };

//???  BigInt CRT(const std::vector<long>& R, const std::vector<long>& M);
  const BigInt& CombinedResidue(const CRTMill& CRT) noexcept;
  const BigInt& CombinedModulus(const CRTMill& CRT) noexcept;
  std::ostream& operator<<(std::ostream& out, const CRTMill& CRT);

} // end of namespace CoCoA

#endif
