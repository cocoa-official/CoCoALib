#ifndef CoCoA_NumTheory_gcd_H
#define CoCoA_NumTheory_gcd_H

//   Copyright (c)  2009,2012,2020  John Abbott  &  Anna M. Bigatti

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


#include <vector>
// using std::vector;

namespace CoCoA
{
  
  class MachineInt; // fwd decl -- defined in MachineInt.H

  long   gcd(const MachineInt& a, const MachineInt& b) noexcept;
  BigInt gcd(const BigInt& A,     const MachineInt& b);
  BigInt gcd(const MachineInt& a, const BigInt& B);
  BigInt gcd(const BigInt& A,     const BigInt& B);
  long   gcd(const std::vector<long>& L) noexcept;
  BigInt gcd(const std::vector<BigInt>& L);
  bool IsCoprime(const MachineInt& a, const MachineInt& b) noexcept;
  bool IsCoprime(const BigInt& A,     const MachineInt& b);  // noexcept (with better impl)
  bool IsCoprime(const MachineInt& a, const BigInt& B);      // noexcept (with better impl)
  bool IsCoprime(const BigInt& A,     const BigInt& B);


  long ExtGcd(long& CofacA, long& CofacB, const MachineInt& a, const MachineInt& b);
//???? signature???  BigInt ExtGcd(BigInt& CoFacA, BigInt& CofacB, const MachineInt& A, const BigInt& B);
//???? signature???  BigInt ExtGcd(BigInt& CoFacA, BigInt& CofacB, const BigInt& A, const MachineInt& B);
  BigInt ExtGcd(BigInt& CoFacA, BigInt& CofacB, const BigInt& A, const BigInt& B);

  struct GcdAndCofacs
  {
    GcdAndCofacs(const BigInt& g, const std::vector<BigInt>& cofacs);
    BigInt myGcd;  // just for this field, I must include BigInt.H... grrr!
    std::vector<BigInt> myCofacs;
  };

  GcdAndCofacs ExtGcd(const std::vector<BigInt>& L);

  long   lcm(const MachineInt& a, const MachineInt& b); // NOT noexcept (arg check)
  BigInt lcm(const BigInt& A,     const MachineInt& b);
  BigInt lcm(const MachineInt& a, const BigInt& B);
  BigInt lcm(const BigInt& a,     const BigInt& b);
  BigInt lcm(const std::vector<long>& L);
  BigInt lcm(const std::vector<BigInt>& L);

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-gcd.H,v 1.6 2022/02/18 14:11:46 abbott Exp $
// $Log: NumTheory-gcd.H,v $
// Revision 1.6  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2022/02/08 20:20:45  abbott
// Summary: Corrected GcdAndCofacs (but had to add include directive).
//
// Revision 1.4  2021/08/04 19:06:00  abbott
// Summary: Removed const (redmine 1606)
//
// Revision 1.3  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.2  2020/02/28 08:57:08  abbott
// Summary: Added gcd & lcm for vector of values; also ExtGcd (prototype)
//
// Revision 1.1  2020/01/26 14:14:31  abbott
// Summary: Finished splitting NumTheory into smaller pieces (redming 1161)
//
//
//

#endif
