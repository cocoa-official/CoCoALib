//   Copyright (c)  1999,2009-2011  John Abbott & Anna M Bigatti

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


#include "CoCoA/NumTheory-CRT.H"

#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
#include <vector>
using std::vector;

namespace CoCoA
{

  CRTMill::CRTMill():
      myR(/*0*/),
      myM(1)
  {}


  void CRTMill::myAddInfo(const MachineInt& res, const MachineInt& mod, CoprimeFlag check)
  {
    if (IsNegative(mod))  CoCoA_THROW_ERROR1(ERR::BadModulus);
    if (!IsSignedLong(mod) || !IsSignedLong(res))
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    const long m = AsSignedLong(mod);
    if (check == CheckCoprimality && gcd(m,myM) != 1)
      CoCoA_THROW_ERROR2(ERR::BadArg, "new modulus not coprime");
    CoCoA_ASSERT(IsCoprime(m,myM));
    long r = uabs(res)%m;
    if (r != 0 && IsNegative(res))  r = m-r;
    const long a = myR%m;
    long k;
    if (m <= MaxSquarableInteger<long>())
      k = SymmRemainder((r-a)*InvMod(myM,m),m); // if both (r-a) & InvMod(..) are SymmRem then can allow m <= 2*MaxSquarableInteger
    else
      k = SymmRemainder(BigInt(r-a)*InvMod(myM,m),m); // use BigInt to avoid possible overflow (don't care about speed)
    myR += k*myM;
    myM *= m;
  }

  void CRTMill::myAddInfo(const BigInt& r, const BigInt& m, CoprimeFlag check)
  {
//???JAA    CoCoA_ASSERT(m > 1);
    if (IsOne(m))  return; //???JAA
    if (IsOne(m))  { if (!IsOne(myM))  CoCoA_THROW_ERROR1(ERR::BadModulus); return; }
    CoCoA_ASSERT((r < m) && (r > -m));
    if (check == CheckCoprimality && gcd(m,myM) != 1)
      CoCoA_THROW_ERROR2(ERR::BadArg, "new modulus not coprime");
    CoCoA_ASSERT(IsCoprime(m,myM));
    const BigInt a = myR%m;
    const BigInt k = SymmRemainder((r-a)*InvMod(myM,m), m); ///???BUG/SLUG??? SymmRemainder(r-a,m) ???
    myR += k*myM;
    myM *= m;
  }



  std::ostream& operator<<(std::ostream& out, const CRTMill& CRT)
  {
    if (!out)  return out;  // short-cut for bad ostreams

    out << "CRT(residue=" << CRT.myR << ", modulus=" << CRT.myM << ")";
    return out;
  }



} // end of namespace CoCoA
