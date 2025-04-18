//   Copyright (c)  2023  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

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


#include "CoCoA/SparsePolyOps-eval.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyRing.H"

#include <iostream>
//#include <vector>
using std::vector;


namespace CoCoA
{

  namespace // anonymous
  {

    // These 2 functions actually do the work.


    std::vector<BigInt> ConvertToVectorBigInt(ConstRefRingElem f)
    {
      if (!IsSparsePolyRing(owner(f)) || !IsZero(characteristic(owner(f))))
        CoCoA_THROW_ERROR("Expected univariate poly with coeffs in ZZ", "EvalUPoly ctor");
        
      if (UnivariateIndetIndex(f) < 0)
        CoCoA_THROW_ERROR(ERR::ReqUnivariate, "EvalUPoly ctor");
      vector<BigInt> C(1+deg(f));
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
        C[deg(PP(it))] = ConvertTo<BigInt>(coeff(it));
      return C;
    }
    

    // --------------------------------------------
    // Evaluation of f in ZZ[x] at integer point

    // Qn:  Separate fns for eval at 1 and -1?

    // "binary" method for evaluation of poly in ZZ[x] at integer pt x=a
    BigInt EvalAt_bin(std::vector<BigInt> C, BigInt a)
    {
      if (IsZero(a))  return C[0]; // but we did make some wasteful copies...
      while (len(C) > 1)
      {
        if (IsOdd(len(C)))  C.push_back(BigInt(0));
        const long n = len(C);
        // HINT: allocate workspace for the product & use mpz_add, mpz_mul
        for (long i=0; i < n/2; ++i)
          C[i] = C[2*i]+a*C[2*i+1];  // maybe use mpz_add_mul  ??? worth it ???
        C.resize(n/2);
        a *= a;
      }
      return C[0];
    }



    // --------------------------------------------
    // Evaluation of f in ZZ[x] at rational point a/b -- returns numerator

    // "binary" method for evaluation of poly in ZZ[x] at rational pt x=a/b
    // !!RETURNS JUST THE NUMERATOR!!
    BigInt EvalAt_bin(std::vector<BigInt> C, BigInt a, BigInt b)
    {
      if (IsZero(b))  CoCoA_THROW_ERROR(ERR::DivByZero, "EvalUPoly/EvalAt");
      if (IsOne(b))  return EvalAt_bin(C, a);
      // HINT: allocate workspace for the products & use mpz_add, mpz_mul
      BigInt ExcessFactor(1);
      while (len(C) > 1)
      {
        if (IsOdd(len(C)))  { C.push_back(BigInt(0));  ExcessFactor *= b; }
        const long n = len(C);
        for (long i=0; i < n/2; ++i)
          C[i] = b*C[2*i]+a*C[2*i+1];
        C.resize(n/2);
        a *= a;
        b *= b;
      }
      return C[0]/ExcessFactor;
    }


  } // end of namespace anonymous



  EvalUPoly::EvalUPoly(ConstRefRingElem f):
      myCoeffs(ConvertToVectorBigInt(f))  // checks univariate, char=0, and each coeff is integer
  {}


  BigInt EvalUPoly::operator()(long a) const
  {
    return EvalAt_bin(myCoeffs, BigInt(a));
  }

  BigInt EvalUPoly::operator()(const BigInt& a) const
  {
    return EvalAt_bin(myCoeffs, a);
  }


  BigInt EvalUPoly::operator()(long n, long d) const
  {
    return EvalAt_bin(myCoeffs, BigInt(n),BigInt(d));
  }


  BigInt EvalUPoly::operator()(const BigInt& n, const BigInt& d) const
  {
    return EvalAt_bin(myCoeffs, n,d);
  }


  std::ostream& EvalUPoly::myOutput(std::ostream& out) const
  {
    if (!out)  return out;
    out << "EvalUPoly(deg=" << len(myCoeffs)-1 << ")";
    return out;
  }

  // Eval f(n)  for f in ZZ[x]
  BigInt EvalAt(ConstRefRingElem f, long n)
  {
    if (n == 0 || IsConstant(f))
      return ConvertTo<BigInt>(ConstantCoeff(f));
    return EvalUPoly(f)(n);
  }

  // Eval f(n)  for f in ZZ[x]
  BigInt EvalAt(ConstRefRingElem f, const BigInt& n)
  {
    if (IsZero(n) || IsConstant(f))
      return ConvertTo<BigInt>(ConstantCoeff(f));
    return EvalUPoly(f)(n);
  }


  // Evaluation of poly in ZZ[x] at rational pt x=n/d
  // !!RETURNS JUST THE NUMERATOR!!
  BigInt EvalAt(ConstRefRingElem f, long n, long d)
  {
    if (n == 0 || IsConstant(f))
      return ConvertTo<BigInt>(ConstantCoeff(f));
    if (d == 1)  return EvalUPoly(f)(n);
    return EvalUPoly(f)(n,d);
  }

  // Evaluation of poly in ZZ[x] at rational pt x=n/d
  // RETURNS NUMERATOR!!
  BigInt EvalAt(ConstRefRingElem f, const BigInt& n, const BigInt& d)
  {
    if (IsZero(n) || IsConstant(f))
      return ConvertTo<BigInt>(ConstantCoeff(f));
    if (IsOne(d))  return EvalUPoly(f)(n);
    return EvalUPoly(f)(n,d);
  }


  std::ostream& operator<<(std::ostream& out, const EvalUPoly& eval)
  {
    return eval.myOutput(out);
  }
  

} // end of namespace CoCoA
