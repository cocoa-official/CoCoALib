//   Copyright (c)  2017-2018, 2021  John Abbott,  Anna M. Bigatti
//   Original authors: Marvin Brandenstein, Alice Moallemy and Carsten Dettmar

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


#include "CoCoA/SparsePolyOps-ideal-RadicalMembership.H"

#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-ideal.H"
//#include "CoCoA/SparsePolyOps-hilbert.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/factor.H"
#include "CoCoA/ideal.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/verbose.H"

#include <vector>
using std::vector;


namespace CoCoA
{

  namespace // anonymous
  {

    vector<RingElem> GensForRabinovich(const ideal& I)
    {
      vector<RingElem> ans;
      const ring& P = RingOf(I);
      const PolyRing RabinovichPolyRing = NewPolyRing(CoeffRing(P), NewSymbols(NumIndets(P)+1));
      const vector<RingElem>& image = indets(RabinovichPolyRing);
      const RingHom Phi = PolyAlgebraHom(P, RabinovichPolyRing, vector<RingElem>(image.begin(),image.begin()+NumIndets(P)));
//      for (const RingElem& g: gens(I))  if (!IsZero(g)) ans.push_back(Phi(radical(g)));
//      for (const RingElem& g: GBasis(I))  ans.push_back(Phi(radical(g)));
      // radical(f) uses gcd and that is currently costly.  So I disabled it.
      // Try if it make a difference with QQ coeffs
      for (const RingElem& g: gens(I))  if (!IsZero(g)) ans.push_back(Phi((g)));
      for (const RingElem& g: GBasis(I))  ans.push_back(Phi((g)));
      return ans;
    }


    bool TryPowersModI(ConstRefRingElem f, const ideal& I, long MaxPower)
    {
      VerboseLog VERBOSE("TryPowersModI");
      RingElem g = NF(f,I);
      long d=1;
      for (; d<MaxPower || !IsZero(g); d*=2)
        g = NF(power(g,2), I);
      VERBOSE(40) << "power = " << d << std::endl;
      return IsZero(g);
    }
    

//    bool IsInRadical0Dim(ConstRefRingElem f, const ideal& I, long& LenQB)


    bool RabinovichTrick(ConstRefRingElem f, const ideal& I, vector<RingElem>& G)
    {
      if (IsElem(f,I)) return true; // simple special case
      long L;
      if (IsZeroDim(I) && (L=len(QuotientBasis(I)))<513) return TryPowersModI(f,I,L); //ConvertTo<long>(MultiplicityQuot(I)));
      const PolyRing& P = RingOf(I);
      if (G.empty())    G = GensForRabinovich(I);
      PolyRing RabinovichPolyRing = owner(G[0]);
      const vector<RingElem>& image = indets(RabinovichPolyRing);
      const RingHom Phi = PolyAlgebraHom(P, RabinovichPolyRing, vector<RingElem>(image.begin(),image.begin()+NumIndets(P)));
      vector<RingElem> newGens = G;
      newGens.push_back(1 - Phi(f) * indet(RabinovichPolyRing, NumIndets(P)));
      return IsOne(ideal(newGens)); 
    }


  } // end of namespace anonymous
  

  bool IsInRadical(ConstRefRingElem f, const ideal& I)
  {
    const char* const FnName = "IsInRadical";
    if (owner(f) != RingOf(I))
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (!IsPolyRing(RingOf(I)) || !IsField(CoeffRing(RingOf(I))))
      CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName); // BUG BUG err should say "polyring over field"
    if (IsZero(I))  return IsZero(f);  // SparsePolyRing

    vector<RingElem> GensI_Rabinovich; // filled, only if needed, by RabinovichTrick
    if (GradingDim(RingOf(I))==0 || !IsHomog(I))
      return RabinovichTrick(f, I, GensI_Rabinovich);
// should the following be a different, homgeneous, algorithm?
    for (RingElem g=f; !IsZero(g); /*CutLF(g) inside loop body*/)
      if (!RabinovichTrick(CutLF(g), I, GensI_Rabinovich))
        return false;
    return true;
  }

  
  bool IsInRadical(const ideal& I, const ideal& J)
  {
    vector<RingElem> GensJ_Rabinovich; // filled, only if needed, by RabinovichTrick
    for (const RingElem& g: gens(I))
      if (!RabinovichTrick(g, J, GensJ_Rabinovich))
        return false;
    return true;
  }



// // MinPowerInIdeal: Determines the minimal power m such that f^m is in the Ideal J
//   long MinPowerInIdeal_naive(ConstRefRingElem f, const ideal& I)
//   {
//     const char* const FnName = "MinPowerInIdeal";
//     if (owner(f) != RingOf(I))
//       CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
//     if (!IsPolyRing(RingOf(I)) || !IsField(CoeffRing(RingOf(I))))
//       CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName); // BUG BUG err should say "polyring over field"

//     if (!IsInRadical(f, I)) return -1;
//     long D = 1;
//     RingElem FpowD = f;
//     while (!IsElem(FpowD, I))
//     {
//       CheckForInterrupt(FnName);
//       FpowD *= f;
//       ++D;
//     }
//     return D;
//   }
  

// not sure if this is really any better than the sequential version below
  long MinPowerInIdeal_binary(RingElem f, const ideal& I)
  {
    const char* const FnName = "MinPowerInIdeal_binary";
    if (owner(f) != RingOf(I))
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (!IsPolyRing(RingOf(I)) || !IsField(CoeffRing(RingOf(I))))
      CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName); // BUG BUG err should say "polyring over field"
//JAA: suggest clearing denoms (before & after NF?)
    f = NF(f,I);
    if (IsZero(f)) return 1;
    if (!IsInRadical(f, I)) return -1;  // WARNING: test can be slow/costly

    long D = 1;
    vector<RingElem> powers;
    RingElem FpowD = f;
    while (!IsZero(FpowD))
    {
      CheckForInterrupt(FnName);
      powers.push_back(FpowD);
      FpowD = NF(FpowD*FpowD, I);
      D *= 2;
    }
    D /= 2;
    long CurrPower = D;
    FpowD = powers.back();
    int n = len(powers)-1;
    while (D > 1)
    {
      CoCoA_ASSERT(n >= 1);
      --n; D /= 2;
      RingElem tmp = NF(FpowD*powers[n], I);
      if (IsZero(tmp)) continue;
      swap(FpowD, tmp);
      CurrPower += D;
    }
    return 1+CurrPower;
  }



// Exercise 5b): Alternative for MinPowerInIdeal using normal form
  // Just compute powers sequentially until we get NF=0.
  long MinPowerInIdeal(RingElem f, const ideal& I)
  {
    const char* const FnName = "MinPowerInIdeal";
    if (owner(f) != RingOf(I))
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (!IsPolyRing(RingOf(I)) || !IsField(CoeffRing(RingOf(I))))
      CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName); // BUG BUG err should say "polyring over field"
//JAA: suggest clearing denoms (before & after NF?)
    f = NF(f,I);
    if (!IsInRadical(f, I)) return -1;  // WARNING: test can be slow/costly
    // simple sequential loop
    long D = 1;
    RingElem FpowD = f;
    while (!IsZero(FpowD))
    {
      CheckForInterrupt(FnName);
      FpowD = NF(FpowD*f, I);
      ++D;
    }
    return D;
  }



  
  // // MinPowerInIdealH (for homogeneous ideals): For every homogeneous
  // // component f_i, it determines the minimal exponent m_i, such that
  // // (f_i)^(m_i) is in the ideal J.
  // // Question: Knowing that, can we derive the minimal exponent m, s.t. f^m is in J?
  // vector<long> MinPowerInIdealH(ConstRefRingElem f, const ideal& J)
  // {
  //   vector<long> presumedExponent;
  //   vector<long> homogExponents; // alternative
  //   if (!IsInRadical(f, J))
  //   {
  //     presumedExponent.push_back(-1);
  //   }
  //   else if (!IsHomog(J))
  //   {
  //     presumedExponent.push_back(MinPowerInIdeal(f, J));
  //   }
  //   // f is in the radical of J, J is a radical ideal
  //   else {
  //     vector<RingElem> homogComp = HomogComp(f);	 

  //     for (int i = 0; i < len(homogComp); i++)
  //     {     
  //       homogExponents.push_back(MinPowerInIdeal(homogComp[i], J)); // alternative
  //     }

  //     // first conjecture: the minimal power m such that the polynomial f^m is in the Ideal J
  //     // equals the minimal power such that every homogenous component of f is in the Ideal J
  //     // second conjecture: the minimal power m such that the polynomial f^m is in the Ideal J
  //     // equals the product of the minimal powers for the homogenous components of f
  //     // such that every homogenous component of f is in the Ideal J
  //     // third conjecture: the minimal power m such that the polynomial f^m is in the Ideal J
  //     // equals the lcm of the minimal powers of the homogenous components 
  //     // such that every homogenous component of f is in the Ideal J
  //     // turns out: conjectures are wrong!
  //     // we could not find a simple connection between those powers
  //     presumedExponent.push_back(MaxElem(homogExponents));
  //     presumedExponent.push_back(Prod(homogExponents));
  //     presumedExponent.push_back(Lcm(homogExponents));
  //     cout << "Here the vector of powers for the homogenous components of the polynomial: " << homogExponents << endl; // alternative     
  //   }
  //   return presumedExponent;
  // }


} // end of namespace CoCoA
