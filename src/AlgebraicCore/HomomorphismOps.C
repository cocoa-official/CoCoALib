//   Copyright (c)  2017  John Abbott, Anna Bigatti
//   Original authors: Thomas Izgin and Filip Skrentny

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

#include "CoCoA/HomomorphismOps.H"

#include "CoCoA/PolyRing.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/CanonicalHom.H"
///#include "CoCoA/apply.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixView.H" // for RowMat
#include "CoCoA/VectorOps.H"
#include "CoCoA/RingZZ.H" // for RingZZ()
#include "CoCoA/verbose.H"

#include <iostream>
using std::ostream;
#include <vector>
using std::vector;

#include <string>
using std::string;

namespace CoCoA
{

  namespace // anonymous
  {

    // This fn checks whether R is a PolyRing or QuotientRing (with a polynomial BaseRing)
    bool IsPolyRingOrQuotient(const ring& R)
    {
      return ( (IsPolyRing(R) && IsField(CoeffRing(R))) ||
               (IsQuotientRing(R) && IsPolyRing(BaseRing(R)) && IsField(CoeffRing(BaseRing(R)))) );
    }


    vector<symbol> concat(const vector<symbol>& L1, const vector<symbol>& L2) // sticks two lists L1 and L2 together to one bigger list
    {
      vector<symbol> ans = L1; 
      ans.insert(ans.end(), L2.begin(), L2.end());
      // for (int i=0; i < len(L2); ++i)
      // {
      //   ans.push_back(L2[i]);
      // }
      return ans;
    }


    // SOONER OR LATER.... following struct is to be moved into RingHom (for SparsePolyRing)

    struct RichHom // this function will help us return something in MakeRichHom
    {
    public:
      RichHom(const RingHom& phi, const ring& Ri, const ring& Si, const vector<RingElem>& SIndetsRSi, const PPMonoidElem& MinSIndetRSi, const ideal& idealRSi, const RingHom& HomS_RSi, const RingHom& HomRS_Ri);
    public:
//  private: // data members
      RingHom Hom;
      ring R; // const ref???
      ring S; // const ref???
      vector<RingElem> SIndetsRS;
      PPMonoidElem MinSIndetRS;
      ideal idealRS;
      RingHom HomS_RS;
      RingHom HomRS_R;
    };

    inline RichHom::RichHom(const RingHom& phi, const ring& Ri, const ring& Si, const vector<RingElem>& SIndetsRSi, const PPMonoidElem& MinSIndetRSi, const ideal& idealRSi, const RingHom& HomS_RSi, const RingHom& HomRS_Ri):
        Hom(phi),
        R(Ri),
        S(Si),
        SIndetsRS(SIndetsRSi),
        MinSIndetRS(MinSIndetRSi),
        idealRS(idealRSi),
        HomS_RS(HomS_RSi),
        HomRS_R(HomRS_Ri)
    {}


    // pseudo-ctor
    RichHom MakeRichHom(const RingHom& Phi) // returns many things (look for RichHom to see in detail, what will be returned)
    {
///      using ::CoCoA::apply; // template fn
      VerboseLog VERBOSE("MakeRichHom");
//      VERBOSE(90) << "starting..." << std::endl;
      CoCoA_ASSERT(IsPolyRingOrQuotient(domain(Phi)));
      CoCoA_ASSERT(IsPolyRingOrQuotient(codomain(Phi)));
      ring R;
      if(IsPolyRing(domain(Phi)))
        R = domain(Phi);
      else
        R = BaseRing(domain(Phi));
//      VERBOSE(90) << "inside.." << std::endl;
      ring S;
      ideal I2 = ideal(zero(S));
      if (IsPolyRing(codomain(Phi)))
      {
        S = codomain(Phi);
        I2 = ideal(zero(S));
      } 
      else
      {
        S = BaseRing(codomain(Phi));
        I2 = DefiningIdeal(codomain(Phi));
      } 
      const int NumS = NumIndets(S);
      const int NumR = NumIndets(R);

      // RS = R[ SX[0],...,SX[NumS-1],RX[0],..., RX[NumR-1] ]
      // with an elim.ordering for the SX[i], i=0,..., NumS-1
      PolyRing RS = NewPolyRing(CoeffRing(R),
                                concat(SymbolRange("SX",0,NumS-1), SymbolRange("RX",0,NumR-1)),
                                NewMatrixOrdering(ElimMat(LongRange(0,NumS-1),NumR+NumS),0));

      const vector<RingElem>  SIndetsRS = indets(RS, "SX");
      PPMonoidElem MinSIndetRS = LPP(SIndetsRS[0]);
      for (const auto& x: SIndetsRS)
        if (LPP(x) < MinSIndetRS)  MinSIndetRS = LPP(x);
      const vector<RingElem> RIndetsRS = indets(RS, "RX");   
      const RingHom HomS_RS = PolyAlgebraHom(S, RS, SIndetsRS); //S->RS   
      vector<RingElem> imgs;
      for (int i=0; i< NumS; ++i) { imgs.push_back(zero(R)); }
      for (int i=0; i < NumIndets(R); ++i) { imgs.push_back(indet(R,i)); }
      const RingHom HomRS_R = PolyAlgebraHom(RS, R, imgs); //RS->R
      vector<RingElem> X = indets(R);  // X in R
      if (IsQuotientRing(domain(Phi)))
      {
        X = CanonicalHom(R, domain(Phi))(X); //X in domain(Phi)
      }
      vector<RingElem> PhiX = Phi(X); // PhiX in codomain(Phi)
    
      if (IsQuotientRing(codomain(Phi)))
      {
        for (int i=0; i<len(PhiX); ++i)
          PhiX[i] = CanonicalRepr(PhiX[i]);
      }
    
      vector<RingElem> PhiXRS = HomS_RS(PhiX); // PhiXRS in RS

      vector<RingElem> Gens;
      for (int i=0; i<NumR; ++i)
        Gens.push_back(RIndetsRS[i]-PhiXRS[i]);

      vector<RingElem> geni2rs;
      for (int i=0; i<len(gens(I2)); ++i)
        geni2rs.push_back(HomS_RS(gens(I2)[i]));

      const ideal J = ideal(RS,geni2rs) + ideal(RS,Gens);
       
      // storing R and S helps with ref-counting?
      VERBOSE(90) << "returning: ideal J is " << J << std::endl;
      return  RichHom(Phi, R,S,SIndetsRS,MinSIndetRS, J, HomS_RS, HomRS_R);
    }

  
    RichHom MakeRichHom_H(const RingHom& Phi) // returns many things (look for RichHom to see in detail, what will be returned)
    {
      VerboseLog VERBOSE("MakeRichHom_H");
      CoCoA_ASSERT(IsPolyRingOrQuotient(domain(Phi)));
      CoCoA_ASSERT(IsPolyRingOrQuotient(codomain(Phi)));
      ring R;
      if(IsPolyRing(domain(Phi)))
        R = domain(Phi);
      else
        R = BaseRing(domain(Phi));
      for (long i=0; i<NumIndets(R); ++i)
        if (deg(Phi(indet(R,i)))==0)
          CoCoA_THROW_ERROR2(ERR::NYI, "ker_H not yet implemented for constant images");
 
      ring S;
      ideal I2 = ideal(zero(S));
      if (IsPolyRing(codomain(Phi)))
      {
        S = codomain(Phi);
        I2 = ideal(zero(S));
      } 
      else
      {
        S = BaseRing(codomain(Phi));
        I2 = DefiningIdeal(codomain(Phi));
      } 
      const int NumS = NumIndets(S);
      const int NumR = NumIndets(R);
      std::vector<RingElem> W;
      RingElem OneZZ(one(RingZZ()));
      for (long i=0; i<NumIndets(S); ++i) W.push_back(OneZZ);
      for (long i=0; i<NumIndets(R); ++i) W.push_back(OneZZ*deg(Phi(indet(R,i))));
      W.push_back(OneZZ); // HomogIndet

      // SRh = R[ SX[0],...,SX[NumS-1], RX[0],..., RX[NumR-1] , h]
      // with an elim.ordering for the SX[i], i=0,..., NumS-1

      PolyRing SRh = NewPolyRing(CoeffRing(R),
                                 concat(concat(SymbolRange("t",0,NumS-1),
                                               SymbolRange("x",0,NumR-1)),
                                        symbols("HomogIndet")),
                                 NewMatrixOrdering(ElimHomogMat(LongRange(0,NumS-1), RowMat(W)),1));
      VERBOSE(90) << "Ring for hom: OrdMat is "
                  << OrdMat(ordering(PPM(SRh))) << std::endl;
      const vector<RingElem>  SIndetsSRh = indets(SRh, "t");
      PPMonoidElem MinSIndetSRh = LPP(SIndetsSRh[0]);
      for (int i=1; i < len(SIndetsSRh); ++i)
      {
        if (LPP(SIndetsSRh[i]) < MinSIndetSRh)
          MinSIndetSRh = LPP(SIndetsSRh[i]);
      }
      const vector<RingElem> RIndetsSRh = indets(SRh, "x");   
      const RingHom HomS_SRh = PolyAlgebraHom(S, SRh, SIndetsSRh); //S->SRh   
      vector<RingElem> imgs;
      for (int i=0; i < NumS; ++i)  imgs.push_back(zero(R));
      for (int i=0; i < NumR; ++i)  imgs.push_back(indet(R,i));
      imgs.push_back(one(R));  // HomogIndet
      const RingHom HomSRh_R = PolyAlgebraHom(SRh, R, imgs); //SRh->R
      vector<RingElem> X = indets(R);  // X in R
      if (IsQuotientRing(domain(Phi)))
        X = CanonicalHom(R, domain(Phi))(X); // X in domain(Phi)
      vector<RingElem> PhiX = Phi(X); // Phi(X) in codomain(Phi)
    
      if (IsQuotientRing(codomain(Phi)))
        for (int i=0; i<len(PhiX); ++i) PhiX[i] = CanonicalRepr(PhiX[i]);
    
      vector<RingElem> PhiXSRh = HomS_SRh(PhiX); // PhiXSRh in SRh

      vector<RingElem> Gens;
      for (int i=0; i<NumR; ++i) Gens.push_back(RIndetsSRh[i]-PhiXSRh[i]);

      vector<RingElem> gensI2SRh = HomS_SRh(gens(I2));
         
      const ideal J = ideal(SRh,gensI2SRh) + ideal(SRh,Gens);
       
      // storing R and S helps with ref-counting?
      return  RichHom(Phi, R,S,SIndetsSRh,MinSIndetSRh, J, HomS_SRh, HomSRh_R);
    }

  
    // This fn checks the input: it must be a RingHom between poly rings or their quotients.
    // FnName represents the function,that is calling CheckInputAffAlgebraHom;
    // look below for examples producing errors
    void CheckInputAffAlgebraHom(const RingHom& Phi, const ErrorContext& ErrCtx)
    {
      ring R = domain(Phi);
      ring S = codomain(Phi);
      if (!IsPolyRingOrQuotient(R) || !IsPolyRingOrQuotient(S))
        CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadArg, "RINGHOM must be between poly rings or their quotients", ErrCtx);
      if (!IsPolyRing(R))  R = BaseRing(R);
      if (!IsPolyRing(S))  S = BaseRing(S);
      if (!IsField(CoeffRing(R)) ||  (CoeffRing(R) != CoeffRing(S)))
        CoCoA_THROW_ERROR_WITH_CONTEXT3(ERR::BadArg, "RINGHOM must be between quotients of K-algebras", ErrCtx);
    }


    // This fn calculates the kernel of RichPhi (  ideal("reduced gröbner basis")  ),
    // in order to understand the last step of this function you might uncomment the
    // couts and look what the function does in our examples (e.g. psi2)
    ideal KerRichHom(const RichHom& RichPhi)
    {
      VerboseLog VERBOSE("KerRichHom");
//      PPMonoidElem ElimX = LPP(product(RichPhi.SIndetsRS));
      const vector<RingElem> RGB = ReducedGBasis(RichPhi.idealRS);
      vector<RingElem> GensKer; GensKer.reserve(len(RGB)); // may be needlessly big
      for (const auto& g: RGB)
        if (LPP(g) < RichPhi.MinSIndetRS) // tests whether LPP "is in" R (elimord)
          GensKer.push_back(RichPhi.HomRS_R(g));

      const ring& DomainPhi = domain(RichPhi.Hom);
      VERBOSE(90) << "about to make ideal(GensKer)" << std::endl;
      if ( GensKer.empty() || (!IsQuotientRing(DomainPhi)) )
        return ideal(DomainPhi, GensKer);

      // Next block of code is to eliminate some redundant generators from the answer.
      // Reduce each gen to NF modulo DefiningIdeal, then compute RGB of these polys.
      const int LenGensKer = len(GensKer);
      // Reduce each elem of GensKer to NF modulo DefiningIdeal
      for (int i=0; i < LenGensKer; ++i)
        GensKer[i] = NF(GensKer[i], DefiningIdeal(DomainPhi));
      vector<RingElem> RGBJ = ReducedGBasis(ideal(GensKer));
      GensKer.clear();  // Re-use GensKer to contain the final answer
      RingHom ImageInQuot = QuotientingHom(DomainPhi);
      for (const auto& g: RGBJ)
        GensKer.push_back(ImageInQuot(g));
      return ideal(DomainPhi, GensKer);
    }  


    ideal KerRichHom_H(const RichHom& RichPhi)
    {
      VerboseLog VERBOSE("KerRichHom_H");
//      PPMonoidElem ElimX = LPP(product(RichPhi.SIndetsRS));
      const ideal I = RichPhi.idealRS;
      //      std::cout << I << std::endl;
      
      const RingElem h(RingOf(I), "HomogIndet");
      const ring& P_h = owner(h); // same as RingOf(I)
      std::vector<RingElem> gensH;
      for (const auto& f: gens(I))
        if (!IsZero(f))
          gensH.push_back(homog(f, h));
      //      std::cout << gensH << std::endl;
      VERBOSE(90) << "about to make RGBh" << std::endl;
      const vector<RingElem> RGBh = ReducedGBasis(ideal(P_h,gensH));
      //      std::cout << RGBh << std::endl;

      int LenRGBh = len(RGBh);
      vector<RingElem> GensKer; GensKer.reserve(LenRGBh); // may be needlessly big

      PPMonoidElem ProdElim=LPP(product(RichPhi.SIndetsRS));
      
      for (int i=0; i < LenRGBh; ++i)
        if (IsCoprime(ProdElim,LPP(RGBh[i])))
          GensKer.push_back(RichPhi.HomRS_R(RGBh[i]));
    
      const ring& DomainPhi = domain(RichPhi.Hom);
      if (!IsQuotientRing(DomainPhi))
        return ideal(DomainPhi, GensKer);

      // Next block of code is to eliminate some redundant generators from the answer.
      // Reduce each gen to NF modulo DefiningIdeal, then compute RGBh of these polys.
      const int LenGensKer = len(GensKer);
      // Reduce each elem of GensKer to NF modulo DefiningIdeal
      for (int i=0; i < LenGensKer; ++i)
        GensKer[i] = NF(GensKer[i], DefiningIdeal(DomainPhi));
      VERBOSE(90) << "about to make RGBhJ" << std::endl;
      vector<RingElem> RGBhJ = ReducedGBasis(ideal(DomainPhi,GensKer));
      LenRGBh = len(RGBhJ);
      RingHom ImageInQuot = QuotientingHom(DomainPhi);
      return ideal(DomainPhi, ImageInQuot(RGBhJ));
    }  


    // This fn checks whether y is in the image of RichPhi,
    // and if so, it returns the preimage(y)
    RingElem Preimage0RichHom(const RichHom& RichPhi, RingElem y) // wasteful copy of y!
    {
      const RingHom& Phi = RichPhi.Hom;
  
      if (IsQuotientRing(codomain(Phi)))
        y = CanonicalRepr(y); // in S

      const RingElem candidate = NF(RichPhi.HomS_RS(y), RichPhi.idealRS);
      if (!IsZero(candidate)  &&  LPP(candidate) >= RichPhi.MinSIndetRS)
        return zero(domain(Phi)); // --> no preimage

      RingElem preim = RichPhi.HomRS_R(candidate);
      if (IsQuotientRing(domain(Phi)))
        preim = CanonicalHom(RichPhi.R,domain(Phi))(preim);

      return preim;
    }

  } // end of namespace anonymous


  ideal ker(const RingHom& phi) //calculates the kernel of Ringhom Phi
  {
    VerboseLog VERBOSE("ker");
    CheckInputAffAlgebraHom(phi, CoCoA_ERROR_CONTEXT);
    VERBOSE(90) << "starting..." << std::endl;
    return KerRichHom(MakeRichHom(phi));
  }


  ideal ker_H(const RingHom& phi)
  {
    CheckInputAffAlgebraHom(phi, CoCoA_ERROR_CONTEXT);
    return KerRichHom_H(MakeRichHom_H(phi));
  }


  // Is phi injective?
  bool IsInjective(const RingHom& phi)
  {
    CheckInputAffAlgebraHom(phi, CoCoA_ERROR_CONTEXT);
    return IsZero(KerRichHom(MakeRichHom(phi)));  
  }


  // Is phi surjective?
  bool IsSurjective(const RingHom& phi)
  {
    CheckInputAffAlgebraHom(phi, CoCoA_ERROR_CONTEXT);
    const RichHom RichPhi = MakeRichHom(phi);
    const ideal LTI = LT(RichPhi.idealRS);
    // SLEDGEHAMMER below!!!  How to do it better?
    const vector<RingElem>& y = RichPhi.SIndetsRS;
    const int nvars = len(y);
    for (int i=0; i < nvars; ++i)
      if (!IsZero(NF(y[i],LTI)))  return false;
    return true;
  }


  // Is y in Im(phi)?
  bool IsInImage(const RingHom& phi, ConstRefRingElem y)
  {
    CheckInputAffAlgebraHom(phi, CoCoA_ERROR_CONTEXT);
    if (owner(y) != codomain(phi))
      CoCoA_THROW_ERROR2(ERR::BadArg, "arg 2 must be in codomain");
    if (IsZero(y)) return true;
    return !IsZero(Preimage0RichHom(MakeRichHom(phi), y));
  } 

  
  //  phi^(-1)(y) if it exists, error otherwise
  RingElem preimage(const RingHom& phi, ConstRefRingElem y)
  {
    CheckInputAffAlgebraHom(phi, CoCoA_ERROR_CONTEXT);
    if (owner(y) != codomain(phi)) 
      CoCoA_THROW_ERROR2(ERR::BadArg, "arg 2 must be in codomain");

    if (IsZero(y)) return zero(domain(phi));
    const RingElem ans = Preimage0RichHom(MakeRichHom(phi), y);
    if (IsZero(ans)) CoCoA_THROW_ERROR2(ERR::BadArg, "Not in image of RingHom");
    return ans;
  }

  //  phi^(-1)(y) if it exists, zero(domain(phi)) otherwise
  RingElem preimage0(const RingHom& phi, ConstRefRingElem y)
  {
    CheckInputAffAlgebraHom(phi, CoCoA_ERROR_CONTEXT);
    if (owner(y) != codomain(phi)) 
      CoCoA_THROW_ERROR2(ERR::BadArg, "arg 2 must be in codomain");
    return Preimage0RichHom(MakeRichHom(phi), y);
  }


} // end of namespace CoCoA
