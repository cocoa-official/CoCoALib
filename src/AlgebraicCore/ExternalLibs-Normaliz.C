//   Copyright (c)  2010-2021  John Abbott and Anna M. Bigatti
//   Authors:  2010-2015 Anna M. Bigatti, Christof Soeger

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


// Source code for Normaliz integration

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#ifdef CoCoA_WITH_NORMALIZ

#include "CoCoA/ExternalLibs-Normaliz.H"
#include "libnormaliz/libnormaliz.h"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/TmpPPVector.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/symbol.H"
#include "CoCoA/QuasiPoly.H"

#include "CoCoA/VectorOps.H"  // just for debugging

//#include "libnormaliz/libnormaliz.h"
using libnormaliz::Cone;
using libnormaliz::ConeProperties;



// #include <list>
// using std::list;
#include <map>
using std::map;
#include <vector>
using std::vector;
#include <ostream>
using std::ostream;
using std::endl;
#include <sstream>  // for ErrorMessage


namespace CoCoA
{

  namespace Normaliz
  {

    libnormaliz::InputType ToInputType(const std::string& TypeString)
    {
      return libnormaliz::to_type(TypeString);
    }

    class ConeImpl: protected IntrusiveReferenceCount
    {
      friend class SmartPtrIRC<const ConeImpl>; // Morally "friend Cone", so it can alter reference count.

      public:
        ConeImpl(libnormaliz::Type::InputType InputType, const std::vector<std::vector<BigInt> >& v);
        ConeImpl(const std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > >& m);
        ConeImpl(libnormaliz::Type::InputType InputType, const PPVector& v);
        friend void SetVerbosityLevel(cone& c, long n); // 0 off/ >1 on
        friend long VerbosityLevel(cone& c);
        friend std::vector<std::vector<BigInt> > HilbertBasis(const cone& c);
        friend std::vector<std::vector<BigInt> > ModuleGenerators(const cone& c);
        friend std::vector<std::vector<BigInt> > TriangulationGenerators(const cone& c);
        friend std::vector<std::vector<BigInt> > ExtremeRays(const cone& c);
        friend std::vector<std::vector<BigInt> > VerticesOfPolyhedron(const cone& c);
        friend std::vector<std::vector<BigInt> > Deg1Elements(const cone& c);
        friend std::vector<std::vector<BigInt> > GeneratorsOfToricRing(const cone& c);
        friend std::vector<std::vector<BigInt> > SupportHyperplanes(const cone& c);
        friend std::vector<std::vector<BigInt> > Equations(const cone& c);
        friend std::vector<std::vector<BigInt> > Congruences(const cone& c);
        friend std::vector<std::vector<BigInt> > ExcludedFaces(const cone& c);
        //friend std::vector<BigInt> HVector(const cone& c);
        friend HPSeries HilbertSeries(const cone& c);
        friend RingElem HilbertPoly(const cone& c);
        friend QuasiPoly HilbertQuasiPoly(const cone& c);
        friend BigRat multiplicity(const cone& c);
        friend std::vector<BigRat> grading(const cone& c);
        friend bool IsPointed(const cone& c);
        friend bool IsInhomogeneous(const cone& c);
        friend bool IsIntegrallyClosed(const cone& c);
        friend bool IsDeg1HilbertBasis(const cone& c);
        friend long EmbeddingDim(const cone& c);
        friend long rank(const cone& c);
        friend long RecessionRank(const cone& c);
        friend long AffineDim(const cone& c);
        friend long ModuleRank(const cone& c);
        friend std::vector<BigInt> dehomogenization(const cone& c);

      public:
        void myComputation(const libnormaliz::ConeProperties& CPs) const;
        void myComputation(libnormaliz::ConeProperty::Enum CP) const;    
        void myComputation(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const;    
        void myComputation() const;   //default: compute everything possible 
        bool isComputed(libnormaliz::ConeProperty::Enum CP) const;    
        bool isComputed(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const;    

      private:
        mutable libnormaliz::Cone<mpz_class> myConeMPZ; // need "mutable" to cater for Normaliz interpretation of "const" (namely, they take it to mean that the repr is const, and not just the object being represented).
    };


    // implementation of cone

    cone::cone(libnormaliz::Type::InputType InputType, const std::vector<std::vector<BigInt> >& v) : mySmartPtr(new ConeImpl(InputType, v)) {}

    cone::cone(const std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > >& m) : mySmartPtr(new ConeImpl(m)) {}

    cone::cone(libnormaliz::Type::InputType InputType, ConstMatrixView M) : mySmartPtr(new ConeImpl(InputType, MatrixToVecVecBigInt(M))) {}

    cone::cone(const ConeImpl* ptr): mySmartPtr(ptr) {}
    cone::cone(const cone& c) : mySmartPtr(c.mySmartPtr) {}
    cone::~cone() {}


    // printing
    ostream& operator<< (ostream& out, const cone& C)
    {
      if (!out) return out;  // short-cut for bad ostreams
      using namespace libnormaliz;
      out << "cone(EmbeddingDim = " << EmbeddingDim(C);
      if (C.isComputed(ConeProperty::ExtremeRays)) {
        out << ", rank = " << rank(C);
        out << ", pointed = " << IsPointed(C);
      }
      if (C.isComputed(ConeProperty::Generators))
        out << ", NumGenerators = " << TriangulationGenerators(C).size();
      if (C.isComputed(ConeProperty::SupportHyperplanes))
        out << ", NumSupportHyperplanes = " << SupportHyperplanes(C).size();
      out << ")";
      return out;
    }
    
    // letVerbose(Default) sets the verbosity, and returns the old value
    void SetVerbosityLevel(long n) {
      libnormaliz::setVerboseDefault(n!=0);
    }

    long VerbosityLevel() {
      const bool b = libnormaliz::setVerboseDefault(0);
      libnormaliz::setVerboseDefault(b);
      return b;
    }

    void SetVerbosityLevel(cone& c, long n) {
      c->myConeMPZ.setVerbose(n!=0);
    }

    long VerbosityLevel(cone& c) {
      const bool b = c->myConeMPZ.setVerbose(0);
      c->myConeMPZ.setVerbose(b);
      return b;
    }

    // computations
    void cone::myComputation(const libnormaliz::ConeProperties& CPs) const { mySmartPtr->myComputation(CPs); }
    void cone::myComputation(libnormaliz::ConeProperty::Enum CP) const { mySmartPtr->myComputation(CP); }
    void cone::myComputation(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const { mySmartPtr->myComputation(CP1, CP2); }
    void cone::myComputation() const { mySmartPtr->myComputation(); }
    bool cone::isComputed(libnormaliz::ConeProperty::Enum CP) const { return mySmartPtr->isComputed(CP); }
    bool cone::isComputed(libnormaliz::ConeProperty::Enum CP1, libnormaliz::ConeProperty::Enum CP2) const { return mySmartPtr->isComputed(CP1, CP2); }
    const ConeImpl* cone::operator->() const { return mySmartPtr.operator->(); }

    namespace  // conversion functions, implementation at the end of the file
    {
//      std::vector<BigInt> LongVToBigIntV(const std::vector<long>& VIn);
      std::vector<BigInt> MPZ_classVToBigIntV(const std::vector<mpz_class>& VIn);
//      std::vector<BigRat> LongVToBigRatV(const std::vector<long>& VIn, const BigInt& denom);
      std::vector<BigRat> MPZ_classVToBigRatV(const std::vector<mpz_class>& VIn, const BigInt& denom);
//      std::vector<long> BigIntVToLongV(const std::vector<BigInt>& VIn);
      std::vector<mpz_class> BigIntVToMPZ_classV(const std::vector<BigInt>& VIn);
//      void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<long> >& l);
      void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<mpz_class> >& l);
//      std::vector<std::vector<long> >  ReturnVecVecLong(const std::vector<std::vector<BigInt> >& v);
//      map< libnormaliz::InputType, vector<vector<long> > >  ReturnMapVecVecLong(const map< libnormaliz::InputType, vector<vector<BigInt> > >& m);
      std::vector<std::vector<mpz_class> >  ReturnVecVecMPZ_class(const std::vector<std::vector<BigInt> >& v);
      map< libnormaliz::InputType, vector<vector<mpz_class> > >  ReturnMapVecVecMPZ_class(const map< libnormaliz::InputType, vector<vector<BigInt> > >& m);
      std::vector<std::vector<BigInt> > PPVectorToVecVecBigInt(const PPVector& ppv);
      void VecVecBigIntToPPVector(PPVector& ppv, const std::vector<std::vector<BigInt> >& M);
    } // end of anonymous namespace


    // implementation of ConeImpl
    ConeImpl::ConeImpl(libnormaliz::Type::InputType InputType, const std::vector<std::vector<BigInt> >& v):
          myConeMPZ(InputType, ReturnVecVecMPZ_class(v))
    {
    }

    ConeImpl::ConeImpl(const std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > >& m):
      myConeMPZ(ReturnMapVecVecMPZ_class(m))
    {
    }

    void ConeImpl::myComputation(const libnormaliz::ConeProperties& CPs) const
    {
        libnormaliz::ConeProperties missing = myConeMPZ.compute(CPs);
        if (missing.any())
        {
          std::ostringstream os;
          os << missing;
          CoCoA_THROW_ERROR2(ERR::ShouldNeverGetHere, "Normaliz cannot compute "+ os.str());
        }
    }

    void ConeImpl::myComputation(libnormaliz::ConeProperty::Enum CP) const
    {
      myComputation(ConeProperties(CP));
    }

    void ConeImpl::myComputation(libnormaliz::ConeProperty::Enum CP1,
                                 libnormaliz::ConeProperty::Enum CP2) const
    {
      myComputation(ConeProperties(CP1, CP2));
    }

    void ConeImpl::myComputation() const
    {
      myComputation(ConeProperties(libnormaliz::ConeProperty::DefaultMode));
    }

    bool ConeImpl::isComputed(libnormaliz::ConeProperty::Enum CP) const
    {
      return myConeMPZ.isComputed(CP);
    }

    bool ConeImpl::isComputed(libnormaliz::ConeProperty::Enum CP1,
                              libnormaliz::ConeProperty::Enum CP2) const
    {
        return isComputed(CP1) && isComputed(CP2);
    }

    // friend functions which are the interface to the user
    std::vector<std::vector<BigInt> > HilbertBasis(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::HilbertBasis);
      ConvertFromNormaliz(v, c->myConeMPZ.getHilbertBasis());
      return v;
    }

    std::vector<std::vector<BigInt> > ModuleGenerators(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::ModuleGenerators);
      ConvertFromNormaliz(v, c->myConeMPZ.getModuleGenerators());
      return v;
    }

    std::vector<std::vector<BigInt> > Deg1Elements(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::Deg1Elements);
      ConvertFromNormaliz(v, c->myConeMPZ.getDeg1Elements());
      return v;
    }

    std::vector<std::vector<BigInt> > TriangulationGenerators(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::Generators);
      //ConvertFromNormaliz(v, c->myConeMPZ.getTriangulationGenerators()); // v 3.8.9
      ConvertFromNormaliz(v, (c->myConeMPZ.getTriangulation()).second.get_elements()); // v 3.8.10      
      return v;
    }

    std::vector<std::vector<BigInt> > ExtremeRays(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      ConvertFromNormaliz(v, c->myConeMPZ.getExtremeRays());
      return v;
    }

    std::vector<std::vector<BigInt> > VerticesOfPolyhedron(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::VerticesOfPolyhedron);
      ConvertFromNormaliz(v, c->myConeMPZ.getVerticesOfPolyhedron());
      return v;
    }

    std::vector<std::vector<BigInt> > GeneratorsOfToricRing(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::OriginalMonoidGenerators);
      ConvertFromNormaliz(v, c->myConeMPZ.getOriginalMonoidGenerators());
      return v;
    }


    // The following constraints depend all on ConeProperty::SupportHyperplanes
    std::vector<std::vector<BigInt> > SupportHyperplanes(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::SupportHyperplanes);
      ConvertFromNormaliz(v, c->myConeMPZ.getSupportHyperplanes());
      return v;
    }

    std::vector<std::vector<BigInt> > Equations(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::Sublattice);
      ConvertFromNormaliz(v, c->myConeMPZ.getSublattice().getEquations());
      return v;
    }

    std::vector<std::vector<BigInt> > Congruences(const cone& c)
    { 
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::Sublattice);
      ConvertFromNormaliz(v, c->myConeMPZ.getSublattice().getCongruences());
      return v;
    }

    // excluded faces are only available if they were explicitly given
    std::vector<std::vector<BigInt> > ExcludedFaces(const cone& c)
    {
      std::vector<std::vector<BigInt> > v;
      c->myComputation(libnormaliz::ConeProperty::ExcludedFaces);
      ConvertFromNormaliz(v, c->myConeMPZ.getExcludedFaces());
      return v;
    }

    HPSeries HilbertSeries(const cone& c)
    { 
      c->myComputation(libnormaliz::ConeProperty::HilbertSeries);
      const libnormaliz::HilbertSeries& HS = c->myConeMPZ.getHilbertSeries();

      return HPSeries(MPZ_classVToBigIntV(HS.getNum()),
                      libnormaliz::to_vector(HS.getDenom()), HS.getShift());
    }

    RingElem HilbertPoly(const cone& c)
    { 
      c->myComputation(libnormaliz::ConeProperty::HilbertSeries);
      const libnormaliz::HilbertSeries& HS = c->myConeMPZ.getHilbertSeries();
      if (HS.getPeriod() != 1)
      {
        CoCoA_THROW_ERROR2(ERR::BadArg, "Hilbert function is a quasi-polynomial of period > 1.  This function works for regular polynomials only.");
      }
      const vector<BigRat> coeffs = MPZ_classVToBigRatV( HS.getHilbertQuasiPolynomial()[0],
                                                         BigIntFromMPZ(HS.getHilbertQuasiPolynomialDenom().get_mpz_t()));
      PolyRing QQt = RingQQt(1);
      const RingElem t = indet(QQt,0);
      RingElem tpower = one(QQt);
      RingElem hpoly(QQt);
      for (long i=0; i<len(coeffs); ++i)
      {
          hpoly += coeffs[i] * tpower;   // NB tpower = t^i
          tpower *= t;
      }

      return hpoly;
    }

    QuasiPoly HilbertQuasiPoly(const cone& c)
    { 
      c->myComputation(libnormaliz::ConeProperty::HilbertSeries);
      const libnormaliz::HilbertSeries& HS = c->myConeMPZ.getHilbertSeries();
      const long period = HS.getPeriod();
      if (period < 1)
        CoCoA_THROW_ERROR2(ERR::BadArg, "Hilbert function not computed.");

      const PolyRing QQt = RingQQt(1);
      const RingElem t = indet(QQt,0);
      RingElem tpower;
      vector<BigRat> coeffs;
      vector<RingElem> qp = vector<RingElem>(period, RingElem(QQt));
      for (long j=0; j<period; ++j)
      {
        coeffs = MPZ_classVToBigRatV( HS.getHilbertQuasiPolynomial()[j],
                                      BigIntFromMPZ(HS.getHilbertQuasiPolynomialDenom().get_mpz_t()));
        tpower = one(QQt);
        for (long i=0; i<len(coeffs); ++i)
        {
            qp[j]  += coeffs[i] * tpower;   // NB tpower = t^i
            tpower *= t;
        }
      }

      return QuasiPoly(qp);
    }

    BigRat multiplicity(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::Multiplicity);
      return BigRatFromMPQ(c->myConeMPZ.getMultiplicity().get_mpq_t());
    }



    bool IsPointed(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      return c->myConeMPZ.isPointed();
    }

    bool IsInhomogeneous(const cone& c)
    {
      // is always known
      return c->myConeMPZ.isInhomogeneous();
    }

    bool IsIntegrallyClosed(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::HilbertBasis);
      return c->myConeMPZ.isIntegrallyClosed();
    }

    bool IsDeg1HilbertBasis(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::HilbertBasis, libnormaliz::ConeProperty::Grading);
      return c->myConeMPZ.isDeg1HilbertBasis();
    }

   

    // dimension and rank invariants
    long EmbeddingDim(const cone& c)
    {
      // is always known
      return c->myConeMPZ.getEmbeddingDim();
    }

    long rank(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      return c->myConeMPZ.getRank();
    }

    std::vector<BigRat> grading(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::Grading);
      return MPZ_classVToBigRatV(c->myConeMPZ.getGrading(), BigIntFromMPZ(c->myConeMPZ.getGradingDenom().get_mpz_t()));
    }

    // only for inhomogeneous case:
    long RecessionRank(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::RecessionRank);
      return c->myConeMPZ.getRecessionRank();
    }

    long AffineDim(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ExtremeRays);
      return c->myConeMPZ.getRank();
    }

    long ModuleRank(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::ModuleRank);
      return c->myConeMPZ.getModuleRank();
    }

    std::vector<BigInt> dehomogenization(const cone& c)
    {
      c->myComputation(libnormaliz::ConeProperty::Dehomogenization);
      return MPZ_classVToBigIntV(c->myConeMPZ.getDehomogenization());
    }

    // not so important in this library
    // size_t getTriangulationSize() const;
    // Integer getTriangulationDetSum() const;


/**************  applications to monomials  **************/

    namespace // anonymous for file local fns
    {
      PPVector HilbertBasis_PPVector(libnormaliz::InputType t, const PPVector& ppv)
      {
        cone c (t, PPVectorToVecVecBigInt(ppv));
        std::vector<std::vector<BigInt> > res = HilbertBasis(c);
        //create an emptyPPVector with the same PPMonoid
        PPVector ppv_res(PPM(ppv), DMR(ppv));
        VecVecBigIntToPPVector(ppv_res,res);
        return ppv_res;
      }

      
      //only for modes where normaliz adds a component, like polytope or rees_algebra
      PPVector HilbertBasis_PPVector(libnormaliz::Type::InputType t, const PPVector& ppv, long sym_pos)
      { 
        const long dim = NumIndets(PPM(ppv));
        if (sym_pos < 0 || sym_pos >= dim)
          CoCoA_THROW_ERROR2(ERR::OutOfRange, "Symbol position");

        vector< vector<BigInt> > vv_input = PPVectorToVecVecBigInt(ppv);
        if (sym_pos != dim-1)
          CoCoA_THROW_ERROR2(ERR::NYI, "Symbol needs to be the last symbol of the PPMonoid (other implementations are missing)");
        for (long i=0; i < len(vv_input); ++i) 
        {
          if (vv_input[i][sym_pos] != 0)
            CoCoA_THROW_ERROR2(ERR::BadArg, "Symbol must not be used in the power products");
          vv_input[i].erase(vv_input[i].begin()+sym_pos);
        }
        if (t != libnormaliz::Type::polytope && t != libnormaliz::Type::rees_algebra)
          CoCoA_THROW_ERROR2(ERR::BadArg, "Invalid InputType for this method");
        
        cone c (t, vv_input);
        vector<vector<BigInt> > res = HilbertBasis(c);
        if (sym_pos != dim-1)
          for (long i=0; i < dim; ++i)  swap(res[i][dim],res[i][sym_pos]);
        //create an emptyPPVector with the same PPMonoid
        PPVector ppv_res(PPM(ppv), DMR(ppv));
        VecVecBigIntToPPVector(ppv_res,res);

        return ppv_res;
      }

      
      PPVector HilbertBasis_PPVector(libnormaliz::InputType t, const std::vector<std::vector<BigInt> >& Mat, PPMonoid ppm)
      {
        cone c (t, Mat);
        std::vector<std::vector<BigInt> > res = HilbertBasis(c);
        //create an empty PPVector with the same PPMonoid
        PPVector ppv_res(ppm,NewDivMaskNull());
        VecVecBigIntToPPVector(ppv_res,res);
        return ppv_res;
      }

      
      PPVector HilbertBasis_PPVector(const cone& C, PPMonoid ppm)
      {
        std::vector<std::vector<BigInt> > res = HilbertBasis(C);
        //create an empty PPVector with the same PPMonoid
        PPVector ppv_res(ppm,NewDivMaskNull());
        VecVecBigIntToPPVector(ppv_res,res);
        return ppv_res;
      }
      
    } // end anonymous namespace

    
    PPVector NormalToricRing(const PPVector& ppv)
    { return HilbertBasis_PPVector(libnormaliz::Type::normalization, ppv); }

    
    PPVector IntClosureToricRing(const PPVector& ppv)
    { return HilbertBasis_PPVector(libnormaliz::Type::integral_closure, ppv); }

    
    PPVector IntClosureMonIdeal(const PPVector& ppv)
    {
      cone c (libnormaliz::Type::rees_algebra, PPVectorToVecVecBigInt(ppv));
      vector< vector<BigInt> > res1;
      for (const auto& v: HilbertBasis(c))
        if (v.back() == 1 )
        {
          res1.push_back(v);
          res1.back().pop_back();
        }

      //create an empty PPVector with the same PPMonoid
      PPVector ppv_res1(PPM(ppv), DMR(ppv));
      VecVecBigIntToPPVector(ppv_res1,res1);
      return ppv_res1;
    }

    
    /* If you want also the whole normalization of the rees algebra then the
     * ring must have an unused symbol that we can use for homogenization.
     * You have to specify its position.
     */
    PPVector IntClosureMonIdeal(const PPVector& ppv, long sym_pos)
    { return HilbertBasis_PPVector(libnormaliz::Type::rees_algebra, ppv, sym_pos); }

    
    PPVector EhrhartRing(const PPVector& ppv, long sym_pos)
    { return HilbertBasis_PPVector(libnormaliz::Type::polytope, ppv, sym_pos); }

/**************  torus invariants and valuation rings  **************/
    PPVector TorusInvariants(const vector< vector<BigInt> >& T, const PPMonoid& ppm)
    {
      if (T.empty())
        CoCoA_THROW_ERROR2(ERR::BadMatrixSize, "matrix must be non-empty");
      if (NumIndets(ppm) != len(T[0]))
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "NumCols and NumIndets must be equal");
      return HilbertBasis_PPVector(libnormaliz::Type::equations, T, ppm);
    }

    
    PPVector FiniteDiagInvariants(const vector< vector<BigInt> >& Cong, const PPMonoid& ppm)
    {
      if (Cong.empty())
        CoCoA_THROW_ERROR2(ERR::BadMatrixSize, "matrix must be non-empty");
      if (NumIndets(ppm) != len(Cong[0])-1)
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "NumCols-1 and NumIndets must be equal");
      return HilbertBasis_PPVector(libnormaliz::Type::congruences, Cong, ppm);
    }

    
    PPVector DiagInvariants(const vector< vector<BigInt> >& T, const vector< vector<BigInt> >& Cong, const PPMonoid& ppm)
    {
      if (T.empty() && Cong.empty())
        CoCoA_THROW_ERROR2(ERR::BadMatrixSize, "at leasst one matrix must be non-empty");
      if (T.empty())
        return FiniteDiagInvariants(Cong, ppm);
      if (Cong.empty())
        return TorusInvariants(T, ppm);

      if (NumIndets(ppm) != len(T[0]))
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "NumCols(T) and NumIndets must be equal");
      if (NumIndets(ppm) != len(Cong[0])-1)
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "NumCols(Cong)-1 and NumIndets must be equal");
      std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > cone_input;
      cone_input[libnormaliz::Type::equations] = T;
      cone_input[libnormaliz::Type::congruences] = Cong;
      cone C(cone_input);
      return HilbertBasis_PPVector(C, ppm);
    }


    PPVector IntersectionValRings (const vector< vector<BigInt> >& V, const PPMonoid& ppm)
    {
      const long dim = len(V[0]);
      std::map< libnormaliz::InputType, std::vector<std::vector<BigInt> > > cone_input;
      cone_input[libnormaliz::Type::inequalities] = V;
      const vector< vector<BigInt> > positive_signs = vector< vector<BigInt> >(1,vector<BigInt>(dim,BigInt(1)));
      cone_input[libnormaliz::Type::signs] = positive_signs;
      cone C(cone_input);
      return HilbertBasis_PPVector(C, ppm);
    }


/********************************************************
 ***               conversion functions               ***
 ********************************************************/

    namespace  // conversion functions
    {

// Compiler says this is never used
      // std::vector<BigInt> LongVToBigIntV(const std::vector<long>& VIn)
      // {
      //   std::vector<BigInt> v;
      //   v.reserve(VIn.size());
      //   for (vector<long>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
      //     v.push_back(BigInt(*it));
      //   return v;
      // }

      std::vector<BigInt> MPZ_classVToBigIntV(const std::vector<mpz_class>& VIn)
      {
        std::vector<BigInt> v; v.reserve(VIn.size());
        for (const mpz_class& N: VIn)
          v.push_back(BigIntFromMPZ(N.get_mpz_t()));
        return v;
      }

// Compiler says this is never used
      // std::vector<BigRat> LongVToBigRatV(const std::vector<long>& VIn, const BigInt& denom)
      // {
      //   std::vector<BigRat> v;
      //   v.reserve(VIn.size());
      //   for (vector<long>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
      //     v.push_back(BigRat(BigInt(*it),denom));
      //   return v;
      // }

      std::vector<BigRat> MPZ_classVToBigRatV(const std::vector<mpz_class>& VIn, const BigInt& denom)
      {
        std::vector<BigRat> v;
        v.reserve(VIn.size());
        for (const mpz_class& N: VIn)
          v.push_back(BigRat(BigIntFromMPZ(N.get_mpz_t()), denom));
        return v;
      }


// Compiler says this is never used
      // std::vector<long> BigIntVToLongV(const std::vector<BigInt>& VIn)
      // {
      //   std::vector<long> v;
      //   long n;
      //   v.reserve(VIn.size());
      //   for (vector<BigInt>::const_iterator it=VIn.begin(); it!=VIn.end(); ++it)
      //     if (IsConvertible(n, *it))  v.push_back(n);
      //     else CoCoA_THROW_ERROR(ERR::BadConvert, "normaliz: BigIntVToLongV");
      //   return v;
      // }


      std::vector<mpz_class> BigIntVToMPZ_classV(const std::vector<BigInt>& VIn)
      {
        std::vector<mpz_class> v;
        v.reserve(VIn.size());
        for (const BigInt& N: VIn)
          v.push_back(mpz_class(mpzref(N)));
        return v;
      }


// Compiler says this is never used
      // void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<long> >& l)
      // {
      //   v.clear();      // not exception safe
      //   v.reserve(l.size());
      //   //      transform(l.begin(), l.end(), v.begin(), LongVToBigIntV);
      //   for (vector<vector<long> >::const_iterator it=l.begin(); it!=l.end(); ++it)
      //     v.push_back(LongVToBigIntV(*it));
      // }    


      void ConvertFromNormaliz(std::vector<std::vector<BigInt> >& v, const std::vector<std::vector<mpz_class> >& L)
      {
        v.clear();      // not exception safe
        v.reserve(L.size());
        //      transform(L.begin(), L.end(), v.begin(), LongVToBigIntV);
        for (const auto& vec: L)
          v.push_back(MPZ_classVToBigIntV(vec));
      }    


// Compiler says this is never used
      // std::vector<std::vector<long> >  ReturnVecVecLong(const std::vector<std::vector<BigInt> >& v)
      // {
      //   std::vector<std::vector<long> > l;
      //   //      transform(v.begin(), v.end(), l.begin(), BigIntVToLongV);
      //   try
      //   {
      //       for (vector<vector<BigInt> >::const_iterator it=v.begin(); it!=v.end(); ++it)
      //         l.push_back(BigIntVToLongV(*it));
      //   } catch (...) {l.clear();}
      //   return l;
      // }


// Compiler says this is never used
      // map< libnormaliz::InputType, vector<vector<long> > >  ReturnMapVecVecLong(const map< libnormaliz::InputType, vector<vector<BigInt> > >& m)
      // {
      //   map< libnormaliz::InputType, vector<vector<long> > > ret;
      //   for (map<libnormaliz::InputType, vector<vector<BigInt> > >::const_iterator it=m.begin(); it!=m.end(); ++it)
      //     ret.insert(make_pair((*it).first, ReturnVecVecLong((*it).second)));
      //   return ret;
      // }
      
      
      std::vector<std::vector<mpz_class> >  ReturnVecVecMPZ_class(const std::vector<std::vector<BigInt> >& vlist)
      {
        std::vector<std::vector<mpz_class> > L; L.reserve(vlist.size());
        //      transform(v.begin(), v.end(), L.begin(), BigIntVToLongV);
        for (const auto& vec: vlist)
          L.push_back(BigIntVToMPZ_classV(vec));
        return L;
      }


      map< libnormaliz::InputType, vector<vector<mpz_class> > >  ReturnMapVecVecMPZ_class(const map< libnormaliz::InputType, vector<vector<BigInt> > >& mlist)
      {
        map< libnormaliz::InputType, vector<vector<mpz_class> > > ret;
//        for (map<libnormaliz::InputType, vector<vector<BigInt> > >::const_iterator it=m.begin(); it!=m.end(); ++it)
        for (const auto& m: mlist)
          ret.insert(make_pair(m.first, ReturnVecVecMPZ_class(m.second)));
        return ret;
      }


      std::vector<std::vector<BigInt> > PPVectorToVecVecBigInt(const PPVector& ppv)
      {
        vector<BigInt> tmp;
        const long n =  len(ppv);
        vector<vector<BigInt> > v(n);
        for (long i=0; i<n; ++i)
        {
          BigExponents(tmp,PP(ppv[i]));
          v[i]=tmp;
        }
        return v;
      }

      void VecVecBigIntToPPVector(PPVector& ppv, const std::vector<std::vector<BigInt> >& M)
      {
        PPMonoid ppm = PPM(ppv);
        const long n =  len(M);
        for (long i=0; i < n; ++i)
        {
          ppv.myPushBack(PPMonoidElem(ppm, M[i]));
        }
      }

    } // end of anonymous namespace

    std::vector<std::vector<BigInt> > MatrixToVecVecBigInt(ConstMatrixView M)
    {
      const ErrorInfo ErrMesg("Matrix entries must be integer", "MatrixToVecVecBigInt (cone ctor)");

      vector<vector<BigInt> > v;
      for (long i=0; i<NumRows(M); ++i)
      {
        v.push_back(vector<BigInt>());
        for (long j=0; j<NumCols(M); ++j)
          v[i].push_back(ConvertTo<BigInt>(M(i,j), ErrMesg));
      }
      return v;
    }

    PPVector MonomialsToPPV(const std::vector<RingElem>& v)
    {
      if (!AreMonomials(v))  CoCoA_THROW_ERROR2(ERR::BadArg, "vector must contain monomials");
      if (v.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
      //convert it to a PPVector
      PPVector ppv(PPM(owner(v[0])), NewDivMaskNull());
      convert(ppv,v);
      return ppv;
    }


  } // namespace Normaliz
} // namespace CoCoA

#endif
