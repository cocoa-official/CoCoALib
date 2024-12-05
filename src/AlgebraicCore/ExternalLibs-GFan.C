//   Copyright (c)  2015 Anna M. Bigatti, Anders Nedergaard Jensen

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


// Source code for GFan integration

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#ifdef CoCoA_WITH_GFAN

#include "CoCoA/ExternalLibs-GFan.H"
#include "gfanlib/gfanlib.h"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/TmpPPVector.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/matrix.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/symbol.H"
#include "CoCoA/RingZZ.H"

#include "CoCoA/VectorOps.H"  // just for debugging

#include <vector>
using std::vector;
#include <ostream>
using std::ostream;
using std::endl;
#include <sstream>  // for ErrorMessage


namespace CoCoA
{

  namespace GFan  // "gfan" is its own namespace
  {

    namespace  // conversion functions, implementation at the end of the file
    {

      gfan::ZMatrix ConvertMat(ConstMatrixView CM)
      {
        const long nrows = NumRows(CM);
        const long ncols = NumCols(CM);
        gfan::ZMatrix ans(nrows, ncols);
        for (long i=0; i<nrows; ++i)
          for (long j=0; j<ncols; ++j)
            ans[i][j] = gfan::Integer(mpzref(ConvertTo<BigInt>(CM(i,j))));
        return ans;
      }
      
      matrix ConvertMat(const gfan::ZMatrix& M)
      {
        const long nrows = M.getHeight();
        const long ncols = M.getWidth();
        matrix CM(NewDenseMat(RingZZ(), nrows, ncols));
        BigInt tmp;
        for (long i=0; i<nrows; ++i)
          for (long j=0; j<ncols; ++j)
          {
            M[i][j].setGmp(mpzref(tmp));
            SetEntry(CM,i,j, tmp);
          }
        return CM;
      }
      

      matrix ConvertMat(const gfan::ZVector& V)
      {
        const long nrows = len(V);
        matrix CM(NewDenseMat(RingZZ(), nrows, 1));
        BigInt tmp;
        for (long i=0; i<nrows; ++i)
        {
          V[i].setGmp(mpzref(tmp));
          SetEntry(CM,i,0, tmp);
        }
        return CM;
      }
      

    } // end of anonymous namespace

    class ConeImpl: protected IntrusiveReferenceCount
    {
      friend class SmartPtrIRC<const ConeImpl>; // Morally "friend Cone", so it can alter reference count.
      
    public:
      ConeImpl(ConstMatrixView IneqMat, ConstMatrixView EqMat);
      
    public:
      friend  ostream& operator<< (ostream& out, const cone& C);
      friend  matrix equations(const cone& C);
      friend  matrix inequalities(const cone& C);
      friend  matrix RelativeInteriorPoint(const cone& C);
      friend  matrix GeneratorsOfSpan(const cone& C);
      friend  matrix GeneratorsOfLinealitySpace(const cone& C);
      friend  matrix GetFacets(const cone& C);
      friend  matrix GetImpliedEquations(const cone& C);
      friend  matrix GetUniquePoint(const cone& C);
      friend  long GetAmbientDimension(const cone& C);
      friend  long GetDimension(const cone& C);
      friend  long GetCodimension(const cone& C);
      friend  long GetDimensionOfLinealitySpace(const cone& C);
      friend  cone GetDimensionOfLinealitySpace(const cone& C, const cone & D);
      friend  bool ContainsPositiveVector(const cone& C);
      

    private: // data members
      gfan::ZCone myZCone;
    };
    


    ConeImpl::ConeImpl(ConstMatrixView IneqMat, ConstMatrixView EqMat):
      myZCone(ConvertMat(IneqMat), ConvertMat(EqMat)) {}
    
    const ConeImpl* cone::operator->() const { return mySmartPtr.operator->(); }


    // implementation of cone

    cone::cone(ConstMatrixView IneqMat, ConstMatrixView EqMat):
      mySmartPtr(new ConeImpl(IneqMat, EqMat)) {}

    cone::~cone() {}


    // printing
    ostream& operator<< (ostream& out, const cone& C)
    {
      if (!out) return out;  // short-cut for bad ostreams
      using namespace gfan;
      out << C->myZCone.toString();
      return out;
    }
    
    matrix equations(const cone& C)
    { return ConvertMat(C->myZCone.getEquations()); }
    
    matrix inequalities(const cone& C)
    { return ConvertMat(C->myZCone.getInequalities()); }
  
    matrix RelativeInteriorPoint(const cone& C)
    { return ConvertMat(C->myZCone.getRelativeInteriorPoint()); }
  
    matrix GeneratorsOfSpan(const cone& C)
    { return ConvertMat(C->myZCone.generatorsOfSpan()); }

    matrix GeneratorsOfLinealitySpace(const cone& C)
    { return ConvertMat(C->myZCone.generatorsOfLinealitySpace()); }

    matrix GetFacets(const cone& C)
    { return ConvertMat(C->myZCone.getFacets()); }

    matrix GetImpliedEquations(const cone& C)
    { return ConvertMat(C->myZCone.getImpliedEquations()); }

    matrix GetUniquePoint(const cone& C)
    { return ConvertMat(C->myZCone.getUniquePoint()); }

    long GetAmbientDimension(const cone& C)
    { return (C->myZCone.ambientDimension()); }

    long GetDimension(const cone& C)
    { return (C->myZCone.dimension()); }

    long GetCodimension(const cone& C)
    { return (C->myZCone.codimension()); }

    long GetDimensionOfLinealitySpace(const cone& C)
    { return (C->myZCone.dimensionOfLinealitySpace()); }

    cone GetDimensionOfLinealitySpace(const cone& C, const cone & D)
    { 
      gfan::ZCone T=intersection(C->myZCone,D->myZCone);
      return cone(ConvertMat(T.getInequalities()),ConvertMat(T.getEquations()));
    }

    bool ContainsPositiveVector(const cone& C)
    { return C->myZCone.containsPositiveVector(); }

  } // namespace gfanlib
} // namespace CoCoA

#endif
