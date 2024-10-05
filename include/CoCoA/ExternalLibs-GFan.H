#ifndef CoCoA_ExternalLibs_GFan_H
#define CoCoA_ExternalLibs_GFan_H

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

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#ifdef CoCoA_WITH_GFAN

#include "CoCoA/SmartPtrIRC.H"

#include <iosfwd>

namespace CoCoA
{
  class BigInt;          // forward declaration -- defined in BigInt.H 
  class ConstMatrixView; // forward declaration -- defined in matrix.H
  class matrix;          // forward declaration -- defined in matrix.H

  namespace GFan
  {

    class ConeImpl; // fwd decl for SmartPtrIRC
    class cone
    {
    public:

      cone(ConstMatrixView IneqMat, ConstMatrixView EqMat);

      explicit cone(const ConeImpl* ptr); //: mySmartPtr(ptr) {}
      // Default copy ctor and dtor works fine, but we implement them to ensure they are not made when including this .H file which will lead to errors because the declaration of ConeImpl is unknown
      ~cone();
      cone(const cone&);
      // Assignment disabled because SmartPtrIRC has no assignment.

      //    void myComputation(const gfan::ConeProperties& CPs) const;

      const ConeImpl* operator->() const; // { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
      //    const ConeImpl* myRawPtr() const { return mySmartPtr.myRawPtr(); } ///< Used by "downcasting" functions
    private: // data members
      SmartPtrIRC<const ConeImpl> mySmartPtr;
    };

    std::ostream& operator<< (std::ostream& out, const cone& C);
    matrix equations(const cone& c);
    matrix inequalities(const cone& c);
    matrix RelativeInteriorPoint(const cone& C);
    matrix GeneratorsOfSpan(const cone& C);
    matrix GeneratorsOfLinealitySpace(const cone& C);
    matrix GetFacets(const cone& C);
    matrix GetImpliedEquations(const cone& C);
    matrix GetUniquePoint(const cone& C);
    long GetAmbientDimension(const cone& C);
    long GetDimension(const cone& C);
    long GetCodimension(const cone& C);
    long GetDimensionOfLinealitySpace(const cone& C);
    bool ContainsPositiveVector(const cone& C);

  } // end of namespace gfanlib

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ExternalLibs-GFan.H,v 1.9 2022/02/25 10:36:37 abbott Exp $
// $Log: ExternalLibs-GFan.H,v $
// Revision 1.9  2022/02/25 10:36:37  abbott
// Summary: Removed unnecessary includes from header filed
//
// Revision 1.8  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.7  2020/01/24 21:37:37  abbott
// Summary: Changed name of include guard
//
// Revision 1.6  2017/03/13 12:17:08  abbott
// Summary: Include PREPROCESSOR_DEFNS.H (in header files which need the info)
//
// Revision 1.5  2015/09/10 16:03:24  abbott
// Summary: Removed ctrl-M before each newline
//
// Revision 1.4  2015/09/04 09:47:12  bigatti
// -- lost new functions with Gfan! (by Anders Nedergaard Jensen)
//
// Revision 1.3  2015/09/02 16:44:25  bigatti
// -- added RelativeInteriorPoint
// -- added conversion for ZVector
//
// Revision 1.2  2015/09/02 16:24:52  bigatti
// -- first real functions
//
// Revision 1.1  2015/09/02 09:29:08  bigatti
// -- first import (in Aarhus)
//

#endif // CoCoA_WITH_GFAN
#endif
