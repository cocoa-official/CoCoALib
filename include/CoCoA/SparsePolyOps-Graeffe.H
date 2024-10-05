#ifndef CoCoA_NumTheory_Graeffe_H
#define CoCoA_NumTheory_Graeffe_H

//   Copyright (c)  2022  John Abbott,  Anna M. Bigatti
//   Original author: Nico Mexis (2022)

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


#include "CoCoA/ring.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyRing.H"

#include <ostream>
#include <vector>

namespace CoCoA
{

    /**
     * Provides a way to sequentially calculate Graeffe iterations of polynomials
     */
    class GraeffeSeq
    {
    public:
        GraeffeSeq() = delete;                    // We don't need a default constructor
        GraeffeSeq(ConstRefRingElem, const long); // Construct a default GraeffeSeq object
        RingElem operator*() const;               // Retrieve the current Graeffe iteration
        GraeffeSeq &operator++();                 // Calculate Newton coeffcients for next iteration
        GraeffeSeq operator++(int);               // Calculate Newton coeffcients for next iteration

    private: // methods
        /**
         * Checks the given polynomial for validity in the GraeffeSeq class
         */
        static const RingElem &CheckElementValidity(const RingElem &f);
        /**
         * Retrieves the coefficients of the given polynomial mapped to the given polynomial ring
         */
        static const std::vector<RingElem> GetCoefficients(const RingElem &f, const ring &Q, const RingElem &xQ);

    private:                                  // data members
        const SparsePolyRing Q;               // ring of temporary poly
        const RingElem xQ;                    // indeterminate of temporary ring
        const long myDeg;                     // degree of original poly
        long myCurrPow = 0;                   // degree to which newton values have been computed already
        mutable long myResultPow = 0;         // degree to which the resulting poly has been calculated
        const std::vector<RingElem> myCoeffs; // coefficients of original poly
        mutable std::vector<RingElem> myA;    // cached coefficients of resulting poly (A_i values from dissertation)
        mutable std::vector<RingElem> myS;    // values extended after each iteration (S_i values from dissertation)
        mutable RingElem myResultCache;       // cached resulting poly
        const RingHom phi;                    // map for mapping poly back to original ring
    };

    /**
     * Provides a way to format-print the given GraeffeSeq object
     */
    std::ostream &operator<<(std::ostream &out, const GraeffeSeq &GSeq);

    /**
     * Returns the n-th order Graeffe transform of f, i.e.
     * a univariate polynomial whose roots are the n-th
     * powers of the roots of f.
     */
    RingElem GraeffeN(RingElem f, long n);

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-Graeffe.H,v 1.2 2022/11/24 21:11:29 abbott Exp $
// $Log: SparsePolyOps-Graeffe.H,v $
// Revision 1.2  2022/11/24 21:11:29  abbott
// Summary: Added copyright notice
//


#endif
