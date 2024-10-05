//   Copyright (c)  2022  John Abbott and Anna M. Bigatti
//   Original author:  Nico Mexis  (2022)

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


#include "CoCoA/SparsePolyOps-Graeffe.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-resultant.H"
#include "CoCoA/SparsePolyOps-RingElem.H"


namespace CoCoA
{

    /**
     * Checks the given polynomial for validity in the GraeffeSeq class
     */
    const RingElem &GraeffeSeq::CheckElementValidity(const RingElem &f)
    {
        const ring &Px = owner(f);
        const ring &P = CoeffRing(Px);
        if (!IsSparsePolyRing(Px))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
        if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
        const long IndetIndex = UnivariateIndetIndex(f);
        if (IndetIndex < 0)
            CoCoA_THROW_ERROR("Expected univariate poly", "GraeffeSeq ctor");
        if (!IsZZ(P) && !IsQQ(P))
            CoCoA_THROW_ERROR(ERR::BadRing, "GraeffeSeq ctor");
        return f;
    }

    /**
     * Retrieves the coefficients of the given polynomial mapped to the given polynomial ring
     */
    const std::vector<RingElem> GraeffeSeq::GetCoefficients(const RingElem &f, const ring &Q, const RingElem &xQ)
    {
        const ring &Px = owner(f);
        const ring &P = CoeffRing(Px);

        if (IsQQ(P))
        {
            const RingElem fQ = PolyRingHom(owner(f), Q, QQEmbeddingHom(CoeffRing(Q)), std::vector<RingElem>(NumIndets(owner(f)), xQ))(f);
            return CoeffVecWRT(monic(fQ), xQ);
        }
        if (IsZZ(P))
        {
            const RingElem fQ = PolyRingHom(owner(f), Q, ZZEmbeddingHom(CoeffRing(Q)), std::vector<RingElem>(NumIndets(owner(f)), xQ))(f);
            return CoeffVecWRT(monic(fQ), xQ);
        }
        CoCoA_THROW_ERROR(ERR::BadRing, "GraeffeSeq ctor");
        return {}; // Shut up compiler warnings
    }

    /**
     * Construct a default GraeffeSeq object
     */
    GraeffeSeq::GraeffeSeq(ConstRefRingElem f, const long N)
        : Q(NewPolyRing(RingQQ(), NewSymbols(NumIndets(owner(f))))),
          xQ(indet(Q, UnivariateIndetIndex(f))),
          myDeg(deg(CheckElementValidity(f))),
          myCoeffs(GetCoefficients(f, Q, xQ)),
          myA(std::vector<RingElem>(myDeg + 1)),
          myS(std::vector<RingElem>(myDeg * N + 1)),
          myResultCache(zero(owner(f))),
          phi(PolyRingHom(Q, owner(f), QQEmbeddingHom(CoeffRing(owner(f))), indets(owner(f))))
    {
    }

    /**
     * Retrieve the current Graeffe iteration
     */
    RingElem GraeffeSeq::operator*() const
    {
        if (myResultPow == myCurrPow)
            return myResultCache;

        RingElem Pn = power(xQ, myDeg);
        for (long k = 1; k <= myDeg; ++k)
        {
            myA[k] = -myS[k * myCurrPow];
            for (long j = 1; j < k; ++j)
            {
                CheckForInterrupt("GraeffeSeq::operator*");
                myA[k] -= myA[j] * myS[(k - j) * myCurrPow];
            }
            myA[k] /= k;
            Pn += myA[k] * power(xQ, myDeg - k);
        }
        myResultCache = phi(Pn / content(Pn));

        myResultPow = myCurrPow;
        return myResultCache;
    }

    /**
     * Calculate Newton coeffcients for next iteration
     */
    GraeffeSeq &GraeffeSeq::operator++()
    {
        const long bound = myDeg * (myCurrPow + 1);
        for (long k = myDeg * myCurrPow + 1; k <= bound; ++k)
        {
            const long index = myDeg - k;
            myS[k] = (index >= 0) ? -k * myCoeffs[index] : zero(Q);
            for (long j = std::max(1L, k - myDeg); j < k; ++j)
            {
                CheckForInterrupt("GraeffeSeq::operator++");
                myS[k] -= myS[j] * myCoeffs[index + j];
            }
        }
        ++myCurrPow;
        return *this;
    }

    /**
     * Calculate Newton coeffcients for next iteration
     */
    GraeffeSeq GraeffeSeq::operator++(int)
    {
        GraeffeSeq copy(*this);
        operator++();
        return copy;
    }

    /**
     * Provides a way to format-print the given GraeffeSeq object
     */
    std::ostream &operator<<(std::ostream &out, const GraeffeSeq &GSeq)
    {
        if (!out)
            return out; // short-cut for bad ostreams
        out << "GraeffeSeq(curr=" << *GSeq << ")";
        return out;
    }

    /**
     * Returns the n-th order Graeffe transform of f, i.e.
     * a univariate polynomial whose roots are the n-th
     * powers of the roots of f.
     */
    RingElem GraeffeN(RingElem f, long n)
    {
        if (n < 1)  CoCoA_THROW_ERROR2(ERR::ReqPositive, "arg2");
        const ring &Px = owner(f);
        if (!IsSparsePolyRing(Px))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
        const ring &P = CoeffRing(Px);
        const long IndetIndex = UnivariateIndetIndex(f);
        if (IndetIndex < 0)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
        if (!IsZZ(P) && !IsQQ(P))  CoCoA_THROW_ERROR1(ERR::BadRing);

        // graeffe and graeffe3 are way faster; apply them as often as possible
        while (n % 2 == 0)
        {
            f = graeffe(f);
            n /= 2;
        }
        while (n % 3 == 0)
        {
            f = graeffe3(f);
            n /= 3;
        }

        if (n == 1)
            return f;

        // For the "rest", use the resultant approach (in bivariate ring)
        const SparsePolyRing Pxy = NewPolyRing(P, NewSymbols(2));
        const RingElem &x = indet(Pxy, 0);
        const RingElem &y = indet(Pxy, 1);
        RingElem F = PolyAlgebraHom(Px, Pxy, std::vector<RingElem>(NumIndets(Px), x))(f); // f in Pxy
        RingHom MapXtoY = PolyAlgebraHom(Pxy,Pxy, std::vector<RingElem>{y,zero(Pxy)});
        // loop below does res(f(y), y^n-x, y) one prime at a time -- seems to be faster
        const auto facs = factor(n);
        const std::vector<long>& plist = facs.myFactors(); // assume in incr order
        const std::vector<long>& mult = facs.myMultiplicities();
        for (int i=len(plist)-1; i >= 0; --i)
        {
          const long p = plist[i];
          for (long j=0; j < mult[i]; ++j)
            F = resultant(MapXtoY(F), power(y,p)-x, 1); // last arg refers to "y" == indet(Pxy,1)
//        RingElem res = resultant(PolyAlgebraHom(Px, Pxy, std::vector<RingElem>(NumIndets(Px), y))(f), power(y, n) - x, 1);
        }
        return PolyAlgebraHom(Pxy, Px, {indet(Px, IndetIndex), one(Px)})(F);
//        return PolyAlgebraHom(Q, Px, {indet(Px, IndetIndex), one(Px)})(res);
    }

}


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-Graeffe.C,v 1.8 2024/08/01 16:52:44 bigatti Exp $
// $Log: SparsePolyOps-Graeffe.C,v $
// Revision 1.8  2024/08/01 16:52:44  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;
// now using CoCoA_THROW_ERROR1/2
//
// Revision 1.7  2024/07/02 15:36:15  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.6  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.5  2024/01/21 19:57:24  abbott
// Summary: Improved err mesg
//
// Revision 1.4  2023/11/11 20:28:48  abbott
// Summary: Improved GraeffeN when order has several primes greater than 3
//
// Revision 1.3  2023/08/06 10:12:56  abbott
// Summary: Remove incorrect "sign trick" (redmine 1758)
//
// Revision 1.2  2022/11/24 21:11:29  abbott
// Summary: Added copyright notice
//
//
