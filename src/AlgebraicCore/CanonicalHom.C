//   Copyright (c)  2007,2009  John Abbott and Anna M. Bigatti

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


#include "CoCoA/CanonicalHom.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/FractionField.H"

// #include <iostream>  // for debugging only

namespace CoCoA
{

  RingHom CanonicalHom(const ring& domain, const ring& codomain)
  {
    if (domain == codomain) return IdentityHom(domain);

    // Check codomain first, as this makes it possible to exploit certain "shortcuts"
    if (IsFractionField(codomain))
    {
///      return EmbeddingHom(codomain)(CanonicalHom(domain, BaseRing(codomain)));
      if (domain == BaseRing(codomain))
        return EmbeddingHom(codomain);
      goto CheckDomain;
    }
    if (IsPolyRing(codomain))
    {
///      return CoeffEmbeddingHom(P)(CanonicalHom(domain, CoeffRing(P)));
      if (domain == CoeffRing(codomain))
        return CoeffEmbeddingHom(codomain);
      goto CheckDomain;
    }
    if (IsQuotientRing(codomain))
    {
      const QuotientRing QR = codomain;
///      return QuotientingHom(QR)(CanonicalHom(domain, BaseRing(QR)));
      if (domain == BaseRing(QR))
        return QuotientingHom(QR);
      goto CheckDomain;
    }
  CheckDomain:
    // Two easy cases:
    if (IsZZ(domain)) return ZZEmbeddingHom(codomain);
    if (IsQQ(domain)) return QQEmbeddingHom(codomain); // NB result is only a partial hom!!

    CoCoA_THROW_ERROR(ERR::CanonicalHomFail, "CanonicalHom(R1,R2)");
    return IdentityHom(codomain); // Never executed; just to keep the compiler quiet.
  }


  RingHom ChainCanonicalHom(const ring& domain, const ring& codomain)
  {
    try { return CanonicalHom(domain, codomain); }
    catch (const CoCoA::ErrorInfo& err) {if (err!=ERR::CanonicalHomFail) throw;}
    if (IsFractionField(codomain))
    {
      return EmbeddingHom(codomain) (ChainCanonicalHom(domain, BaseRing(codomain)));
    }
    if (IsPolyRing(codomain))
    {
      return CoeffEmbeddingHom(codomain) (ChainCanonicalHom(domain, CoeffRing(codomain)));
    }
    if (IsQuotientRing(codomain))
    {
      const QuotientRing QR = codomain;
      return QuotientingHom(QR) (ChainCanonicalHom(domain, BaseRing(QR)));
    }
    CoCoA_THROW_ERROR(ERR::CanonicalHomFail, "ChainCanonicalHom(R1,R2)");
    return IdentityHom(codomain); // Never executed; just to keep the compiler quiet.
  }


} // end of namespace CoCoA
