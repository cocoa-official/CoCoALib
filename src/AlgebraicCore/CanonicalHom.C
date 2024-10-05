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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/CanonicalHom.C,v 1.15 2022/02/18 14:11:53 abbott Exp $
// $Log: CanonicalHom.C,v $
// Revision 1.15  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.14  2021/01/07 15:07:02  abbott
// Summary: Corrected copyright
//
// Revision 1.13  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.12  2019/12/11 14:50:48  abbott
// Summary: Removed some cruft
//
// Revision 1.11  2018/05/17 16:05:43  bigatti
// -- renamed TmpChainCanonicalHom --> ChainCanonicalHom
//
// Revision 1.10  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.9  2014/07/08 08:33:18  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.8  2014/07/07 12:12:17  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.7  2012/02/10 10:26:40  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.6  2012/02/08 15:07:08  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.5  2011/02/18 12:56:08  bigatti
// -- added TmpChainCanonicalHom
//
// Revision 1.4  2009/07/24 14:21:03  abbott
// Cleaned up include directives, and added some fwd decls.
//
// Revision 1.3  2009/07/24 12:27:46  abbott
// Added some comments.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2007/03/07 11:18:01  bigatti
// -- fixed CanonicalHom goto flag
//
// Revision 1.1  2007/03/05 21:25:02  cocoa
// New CanonicalHom pseudo-ctor.
//
//
