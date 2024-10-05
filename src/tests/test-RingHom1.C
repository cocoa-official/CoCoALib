//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BuildInfo.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ring.H"
#include "CoCoA/VectorOps.H"  // for printing indets(QQxyzt)


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;



//----------------------------------------------------------------------
// Test for RingHom
// Z->Q, Q->Q[x,y,z,t], Z->Fp, Q->Fpx, Q[x,y,z,t]->Fp[x], Fp[x]->Fp
// CompositeHom
//----------------------------------------------------------------------
namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations(UseNonNegResidues);

    // Create the rings we will be using.
    const FractionField QQ = RingQQ(); // used to fail for non-uniqueness of ZZ
    const QuotientRing Fp = NewZZmod(7);

    const PolyRing QQxyzt = NewPolyRing(QQ, symbols("x,y,z,t"));
    const PolyRing Fpx = NewPolyRing(Fp, symbols("x"));

//   const ring R1 = NewPolyRing(Qxyzt, symbols("alpha"));
//   const ring R2 = NewPolyRing(NewFractionField(R1), symbols("a"));
//   const ring R3 = NewPolyRing(NewFractionField(R2), symbols("b"));

    // Some simple homomorphisms...
//   RingHom Z2QQ = EmbeddingHom(QQ);
//   RingHom QQ2QQxyzt = CoeffEmbeddingHom(QQxyzt);
//   RingHom Z2Fp = QuotientingHom(Fp);
    RingHom ZZ2QQ = CanonicalHom(RingZZ(), QQ);
    RingHom QQ2QQxyzt = CanonicalHom(QQ, QQxyzt);
    RingHom ZZ2Fp = CanonicalHom(RingZZ(), Fp);
    RingHom QQ2Fpx = InducedHom(QQ, ZZEmbeddingHom(Fpx)); // only partial!
//   RingHom QQxyzt2R3 = ChainCanonicalHom(QQxyzt, R3);

//   CoCoA_ASSERT_ALWAYS( QQxyzt2R3(3) == RingElem(R3, 3) );

    // A homomorphism from the polynomial ring QQxyzt to Fpx.
//     vector<RingElem> v1;
//     v1.push_back(zero(Fpx));       // x |-> 0
//     v1.push_back(indet(Fpx,0));    // y |-> x in Fpx
//     v1.push_back(RingElem(Fpx,1)); // z |-> 1 in Fpx
//     v1.push_back(zero(Fpx));       // t |-> 0
    RingHom QQxyzt2Fpx = PolyRingHom(QQxyzt, Fpx, QQ2Fpx, "0,x,1,0");
    
    // A homomorphism from Fpx to Fp.
//     vector<RingElem> v2;
//     v2.push_back(RingElem(Fp,-1)); // x |-> -1
    RingHom Fpx2Fp = PolyRingHom(Fpx, Fp, IdentityHom(Fp), "-1");

    // testing phi(vector):
    cout << "QQxyzt2Fpx: " << indets(QQxyzt) << " map to " << QQxyzt2Fpx(indets(QQxyzt)) << endl;
    cout << "Fpx2Fp:     " << indets(Fpx) << " map to " << Fpx2Fp(indets(Fpx)) << endl;

    // Compose four homomorphisms
    RingHom CompositeHom = Fpx2Fp(QQxyzt2Fpx(QQ2QQxyzt(ZZ2QQ)));
    cout << "Composite homomorphism is " << CompositeHom << endl << endl;

    for (int i=1; i <= 10; ++i)
    {
      RingElem image = CompositeHom(i); // i is a machine int, and is mapped into ZZ automatically.
      cout << "Using composite hom: " << i << " maps to "
           << image << " in " << owner(image) << endl;
    }

    cout << "Image of  3  under ZZ2QQ      is " << ZZ2QQ(3) << endl;
    cout << "Image of 3/5 under QQ2QQxyzt  is " << QQ2QQxyzt(RingElem(QQ,3)/5) << endl;
    cout << "Image of 3/5 under QQxyzt2Fpx is " << QQxyzt2Fpx(RingElem(QQxyzt,3)/5) << endl;
    cout << "Image of 3/5 under QQ2Fpx     is " << QQ2Fpx(3/RingElem(QQ,5)) << endl; // works fine
    try
    {
      cout << "Image of 1/7 under QQ2Fpx is " << std::flush;
      cout << QQ2Fpx(1/RingElem(QQ,7)) << endl; // generates an error
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      if (err != ERR::BadPartialRingHomArg) throw; // rethrow if err is not the expected error.
      cout << "***" << endl
           << "Program correctly detected improper use of a partial ring hom." << endl;
    }
  }

} // end of namespace CoCoA


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
