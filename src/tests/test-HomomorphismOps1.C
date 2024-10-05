//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti
//   Original authors: Thomas Izgin and Filip Skrentny

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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/ideal.H"
#include "CoCoA/HomomorphismOps.H"

#include <vector>
using std::vector;

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test does virtually nothing, but is a handy template if you want
// to create your own test code: just copy this file and add your code
// after the line with "PUT YOUR CODE HERE"
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;
    const ring QQ = RingQQ();


    //------------------------ Example 1: phi,psi: PolyRing-> PolyRing ---------------------------------------------------
    PolyRing R = NewPolyRing(QQ, symbols("x,y")); //  R=QQ[x,y]
    PolyRing S = NewPolyRing(QQ, symbols("x,y")); //  S=QQ[x,y]

//     vector<RingElem> RImage1;
//     RImage1.push_back(RingElem(S,"y"));
//     RImage1.push_back(RingElem(S,"y"));
    
    RingHom phi = PolyAlgebraHom(R,S, "y,y"); 
    // phi: R --> S sending  x |--> y,  y |--> y
  
    CoCoA_ASSERT_ALWAYS(!IsInjective(phi));
    CoCoA_ASSERT_ALWAYS(ker(phi) == ideal(RingElem(R,"x-y")));

    CoCoA_ASSERT_ALWAYS(!IsSurjective(phi));
    CoCoA_ASSERT_ALWAYS(IsInImage(phi, one(S)));
    CoCoA_ASSERT_ALWAYS(IsOne(phi(preimage(phi, one(S)))));

    CoCoA_ASSERT_ALWAYS(IsInImage(phi, RingElem(S, "y")));
    CoCoA_ASSERT_ALWAYS(phi(preimage(phi, RingElem(S, "y"))) == RingElem(S, "y"));

    //-------------------------------------------------------
    // Now a function which is injective

//     vector<RingElem> RImage2;
//     RImage2.push_back(RingElem(S,"y"));
//     RImage2.push_back(RingElem(S,"x"));

    RingHom psi = PolyAlgebraHom(R,S, "y,x");
    // psi: R --> S sending  x |--> y,  y |--> x
  
    CoCoA_ASSERT_ALWAYS(IsInjective(psi));
    CoCoA_ASSERT_ALWAYS(IsZero(ker(psi)));
    CoCoA_ASSERT_ALWAYS(IsSurjective(psi));
    CoCoA_ASSERT_ALWAYS(IsInImage(psi,RingElem(S,"y^2")));
    CoCoA_ASSERT_ALWAYS(phi(preimage(psi, RingElem(S, "y^2"))) == RingElem(S, "y^2"));
    

    //----------------------- Example 2: phi2: QuotientRing -> PolyRing -----------------------------------------------

    // R=QQ[x,y] and S=QQ[x,y], phi as in Example 1
    
    ideal I1(RingElem(R,"x^2-y^2"));
    QuotientRing RmodI1 = NewQuotientRing(R, I1); // RmodI1= R/I1
  
    RingHom phi2 = InducedHom(RmodI1, phi);
    // phi2:  QQ[x,y]/<x^2-y^2> --> QQ[x,y] sending  (x) |--> y,  (y) |--> y
    CoCoA_ASSERT_ALWAYS(!IsInjective(phi2));
    CoCoA_ASSERT_ALWAYS(ker(phi2) == ideal(RingElem(RmodI1, "x-y")));

    CoCoA_ASSERT_ALWAYS(!IsSurjective(phi2));
    CoCoA_ASSERT_ALWAYS(!IsInImage(phi2,RingElem(S,"x*y^2")));
    CoCoA_ASSERT_ALWAYS(IsZero(preimage0(phi2, RingElem(S, "x*y^2"))));
    CoCoA_ASSERT_ALWAYS(IsInImage(phi2,RingElem(S,"y")));
    CoCoA_ASSERT_ALWAYS(phi2(preimage(phi2, RingElem(S, "y"))) == RingElem(S, "y"));
    

  //----------------------- Example 3.1: phi3: PolyRing -> QuotientRing -------------------------------------------------

     ideal I2(RingElem(S,"y+1")); // I2=<y+1> in S
     QuotientRing SmodI2 = NewQuotientRing(S, I2);
     RingHom phi3 = QuotientingHom(SmodI2);
//   phi3: QQ[x,y] --> QQ[x,y]/<y+1>,  x |--> (x), y |--> (y)
     CoCoA_ASSERT_ALWAYS(!IsInjective(phi3));
     CoCoA_ASSERT_ALWAYS(ker(phi3) == ideal(RingElem(S, "y+1")));
     CoCoA_ASSERT_ALWAYS(IsSurjective(phi3));
     CoCoA_ASSERT_ALWAYS(phi3(preimage(phi3, RingElem(SmodI2,"y^3+x"))) == RingElem(SmodI2,"y^3+x"));


     // Now a few more complicated examples
     //----------------------- Example 3.2: phi4: PolyRing -> QuotientRing -------------------------------------------------
   
     
     ring P1 = NewPolyRing(QQ, symbols("x,y,z"));

     ring P2 = NewPolyRing(QQ, symbols("a,b"));
     I2 = ideal(RingElem(P2, "a^2+4*b^2-1"));
     QuotientRing Q2 = NewQuotientRing(P2,I2);

//      vector<RingElem> images;
//      images.push_back(RingElem(Q2, "2*b"));
//      images.push_back(RingElem(Q2, "0"));
//      images.push_back(RingElem(Q2, "a"));

     RingHom phi4 = PolyAlgebraHom(P1, Q2, "2*b,0,a");

     // phi4: QQ[x,y,z] --> QQ[a,b]/<a^2+4b^2-1> sending  x |--> (2b)   y |--> (0)   z |--> (a)
     CoCoA_ASSERT_ALWAYS(!IsInjective(phi4));
     CoCoA_ASSERT_ALWAYS(ker(phi4) == ideal(RingElems(P1, "y, x^2+z^2-1")));
     CoCoA_ASSERT_ALWAYS(IsSurjective(phi4));
     CoCoA_ASSERT_ALWAYS(phi4(preimage(phi4, RingElem(Q2,"a"))) == RingElem(Q2,"a"));
    
    

 //----------------------- Example 4: psi2: QuotientRing -> QuotientRing ------------------------

     
     I1 = ideal(RingElem(P1, "x^2+y^2+z^2-1"));
     QuotientRing Q1 = NewQuotientRing(P1, I1);
     RingHom psi2 = InducedHom(Q1, phi4);
     // psi2: QQ[x,y,z]/<x^2+y^2+z^2-1> --> QQ[a,b]/<a^2+4b^2-1> sending
     // (x) |--> (2b)    (y) |--> (0)   (z) |--> (a)
     CoCoA_ASSERT_ALWAYS(!IsInjective(psi2));
     CoCoA_ASSERT_ALWAYS(ker(psi2) == ideal(RingElem(Q1, "y")));
     CoCoA_ASSERT_ALWAYS(IsSurjective(psi2));
     CoCoA_ASSERT_ALWAYS(psi2(preimage(psi2, RingElem(Q2,"a"))) == RingElem(Q2,"a"));

  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
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
