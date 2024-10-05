//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti
//   Orig author:  2007 Massimo Caboara

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
#include "CoCoA/TmpIsTree.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/VectorOps.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>



namespace CoCoA
{

  void TestCycle()
  {
    SparsePolyRing P(NewPolyRing(RingQQ(), symbols("x,y,z,t,u,v")));
    RingElem x = indet(P,0),  y = indet(P,1),  z = indet(P,2);
    RingElem t = indet(P,3),  u = indet(P,4),  v = indet(P,5);

    PolyList PL;
    PL.push_back(x*y);
    PL.push_back(y*z);
    PL.push_back(z*t);
    PL.push_back(t*u);
    PL.push_back(u*v);
    PL.push_back(v*x);
    FacetComplex Cycle(P,PL);
    std::list<facet> Result=Cycle.myIsTreeNoOpt();
    cout << "myIsTreeNoOpt: Cycle in Complex = " << Result << endl;
    Result=Cycle.myIsTreeNoOpt();
    cout << "myIsTreeNoOpt: Cycle in Complex = " << Result << endl;
    Result=Cycle.myIsTreeCBNoOpt();
    cout << "myIsTreeCBNoOpt: Cycle in Complex = " << Result << endl;
    Result=Cycle.myIsTreeCBOpt();
    cout << "myIsTreeCBOpt: Cycle in Complex = " << Result << endl;
  }

  void TestLine()
  {
    SparsePolyRing P(NewPolyRing(RingQQ(), symbols("x,y,z,t,u,v")));
    RingElem x = indet(P,0),  y = indet(P,1),  z = indet(P,2);
    RingElem t = indet(P,3),  u = indet(P,4),  v = indet(P,5);

    PolyList PL;
    PL.push_back(x*y);
    PL.push_back(y*z);
    PL.push_back(z*t);
    PL.push_back(t*u);
    PL.push_back(u*v);
    FacetComplex Cycle(P,PL);
    CoCoA_ASSERT_ALWAYS(Cycle.myIsTreeNoOpt().empty());
    CoCoA_ASSERT_ALWAYS(Cycle.myIsTreeCBNoOpt().empty());
    CoCoA_ASSERT_ALWAYS(Cycle.myIsTreeCBOpt().empty());
    //  std::list<facet> Result=Cycle.myIsTreeNoOpt();
    //  cout << "myIsTreeNoOpt: Cycle in Complex = " << Result << endl;
    //  Result=Cycle.myIsTreeNoOpt(); 
    //  cout << "myIsTreeNoOpt: Cycle in Complex = " << Result << endl;
    //  Result=Cycle.myIsTreeCBNoOpt();
    //  cout << "myIsTreeCBNoOpt: Cycle in Complex = " << Result << endl;
    //  Result=Cycle.myIsTreeCBOpt();
    //  cout << "myIsTreeCBOpt: Cycle in Complex = " << Result << endl;
  }

  void TestLine1()
  {
    constexpr unsigned int K=1200;
    constexpr unsigned int W=120;
  
    SparsePolyRing P(NewPolyRing(RingQQ(), NewSymbols(K+W)));
    RingElem T(P);
    PolyList PL;
    for (unsigned int i=0;i!=K;++i)
    {  
      T=one(P);
      for (unsigned int j=i;j!=i+W;++j)
        T=T*indet(P,j);
      PL.push_back(T);
    }
    FacetComplex Cycle(P,PL);
    CoCoA_ASSERT_ALWAYS(Cycle.myIsTreeCBOpt().empty());
  }//TestLine1


  void program()
  {
    GlobalManager CoCoAFoundations;
  
    TestCycle();
    TestLine();
    //  TestLine1();// longer op, ~10s
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
