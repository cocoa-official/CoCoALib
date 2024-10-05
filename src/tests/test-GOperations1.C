//   Copyright (c)  2007-2012  John Abbott, Anna M. Bigatti

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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpGOperations.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/matrix.H"
#include "CoCoA/submodule.H"
#include "CoCoA/symbol.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for ideal/module operations (SAME EXAMPLES AS cocoa5.ts)
// functions: [GBasis], NF, IsElem, IsContained
// environments: DMP (Fp, Q), DMPI (Q), DMPII
//----------------------------------------------------------------------
namespace CoCoA
{

  //[Id = "cocoa5-01", Descr = "ReducedGBasis5x Module posto"];
  void test1()
  {
    //Use Q[x,y,z],PosTo;
    SparsePolyRing P(NewPolyRing(RingQQ(),symbols("x,y,z")));
    RingElem x = indet(P,0),  y = indet(P,1),  z = indet(P,2);
    FreeModule FM = NewFreeModule(P, 2, WDegPosnOrd);
    const vector<ModuleElem>& e = gens(FM);
    vector<ModuleElem> g;
    g.push_back(x*e[0] + z*e[1]);
    g.push_back(z*e[0] + y*e[1]);
    FGModule M = submodule(FM, g);
    cout << "TidyGens(M) = " << TidyGens(M) << endl;
    // [Vector(z, y), Vector(x, z), Vector(0, xy - z^2)];
  }


  //[Id = "cocoa5-02", Descr = "ReducedGBasis5x Module ToPos"];
  void test2()
  {
    //Use Q[x,y,z],ToPos;
    SparsePolyRing P(NewPolyRing(RingQQ(),symbols("x,y,z")));
    RingElem x = indet(P,0),  y = indet(P,1),  z = indet(P,2);
    FreeModule FM = NewFreeModule(P, 2, OrdPosn);
    const vector<ModuleElem>& e = gens(FM);
    vector<ModuleElem> g;
    g.push_back(x*e[0] + z*e[1]);
    g.push_back(z*e[0] + y*e[1]);
    FGModule M = submodule(FM, g);
    cout << "TidyGens(M) = " << TidyGens(M) << endl;
    // [Vector(z, y), Vector(x, z)];
  }


  //Test := Record[Id = "cocoa5-03", Descr = "GBasis5x param"];
  void test3()
  {
    //Use Z/(32003)[a, x,y], Ord(M);
    PolyRing Fpa = NewPolyRing(NewRingFp(32003), symbols("a"));
    ring K = NewFractionField(Fpa);
    PolyRing P(NewPolyRing(K, symbols("x,y")));
    RingElem x = indet(P,0),  y = indet(P,1);
    RingElem a = RingElem(P, symbol("a"));
    //    I := Ideal((a-1)*x+(a^2+a)*y, (a+1)*x + y);
    ideal I = ideal((a-1)*x+(power(a,2)+a)*y, (a+1)*x + y);
    vector<RingElem> GB = TidyGens(I);
    CoCoA_ASSERT_ALWAYS("cocoa5-3" && GB[0]==x && GB[1]==y);
  }

  // Test := Record[Id := "cocoa5-04", Descr := "GBasis5x Module shifts"];

  // Test.Input :="
  // Use QQ[x,y,z], Weights([1,2,4]);
  // M := Module([x^2-y,1],[x^4-z,y^2]);
  // Info5 := Record[];
  // Info5.ModuleShifts := Mat([[0,2]]);
  // X := GBasis5x(M, Info5);
  // X = [ Vector(x^2 - y, 1), Vector(y^2 - z, y^2 - x^2 - y)];
  // ";    CoCoAServerRegister(Test, 1 /* times True */);

  // -------------------------------
  // -- TEST
  // Test := Record[Id := "cocoa5-05", Descr := "ReducedGBasis5 Module OrdMat"];

  // Test.Input :="
  // OrdMat := Mat([[1,1,1],[2,1,1],[1,1,0]]); // Ring Grading (first 2 Rows) 
  //                                           // Plus order (last row)
  // Use ZZ/(101)[x,y,z], Ord(OrdMat), ToPos;
  // M := Module([y-x,0,0], [x,0,z], [0,y^2-z^2,0]);
  // X := ReducedGBasis5(Module(Gens(M)));
  // X= [Vector(x,0,z),Vector(y,0,z),Vector(0, 0, x*z-y*z),Vector(0, y^2-z^2, 0)];
  // ";    CoCoAServerRegister(Test, 1 /* times True */);

  // -------------------------------
  // -- TEST
  // Test := Record[Id := "cocoa5-06", Descr := "ReducedGBasis5x Module shifts gradingdim2"];

  // Test.Input :="
  // OrdMat := Mat([[1,1,1],[2,1,1],[1,1,0]]); // Ring Grading (first 2 Rows) 
  //                                           // Plus order (last row)
  // Use ZZ/(101)[x,y,z], Ord(OrdMat), ToPos;
  // M := Module([y-x,0,0], [x,0,z], [0,y^2-z^2,0]);
  // Info5 := Record[];
  // Info5.OrdMat := OrdMat;
  // Info5.GradingDim := 2;
  // Info5.ModuleShifts := Mat([[3,1,2],[2,2,5]]); // GrDim rows!!
  // X := ReducedGBasis5x(M, Info5);
  // X=[Vector(0, y^2-z^2, 0),Vector(x, 0, z),Vector(y, 0, z),Vector(0, 0, x*z - y*z)];
  // ";    CoCoAServerRegister(Test, 1 /* times True */);

  // Test := Record[Id := "cocoa5-07", Descr := "Elim5"];
  void test7()
  {
    matrix W = NewDenseMat(RowMat(RingElems(RingZZ(), "7, 4, 3, 1, 1, 1")));
    PolyRing P = NewPolyRingWeights(RingQQ(), symbols("x,y,z,w[3],w[4],w[5]"), W);
    ideal I(RingElems(P,
          "x - 7*w[4]^7 - 9*w[3]*w[4]*w[5]^5, y - 6*w[4]^4 - 2*w[3]^2*w[4]*w[5], z - 5*w[4]^3 - 2*w[3]*w[5]^2"));
    ideal E  = LT(elim(I, RingElems(P,"w[3],w[4]")));
    CoCoA_ASSERT_ALWAYS("cocoa5-7" && E == ideal(RingElems(P,"x^7")));
  }
  

  // -------------------------------
  // -- TEST
  // Test := Record[Id := "cocoa5-08", Descr := "Intersection5"];

  // Test.Input :="
  // Use QQ[x,y,z], Weights(1,2,1);
  // I := Ideal(x*y, z^2);
  // J := Ideal(y*z, x-z);
  // I5 := Intersection5(I, J);
  // I4 := Intersection(I, J);
  // I4 = I5;
  // ";    CoCoAServerRegister(Test, 1 /* times True */);

  // -------------------------------
  // -- TEST
  // Test := Record[Id := "cocoa5-09", Descr := "Intersection5 ModuleShifts"];

  // Test.Input :="
  // Use QQ[x,y,z], Weights(1,2,1);
  // M := Module([x^2*y, 0], [x*z^2, 0], [-x*z + x^2, -z + x]);
  // N := Module([x*z^2, 0], [x*y*z, y*z], [-x*z + x^2, -z + x]);
  // I5 := Intersection5x(M, N, Record[ModuleShifts := [[1,2]]]);
  // I4 := Intersection(M, N);
  // I4 = I5;
  // ";    CoCoAServerRegister(Test, 1 /* times True */);

  // -------------------------------
  // -- TEST
  // Test := Record[Id := "cocoa5-10", Descr := "Intersection5 Module"];

  // Test.Input :="
  // Use ZZ/(101)[x,y,z,t],PosTo;
  // M:=Module([(x+y+t)^3,y^3],[(x-y-z)^3,z^3],[x^2-y^2,x^2-z^2],[x*y*z+y*z*t+z*t*x,x^3]);
  // N:=Module([x^2+x*y,y^2]);

  // I4 := Intersection(M,N);// Long
  // I5 := Intersection5(M,N);// Long translation
  // I4 = I5;
  // ";    CoCoAServerRegister(Test, 1 /* times True */);


  //[Id = "cocoa5-11", Descr = "Intersection5x 2 param"];
  void test11()
  {
    SparsePolyRing R = NewPolyRing(NewZZmod(32003), symbols("a"));
    FractionField K = NewFractionField(R);
    SparsePolyRing P(NewPolyRing(K, symbols("x,y")));
    RingElem a = RingElem(P, symbol("a"));
    RingElem x = indet(P,0);
    RingElem y = indet(P,1);
// II := Ideal(x-y);
// I := (a-1) * x * II;
// J := (a+1) * y * II;
    ideal I = ideal((a-1) * x * (x-y));
    ideal J = ideal((a+1) * y * (x-y));
// X := Intersection5x(I,J, Record(NumParams=1));
    ideal X = intersect(I, J);
// X = x*y*II;
    CoCoA_ASSERT_ALWAYS("cocoa5-11" && X == ideal(x*y*(x-y)));
  }



  //[Id = "cocoa5-26",Descr = "Homogenized5x "];
  void test26()
  {
// M1:=Mat([[2, 1, 1,0],
//          [1, 2, 0,1],
//          [0, 1, 0,1],
//          [1,0,0,0]]);
// Use Q[x,y,a,b],Ord(M1),PosTo;
// Info5 := Record();
// Info5.GradingDim := 2;
// Info5.OrdMat :=M1;
// I := Ideal(x+1,y^2+x);
// X := Homogenized5x([a,b],I,Info5);
// X = Ideal(a^2b + x, xb^3 + y^2, x^2b^2 - y^2a^2, y^2a^4 + x^3b);
  }


  //[Id = "cocoa5-28",Descr = "Homogenized5 "];
  void test28()
  {
// Use Q[x,y,z,h];
// I := Ideal(x^2y+xy+1, xy^2+y^2+1);
// X := Homogenized5([h],I);
// X = Ideal(x - y, y^3 + y^2h + h^3);
    SparsePolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z,h"));
    RingElem x = RingElem(P,2);
    RingElem y = RingElem(P,3);
    RingElem z = RingElem(P,5);
    RingElem h = RingElem(P,7);
    RingElem y2 = power(y,2);
    ideal I(power(x,2)*y+x*y+1, x*y2+y2+1);
    ideal J(homog(I, h));
    CoCoA_ASSERT_ALWAYS("cocoa5-28" && J == ideal(x - y, power(y,3) + y2*h + power(h,3)));
  }


  //[Id = "cocoa5-30", Descr = "Intersection5x multi-param"];
  void test30()
  {
// II := Ideal(x-y);
// I := (a+b-1) * x * II;
// J := (a-b^2+1) * y * II;
// X := Intersection5x(I,J, Record(NumParams=2));
// X=x*y*II;
    SparsePolyRing R(NewPolyRing(NewZZmod(32003),symbols("a,b")));
    SparsePolyRing P(NewPolyRing(NewFractionField(R), symbols("x,y")));

    ideal II = ideal(RingElem(P, "x-y"));
    ideal I = ideal(RingElem(P, "x * (a+b-1)")) * II;
    ideal J = ideal(RingElem(P, "y * (a-b^2+1)")) * II;
    ideal X = intersect(I,J);
    CoCoA_ASSERT_ALWAYS("cocoa5-30" && X == ideal(RingElem(P, "x*y"))*II);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;
  
    test1();
    test2();
    test3();
    test7();  // elim
    test11();
    test30();
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
