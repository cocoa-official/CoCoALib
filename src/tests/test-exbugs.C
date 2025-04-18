//   Copyright (c)  2005,2019  John Abbott, Anna M. Bigatti

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


#include "CoCoA/library.H"


// Includes from the standard C++ library
#include <algorithm>
using std::max;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

// Many exbugs: most checked in own procedure.

namespace CoCoA
{

  // Minor bugs that were never put into redmine:
  void TestSundry()
  {
    ring Pxy = NewPolyRing(RingQQ(), symbols("x,y"));
    CoCoA_ASSERT_ALWAYS(!IsDenseUPolyRing(Pxy)); // used to give true! 8-O
    ring Px = NewPolyRing_DUP(RingQQ(), symbol("x"));
    CoCoA_ASSERT_ALWAYS(IsDenseUPolyRing(Px));
  }


  // This bug originally reported by Dag Arneson 2005/04/01
  // It used to provoke SEGV.
  void OldBug1()
  {
    // 2019-09-11 JAA: presumably there was a bug in NewQuotientRing (redmine ???)
    PolyRing Qp = NewPolyRing(RingQQ(), SymbolRange("y", 0, 2)); // Q[y1..yn]

    RingElem foo(Qp);
    foo = 5*indet(Qp,0)*indet(Qp,1)+2*indet(Qp,0)*indet(Qp,2);
    ideal I(foo);
    QuotientRing R = NewQuotientRing(Qp,I);
  }


  void Bug_FactorMultiplicity()
  {
    const BigInt N = factorial(21);
    CoCoA_ASSERT_ALWAYS(FactorMultiplicity(9,N) == 4); // used to be infinite loop
  }


  void redmine858()
  {
    // Check that floor for TwinFloats does not throw ERR::ShouldNeverGetHere
    const ring RR = NewRingTwinFloat(64);
    const RingElem x = one(RR);
    for (int k = 1; k < 200; ++k)
    {
      const RingElem smaller = x - power(BigRat(1,2), k); // k=145 used to trigger ERR::ShouldNeverGetHere in floor
      try
      {
        const BigInt N = floor(smaller);
        CoCoA_ASSERT_ALWAYS(IsZero(N) || IsOne(N));
      }
      catch (const RingTwinFloat::InsufficientPrecision& err)
      {}
    }
  }


  void redmine870()
  {
    using std::vector;
    const ring Qxy = NewPolyRing(RingQQ(), symbols("x,y"));
    const RingElem x = indet(Qxy,0);
    const RingElem y = indet(Qxy,1);
    vector<RingElem> v;
    v.push_back(x); v.push_back(y);
    const ideal I1 = ideal(v);
    v.clear();
    v.push_back(x-1); v.push_back(y-1);
    const ideal I2 = ideal(v);
    // These two ideals should be equal:
    const ideal I3 = I1*I2;
    const ideal I4 = I2*I1;
    // Check explicitly that GBasis elems reduce to zero
    const vector<RingElem>& GB3 = GBasis(I3);
    const vector<RingElem>& GB4 = GBasis(I4);
    CoCoA_ASSERT_ALWAYS(len(GB3) >= 2);
    CoCoA_ASSERT_ALWAYS(len(GB4) >= 2);
    for (int i=0; i < len(GB3); ++i)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I3)));
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I4)));
    }
    for (int i=0; i < len(GB4); ++i)
    {
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I3)));
      CoCoA_ASSERT_ALWAYS(IsZero(NF(GB3[i],I4)));
    }
  }



  void redmine951comment7()
  {
    // gcd call in last line used to SEGV.
    const ring F101 = NewZZmod(101);
    const ring P = NewPolyRing(F101, symbols("z"));
    const RingElem z = indet(P,0);
    const RingElem f = RingElem(P, "-19*z^111 -26*z^110 +23*z^109 -28*z^108 -36*z^107 +10*z^106 -47*z^105 +16*z^104 -19*z^103 -20*z^102 +20*z^101 -39*z^100 +45*z^99 +42*z^98 +15*z^97 -36*z^96 +40*z^95 -11*z^94 -24*z^93 -21*z^92 -8*z^91 +38*z^90 -44*z^89 -50*z^88 +13*z^87 +12*z^86 +37*z^85 +46*z^84 +37*z^83 -45*z^82 +35*z^81 -3*z^80 +39*z^79 -12*z^78 +45*z^77 +10*z^76 -22*z^75 -25*z^74 +20*z^73 -2*z^72 +18*z^71 +34*z^70 +41*z^69 -38*z^68 +19*z^67 -49*z^66 -48*z^65 +47*z^64 -20*z^63 -31*z^62 -8*z^61 -7*z^60 -4*z^59 -19*z^58 -18*z^57 +28*z^56 +9*z^55 +49*z^54 -13*z^53 +21*z^52 -14*z^51 +8*z^50 +22*z^49 -9*z^48 -18*z^47 -43*z^46 +48*z^45 -11*z^44 +21*z^43 +45*z^42 +20*z^41 -10*z^40 +37*z^39 -49*z^38 +5*z^37 -23*z^36 +38*z^35 -13*z^34 +12*z^33 +32*z^32 -44*z^31 +19*z^30 -25*z^29 -21*z^28 -36*z^27 -47*z^26 +19*z^25 +17*z^24 -30*z^23 +34*z^22 +43*z^21 +5*z^20 +3*z^19 +15*z^18 -36*z^17 +28*z^16 -49*z^15 +23*z^14 -32*z^13 -10*z^12 -50*z^11 -29*z^10 +30*z^9 -30*z^8 +5*z^7 +34*z^6 +44*z^5 -50*z^4 -44*z^3 +44*z^2");

    const RingElem Df = deriv(f, z);
    CoCoA_ASSERT_ALWAYS(monic(gcd(f,Df)) == z);
  }


  void redmine1185()
  {
    // See redmine #1185: RemainingFactor used to have wrong sign
    const ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    RingElem f(P, "(y+z)*(y^2-x)");
    factorization<RingElem> facs = factor(f);
    const std::vector<RingElem>& fac = facs.myFactors();
    const std::vector<long>& mult = facs.myMultiplicities();
    CoCoA_ASSERT_ALWAYS(len(fac) == 2);
    CoCoA_ASSERT_ALWAYS(mult[0] == 1);
    CoCoA_ASSERT_ALWAYS(mult[1] == 1);

    const RingElem g = facs.myRemainingFactor() * fac[0] * fac[1];
    CoCoA_ASSERT_ALWAYS(f == g);
  }


  void redmine1300()
  {
    vector<BigInt> v;
    v.push_back(BigInt(2));
    v.push_back(BigInt(-3));
    CoprimeFactorBasis_BigInt CFB;
    CFB.myAddInfo(v);
    const vector<BigInt> B = FactorBase(CFB);
    CoCoA_ASSERT_ALWAYS(len(B) == 2);
    CoCoA_ASSERT_ALWAYS(B[0] == 2);
    CoCoA_ASSERT_ALWAYS(B[1] == 3);
  }


  void redmine1310()
  {
    const PolyRing QQx = NewPolyRing(RingQQ(), symbols("x"));
    const RingElem x = indet(QQx,0);
    const RingElem f = RingElem(QQx, "x^6 +4*x^5 -x^3 -4*x^2 -1");
    const BigRat RBf = RootBound(f);
    CoCoA_ASSERT_ALWAYS(RBf > 4); // should be greater than 4.0009
    const RingElem g = RingElem(QQx, "x^4 -288*x^3 -593*x +256");
    const BigRat RBg = RootBound(g);
    CoCoA_ASSERT_ALWAYS(RBg > 288); // should be greater than about 288.01
  }


  void redmine1322()
  {
    PolyRing QQx = NewPolyRing(RingQQ(), symbols("x"));
    ideal I(zero(QQx));
    const auto RGB = GBasis(I);
  }


  void redmine1379()
  {
    PolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    ideal I1(RingElems(P, "x,y"));
    CoCoA_ASSERT_ALWAYS(!IsZeroDim(I1));
    ideal I2(RingElems(P, "x,z"));
    CoCoA_ASSERT_ALWAYS(!IsZeroDim(I2));
    ideal I3(RingElems(P, "y,z"));
    CoCoA_ASSERT_ALWAYS(!IsZeroDim(I3));
    ideal I4(RingElems(P, "x,y,z"));
    CoCoA_ASSERT_ALWAYS(IsZeroDim(I4));

    ideal I5(RingElems(P,"y*z -2*z^5 -2*z^4 +2*z^3,"
	   "y^2 -16*z^5 -12*z^4 +12*z^3 +12*z^2,"
	   "x*y -2*x*z -5*y +6*z^5 +2*z^3 -8*z^2 +10*z,"
	   "x^2 +10*x*z -10*x -51*z^5 +74*z^4 +16*z^3 -14*z^2 -50*z +25,"
	   "x*z^2 -x*z +5*z^5 -6*z^4 -4*z^3 +5*z,"
	   "z^6 -z^5 -z^4 +z^2,  z^3 -z -1"));
    CoCoA_ASSERT_ALWAYS(IsZeroDim(I5));
  }


  void redmine1382()
  {
    PolyRing QQx = NewPolyRing(RingQQ(), symbols("x"));
    vector<RingElem> v; // empty vector
    try
    {
      RingElem g = CommonDenom(v);
      CoCoA_ASSERT_ALWAYS(!"redmine1382 is back!");
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      if (err != ERR::ReqNonEmpty) throw;
    }
  }


  void redmine1449()
  {
    PolyRing QQx = NewPolyRing(RingQQ(), symbols("x"));
    const RingElem& x = indet(QQx,0);
    vector<long> L = {21458,21469,21470,21484,21496,21497};
    const RingElem g = (x+11)*(x*x+1);
    for (auto n: L)
    {
      auto facs = factor(g*(power(x,4)+n));
      CoCoA_ASSERT_ALWAYS(len(facs.myMultiplicities()) == 3);
    }
  }


  void redmine1331()
  {
    ring ZZ25 = NewZZmod(25);
    ConstMatrix I7 = IdentityMat(ZZ25, 7);
    matrix adj7 = adj(I7);
    CoCoA_ASSERT_ALWAYS(adj7 == I7);
  }


    void TestBigIntProxy()
  {
    // The call to max used to cause a SEGV -- an undesired "feature" of the proxy class used to reduce wasteful copying upon returning BigInt values.
    // This bug was originally reported by Max Caboara.

    const BigInt a(6);
    const BigInt b(10);
    const BigInt c(15);

    const BigInt d = max(gcd(a,b), gcd(b,c));
    CoCoA_ASSERT_ALWAYS(d == 5);
  }


  void TestPolyPrinting()
  {
    constexpr int p = 13; // a small prime number
    PolyRing P1 = NewPolyRing(NewZZmod(p), symbols("x"));
    const RingElem x = indet(P1,0);
    cout << x-1 << endl;

    constexpr int n = 20; // a small composite number
    PolyRing P2 = NewPolyRing(NewZZmod(n), symbols("y"));
    const RingElem y = indet(P2,0);
    cout << y-1 << endl;

    const BigInt P = NextProbPrime(power(10,20)); // a large prime number
    PolyRing P3 = NewPolyRing(NewZZmod(P), symbols("z"));
    const RingElem z = indet(P3, 0);
    cout << z-1 << endl;

    const BigInt N = power(10,20); // a large composite number
    PolyRing P4 = NewPolyRing(NewZZmod(N), symbols("t"));
    const RingElem t = indet(P4, 0);
    cout << t-1 << endl;
  }


  void TestIsIndetPosPower()
  {
    PolyRing P1 = NewPolyRing(NewZZmod(3), symbols("x"));
    CoCoA_ASSERT_ALWAYS(!IsIndetPosPower(one(P1)));
    CoCoA_ASSERT_ALWAYS(IsIndetPosPower(RingElem(P1,"x")));
    CoCoA_ASSERT_ALWAYS(IsIndetPosPower(RingElem(P1,"x^2")));
    CoCoA_ASSERT_ALWAYS(!IsIndetPosPower(RingElem(P1,"x^2+1")));
    PolyRing P2 = NewPolyRing(NewZZmod(3), symbols("x,y"));
    CoCoA_ASSERT_ALWAYS(!IsIndetPosPower(one(P2)));
    CoCoA_ASSERT_ALWAYS(IsIndetPosPower(RingElem(P2, "y")));
    CoCoA_ASSERT_ALWAYS(IsIndetPosPower(RingElem(P2, "y^2")));
    CoCoA_ASSERT_ALWAYS(IsIndetPosPower(RingElem(P2, "x")));
    CoCoA_ASSERT_ALWAYS(!IsIndetPosPower(RingElem(P2, "x*y")));
  }


  void TestCoeffEmbeddingHom()
  {
    // Checks something to do with CoeffEmbeddingHom...
    const ring QQ = RingQQ();
    SparsePolyRing QQx = NewPolyRing(QQ, symbols("x"));
    const vector<RingElem>& x = indets(QQx);

    const RingElem half(QQ, "1/2");
    const RingHom phi = CoeffEmbeddingHom(QQx);
    RingElem f(QQx);
    f = x[0] + phi(half);
    CoCoA_ASSERT_ALWAYS(f == RingElem(QQx, "x+1/2"));
    const ideal I = ideal(f);
    vector<RingElem> GBasis = TidyGens(I);
    CoCoA_ASSERT_ALWAYS(!GBasis.empty() && GBasis[0] == f);
  }


  void TestNewMatrixOrdering()
  {
    // This test checks that a matrix with non integer entries
    // is rejected by NewMatrixOrdering.

    matrix M = NewDenseMat(RingQQ(),2,2);
    SetEntry(M,0,0,BigRat(1,2));
    SetEntry(M,0,1,BigRat(1));
    SetEntry(M,1,0,BigRat(3,2));
    SetEntry(M,1,1,BigRat(2));

    try
    {
      PPOrdering ord = NewMatrixOrdering(M,1);
      CoCoA_ASSERT_ALWAYS(false && "SHOULD NEVER GET HERE!");
    }
    catch (const ErrorInfo& err)
    {
      // Dispose of expected error; o/w rethrow
      if (err != ERR::BadArg) throw;
    }

    // This test checks that a matrix over a ring of non-zero
    // characteristic is rejected by Newmatrixordering.
    matrix M2 = NewDenseMat(NewZZmod(5),2,2);
    SetEntry(M,0,0,1);
    SetEntry(M,0,1,1);
    SetEntry(M,1,0,1);
    SetEntry(M,1,1,2);

    try
    {
      PPOrdering ord = NewMatrixOrdering(M2,1);
      CoCoA_ASSERT_ALWAYS(false && "SHOULD NEVER GET HERE!");
    }
    catch (const ErrorInfo& err)
    {
      // Dispose of expected error; o/w rethrow
      if (err != ERR::BadRing) throw;
    }
  }


  void redmine1484()
  {
    const ring R = NewPolyRing(RingQQ(), symbols("x,y,z"));
    const ideal I(RingElem(R,"x-y"));
    const QuotientRing QR = NewQuotientRing(R,I);
    const RingHom phi = PolyAlgebraHom(R,QR,"x,y,z");
    const ideal K = ker(InducedHom(QR, phi));
    CoCoA_ASSERT_ALWAYS(IsZero(K));
  }


  void redmine1565()
  {
    // Checks IsInRadical when ideal has a zero gen
    const ring QQ = RingQQ();
    const SparsePolyRing QQx = NewPolyRing(QQ, symbols("x"));
    const RingElem& x = RingElem(QQx,"x");
//    SetVerbosityLevel(99);
    CoCoA_ASSERT_ALWAYS(IsInRadical(x, ideal(RingElems(QQx, "x^2,0"))));
  }

  void redmine1570()
  {
    // Occasional error in FloorLogBase when very close to exact power
    const BigRat q = BigRatFromString("531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219/531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728128");
    for (int pwr = 0; pwr < 10; ++pwr)
      CoCoA_ASSERT_ALWAYS(pwr-10 == FloorLog10(q*power(10,pwr)));
  }

  void redmine1610()
  {
    // Checks IsInRadical when ideal is (0), or GradingDim=0
    const ring QQ = RingQQ();
    const SparsePolyRing QQx = NewPolyRing(QQ, symbols("x,y"), lex);
    const RingElem& x = RingElem(QQx,"x");
    vector<RingElem> v;
    ideal I = ideal(QQx, v);
    CoCoA_ASSERT_ALWAYS(!IsInRadical(x,I));
    v.push_back(x*x);
    I = ideal(QQx, v);
    CoCoA_ASSERT_ALWAYS(!IsInRadical(x,I));
  }

  void redmine1690()
  {
    // resultant when one input is factor of the other
    ring P = NewPolyRing(RingQQ(), symbols("x"));
    RingElem x = indet(P,0);
    RingElem f = x*x*x-2;
    RingElem g = 3*x/x;
    RingElem R = resultant(f,g,0);
    CoCoA_ASSERT_ALWAYS(R == 27);
    vector<RingElem> seq = SubresultantSeq(f,g,0);
    CoCoA_ASSERT_ALWAYS(len(seq) == 1 && seq[0] == 27);

    f = 27*x*x*x-1;
    g = 3*x-1;
    R = resultant(f,g,0);
    CoCoA_ASSERT_ALWAYS(IsZero(R));
    seq = SubresultantSeq(f,g,0);
    CoCoA_ASSERT_ALWAYS(len(seq) == 1 && seq[0] == 0);
  }


  void redmine1710()
  {
    const ring P = NewPolyRing(RingZZ(), symbols("x"));
    const RingElem f(P, "2*x+2");
    CoCoA_ASSERT_ALWAYS(IsSqFree(f));
    CoCoA_ASSERT_ALWAYS(IsIrred(f));
    const auto facs1 = factor(-f);
    CoCoA_ASSERT_ALWAYS(facs1.myRemainingFactor() == -2);
    const auto facs2 = SqFreeFactor(-f);
    CoCoA_ASSERT_ALWAYS(facs2.myRemainingFactor() == -2);
  }

  void bug_20221125()
  {
    // bad bugs in IsPalindromic
    const ring P = NewPolyRing(RingZZ(), symbols("x"));
    const RingElem x = indet(P,0);
    const RingElem f(P, "x+1");
    CoCoA_ASSERT_ALWAYS(IsPalindromic(f));
    CoCoA_ASSERT_ALWAYS(IsPalindromic(2*f));
    CoCoA_ASSERT_ALWAYS(IsPalindromic(power(f,4)));
    CoCoA_ASSERT_ALWAYS(IsPalindromic(power(f,4)-6*x*x));

    CoCoA_ASSERT_ALWAYS(!IsPalindromic(f-2));
  }


  void redmine1714()
  {
    const ring P = NewPolyRing(RingQQ(), symbols("x"));
    const RingElem f(P, "x^4-x");
    const ring PmodI = NewQuotientRing(P,ideal(f));
    const RingElem x(PmodI, "x");
    CoCoA_ASSERT_ALWAYS(IsDivisible(RingElem(PmodI,"x^3-2*x+1"), x*x-1));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(x, x-1));
  }


  void redmine1725()
  {
    const ring Px = NewPolyRing(RingZZ(), symbols("x"));
    const RingElem f(Px, "x^2-3*x+2");
    const RingElem g(Px, "4*x^2-6*x+2");
    const RingElem R = resultant(f,g,0); // 0 = indet index of x
    CoCoA_ASSERT_ALWAYS(IsZero(R));
    CoCoA_ASSERT_ALWAYS(IsZero(resultant(f,g/2,0))); // 0 = indet index of x

    const ring Pxy = NewPolyRing(RingQQ(), symbols("x,y"));
    const RingElem f2(Pxy, "x^4+x^3+2*x^2-x+1");
    const RingElem g2(Pxy, "(x*y)^4+(x*y)^3+2*(x*y)^2-x*y+1");
    const RingElem R2 = resultant(f2,g2,0); // 0 = indet index of x
    CoCoA_ASSERT_ALWAYS(!IsZero(R2) && deg(R2,1) == 16); // 1 = indet index of y
  }


  void bug_20221219_rank()
  {
    // sent by email from Bernhard Andraschko 2022-12-19
    const matrix M1 = NewDenseMat(RingQQ(), 0,1);
    CoCoA_ASSERT_ALWAYS(rk(M1) == 0);
    const matrix M2 = NewDenseMat(RingQQ(), 1,0);
    CoCoA_ASSERT_ALWAYS(rk(M2) == 0);
  }


  void bug_20231209_LRS()
  {
    ring Px = NewPolyRing(RingZZ(), symbols("x"));
    RingElem f(Px,"16*x^8 +56*x^7 +120*x^6 +184*x^5 +209*x^4 +184*x^3 +120*x^2 +56*x +16");
    RingElem g(Px, "16*x^8 -64*x^7 +120*x^6 -146*x^5 +149*x^4 -146*x^3 +120*x^2 -64*x +16");

    CoCoA_ASSERT_ALWAYS(IsLRSDegenerate(f));
    CoCoA_ASSERT_ALWAYS(IsLRSDegenerate(g));
    CoCoA_ASSERT_ALWAYS(IsLRSDegenerateOrder(f*g,3));
  }

  void bug_20231222_cyclotomic()
  {
    ring P = NewPolyRing(RingZZ(), symbols("x"));
    const RingElem x = indet(P,0);
    RingElem c3 = cyclotomic(3,x);
    RingElem c4 = cyclotomic(4,x);
    RingElem c6 = cyclotomic(6,x);
    vector<long> L = CyclotomicFactorIndices(c3*c4*c6);

    CoCoA_ASSERT_ALWAYS(L[0] == 3);
    CoCoA_ASSERT_ALWAYS(L[1] == 4);
    CoCoA_ASSERT_ALWAYS(L[2] == 6);
  }


  void redmine1790()
  {
    PolyRing Px = NewPolyRing(RingQQ(), symbols("x"));
    PolyRing Py = NewPolyRing(RingQQ(), symbols("x"));
    ideal I(zero(Px));
    ideal Ignore1 = saturate(I,I);
    ideal J(zero(Py));
    try { ideal Ignore2 = saturate(I,J); } catch (const CoCoA::ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::MixedRings); }
  }

  //-------------------------------------------------------
  void program()
  {
    GlobalManager CoCoAFoundations;

    OldBug1();
    Bug_FactorMultiplicity();
    redmine858();
    redmine870();
    redmine951comment7();
    redmine1185();
    redmine1300();
    redmine1310();
    redmine1322();
    redmine1331();
    redmine1379();
    redmine1382();
    redmine1449();
    redmine1484();
    redmine1565();
    redmine1570();
    redmine1690();
    redmine1710();
    redmine1714();
    bug_20221125();
    bug_20221219_rank();
    redmine1725();
    bug_20231209_LRS();
    bug_20231222_cyclotomic();
    redmine1790();

    TestSundry();
    TestBigIntProxy();
    TestPolyPrinting(); // produces output!
    TestIsIndetPosPower();
    TestCoeffEmbeddingHom();
    TestNewMatrixOrdering();
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
