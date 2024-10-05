//   Copyright (c)  2013  Anna M. Bigatti

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
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  void test_Qi(bool print)
  {
    if (print) std::cout << "------- Q[i] --------" << std::endl;
    PolyRing Qi = NewPolyRing(RingQQ(), symbols("i"));
    QuotientRing QR = NewQuotientRing(Qi, ideal(1+power(indet(Qi,0),2)));
    RingElem i(QR, symbol("i"));
    RingElem z1 = 2+3*i;
    if (print) std::cout << "z1 = " << z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(IsInvertible(z1));
    CoCoA_ASSERT_ALWAYS(!IsZeroDivisor(z1));
    CoCoA_ASSERT_ALWAYS(1/z1 == (-3*i +2)/13);
    if (print) std::cout << "1/z1 = " << 1/z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(IsIntegralDomain(QR));
    CoCoA_ASSERT_ALWAYS(IsField(QR));
    ideal I(2+i, 3*i);
    if (print) std::cout << "I = " << I << std::endl;
    CoCoA_ASSERT_ALWAYS(IsOne(I));
    CoCoA_ASSERT_ALWAYS(len(TidyGens(I)) == 1);
    // std::cout << "IsPrime(I) = " << IsPrime(I) << std::endl;

    if (print) std::cout << "------- Q(i)[x,y] --------" << std::endl;
    PolyRing QR_xy = NewPolyRing(QR, symbols("x,y"));
    i = RingElem(QR_xy, symbol("i"));
    z1 = 2+3*i;
    if (print) std::cout << "z1 = " << z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(IsInvertible(z1));
    CoCoA_ASSERT_ALWAYS(!IsZeroDivisor(z1));
    CoCoA_ASSERT_ALWAYS(1/z1 == (-3*i +2)/13);
    if (print) std::cout << "1/z1 = " << 1/z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(IsIntegralDomain(QR));
    CoCoA_ASSERT_ALWAYS(IsField(QR));
  }


  void test_xy(bool print)
  {
    if (print) std::cout << "------- QQ[x,y]/(x*y) --------" << std::endl;
    PolyRing Qxy = NewPolyRing(RingQQ(), symbols("x,y"));
    QuotientRing QR = NewQuotientRing(Qxy, ideal(product(indets(Qxy))));
    RingElem x(QR, symbol("x"));
    RingElem y(QR, symbol("y"));
    RingElem z1 = x;
    if (print) std::cout << "z1 = " << z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(!IsInvertible(z1));
    CoCoA_ASSERT_ALWAYS(IsZeroDivisor(z1));
    CoCoA_ASSERT_ALWAYS(!IsIntegralDomain(QR));
    CoCoA_ASSERT_ALWAYS(!IsField(QR));
    z1 = x+y+1;
    CoCoA_ASSERT_ALWAYS(!IsInvertible(z1));
    CoCoA_ASSERT_ALWAYS(!IsZeroDivisor(z1));
    ideal I(2+x);
    if (print) std::cout << "I = " << I << std::endl;
    CoCoA_ASSERT_ALWAYS(!IsOne(z1));
    if (print) std::cout << "TidyGens(I) = " << TidyGens(I) << std::endl;
    CoCoA_ASSERT_ALWAYS(len(TidyGens(I)) == 2);
    // std::cout << "IsPrime(I) = " << IsPrime(I) << std::endl;
  }


  void test_irx(bool print)
  {
    if (print) std::cout << "----- QQ[i,r,x]/(i^2+1, r^2-3) ------" << std::endl;
    PolyRing Qirx = NewPolyRing(RingQQ(), symbols("i,r, x"));
    ideal I(RingElem(Qirx,"i^2+1"), RingElem(Qirx, "r^2-3"), RingElem(Qirx, "x^2"));
    QuotientRing QR = NewQuotientRing(Qirx, I);
    RingElem i(QR, symbol("i"));
    RingElem r(QR, symbol("r"));
    RingElem x(QR, symbol("x"));
    CoCoA_ASSERT_ALWAYS(IsZeroDivisor(x));
    CoCoA_ASSERT_ALWAYS(!IsField(Qirx));
    RingElem z1 = i;
    if (print) std::cout << "z1 = " << z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(IsInvertible(z1));
    CoCoA_ASSERT_ALWAYS(!IsZeroDivisor(z1));
    CoCoA_ASSERT_ALWAYS(1/z1 == -i);
    if (print) std::cout << "1/z1 = " << 1/z1 << std::endl;
    //  CoCoA_ASSERT_ALWAYS(IsIntegralDomain(QR));
    //  CoCoA_ASSERT_ALWAYS(IsField(QR));
    z1 = i+r+1;
    if (print) std::cout << "z1 = " << z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(IsInvertible(z1));
    CoCoA_ASSERT_ALWAYS(!IsZeroDivisor(z1));
    if (print) std::cout << "13/z1 = " << 13/z1 << std::endl;
    CoCoA_ASSERT_ALWAYS(IsOne((13/z1)*(z1)/13));
  }


  void program()
  {
    GlobalManager CoCoAFoundations;
    std::cout << std::boolalpha;

    bool print = false;
    //  print = true;
    test_Qi(print);
    test_xy(print);
    test_irx(print);
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
