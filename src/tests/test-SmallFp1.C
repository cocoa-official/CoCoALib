//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


namespace CoCoA
{

  void TestFermatLittleThm(const SmallFpImpl& Fp)
  {
    const long p = Fp.myModulus();
    for (long A=1; A < p; ++A)
    {
      SmallFpImpl::value a = Fp.myReduce(A);
      CoCoA_ASSERT_ALWAYS(IsOne(Fp.myPower(a, p-1)));
    }
  }

  void TestAddition(const SmallFpImpl& Fp)
  {
    const long p = Fp.myModulus();
    for (long A=0; A < p ; ++A)
      for (long B=0; B < p; ++B)
      {
        const SmallFpImpl::value a = Fp.myReduce(A);
        const SmallFpImpl::value b = Fp.myReduce(B);
        CoCoA_ASSERT_ALWAYS(Fp.myAdd(a,b) == Fp.myAdd(b,a));
        CoCoA_ASSERT_ALWAYS(Fp.myAdd(a,b) == Fp.myReduce(A+B));
      }
  }

  void TestNegation(const SmallFpImpl& Fp)
  {
    const long p = Fp.myModulus();
    for (long A=0; A < p ; ++A)
    {
      const SmallFpImpl::value a = Fp.myReduce(A);
      const SmallFpImpl::value b = Fp.myNegate(a);
      CoCoA_ASSERT_ALWAYS(IsZero(Fp.myAdd(a,b)));
      CoCoA_ASSERT_ALWAYS(b == Fp.myReduce(p-A));
    }
  }
  
  void TestMultiplication(const SmallFpImpl& Fp)
  {
    const long p = Fp.myModulus();
    for (long A=0; A < p ; ++A)
      for (long B=0; B < p; ++B)
      {
        const SmallFpImpl::value a = Fp.myReduce(A);
        const SmallFpImpl::value b = Fp.myReduce(B);
        CoCoA_ASSERT_ALWAYS(Fp.myMul(a,b) == Fp.myMul(b,a));
        CoCoA_ASSERT_ALWAYS(Fp.myMul(a,b) == Fp.myReduce(A*B));
      }
  }

  void TestReciprocal(const SmallFpImpl& Fp)
  {
    const long p = Fp.myModulus();
    for (long A=1; A < p ; ++A)
    {
      const SmallFpImpl::value a = Fp.myReduce(A);
      const SmallFpImpl::value b = Fp.myRecip(a);
      CoCoA_ASSERT_ALWAYS(IsOne(Fp.myMul(a,b)));
      CoCoA_ASSERT_ALWAYS(a == Fp.myRecip(b));
    }
  }


  void TestExport(const SmallFpImpl& Fp)
  {
    const long p = Fp.myModulus();
    const long lwb = -((p-1)/2); // integer division
    for (long A=0; A < p ; ++A)
    {
      const SmallFpImpl::value a = Fp.myReduce(A);
      const long x = Fp.myExport(a);
      CoCoA_ASSERT_ALWAYS(lwb <= x && x < p);
      CoCoA_ASSERT_ALWAYS(x < 0 || x == A);
    }
  }

    

  void program()
  {
    GlobalManager CoCoAFoundations;

    vector<long> plist;
    plist.push_back(2);
    plist.push_back(3);
    plist.push_back(101);
    plist.push_back(1009);
    for (int i=0; i < len(plist); ++i)
    {
      const long p = plist[i];
      SmallFpImpl Fp(p);
      TestFermatLittleThm(Fp);
      TestAddition(Fp);
      TestNegation(Fp);
      TestMultiplication(Fp);
      TestReciprocal(Fp);
      TestExport(Fp);
    }
  }

}  // end of namespace CoCoA

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
