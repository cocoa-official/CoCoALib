// Copyright (c) 2012 John Abbott, Anna Bigatti

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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/PPMonoidEvOv.H"
#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/PPMonoidOv.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"
//#include "CoCoA/DenseMatrix.H"
//#include "CoCoA/matrix.H"

#include<iostream>
using std::cout;
using std::cerr;
using std::endl;
#include<vector>
using std::vector;

//----------------------------------------------------------------------
// First basic test for GeneralHom and RestrictionHom
//----------------------------------------------------------------------
namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    PPMonoid PPM1 = NewPPMonoid(symbols("x,y,z"), StdDegRevLex);
    PPMonoid PPM2 = NewPPMonoid(symbols("alpha,beta,gamma,delta"), lex);

    const vector<PPMonoidElem>& x = indets(PPM1);
    const vector<PPMonoidElem>& alpha = indets(PPM2);

    vector<PPMonoidElem> images;
    images.push_back(alpha[0]*alpha[1]);
    images.push_back(alpha[1]*alpha[2]);
    images.push_back(alpha[2]*power(alpha[3],4));

    PPMonoidHom phi = GeneralHom(PPM1, images);
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
        {
          PPMonoidElem t = power(x[0],i) * power(x[1],j) * power(x[2],k);
          CoCoA_ASSERT_ALWAYS(phi(t) == power(images[0],i) * power(images[1],j) * power(images[2],k));
        }

    PPMonoidHom psi = RestrictionHom(PPM1, vector<long>(1,1));  
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
        {
          PPMonoidElem t = power(x[0],i) * power(x[1],j) * power(x[2],k);
          CoCoA_ASSERT_ALWAYS(psi(t) == power(x[1],j));
        }
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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
