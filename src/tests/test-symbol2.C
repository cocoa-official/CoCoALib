//   Copyright (c)  2014  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/symbol.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;
#include <string>
using std::string;



namespace CoCoA
{

  void test_good_input(const string& str)
  {
    istringstream in(str);
    symbol x("init");
    in >> x;
    CoCoA_ASSERT_ALWAYS(in.good());
    in.peek(); // must set eofbit
    CoCoA_ASSERT_ALWAYS(in.eof());
  }

  void test_good_input_prefix(const string& str)
  {
    istringstream in(str);
    symbol x("init");
    in >> x;
    CoCoA_ASSERT_ALWAYS(in.good());
    in.peek(); // must NOT set eofbit
    CoCoA_ASSERT_ALWAYS(in.good());
  }

  void test_bad_input(const string& str)
  {
    istringstream in(str);
    bool ReadFailed = false;
    try
    {
      symbol x("init");
      in >> x;
    }
    catch (const ErrorInfo&)
    {
      ReadFailed = true;
    }
    CoCoA_ASSERT_ALWAYS(ReadFailed);
  }

  void program()
  {
    // This test does virtually nothing, but is a handy template if you want
    // to create your own test code: just copy this file and add your code
    // after the line below -- remember that you must not use CoCoALib functions
    // without first creating a GlobalManager.
    GlobalManager CoCoAFoundations;

    test_good_input("x");
    test_good_input("xyz");
    test_good_input("x[1]");
    test_good_input("xyz[1]");
    test_good_input("x[0]");
    test_good_input("xyz[0]");
    test_good_input("x[-1]");
    test_good_input("xyz[-1]");
    test_good_input("x[1,1]");
    test_good_input("xyz[1,1]");
    test_good_input("x[1,0]");
    test_good_input("xyz[1,0]");
    test_good_input("x[1,-1]");
    test_good_input("xyz[1,-1]");
    test_good_input("x[1,0,1]");
    test_good_input("xyz[1,0,1]");
    test_good_input("x[1,0,0]");
    test_good_input("xyz[1,0,0]");
    test_good_input("x[1,0,-1]");
    test_good_input("xyz[1,0,-1]");

    test_good_input(" x");
    test_good_input("xyz[ 1,0,-1]");
    test_good_input("xyz[1, 0,-1]");
    test_good_input("xyz[1,0 ,-1]");
    test_good_input("xyz[1,0, -1]");
    test_good_input("xyz[1,0,-1 ]");
    test_good_input("xyz[ 1 , 0 , -1 ]");
    test_good_input("xyz[   1   ,   0   ,   -1   ]");

    test_good_input("x1");
    test_good_input("alpha1");
    test_good_input("alpha999999999999999999999999999999");

    // the input has a good prefix, and a non-empty suffix
    test_good_input_prefix("x ");
    test_good_input_prefix("x y");
    test_good_input_prefix("x1 0");
    test_good_input_prefix("x [1]");
    test_good_input_prefix("x [");
    test_good_input_prefix("x[1] ");
    test_good_input_prefix("x@");

    // Definite bad inputs
    test_bad_input("");
    test_bad_input("@");
    test_bad_input("x[999999999999999999999999999999999]"); //overflow
    test_bad_input("x[-999999999999999999999999999999999]");//overflow
    test_bad_input("x[");
    test_bad_input("x[1");
    test_bad_input("x[1,");
    test_bad_input("x[1,0");
    test_bad_input("x[1,0,-1");
    test_bad_input("x[]");
    test_bad_input("x[,]");
    test_bad_input("x[1 0]");
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
