//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/BigRatOps.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <sstream>
using std::istringstream;
#include <string>
using std::string;


namespace CoCoA
{

  void TestBigRatRead(const string& InputString, const BigRat& ExpectedValue, const string& ExpectedTail)
  {
    istringstream InputStream(InputString);
    BigRat q;
    InputStream >> q;
    CoCoA_ASSERT_ALWAYS(InputStream); // check GoodBit
    if (q != ExpectedValue)
    {
      cerr << "InputString: `" << InputString << "'" << endl;
      cerr << "ExpectedValue: " << ExpectedValue << endl;
      cerr << "Read value: " << q << endl;
      CoCoA_ASSERT_ALWAYS(!"Read wrong value!");
    }

    const int TailLen = len(ExpectedTail);
    for (int i=0; i < TailLen; ++i)
    {
      CoCoA_ASSERT_ALWAYS(InputStream); // check GoodBit
      CoCoA_ASSERT_ALWAYS(InputStream.get() == ExpectedTail[i]);
    }
    // Check that InputStream is now at EOF...
    CoCoA_ASSERT_ALWAYS(InputStream);
    InputStream.get(); // should trigger EOF
    CoCoA_ASSERT_ALWAYS(InputStream.eof());
  }

  void TestBigRatReadError(const string& InputString, const string& ExpectedTail)
  {
    istringstream InputStream(InputString);
    bool ReadFailed = false;
    try
    {
      BigRat q;
      InputStream >> q;
    }
    catch (const ErrorInfo&)
    {
      ReadFailed = true;
    }
    CoCoA_ASSERT_ALWAYS(ReadFailed);

    const int TailLen = len(ExpectedTail);
    for (int i=0; i < TailLen; ++i)
    {
      CoCoA_ASSERT_ALWAYS(InputStream); // check GoodBit
      CoCoA_ASSERT_ALWAYS(InputStream.get() == ExpectedTail[i]);
    }
    // Check that InputStream is now at EOF...
    CoCoA_ASSERT_ALWAYS(InputStream);
    InputStream.get(); // should trigger EOF
    CoCoA_ASSERT_ALWAYS(InputStream.eof());
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    TestBigRatRead("1/2", BigRat(1,2), "");
    TestBigRatRead(" 1/2", BigRat(1,2), "");
    TestBigRatRead("1/2 ", BigRat(1,2), " ");
    TestBigRatRead(" 1/2 ", BigRat(1,2), " ");

    TestBigRatRead("1", BigRat(1), "");
    TestBigRatRead("1 ", BigRat(1), " ");
    TestBigRatRead(" 1", BigRat(1), "");
    TestBigRatRead(" 1 ", BigRat(1), " ");

    TestBigRatReadError("1/", "");
    TestBigRatReadError("1/ ", " ");
    TestBigRatReadError("1/ 2", " 2");
    TestBigRatRead("1 /", BigRat(1), " /");
    TestBigRatRead("1 /2", BigRat(1), " /2");
    TestBigRatRead("1/2/3", BigRat(1,2), "/3");

    TestBigRatReadError("1/x", "x");
    TestBigRatReadError(" 1/x", "x");
    TestBigRatReadError("1/x ", "x ");
    TestBigRatReadError(" 1/x ", "x ");

    TestBigRatRead("0.5", BigRat(1,2), "");
    TestBigRatRead(" 0.5", BigRat(1,2), "");
    TestBigRatRead("0.5 ", BigRat(1,2), " ");
    TestBigRatRead(" 0.5 ", BigRat(1,2), " ");
    TestBigRatReadError(".5", ".5");
    TestBigRatReadError(" .5", ".5");

    TestBigRatRead("1.", BigRat(1), "");
    TestBigRatRead("1. ", BigRat(1), " ");
    TestBigRatRead(" 1.", BigRat(1), "");
    TestBigRatRead(" 1. ", BigRat(1), " ");

    TestBigRatRead("1. 2", BigRat(1), " 2");
    TestBigRatRead("1 .", BigRat(1), " .");
    TestBigRatRead("1 .2", BigRat(1), " .2");
    TestBigRatRead("0.5.3", BigRat(1,2), ".3");

    TestBigRatRead("1./2", BigRat(1), "/2");
    TestBigRatRead("1/2.3", BigRat(1,2), ".3");
    TestBigRatRead("0.5/2", BigRat(1,2), "/2");
    TestBigRatRead("0.5.3", BigRat(1,2), ".3");
  }

}

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
