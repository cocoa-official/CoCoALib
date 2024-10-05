// Copyright (c) 2013  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates use of some of the \"conversion\" functions. \n";

const string LongDescription =
  "This program illustrates use of some of the \"conversion\" functions.  \n"
  "These functions convert integer/rational values from one type into     \n"
  "another preserving the value -- if the conversion cannot be safely made\n"
  "then this is indicated (via a boolean or an error).                    \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // The function IsConvertible checks whether a value may be faithfully
    // converted into another type, and if so, it returns true and puts the
    // converted value into the first parameter.

    cout << "Examples of use of IsConvertible..." << endl;
    int n_int = 0;
    long n_long = 0;
    double d = 0;

    if (IsConvertible(n_int, power(2,50)))
      cout << "2^50 fits into an int: " << n_int << endl;
    else
      cout << "2^50 does not fit into an int" << endl;

    if (IsConvertible(n_long, power(2,50)))
      cout << "2^50 fits into a long: " << n_long << endl;
    else
      cout << "2^50 does not fit into a long" << endl;

    // WARNING!  Converting to a double may lose some precision!
    if (IsConvertible(d, power(2,50)))
      cout << "2^50 fits into a double: " << d << endl;
    else
      cout << "2^50 does not fit into a double" << endl;


    //////////////////////////////////////////////////////////////////
    // If you know that the value will surely fit then ConvertTo is easier:
    cout << endl << "Examples of use of ConvertTo..." << endl;

    BigInt N = power(2,13);
    cout << "2^13 fits into an int: " << ConvertTo<int>(N) << endl;

    cout << "2^13 fits into a long: " << ConvertTo<long>(N) << endl;

    cout << "2^13 fits into a double: " << ConvertTo<double>(N) << endl;

    //////////////////////////////////////////////////////////////////
    // If you want to use ConvertTo with a custom error message
    cout << endl << "Examples of use of ConvertTo with error message..." << endl;

    const ErrorInfo ErrMesg("Does not fit!", "Quart into pint pot");
    // or   const ErrorInfo ErrMesg(ERR::ArgTooBig, "Quart into pint pot");
    cout << "2^13 fits into an int: " << ConvertTo<int>(N, ErrMesg) << endl;
    // example of error:
    try
    {
      cout << "*** Can we convert " << power(2,100) << " into an int?" << endl;
      cout << ConvertTo<int>(power(2,100), ErrMesg) << endl;
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      cout << "*** The program correctly says it cannot with this message:";
      ANNOUNCE(cout, err);
    }
    //////////////////////////////////////////////////////////////////
    // The conversion fns can convert FROM a RingElem but not to one;
    // to create a new RingElem you must use a RingElem ctor.

    cout << endl << "We can convert from a ring elem to an integer or rational" << endl;
    RingElem a(RingZZ(), 3);
    cout << "We convert the ring elem " << a << " into an integer " << ConvertTo<long>(a) << endl;

    // HOWEVER to make a ring elem from a double you must first convert
    // the double into a BigRat:
    constexpr double pi = 3.14159265358979323846264338;
    RingElem q(RingQQ(), ConvertTo<BigRat>(pi));
    cout << "We have converted the double " << pi << " into the ring elem " << q << endl;
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-convert1.C,v 1.6 2024/03/08 19:59:58 abbott Exp $
// $Log: ex-convert1.C,v $
// Revision 1.6  2024/03/08 19:59:58  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.5  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2015/06/25 16:40:28  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.3  2014/03/19 15:44:12  abbott
// Summary: Modified so that expected output appears on cout (and not cerr)
// Author: JAA
//
// Revision 1.2  2014/02/28 16:08:58  bigatti
// -- added example of ConvertTo with custom error message
//
// Revision 1.1  2013/03/26 14:54:11  abbott
// New example for conversion fns.
//
//
