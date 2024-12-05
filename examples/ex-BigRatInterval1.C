// Copyright (c) 2018  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to create objects of type BigRatInterval,  \n"
  "and how to do basic arithmetic with them.    \n";

const string LongDescription =
  "This example shows how to create objects of type BigRatInterval,  \n"
  "and how to do basic arithmetic with them.    \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    BigRatInterval I1(BigRat(-1,2), BigRat(4,5));
    BigRatInterval I2(BigRat(2,3), BigRat(5,2));

    cout << "I1 = " << I1 << endl;
    cout << "I2 = " << I2 << endl;

    cout << "I1+I2 = " << I1+I2 << endl;
    cout << "I1-I2 = " << I1-I2 << endl;
    cout << "I1*I2 = " << I1*I2 << endl;
    cout << "I1/I2 = " << I1/I2 << endl;

    cout << "square(I1) = " << square(I1) << endl;
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
  // catch (const CoCoA::InterruptReceived& intr)
  // {
  //   cerr << endl
  //        << "------------------------------" << endl
  //        << ">>>  CoCoALib interrupted  <<<" << endl
  //        << "------------------------------" << endl
  //        << "-->>  " << intr << "  <<--" << endl;
  //   return 2;
  // }
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
