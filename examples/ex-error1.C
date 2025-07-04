// Copyright (c) 2007  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to throw, catch and create your own errors. \n";

const string LongDescription =
  "CoCoALib uses C++ exceptions to inform the user that a problem has occurred.    \n"
  "This simple example illustrates the CoCoALib convention for reporting an error: \n"
  "call the macro CoCoA_THROW_ERROR1 or CoCoA_THROW_ERROR2 both of which construct \n"
  "the \"ErrorInfo\" object and throw it.  It also shows how to catch these        \n"
  "CoCoALib exceptions, which are always objects of type \"ErrorInfo\".            \n"
  "The most useful part of a CoCoALib exception is its error code: after catching  \n"
  "a CoCoALib exception we can check whether its code is one we expected (and      \n"
  "can deal with).  The full list of possible error codes can be found in the      \n"
  "header file CoCoA/error.H.             \n";

//----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;

  // -------------------------------------------------------
  // Part 1

  // Here we make a CoCoALib function throw an exception.
  BigInt zero;
  BigInt one(1);
  // Since we suspect that the next operation may cause an exception to
  // be thrown, we put it inside a try...catch construct.
  try   ///// ??? a better example ???
  {
    const BigInt infinity = one/zero; // Of course, this will fail.
    cout << "THIS SHOULD NEVER BE EXECUTED" << endl;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    if (err != ERR::DivByZero)  throw; // rethrow if error was not DivByZero
    // Now we handle the expected error; in this case we'll simply print it out.
    cout << "Caught and handled this object: " << err << endl << endl;
  }

  //------------------------------------------------------------------
  // Part 2

  // We use the macro CoCoA_THROW_ERROR1 to signal a CoCoALib error: we must
  // specify the error code (see CoCoA/error.H).
  // We put the call to CoCoA_THROW_ERROR1 in a try...catch construct so that we can simply
  // discard the exception (otherwise it would seem as though an error really did occur).
  try
  {
    CoCoA_THROW_ERROR1(ERR::DivByZero);
  }
  catch (const CoCoA::ErrorInfo&) { } // Simply discard the exception

  //------------------------------------------------------------------
  // Part 3

  // If you can find no suitable error code, you can create a "nonstandard" error;
  // note that all errors created this way have the same code, viz. ERR::nonstandard.
  // Just use the macro CoCoA_THROW_ERROR with an explanatory string as the first arg.
  try
  {
    CoCoA_THROW_ERROR("Funky condition not satisfied", "program");
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    // Check we've caught the right sort of nonstandard error, err.what() is the
    // error message as a C string.
    if (err != ERR::nonstandard ||
        err.what() != string("Funky condition not satisfied"))
      throw;
    // We can also print out nonstandard errors.
    cout << "Caught an ad hoc error: " << err << endl;
  }

  // Finally, if you want to print out an ErrorInfo object in an eye-catching way,
  // you can use the function ANNOUNCE -- see its use in the procedure main below.
}

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
