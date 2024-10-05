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
  "call the macro CoCoA_THROW_ERROR which both constructs the \"ErrorInfo\" value  \n"
  "and throws it.  It also shows how to catch these CoCoALib exceptions, which are \n"
  "always objects of type \"ErrorInfo\".  The most useful part of a CoCoALib       \n"
  "exception is its error code: after catching a CoCoALib exception we can check   \n"
  "whether its code is one we expected (and can deal with).  The full list of      \n"
  "possible error codes can be found in the header file CoCoA/error.H.             \n";

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
    BigInt infinity = one/zero; // Of course, this will fail.
    cout << "THIS SHOULD NEVER BE EXECUTED" << endl;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    if (err != ERR::DivByZero) throw; // rethrow if error was not DivByZero
    // Now we handle the expected error; in this case we'll simply print it out.
    cout << "Caught and handled this object: " << err << endl << endl;
  }

  //------------------------------------------------------------------
  // Part 2

  // We use the macro CoCoA_THROW_ERROR to signal a CoCoALib error: we must
  // specify the error code (see CoCoA/error.H) and the function reporting
  // the problem.
  // We put the call to CoCoA_THROW_ERROR in a try...catch construct so that we can simply
  // discard the exception (otherwise it would seem as though an error really did occur).
  try
  {
    CoCoA_THROW_ERROR(ERR::DivByZero, "program");
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-error1.C,v 1.10 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-error1.C,v $
// Revision 1.10  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.9  2020/06/17 18:59:52  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.8  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.7  2015/06/29 15:52:45  bigatti
// *** empty log message ***
//
// Revision 1.6  2012/09/21 14:03:17  bigatti
// -- minor change
//
// Revision 1.5  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.2  2008/07/21 07:58:51  abbott
// Several cosmetic changes to examples.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.1  2007/03/07 18:06:24  cocoa
// Example about dealing with errors/exceptions in CoCoALib.
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
