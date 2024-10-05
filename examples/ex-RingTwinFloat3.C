// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program showing how some RingTwinFloat values can have more precision   \n"
  "than that requested.                                                    \n";

const string LongDescription =
  "Example showing that certain RingTwinFloat values may have a precision  \n"
  "higher than that requested.  In this case we request 64 bits (i.e.      \n"
  "about 19 decimal digits), but in fact we can remove the nineteen most   \n"
  "significant digits and still have a result with the requested precision.\n"
  "So the original value of the variable third did in fact have at         \n"
  "least 127 bits correct (i.e. about 38 decimal digits).                   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring RR = NewRingTwinFloat(64); // Request 64 bits, about 19 decimal digits of precision.
    RingElem third = one(RR)/3;
    for (size_t i=0; i < 19; ++i)
    {
      cout << "Variable third=" << third << endl;
      cout << "Removing a leading significant digit..." << endl;
      third = 10*third-3;
    }
    cout << endl
         << "We conclude that the variable `third' initially had at least 93 bits of precision." << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingTwinFloat3.C,v 1.6 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingTwinFloat3.C,v $
// Revision 1.6  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2016/03/25 19:57:27  abbott
// Summary: Minor improvement to a printed string
//
// Revision 1.4  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.3  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 15:59:00  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/05/22 15:52:16  cocoa
// Added preprocess-disg algorithm to ApproxPts.
// Sundry minor improvements.
//
