// Copyright (c) 2015  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example illustrates basic use of SmallFpImpl for arithmetic in a\n"
  "small prime finite field.  ex-SmallFp3 shows how to compute more     \n"
  "efficiently, but it is not for the faint-hearted!                    \n";

const string LongDescription =
  "This example illustrates basic use of SmallFpImpl for arithmetic in a\n"
  "small prime finite field.  ex-SmallFp3 shows how to compute more     \n"
  "efficiently, but it is not for the faint-hearted!                    \n"
  "Here we see how to create a SmallFpImpl object, and its use for basic\n"
  "arithmetic (add, sub, mul, div, power) on finite field elements.     \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    const int p = NextPrime(4096);
    SmallFpImpl ModP(p);
    SmallFpImpl::value a,b,c; // initially 0
    a = ModP.myReduce(99);    // a = 99; mod p
    b = ModP.myReduce(28);    // b = 28; mod p
    c = ModP.myReduce(-3);    // c = -3; mod p
    a = ModP.myAdd(b,c);      // a = b+c;  also mySub, myMul, myDiv
    a = ModP.myNegate(a);     // a = -a;
    a = ModP.myRecip(a);      // a = 1/a;
    // Now verify Fermat's Little theorem
    b = ModP.myPower(a,p-1);  // b = a^(p-1);  where ^ means "power"
    if (!IsZero(a) && !IsOne(b)) cerr << "Fermat's Little Theorem failed!" << endl;

    a = ModP.myAddMul(a,b,c); // a = a+b*c;

    const long B = ModP.myExport(b); // deliver value of b as a long (see ex-SmallFp2)
    if (B != 1) cerr << "Exported value " << B << " should have been 1." << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-SmallFp1.C,v 1.7 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-SmallFp1.C,v $
// Revision 1.7  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2015/11/04 10:11:30  abbott
// Summary: Changed to a simpler example (old example is now ex-SmallFp3)
//
// Revision 1.5  2015/06/29 15:49:01  bigatti
// *** empty log message ***
//
// Revision 1.4  2015/06/29 13:17:29  bigatti
// -- code inside namespace CoCoA
//
// Revision 1.3  2014/04/03 15:43:07  abbott
// Summary: Reduced size of prime p (o/w too slow with debugging on some machines)
// Author: JAA
//
// Revision 1.2  2013/05/27 14:48:18  abbott
// Added typedef for FpElem to make code more readable.
//
// Revision 1.1  2013/05/27 12:55:04  abbott
// Added new example ex-SmallFp1.C
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
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
