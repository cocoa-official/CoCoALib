// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows some operations on ideals. \n";

const string LongDescription =
  "This short example shows how to perform various operations  \n"
  "involving ideals.  Certain operations are valid only for    \n"
  "ideals in a polynomial ring.                                \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void IdealOperations(const ideal& I)
  {
    cout << endl << "--------------------------------------------" << endl;
    cout << "Information about " << I << endl;

    // ring and generators
    const ring& R = RingOf(I);
    cout << "The ideal is in the ring  " << R << endl;
    cout << "How many generators does I have?  " << NumGens(I) << endl;
    cout << "The generators of I are  " << gens(I) << endl;

    // Various properties
    cout << "Is the ideal zero?  " << IsZero(I) << endl;
    cout << "Is the ideal one (the whole ring)?  " << IsOne(I) << endl;
//NYI    cout << "Is the ideal prime: " << IsPrime(I) << endl;
//NYI    cout << "Is the ideal maximal: " << IsMaximal(I) << endl;

    // Testing membership/containment
    cout << "Is 3 in I?  " << IsElem(3, I) << endl;
    cout << "Is I contained in ideal(0)?  " << IsContained(I, ideal(zero(R))) << endl;

    if (!IsPolyRing(R)) return;
    // Specific operations for ideals in a polynomial ring:
    cout << endl;
    cout << "Groebner basis for I is  " << GBasis(I) << endl;
    cout << "`leading term ideal' of I is  " << LT(I) << endl;
    cout << "Is I zero-dimensional?  " << IsZeroDim(I) << endl;
    cout << "Is I homogeneous?  " << IsHomog(I) << endl;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // First we create a PolyRing in which to compute -- see also ex-PolyRing1
    PolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    // Two polynomials...
    RingElem f = RingElem(P, "x*y-1");
    RingElem g = RingElem(P, "x^2+y^2-1");

    // Now create some ideals
    ideal I1 = ideal(f);
    IdealOperations(I1);
    
    ideal I2 = ideal(f,g);
    IdealOperations(I2);

    IdealOperations(ideal(one(P)));
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-ideal2.C,v 1.4 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-ideal2.C,v $
// Revision 1.4  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.3  2017/07/08 19:07:02  abbott
// Summary: Removed comment out (dodgy) code for reporting unhandled interrupts
//
// Revision 1.2  2017/04/27 15:24:51  bigatti
// -- changed ReadExpr --> RingElem
//
// Revision 1.1  2017/02/08 16:59:21  abbott
// Summary: New examples for ideals
//
// * examples/ex-ideal1.C:
//
// * examples/ex-ideal2.C:
//
// Revision 1.11  2016/11/18 18:05:15  abbott
// Summary: Added commented out code to catch InterruptReceived
//
// Revision 1.10  2015/06/29 14:23:19  abbott
// Summary: added missing CoCoA:: prefix
// Author: JAA
//
// Revision 1.9  2015/06/29 13:25:54  bigatti
// -- code in namespace CoCoA
//
// Revision 1.8  2015/06/25 14:19:02  abbott
// Summary: Added call to CoCoA::BuildInfo::Printall
// Author: JAA
//
// Revision 1.7  2013/05/28 07:07:04  bigatti
// -- added "cout << boolalpha": useful for testing
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
