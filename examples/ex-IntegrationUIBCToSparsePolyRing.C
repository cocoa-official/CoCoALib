// Copyright (c) 2015  John Abbott,  Anna M Bigatti
// Orig author: 2015  Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a template file for example programs.  \n"
  "The program itself does nothing whatsoever.    \n";

const string LongDescription =
  "Make a copy of this file (called \"foo.C\", say) and put your code \n"
  "inside the procedure \"program\".                                  \n"
  "To compile your file in the examples directory just do this:       \n"
  "  make foo                                                         \n";
//----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  cout << boolalpha; // so that bools print out as true/false
  ring Q = RingQQ();

  SparsePolyRing P = NewPolyRing(Q, symbols("x[0],x[1]"));
  const vector<RingElem>& x = indets(P);
  ideal I = ideal(x[0]*x[0], x[0] * x[1]);
  std::cout <<  "JanetBasis(I)    = " << Involutive::JanetBasis(I) << std::endl;
  std::cout << "IsDeltaRegular(I) = " << Involutive::IsDeltaRegular(I) << std::endl;
  std::cout << "IsMonomial(I)     = " << Involutive::IsMonomial(I) << std::endl;
  std::cout << "IsHomogeneous(I)  = " << Involutive::IsHomogeneous(I) << std::endl;
  std::map<PPMonoidElem, std::vector<bool> > multVars(Involutive::MultVars(I));

  std::cout << "MultVars" << std::endl;
  for (std::map<PPMonoidElem, std::vector<bool> >::iterator i = multVars.begin(); i != multVars.end(); ++i)
  {
    std::cout << i->first << ": " << i->second<< std::endl;
  }
  std::cout << std::endl;
  std::map<PPMonoidElem, std::vector<bool> > nonMultVars(Involutive::NonMultVars(I));

  std::cout << "NonMultVars" << std::endl;
  for (std::map<PPMonoidElem, std::vector<bool> >::iterator i = nonMultVars.begin(); i != nonMultVars.end(); ++i)
  {
    std::cout << i->first << ": " << i->second<< std::endl;
  }
  std::cout << std::endl;
}

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
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
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-IntegrationUIBCToSparsePolyRing.C,v 1.5 2022/03/04 11:29:17 bigatti Exp $
// $Log: ex-IntegrationUIBCToSparsePolyRing.C,v $
// Revision 1.5  2022/03/04 11:29:17  bigatti
// Summary: added dates for Orig author
//
// Revision 1.4  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.3  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.2  2016/05/03 13:29:35  abbott
// Summary: Major update to Mario's Janet/Pommaret code
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
