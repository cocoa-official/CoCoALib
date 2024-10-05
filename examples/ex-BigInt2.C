// Copyright (c) 2004-2007,2009  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program to calculate the modular inverse of a number.  \n"
  "You must give as input the modulus and the residue;    \n"
  "the residue need not be reduced.                       \n";

const string LongDescription =
  "This program illustrates that BigInt values can be used much like normal C++ \n"
  "ints except that there is a almost no limit on the magnitude of the values.  \n"
  "NB If you need extreme efficiency then use the GMP library directly.         \n"
  "Contrast this example with ex-RingZZ1.                                       \n";
//-----------------------------------------------------------------------------


namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Input the modulus: ";
    BigInt modulus;
    cin >> modulus;
    if (!cin || modulus < 2) { cerr << "The modulus must be an integer >= 2\n"; return; }

    cout << "Input the residue: ";
    BigInt residue;
    cin >> residue;
    if (!cin) { cerr << "The residue must be an integer.\n"; return; }

    if (gcd(residue, modulus) != 1)
    {
      cout << "No inverse exists because the residue and modulus have a" << endl
           << "non-trivial common factor: " << gcd(residue, modulus) << endl;
      return;
    }

    const BigInt inverse = InvMod(residue, modulus); // defined in NumTheory
    cout << "The inverse of " << residue << " modulo " << modulus
         << " is " << inverse << endl;

    // We confirm that the result is right.
    if ((inverse*residue)%modulus != 1)
      cout << "THIS SHOULD NEVER BE PRINTED!" << endl;
  }

} // end of namespace CoCoA


// We write main() like this so we can handle uncaught CoCoA errors in
// a sensible way (i.e. by announcing them).
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error" << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-BigInt2.C,v 1.5 2022/02/13 09:56:56 abbott Exp $
// $Log: ex-BigInt2.C,v $
// Revision 1.5  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.3  2012/02/10 17:20:25  abbott
// Changed RingZ into RingZZ.
//
// Revision 1.2  2012/01/12 11:02:57  bigatti
// -- added comment
//
// Revision 1.1  2011/08/26 10:19:48  bigatti
// -- renamed after ZZ->BigInt, QQ->BigRat
//
// Revision 1.3  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2010/03/22 10:15:56  abbott
// Renamed previous ex-ZZ1 to ex-ZZ2, to make way for a simpler ex-ZZ1.
//
// Revision 1.5  2009/12/09 13:23:52  abbott
// Added better checks on the inputs.
//
// Revision 1.4  2009/07/08 12:27:15  abbott
// Minor cleaning.
//
// Revision 1.3  2009/06/05 12:14:56  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInteger.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
// Revision 1.2  2007/03/22 22:45:31  abbott
// Removed spaces at ends of lines.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 16:20:57  bigatti
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
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/12/09 15:08:42  cocoa
// -- added log info
//
