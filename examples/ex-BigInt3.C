// Copyright (c) 2004  John Abbott,  Anna M. Bigatti
// Orig authors: 2004 Daniele Venzano; modified by John Abbott.

#include "CoCoA/library.H"

using namespace std;

//-----------------------------------------------------------------------------
const string ShortDescription =
  "Program to find a (probable) prime with a specified number of bits.\n";


const string LongDescription =
  "This program shows that BigInts can be used much like normal C++ ints\n"
  "with the advantage that there is almost no limit on the magnitude of \n"
  "the values.  Here we generate random BigInts and test them for being \n"
  "probable primes -- stopping as soon as we find a likely prime.       \n"
  "NB If you need extreme efficiency then use the GMP library directly. \n";
//-----------------------------------------------------------------------------


namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Enter the number of bits the prime should have: ";
    long nbits;
    cin >> nbits;
    if (!cin) CoCoA_THROW_ERROR(ERR::InputFail, "ex-BigInt3: input was not a (small) integer");
    if (nbits < 2)
    {
      cout << "There are no primes with fewer than 2 bits." << endl;
      return;
    }

    cout << endl << "Starting search for a prime with " << nbits << " bits..." << endl;
    if (nbits > 9999) cout << "WARNING: the computation will take a very long time!" << endl;
    else if (nbits > 999) cout << "NOTE: this computation may take some time." << endl;

    const BigInt min = power(2, nbits-1);
    const BigInt max = 2*min-1;

    BigInt candidate = RandomBigInt(min, max);
    while (!IsProbPrime(candidate))
    {
      candidate = RandomBigInt(min, max);
    }

    cout << candidate << " seems to be a prime! (" << nbits << " bits)" << endl;

    cout << "And further test with more iterations ";
    if (IsProbPrime(candidate, 50))
      cout << "confirms";
    else
      cout << "rejects";
    cout << " this result." << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-BigInt3.C,v 1.6 2022/03/04 11:29:17 bigatti Exp $
// $Log: ex-BigInt3.C,v $
// Revision 1.6  2022/03/04 11:29:17  bigatti
// Summary: added dates for Orig author
//
// Revision 1.5  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.3  2015/06/25 16:40:27  abbott
// Summary: Moved code into namespace CoCoA (redmine 739)
// Author: JAA
//
// Revision 1.2  2013/10/30 09:54:35  bigatti
// -- lower case in (c)
//
// Revision 1.1  2012/03/30 10:36:47  abbott
// Renamed ex-BigIntPrime1 to ex-BigInt3.
//
// Revision 1.2  2011/12/23 14:55:01  bigatti
// -- added comment "defined in .."
//
// Revision 1.1  2011/08/26 10:19:48  bigatti
// -- renamed after ZZ->BigInt, QQ->BigRat
//
// Revision 1.6  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.5  2010/12/26 13:02:03  abbott
// Finished changing "GlobalXXXput()" into the corresponding standard C++ stream.
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
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
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/22 17:25:04  bigatti
// -- added printout of number of bits used
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
// Revision 1.1  2005/03/02 18:46:41  cocoa
// Added new types ConstRefMatrix, and RefMatrix following along
// the lines of ConstRefRingElem and RefRingElem.  The semantics
// should be a bit clearer now.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.1  2004/11/29 16:50:28  cocoa
// -- first import
//
