// Copyright (c) 2005,2007  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a basic example about the creation and use of the ring of integers.\n"
  "It illustrates the CoCoALib \"philosophy\" of first creating a ring, and   \n"
  "then computing with values in that ring.                                   \n"
  "The C++ commands for performing arithmetic on RingElems have a natural     \n"
  "syntax (except we cannot use ^ for powers).                                \n"
  "It warns about \"mixed ring arithmetic\", which is forbidden in CoCoALib.  \n";

const string LongDescription =
  "To calculate with elements of a ring we must first create the   \n"
  "ring, then we can create C++ objects of type RingElem which     \n"
  "belong to the ring -- a RingElem can change its value but not   \n"
  "the ring to which it belongs.                                   \n";
//---------------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring ZZ = RingZZ();
    RingElem n(ZZ);    // n is a RingElem belonging to ZZ; initial value is 0.

    // To find the ring to which a RingElem belongs use the function owner:
    cout << "The initial value of n is " << n
         << "; this is an element of " << owner(n) << "." << endl;

    //----------------------------------------------------------------------
    // We can do arithmetic on RingElems.
    n = 6;
    n = n+4;
    RingElem m = 12*n + n/2; // rhs specifies the ring to which m belongs
    cout << m << " = one hundred and twenty five" << endl;

    //----------------------------------------------------------------------
    // In a ring only exact division is allowed (contrast with BigInt values!)
    n = m/5; // OK because m is exactly divisible by 5.
    //  However n = m/7;  would throw a CoCoA::ErrorInfo exception with code ERR::BadQuot

    //----------------------------------------------------------------------
    // The next block contains a caveat.
    {
      // It is no longer possible to create a "different" copy of RingZZ:
      ring ZZ2 = RingZZ();
      if (ZZ != ZZ2)
        cerr << "THIS SHOULD NEVER BE PRINTED!" << endl;
      // in general even if R and S are canonically isomorphic rings, but not identical
      // C++ objects, they are regarded as being different.
      // As an exception, all copies of RingZZ() *are* identical C++ object.

      RingElem n2(ZZ2); // This is a RingElem belonging to ZZ2.
      // Since n and n2 do belong to the same identical ring we can
      // combine them arithmetically...
      try
      {
        n+n2; // This used to cause an exception...
      }
      catch (const CoCoA::ErrorInfo& err)
      {
        cerr << "THIS SHOULD NEVER BE PRINTED!" << endl;
        if (err != ERR::MixedRings) throw;
        // It was the exception we used to have
      }
    }
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingZZ1.C,v 1.5 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingZZ1.C,v $
// Revision 1.5  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.3  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.2  2012/03/16 15:45:34  abbott
// Realigned one line.
//
// Revision 1.1  2012/02/10 13:29:27  bigatti
// -- was ex-RingZ1.C
//
// Revision 1.7  2012/02/08 17:48:04  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.6  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.3  2008/10/07 12:08:55  abbott
// Corrected long description (grammatical error).
//
// Revision 1.2  2008/07/21 07:58:51  abbott
// Several cosmetic changes to examples.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/08 22:26:27  cocoa
// Removed try..catch constructs from some "simple" examples.
//
// Revision 1.6  2007/03/07 14:33:20  bigatti
// -- minor: text alignment
//
// Revision 1.5  2007/03/06 13:05:29  bigatti
// -- fixed comments
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 16:15:37  bigatti
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
// Revision 1.3  2004/12/09 15:08:42  cocoa
// -- added log info
//
