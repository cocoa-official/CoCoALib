// Copyright (c) 2005,2007  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Some simple computations with rational numbers.                       \n"
  "This example illustrates how to create the field of rational numbers. \n"
  "It shows that we can compute 7/3 in QQ but not in ZZ.                 \n"
  "It shows how to map an integer into a rational number.                \n";

const string LongDescription =
  "Familiarize yourself with the example ex-RingZZ1.C before proceeding. \n"
  "As C++ does not natively have any rings, we must construct them from  \n"
  "scratch.                                                              \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Creation of the rings we shall be working in.
    ring ZZ = RingZZ();                // ZZ is the ring of integers.
    ring QQ = RingQQ();                // Create field of rationals directly.
    ring K  = NewFractionField(ZZ);    // Result is identical to QQ.
    cout << "RingQQ() == NewFractionField(ZZ)  is " << (QQ == K) << endl;
  
    // Now we have created the rings we can start using their elements:
    RingElem n(ZZ);       // n is a RingElem belonging to Z; initial value is 0.
    RingElem p(QQ);       // p is a RingElem belonging to Q; initial value is 0.

    //----------------------------------------------------------------------
    // Calculate and print 7/3 as an element of Q...
    p = 7;   // p is the rational number 7
    p = p/3; // now is the rational number 7/3
    cout << "7 divided by 3 is " << p 
         << " (an element of " << QQ << ")" << endl;
    // Why did I not simply write "p = 7/3;"?
    // What value would this assign to p?  And why?
    // I cannot use 7.0/3.0 either as C++ would compute only an approximation
    // to the true value.

    //----------------------------------------------------------------------
    // We cannot compute 7/3 as an element of ZZ as the division is not exact.
    // n = 7;
    // n /= 3; // <-- this would throw a CoCoA::ErrorInfo with code ERR::BadQuot
    cout << "But in " << ZZ << " we cannot compute 7 divided by 3." << endl;

    //----------------------------------------------------------------------
    // We can map n into QQ using a homomorphism (here a canonical ringhom).
    // See ex-RingHom* for more general examples of ringhoms.
    RingHom phi = CanonicalHom(ZZ, QQ);
    RingHom psi = EmbeddingHom(QQ); // same homomorphism
    cout << "n = " << n << " in " << owner(n) << endl;
    cout << "phi(n) = " << phi(n) << " in " << owner(phi(n)) << endl;
    cout << "psi(n) = " << psi(n) << " in " << owner(psi(n)) << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingQQ1.C,v 1.6 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingQQ1.C,v $
// Revision 1.6  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2021/10/07 14:48:19  abbott
// Summary: Revised section about CanonicalHom. EmbeddingHom
//
// Revision 1.4  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.3  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.2  2012/02/10 17:19:23  abbott
// Changed RingZ into RingZZ (in a string).
//
// Revision 1.1  2012/02/10 13:28:57  bigatti
// -- was ex-RingQ1.C
//
// Revision 1.5  2012/02/08 17:48:04  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2008/07/21 08:00:15  abbott
// Delete useless comment.
//
// Revision 1.2  2008/07/21 07:58:51  abbott
// Several cosmetic changes to examples.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 22:26:27  cocoa
// Removed try..catch constructs from some "simple" examples.
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 16:15:37  bigatti
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
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
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
// Revision 1.3  2004/11/22 17:10:04  cocoa
// -- changed  "PolyRing" --> "AsPolyRing",  "FractionField" --> "AsFractionField"
//
