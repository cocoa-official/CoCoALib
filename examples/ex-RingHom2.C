// Copyright (c) 2005  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Operations between elements of different rings are not allowed \n"
  "but we can use homomorphisms to map the elements into the      \n"
  "same ring. \n";

const string LongDescription =
  "The example in this file shows how to create and use some      \n"
  "homomorphisms between rings.  In particular, it gives a simple \n"
  "example of mixed ring arithmetic: the user must map all values \n"
  "into a single ring before combining them arithmetically.       \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Create some coefficient rings
    ring ZZ = RingZZ();
    ring QQ = RingQQ();
    ring Fp = NewZZmod(32003);
    ring P = NewPolyRing(ZZ, symbols("alpha,beta,gamma"));   // P is ZZ[alpha,beta,gamma]

    RingElem c(ZZ), f(P), q(QQ), a(Fp);
    c = 3;
    f = 100;
    q = 5;
    a = -1;
  
    cout << "c = " << c << "   \t in " << owner(c) << endl;
    cout << "f = " << f << "   \t in " << owner(f) <<  endl;
    cout << "q = " << q << "   \t in " << owner(q) <<  endl;
    cout << "a = " << a << "   \t in " << owner(a) <<  endl;
    cout << endl;

    // THE WRONG WAY TO DO IT...
    try
    {
      cout << "c * f = " << c * f;
      // Mixed ring arithmetic triggers an exception ==> WE NEVER REACH HERE.
    }
    catch (const CoCoA::ErrorInfo& err) 
    {
      if (err != ERR::MixedRings) throw; // rethrow any unexpected exceptions
      cout << " !!!  As expected we get this error:" << endl
           << "      " << message(err) << endl
           << "The correct way is to use ring homomorphisms as follows:" << endl;
    }

    // THE RIGHT WAY TO DO IT...
    RingHom emb = CanonicalHom(ZZ, P);
    // or, equivalently,  CoeffEmbeddingHom(AsPolyRing(P));
    cout << "emb(c) * f = " << emb(c) * f << endl << endl;

    RingHom frac = CanonicalHom(ZZ, QQ);
    // or, equivalently,  EmbeddingHom(AsFractionField(QQ));
    cout << "frac(c)/4 + q = " << frac(c)/4 + q << endl << endl;

    RingHom mod = CanonicalHom(ZZ, Fp);
    // or, equivalently,  QuotientingHom(AsQuotientRing(Fp));
    cout << "mod(c) - a = " << mod(c) - a << endl << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingHom2.C,v 1.10 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingHom2.C,v $
// Revision 1.10  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.9  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.8  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.7  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.6  2015/06/26 14:54:55  abbott
// Summary: Revised after update to exceptions/errors
// Author: JAA
//
// Revision 1.5  2012/02/08 17:45:51  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/10/22 09:18:18  abbott
// Changed so that all expected is on cout -- previously one message was on clog.
//
// Revision 1.2  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 10:24:51  bigatti
// -- minor: fixed comment
//
// Revision 1.5  2007/03/07 14:31:46  bigatti
// -- now uses CanonicalHom
// -- minor cleaning
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 16:11:12  bigatti
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
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.2  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
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
