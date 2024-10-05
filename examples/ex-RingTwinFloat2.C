// Copyright (c) 2006  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program exhibiting a way of using ever higher precisions...            \n"
  "This example shows how failure can be a success: this pathological     \n"
  "computation produces the **same wrong result** when using normal       \n"
  "floating point arithmetic at any given finite precision!  However,     \n"
  "since twin floats are self-checking, we detect that there is a problem.\n";

const string LongDescription =
  "Example showing iterative increase of precision using RingTwinFloat    \n"
  "until the answer is found (or some maximum precision is reached).      \n"
  "This program will always fail to find the limit: J-M Muller's sequence \n"
  "actually converges to 6 (rather slowly), however it is pathological    \n"
  "because it converges to 100 using any finite precision arithmetic.     \n"
  "RingTwinFloat detects the onset of pathological convergence to 100, and\n"
  "throws an InsufficientPrecision exception.                             \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // This function computes J-M Muller's sequence in the ring R until two
  // successive values are "practically equal" (see documentation).
  // It returns the value the sequence converged to.
  // Note: we use the initial values 2 and -4 which can be represented exactly;
  // the original sequence had rational initial values.
  RingElem MullerSeq(const RingTwinFloat& RR)
  {
    RingElem Vprev2(RR, 2);  // These starting values can be represented exactly.
    RingElem Vprev1(RR, -4);

    RingElem Vcurr(RR);
    while (!IsPracticallyEqual(Vprev1, Vprev2))
    {
      Vcurr = 111 - 1130/Vprev1 + 3000/(Vprev1*Vprev2); // Muller's recurrence.
      Vprev2 = Vprev1;
      Vprev1 = Vcurr;
    }
    return Vprev1;
  }


  // This loop tries ever higher precisions until the answer is found.
  // Arbitrarily limit precision to 4096 so it cannot run forever.
  RingElem MullerLoop()
  {
    for (int BitPrec = 32; BitPrec <= 4096; BitPrec *= 2)
    {
      try { return MullerSeq(NewRingTwinFloat(BitPrec)); }
      catch (const RingTwinFloat::InsufficientPrecision&)
      {
        // Inform user about the failure...
        cout << "A bit precision of " << BitPrec << " was not sufficient." << endl;
      }
    }
    CoCoA_THROW_ERROR("Required too much precision... giving up!", "MullerLoop");
    return zero(RingZZ()); // Never executed; just to keep compiler quiet.
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    try
    {
      const RingElem limit = MullerLoop();
      cout << "Muller's sequence converges to " << limit << endl;
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      if (message(err) != "Required too much precision... giving up!") throw;
      cout << endl
           << "As expected we caught this error: " << message(err) << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingTwinFloat2.C,v 1.17 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RingTwinFloat2.C,v $
// Revision 1.17  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.16  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.15  2016/07/20 08:42:45  abbott
// Summary: Improved a comment
//
// Revision 1.14  2016/03/25 19:57:00  abbott
// Summary: Added a blank line
//
// Revision 1.13  2015/11/21 19:15:39  abbott
// Summary: Added missing const
//
// Revision 1.12  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.11  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.10  2015/06/26 14:54:55  abbott
// Summary: Revised after update to exceptions/errors
// Author: JAA
//
// Revision 1.9  2014/08/16 13:21:22  abbott
// Summary: Several minor improvements
// Author: JAA
//
// Revision 1.8  2012/02/08 17:48:04  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.7  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.6  2010/10/22 09:19:03  abbott
// Changed so that all expected output is on cout (previously some messages were on clog).
//
// Revision 1.5  2009/01/23 16:11:44  abbott
// Improved string ShortDescription.
//
// Revision 1.4  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.3  2008/10/07 12:08:01  abbott
// Cleaned the long and short descriptions.  Removed usless #include.
// Shortened main loop (so that program runs a bit faster).
//
// Revision 1.2  2007/09/24 14:11:58  abbott
// Changed initial values (to integers), and reduced max precision to make program faster.
//
// Revision 1.1  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/06 13:04:06  bigatti
// -- cleaned: it catches the error it throws and performs 1 fewer operation
//
// Revision 1.6  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.5  2007/02/22 17:24:32  bigatti
// -- added a comment in ShortDescription
//
// Revision 1.4  2007/02/12 15:59:00  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.3  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.2  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//

