// Copyright (c) 2007  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows how to use the rings ZZ and QQ,          \n"
  "and how to perform various operations on ring elements (RingElem).  \n";

const string LongDescription =
  "Use of the fundamental rings ZZ and QQ.                              \n"
  "Creation of ring elements (C++ type RingElem).                       \n"
  "Operations allowed on elements of the same ring:                     \n"
  "  zero(R) and one(R)                                                 \n"
  "  a + b, a - b, a * b                                                \n"
  "  -a                                                                 \n"
  "  a = b   (assignment)                                               \n"
  "  a == b  and  a != b  (comparison)                                  \n"
  "  IsZero(a), IsOne(a), IsMinusOne(a)                                 \n"
  "Moreover other operations might be allowed, for example:             \n"
  "  a > b   if the ring is ordered                                     \n"
  "  a / b   if exact division is possible (and implemented!)           \n"
  "See ex-RingHom*.C for how to move elements from one ring to another. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // First we specify the ring(s) in which we want to compute.
    // The next two lines get the CoCoALib rings ZZ and QQ:
    ring ZZ = RingZZ();
    ring QQ = RingQQ();

    cout << "A reminder: C++ prints out boolean values like this:" << endl
         << " true   prints out as " << true << endl
         << " false  prints out as " << false << endl
         << endl;

    // Once we've created the ring(s) we shall compute in we can create
    // RingElem values belonging to those rings.
    // Whenever you create a RingElem, you must say to which ring it belongs.

    // Variables a & b belong to ZZ.
    RingElem a(ZZ), b(ZZ);
    a = 1234;
    b = one(ZZ);

    cout << " a is " << a << "  and b is " << b << endl;

    cout << " a * (b + 3)    gives  " << (a * (b + 3)  ) << endl;
    cout << " a / 2          gives  " << (a / 2        ) << endl;
    // NB a/4 is not in ZZ; trying to compute it would throw a CoCoA Error (see ex-error1)
    cout << " IsZero(a)      gives  " << (IsZero(a)    ) << endl;
    cout << " a - b == 1233  gives  " << (a - b == 1233) << endl;
    cout << endl;

    // Variables p & q belong to QQ.
    RingElem p(QQ), q(QQ);
    p = 1234;
    q = BigRat(3,2);// NOT simply 3/2, because C++ reads it as an integer division

    cout << " p is " << p << "  and q is " << q << endl;
    cout << " p / q          gives  " << (p / q        ) << endl;
    cout << " IsOne(q/q)     gives  " << (IsOne(q/q)   ) << endl;
    cout << " p > q          gives  " << (p > q        ) << endl;

    cout << " q + b          gives  " << (q + b        ) << endl;
    // The sum q+b added elements of different rings (QQ and ZZ in this case);
    // CoCoALib automatically converts b into an element of QQ, then adds.
    // In general one can use ring homomorphisms to "move" values from one
    // ring into another:  see ex-RingHom*.C
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-ring1.C,v 1.11 2022/02/13 09:57:01 abbott Exp $
// $Log: ex-ring1.C,v $
// Revision 1.11  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.10  2021/10/07 14:39:44  abbott
// Summary: Updated now that CoCoALib "coerce" ringelems
//
// Revision 1.9  2015/06/29 15:52:45  bigatti
// *** empty log message ***
//
// Revision 1.8  2015/06/29 13:17:29  bigatti
// -- code inside namespace CoCoA
//
// Revision 1.7  2012/03/30 09:27:47  bigatti
// -- using BigRat instead of RingElem(QQ ..)
//
// Revision 1.6  2012/02/08 17:52:17  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.5  2011/10/12 15:52:00  abbott
// Changed names of local vars to ZZ & QQ.
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2008/07/21 07:58:51  abbott
// Several cosmetic changes to examples.
//
// Revision 1.2  2007/03/22 22:44:42  abbott
// Removed spaces at ends of lines.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/07 18:09:19  bigatti
// -- first import
//
