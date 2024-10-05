// Copyright (c) 2003-2006  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "An example about RingWeyl, the interface is not quite settled yet.\n";

const string LongDescription =
  "This shows a computation of a Groebner Basis.\n"
  "All these examples about RingWeyl will probably be merged into one.\n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void test(ideal I)
  {
    cout << "gens(I) = " << gens(I) << endl;
    cout << "TidyGens(I) = " << TidyGens(I) << endl;
    cout << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    vector<symbol> names = symbols("u,v,x,y,t");

    vector<long> ElimIndets;
    ElimIndets.push_back(0); // elim u
    ElimIndets.push_back(1); // elim v
    SparsePolyRing WA = NewWeylAlgebra(RingQQ(), names, ElimIndets);

    RingElem u(WA, symbol("u"));
    RingElem v(WA, symbol("v"));
    RingElem x(WA, symbol("x"));
    RingElem y(WA, symbol("y"));
    RingElem t(WA, symbol("t"));
    RingElem dx(WA, symbol("dx"));
    RingElem dy(WA, symbol("dy"));
    RingElem dt(WA, symbol("dt"));

    test(ideal(1-u*v,  t*u-x,  v*dt+dx));
    test(ideal(1-u*v,  t*u -x*x -y*y*y,  2*x*v*dt +dx, 3*y*y*v*dt +dy));
    test(ideal(1-u*v,  t*u -power(x,3) +power(y,4), 3*x*x*v*dt +dx,  -4*y*y*y*v*dt +dy));
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

//output
// gens(I) = [-x[0]*x[1] +1,  x[0]*x[4] -x[2],  x[1]*d[4] +d[2]]
// TidyGens(I) = [x[1]*d[4] +d[2],  x[0]*x[4] -x[2],  x[0]*x[1] -1,  x[0]*d[2] +d[4],  x[2]*d[2] +x[4]*d[4] +1,  x[1]*x[2] -x[4]] 

//output
//gens(I)=[-x[0]*x[1] +(1),x[0]*x[4] +-x[3]^3 +-x[2]^2,(2)*x[1]*x[2]*d[4] +d[2],(3)*x[1]*x[3]^2*d[4] +d[3]]

/* 
 TidyGens := [
 (2)*x[1]*x[2]*d[4] +d[2],
 (3)*x[1]*x[3]^2*d[4] +d[3],
 x[3]^2*d[2] +(-2)/(3) *x[2]*d[3],
 x[0]*x[4] +-x[3]^3 +-x[2]^2,
 (-2)*x[0]*x[1]*x[2] +(2)*x[2],
 x[0]*d[2] +(-2)/(-1) *x[2]*d[4],
 -x[2]^2*d[2] +(-2)/(3) *x[2]*x[3]*d[3] +(2)/(-1) *x[2]*x[4]*d[4] +(-2)*x[2],
 (2)/(-3) *x[2]*x[3]^3*d[3] +(-2)*x[2]*x[3]^2*x[4]*d[4] +(-2)/(3) *x[2]^3*d[3] 
 +(2)/(-1) *x[2]*x[3]^2,
 -x[3]^3*d[3] +(-3)*x[3]^2*x[4]*d[4] +-x[2]^2*d[3] +(3)/(-1) *x[3]^2,
 x[2]*d[2]^2 +(-2)/(-3) *x[3]*d[2]*d[3] +(2)*x[4]*d[2]*d[4] +(3)*d[2],
 -x[2]*d[2] +(-2)/(3) *x[3]*d[3] +(2)/(-1) *x[4]*d[4] +(-2),
 (4)/(-3) *x[1]*x[3]*d[3]*d[4] +(-4)*x[1]*x[4]*d[4]^2 +(-6)*x[1]*d[4] +d[2]^2,
 (-9)*x[1]*x[3]*x[4]*d[4]^2 +(-15)/(2) *x[1]*x[3]*d[4] +(9)/(4) *x[3]*d[2]^2 +d[3]^2,
 (-4)*x[1]*x[4]^2*d[4]^3 +(-12)*x[1]*x[4]*d[4]^2 +(-35)/(9) *x[1]*d[4] +(-1)/(3)
 *x[3]*d[2]^2*d[3] +x[4]*d[2]^2*d[4] +(-4)/(27) *d[3]^3 +(1)/(2) *d[2]^2,
 (2)/(-3) *x[0]*x[3]*d[3] +(-2)*x[3]^3*d[4],
 (-2)/(3) *x[0]*x[2]*d[3] +(2)/(-1) *x[2]*x[3]^2*d[4],
 -x[0]*d[3] +(-3)*x[3]^2*d[4],
 (-2)*x[0]*x[1] +(2),
 -x[1]*x[3]^3 +-x[1]*x[2]^2 +x[4]];


-- ????? -- current code returns this... not quite the same
TidyGens(I) = [x[0]*x[1] -1,  x[0]*x[4] -x[3]^3 -x[2]^2,  x[1]*x[3]^3 +x[1]*x[2]^2 -x[4],  x[1]*x[2]*d[4] +1/2*d[2],  x[0]*d[2] +2*x[2]*d[4],  x[3]^3*d[2] -2/3*x[2]*x[3]*d[3],  x[1]*x[3]^2*d[4] +1/3*d[3],  x[2]*d[2] +2/3*x[3]*d[3] +2*x[4]*d[4] +2,  x[0]*x[3]*d[3] +3*x[3]^3*d[4],  x[1]*x[3]*d[3]*d[4] +3*x[1]*x[4]*d[4]^2 +9/2*x[1]*d[4] -3/4*d[2]^2,  x[3]^4*d[3] +3*x[3]^3*x[4]*d[4] +x[2]^2*x[3]*d[3] +3*x[3]^3,  x[1]*x[3]*x[4]*d[4]^2 +5/6*x[1]*x[3]*d[4] -1/4*x[3]*d[2]^2 -1/9*d[3]^2,  x[3]^2*d[2] -2/3*x[2]*d[3],  x[3]^3*d[3] +3*x[3]^2*x[4]*d[4] +x[2]^2*d[3] +3*x[3]^2,  x[0]*x[2]*d[3] +3*x[2]*x[3]^2*d[4],  x[0]*d[3] +3*x[3]^2*d[4],  x[1]*x[4]^2*d[4]^3 +3*x[1]*x[4]*d[4]^2 +35/36*x[1]*d[4] +1/12*x[3]*d[2]^2*d[3] -1/4*x[4]*d[2]^2*d[4] +1/27*d[3]^3 -1/8*d[2]^2]
*/

//output
//gens(I)=[
//-x[0]*x[1] +(1),
//x[0]*x[4] +x[3]^4 +-x[2]^3,
//(3)*x[1]*x[2]^2*d[4] +d[2],
//(-4)*x[1]*x[3]^3*d[4] +d[3]]
 
//TidyGens(I)=[
// (3)*x[1]*x[2]^2*d[4] +d[2],
// (-4)*x[1]*x[3]^3*d[4] +d[3],
// x[3]^3*d[2] +(-3)/(-4) *x[2]^2*d[3],
// x[0]*x[4] +x[3]^4 +-x[2]^3,
// (-3)*x[0]*x[1]*x[2]^2 +(3)*x[2]^2,
// x[0]*d[2] +(-3)/(-1) *x[2]^2*d[4],
// -x[2]^3*d[2] +(3)/(-4) *x[2]^2*x[3]*d[3] +(3)/(-1) *x[2]^2*x[4]*d[4] +(-3)*x[2]^2,
// (-3)/(4) *x[2]^2*x[3]^4*d[3] +(-3)*x[2]^2*x[3]^3*x[4]*d[4] +(-3)/(-4) *x[2]^5*d[3] 
// +(3)/(-1) *x[2]^2*x[3]^3,
// (2)/(-1) *x[2]*x[3]^4*d[3] +(8)/(-1) *x[2]*x[3]^3*x[4]*d[4] +(-2)/(-1) *x[2]^4*d[3] 
// +(-8)*x[2]*x[3]^3, 
// -x[3]^4*d[3] +(4)/(-1) *x[3]^3*x[4]*d[4] +x[2]^3*d[3] +(-4)*x[3]^3,
// x[2]*d[2]^2 +(3)/(4) *x[3]*d[2]*d[3] +(3)*x[4]*d[2]*d[4] +(4)*d[2],
// (-2)*x[2]^2*d[2] +(3)/(-2) *x[2]*x[3]*d[3] +(6)/(-1) *x[2]*x[4]*d[4] 
// +(-6)*x[2],-x[2]*d[2] +(-3)/(4) *x[3]*d[3] +(-3)*x[4]*d[4] +(3)/(-1) ,
// (-9)/(4) *x[1]*x[2]*x[3]*d[3]*d[4] +(-9)*x[1]*x[2]*x[4]*d[4]^2 
// +(-12)*x[1]*x[2]*d[4] +d[2]^2,
// (-3)/(4) *x[1]*x[3]^2*d[3]^2*d[4] +(-6)*x[1]*x[3]*x[4]*d[3]*d[4]^2 
// +(-12)*x[1]*x[4]^2*d[4]^3 +(39)/(-4) *x[1]*x[3]*d[3]*d[4] 
// +(48)/(-1) *x[1]*x[4]*d[4]^2 +(80)/(-3) *x[1]*d[4] +(-4)/(9) *d[2]^3,
// (-32)/(-1) *x[1]*x[3]^2*x[4]*d[3]*d[4]^2 +(-64)/(-1) *x[1]*x[3]*x[4]^2*d[4]^3 
// +(28)*x[1]*x[3]^2*d[3]*d[4] +(256)*x[1]*x[3]*x[4]*d[4]^2 +(1064)/(9) *x[1]*x[3]*d[4] 
// +(-64)/(-27) *x[3]*d[2]^3 +d[3]^3,
// (-16)/(-1) *x[1]*x[2]*x[3]^2*x[4]*d[4]^2 +(-28)/(-3) *x[1]*x[2]*x[3]^2*d[4] 
// +(16)/(-9) *x[3]^2*d[2]^2 +x[2]*d[3]^2,
// (-9)/(2) *x[1]*x[3]*x[4]^2*d[3]*d[4]^3 +(-12)*x[1]*x[4]^3*d[4]^4 
// +(27)/(-2) *x[1]*x[3]*x[4]*d[3]*d[4]^2 +(81)/(-1) *x[1]*x[4]^2*d[4]^3 
// +(-427)/(96) *x[1]*x[3]*d[3]*d[4] +(-332)/(3) *x[1]*x[4]*d[4]^2 
// +(-329)/(16) *x[1]*d[4] +(1)/(18) *x[3]*d[2]^3*d[3] +(-4)/(9) *x[4]*d[2]^3*d[4] 
// +(-3)/(-128) *d[3]^4 +(-1)/(3) *d[2]^3,
// (8)*x[1]*x[3]^2*x[4]^2*d[4]^3 +(-20)/(-1) *x[1]*x[3]^2*x[4]*d[4]^2 
// +(-77)/(-18)*x[1]*x[3]^2*d[4] +(8)/(27) *x[3]^2*d[2]^3 +(-1)/(-8) *x[3]*d[3]^3 
// +x[4]*d[3]^2*d[4] +(7)/(8) *d[3]^2,
// (-9)/(2) *x[1]*x[2]*x[3]*x[4]^2*d[4]^3 +(87)/(-8) *x[1]*x[2]*x[3]*x[4]*d[4]^2 
// +(35)/(-16) *x[1]*x[2]*x[3]*d[4] +(-1)/(8) *x[3]^2*d[2]^2*d[3] 
// +(1)/(2) *x[3]*x[4]*d[2]^2*d[4] +(-9)/(-128) *x[2]*d[3]^3 +(1)/(24) *x[3]*d[2]^2,
// (64)/(-3) *x[1]*x[3]*x[4]^3*d[4]^4 +(-120)*x[1]*x[3]*x[4]^2*d[4]^3 
// +(7)/(-54) *x[1]*x[3]^2*d[3]*d[4] +(3368)/(-27) *x[1]*x[3]*x[4]*d[4]^2 
// +(119)/(-9) *x[1]*x[3]*d[4] +(32)/(81) *x[3]^2*d[2]^3*d[3] 
// +(64)/(-81) *x[3]*x[4]*d[2]^3*d[4] +(-1)/(-6) *x[3]*d[3]^4 +x[4]*d[3]^3*d[4] 
// +(8)/(27) *x[3]*d[2]^3 +(-9)/(-8) *d[3]^3,
// (3)/(-1) *x[1]*x[2]*x[4]^3*d[4]^4 +(33)/(-2) *x[1]*x[2]*x[4]^2*d[4]^3 
// +(-265)/(16) *x[1]*x[2]*x[4]*d[4]^2 +(455)/(-288) *x[1]*x[2]*d[4] 
// +(-1)/(-48) *x[3]^2*d[2]^2*d[3]^2 +(-1)/(12) *x[3]*x[4]*d[2]^2*d[3]*d[4] 
// +(-1)/(-3) *x[4]^2*d[2]^2*d[4]^2 +(-3)/(256) *x[2]*d[3]^4 +(5)/(144) *x[3]*d[2]^2*d[3] 
// +(-13)/(-18) *x[4]*d[2]^2*d[4] +(-67)/(-432) *d[2]^2,
// (-12)*x[1]*x[4]^4*d[4]^5 +(-120)*x[1]*x[4]^3*d[4]^4 
// +(7)/(-24) *x[1]*x[3]*x[4]*d[3]*d[4]^2 +(-3581)/(12) *x[1]*x[4]^2*d[4]^3 
// +(7)/(-24) *x[1]*x[3]*d[3]*d[4] +(-1055)/(6) *x[1]*x[4]*d[4]^2 
// +(-1589)/(144) *x[1]*d[4] +(-1)/(12) *x[3]^2*d[2]^3*d[3]^2 
// +(-2)/(-9) *x[3]*x[4]*d[2]^3*d[3]*d[4] +(-4)/(9) *x[4]^2*d[2]^3*d[4]^2 
// +(9)/(-256) *x[3]*d[3]^5 +(3)/(-16) *x[4]*d[3]^4*d[4] +(-7)/(36) *x[3]*d[2]^3*d[3] 
// +(8)/(-9) *x[4]*d[2]^3*d[4] +(33)/(-128) *d[3]^4 +(-31)/(108) *d[2]^3,
// (-3)/(4) *x[0]*x[3]*d[3] +(3)*x[3]^4*d[4],
// (-3)/(-4) *x[0]*x[2]^2*d[3] +(3)/(-1) *x[2]^2*x[3]^3*d[4],
// (-2)*x[0]*x[2]*d[3] +(8)*x[2]*x[3]^3*d[4],
// -x[0]*d[3] +(-4)/(-1) *x[3]^3*d[4],
// (-6)*x[0]*x[1]*x[2] +(6)*x[2],
// -x[0]*x[1] +(-1)/(-1) ,
// x[1]*x[3]^4 +-x[1]*x[2]^3 +x[4]] 

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingWeyl4.C,v 1.10 2015/07/27 11:50:50 bigatti Exp $
// $Log: ex-RingWeyl4.C,v $
// Revision 1.10  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.9  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.8  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.7  2012/02/08 17:48:04  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.6  2011/03/08 18:04:05  bigatti
// -- changed size_t into long
//
// Revision 1.5  2010/12/22 13:14:29  bigatti
// -- updated with new RingElem ctor (with symbol)
//
// Revision 1.4  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.3  2010/05/14 09:45:29  bigatti
// -- improved syntax/style
//
// Revision 1.2  2007/09/24 14:12:37  abbott
// Added missing newline in a string.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/26 17:39:46  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
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
// Revision 1.1  2006/08/30 15:29:24  cocoa
// -- from WeylAlgebra 3, 5, 7
//
