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
    ElimIndets.push_back(2); // elim x
    ElimIndets.push_back(3); // elim y
    ElimIndets.push_back(7); // elim dx
    ElimIndets.push_back(8); // elim dy
    SparsePolyRing WA = NewWeylAlgebra(RingQQ(), names, ElimIndets);

    RingElem x(WA, symbol("x"));
    RingElem y(WA, symbol("y"));
    RingElem u(WA, symbol("u"));
    RingElem t(WA, symbol("t"));

    RingElem dx(WA, symbol("dx"));
    RingElem dy(WA, symbol("dy"));
    RingElem dt(WA, symbol("dt"));

    // was in WeylAlgebra8.C
    vector<RingElem> g;
    //1-uv,tu-x,vdt+dx
    g.push_back(power(x,3)-power(y,4));
    g.push_back(4*power(y,3)*dx +3*x*x*dy);
    g.push_back(-4*power(x,3)*dx -3*x*x*y*dy -12*x*x*(-t-1) -12*x*x);
    g.push_back(-3*x*x*power(y,4)*dy-12*x*x*power(y,3)*(-t-1)+3*power(x,5)*dy-12*x*x*power(y,3));
    g.push_back(-2*x*power(y,4)*dy -8*x*power(y,3)*(-t-1)+2*power(x,4)*dy -8*x*power(y,3));
    g.push_back(-power(y,4)*dy -4*power(y,3)*(-t-1)+power(x,3)*dy -4*power(y,3));
    g.push_back(4*x*dx*dx + 3*y*dx*dy +12*dx*(-t-1)+16*dx);
    g.push_back(-2*x*x*dx -3/2*x*y*dy -6*x*(-t-1)-6*x);
    g.push_back(-4*x*dx-3*y*dy-12*(-t-1)-12);
    test(ideal(WA, g));

    // was in WeylAlgebra6.C
    vector<RingElem> f;
    f.push_back(x*x+power(y,3));
    f.push_back(3*y*y*dx -2*x*dy);
    // f.push_back(-x*x*dx +(-2)/(3)*x*y*dy +(2)/(-1)*x*(-t-1) +(-2)*x);
    // f.push_back((2)/(-3)*x*power(y,3)*dy+(-2)*x*y*y*(-t-1)+(-2)/(3)*power(x,3)*dy-2*x*y*y);
    // f.push_back(-power(y,3)*dy -3*y*y*(-t-1)-x*x*dy -3*y*y);
    // f.push_back(x*dx*dx + 2/3*y*dx*dy +2*dx*(-t-1)+3*dx);
    f.push_back(-3*x*dx -2*y*dy -6*(-t-1)-2);
    test(ideal(WA, f));  
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
//gens(I)=[
//-x[3]^4 +x[2]^3,
//(4)*x[3]^3*d[2] +(3)*x[2]^2*d[3],
//(-4)*x[2]^3*d[2] +(-3)*x[2]^2*x[3]*d[3] +(12)*x[2]^2*x[4],
//(-3)*x[2]^2*x[3]^4*d[3] +(3)*x[2]^5*d[3] +(12)*x[2]^2*x[3]^3*x[4],
//(-2)*x[2]*x[3]^4*d[3] +(2)*x[2]^4*d[3] +(8)*x[2]*x[3]^3*x[4],
//-x[3]^4*d[3] +x[2]^3*d[3] +(4)*x[3]^3*x[4],
//(4)*x[2]*d[2]^2 +(3)*x[3]*d[2]*d[3] +(-12)*x[4]*d[2] +(4)*d[2],
//(-2)*x[2]^2*d[2] +-x[2]*x[3]*d[3] +(6)*x[2]*x[4],
//(-4)*x[2]*d[2] +(-3)*x[3]*d[3] +(12)*x[4]]
 
//TidyGens(I)=[
//(-4)*x[2]*d[2] +(4)/(-1) ,
//(-1)/(-2) *x[2]*x[3]*d[3],
//(-3)*x[3]^2*d[3]^2 +(48)*x[4]^2 +(20)*x[4] +(-4)/(-3) ,
//(4)*x[3]^3*d[2],
//(-3)*x[3]^4*d[3] +(3)*x[2]^3*d[3],
//(9)/(4) *x[2]^2*d[3]^3 +(12)*x[3]^2*x[4]*d[2]*d[3] +(-76)/(-1) *x[3]*x[4]*d[2] 
//+(130)/(3) *x[3]*d[2],
//(64)/(3) *x[3]^2*x[4]*d[2]^2*d[3] +(304)/(9) *x[3]^2*d[2]^2*d[3] +(12)*x[2]*x[4]*d[3]^3 +(8)*x
//[2]*d[3]^3 +(32)*x[3]*d[2]^2,
//(-2)/(-1) *x[2]^2*x[4]*d[3]^2,
//(12)*x[3]^2*x[4]^2*d[2]*d[3] +(-35)/(-1) *x[3]^2*x[4]*d[2]*d[3] 
//+(-358)/(-3) *x[3]*x[4]*d[2] +(520)/(9) *x[3]*d[2],
//(12)*x[3]*x[4]^3*d[2]*d[3] +(-289)/(4) *x[3]*x[4]*d[2]*d[3] 
//+(126)*x[4]^2*d[2] +(797)/(-9) *x[4]*d[2] +(5915)/(-108) *d[2],
//(-9)/(2) *x[3]^2*d[2]*d[3] +(18)*x[3]*x[4]*d[2] +(-6)/(-1) *x[3]*d[2],
//(-48)/(-1) *x[3]*x[4]^3*d[2] +(-133)*x[3]*x[4]*d[2] +(-1495)/(18) *x[3]*d[2],
//(-3)/(-1) *x[3]*x[4]^2*d[2]*d[3] +(-35)/(-4) *x[3]*x[4]*d[2]*d[3] 
//+(-9)*x[4]^2*d[2] +(-151)/(-12) *x[4]*d[2] +(-455)/(-72) *d[2],
//(-9)/(-2) *x[3]*d[2]*d[3] +(18)/(-1) *x[4]*d[2] +(6)/(-1) *d[2],
//(32)/(-1) *x[4]^6*d[2] +(216)/(-1) *x[4]^5*d[2] +(-227623)/(216) *x[4]^2*d[2] 
//+(3218893)/(-2592) *x[4]*d[2] +(17689945)/(-46656) *d[2],
//(-3)*x[3]*x[4]^6*d[3] +(-81)/(4) *x[3]*x[4]^5*d[3] +(913)/(-16) *x[3]*x[4]^4*d[3] 
//+(-475451)/(-2304) *x[3]*x[4]^2*d[3] +(-7744249)/(-6912) *x[4]^2 
//+(-214824775)/(-248832) *x[4] +(-279522425)/(-1492992) ,
//(2)*x[2]*x[4]^7 +(27)/(2) *x[2]*x[4]^6 +(-913)/(-24) *x[2]*x[4]^5,
//(4)*x[4]^7 +(27)*x[4]^6 +(-3609797)/(20736) *x[4]^2 +(163334989)/(-746496) *x[4] 
//+(103568465)/(-1492992) ,
//(-3)*x[3]^2*x[4]^3*d[3] +(-39)/(4) *x[3]^2*x[4]^2*d[3] +(-83)/(8) *x[3]^2*x[4]*d[3] 
//+(-65)/(18) *x[3]^2*d[3],
//(-7)/(-1) *x[3]*x[4]^3*d[3] +(-91)/(-4) *x[3]*x[4]^2*d[3] +(-1225)/(-12) *x[4]^2 
//+(-32851)/(-432) *x[4] +(-41405)/(-2592) ,
//(-144)/(-7) *x[4]^4*d[2] +(87)/(-1) *x[4]^2*d[2] +(-1955)/(21) *x[4]*d[2] 
//+(-325)/(12) *d[2],
//(48)/(7) *x[4]^5 +(5225)/(126) *x[4]^2 +(-7139)/(-144) *x[4] +(13195)/(864) ,
//(48)/(7) *x[3]*x[4]^4 +(-6187)/(-252) *x[3]*x[4] +(-1235)/(-72) *x[3],
//(-6)/(-7) *x[2]*x[4]^4,
//(-125)/(32) *x[2]^2*d[3]^2,
//(-12)/(-1) *x[3]*x[4]*d[3] +(48)/(-1) *x[4]^2 +(16)/(-1) *x[4],
//(-455)/(-54) *x[3]*d[3] +(-910)/(27) *x[4] +(910)/(-81) ,
//(36)/(-65) *x[2]*x[4]^3,
//(-65)/(-21) *x[2]*x[4]^2,
//(-81)/(-130) *x[2]*x[4],
//(160)/(9) *x[3]^2*x[4]^2*d[2] +(-1690)/(81) *x[3]^2*d[2],
//(-288)/(-25) *x[3]^2*x[4]*d[2] +(-312)/(-25) *x[3]^2*d[2],
//(18)*x[3]*x[4]^2*d[2] +(69)/(2) *x[3]*x[4]*d[2] +(65)/(4) *x[3]*d[2],
//(12)/(-1) *x[4]^3*d[2] +(30)/(-1) *x[4]^2*d[2] +(-97)/(4) *x[4]*d[2] +(-455)/(72) *d[2],
//(3)*x[4]^4 +(-203)/(16) *x[4]^2 +(1955)/(-144) *x[4] +(-2275)/(576) ,
//(36)/(7) *x[3]*x[4]^3 +(57)/(-4) *x[3]*x[4] +(1495)/(-168) *x[3],
//(13)/(4) *x[2]^2*d[3],
//(8)*x[2]*d[3],
//(-3456)/(-455) *x[4]^3*d[3] +(-1728)/(-91) *x[4]^2*d[3] +(-6984)/(-455) *x[4]*d[3] +(4)*d[3],
//(-24)/(-5) *x[3]^2*x[4]^2 +(-169)/(30) *x[3]^2,
//(48)/(13) *x[3]^3*x[4],
//-x[3]^4 +x[2]^3,
//(-12)/(-1) *x[3]^3,
//(-70)/(3) *x[3]^2*x[4] +(455)/(-18) *x[3]^2,
//(144)/(-5) *x[3]*x[4]^2 +(-276)/(5) *x[3]*x[4] +(26)/(-1) *x[3],
//(144)/(-7) *x[4]^3 +(-360)/(7) *x[4]^2 +(-291)/(7) *x[4] +(-65)/(6) ,
//(-7)/(-24) *x[2]]

//output
//gens(I)=[
//x[3]^3 +x[2]^2,
//(3)*x[3]^2*d[2] +(-2)*x[2]*d[3],
//(-3)*x[2]*d[2] +(-2)*x[3]*d[3] +(6)*x[4] +(4)]
 
//TidyGens(I)=[
// (-3)*x[2]*d[2] +(-2)*x[3]*d[3] +(6)*x[4] +(4),
// (3)*x[3]^2*d[2] +(-2)*x[2]*d[3],
// (-2)*x[3]^3*d[3] +(-2)*x[2]^2*d[3] +(-6)*x[3]^2,
// x[3]^3 +x[2]^2,
// (6)*x[3]^2*x[4] +(10)*x[3]^2,
// (-2)*x[2]*x[4]*d[3] +(-10)/(3) *x[2]*d[3],
// (-2)*x[3]*x[4]*d[3]^2 +(-10)/(3) *x[3]*d[3]^2 +(6)*x[4]^2*d[3] 
// +(-15)/(-1) *x[4]*d[3] +(-25)/(-3) *d[3],
// x[2]^2*x[4] +(5)/(3) *x[2]^2,
// (6)*x[2]*x[4]^2 +(-50)/(3) *x[2],
// (-2)*x[3]*x[4]^2*d[3] +(50)/(9) *x[3]*d[3] +(-23)*x[4]^2 +(-133)/(2) *x[4] +(845)/(-18) ,
// (-18)/(-1) *x[3]*x[4]^3 +(-331)/(2) *x[3]*x[4] +(385)/(-2) *x[3],
// (6)*x[4]^4 +(601)/(-6) *x[4]^2 +(-665)/(3) *x[4] +(-275)/(2) ,
// (9)*x[3]*x[4]*d[3] +(15)*x[3]*d[3] +(-27)*x[4]^2 +(153)/(-2) *x[4] +(-105)/(2) ,
// (12)*x[4]^3*d[3] +(60)*x[4]^2*d[3] +(299)/(3) *x[4]*d[3] +(55)*d[3],
// (-8)/(3) *x[4]^3 +(40)/(-3) *x[4]^2 +(598)/(-27) *x[4] +(-110)/(9) ,
// (3)*x[3]*x[4]^2 +(21)/(2) *x[3]*x[4] +(-55)/(-6) *x[3],
// (-1)/(-6) *x[2]*x[4] +(-5)/(-18) *x[2]]

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingWeyl5.C,v 1.10 2015/07/27 11:50:49 bigatti Exp $
// $Log: ex-RingWeyl5.C,v $
// Revision 1.10  2015/07/27 11:50:49  bigatti
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
// Revision 1.6  2011/03/30 09:11:55  bigatti
// -- removed failing assignment to "z"
//
// Revision 1.5  2011/03/08 18:03:52  bigatti
// -- changed size_t into long
// -- using RingElem(ring, symbol) ctor instead of "derivation"
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
// Revision 1.1  2006/08/30 15:42:21  cocoa
// -- from old examples WeylAlgebra8, 6
//
// Revision 1.1  2006/08/30 15:29:24  cocoa
// -- from WeylAlgebra 3, 5, 7
//
