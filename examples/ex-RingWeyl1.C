// Copyright (c) 2003-2006  John Abbott, Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "An example about RingWeyl, the interface is not quite settled yet.\n";

const string LongDescription =
  "All these examples about RingWeyl will probably be merged into one.\n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    vector<symbol> names = symbols("u,v,x,y,z,t");

    vector<long> ElimIndets;
    ElimIndets.push_back(0); // elim u
    ElimIndets.push_back(1); // elim v
    SparsePolyRing WA = NewWeylAlgebra(RingQQ(), names, ElimIndets);

    RingElem u = indet(WA, 0);
    RingElem v = indet(WA, 1);
    RingElem x = indet(WA, 2);
    RingElem y = indet(WA, 3);
    RingElem z = indet(WA, 4);
    RingElem t = indet(WA, 5);

    RingElem dx = RingElem(WA, "dx");
    RingElem dy = RingElem(WA, "dy");
    RingElem dz = RingElem(WA, "dz");
    RingElem dt = RingElem(WA, "dt");

    cout << "indets(WA) = " << indets(WA) << endl;
    cout << "x*dx = " << x*dx << "    dx*x = " << dx*x << endl;
    cout << "y*dy = " << y*dy << "    dy*y = " << dy*y << endl;
    cout << "z*dz = " << z*dz << "    dz*z = " << dz*z << endl;
    cout << "t*dt = " << t*dt << "    dt*t = " << dt*t << endl;
    cout << "wdeg(x) = " << wdeg(x) <<"  wdeg(dx) = " <<wdeg(dx) <<endl;

    ideal K = ideal(x, dx);
    cout << "gens(K) = " << gens(K) << endl;
    cout << "TidyGens(K) = " << TidyGens(K) << endl;

    //-------
    vector<RingElem> g;
    g.push_back(1-u*v);
    g.push_back(t*u - (x*z+y)*(power(x,6)-power(y,7)));
    g.push_back(v*(z*(power(x,6)-power(y,7)) + (x*z+y)*6*power(x,5))*dt + dx);
    g.push_back(v*((power(x,6)-power(y,7)) - 7*power(y,6)*(x*z+y) )*dt + dy);
    g.push_back(v*(x*(power(x,6)-power(y,7)))*dt + dz);
  
    ideal I = ideal(g);
    cout << "gens(I) = " << gens(I) << endl;
    // ----CoCoA can't manage this example (yet)----
    // cout << "TidyGens(I) = " << TidyGens(I) << endl;

    cout << endl;

    //-------
    g.clear();
    //1-uv, tu-x, vdt+dx
    g.push_back(1 -u*v);
    g.push_back(t*u -x*x*x -y*y*y -z*z*z);
    g.push_back(3*x*x*v*dt +dx);
    g.push_back(3*y*y*v*dt +dy);
    g.push_back(3*z*z*v*dt +dz); 

    // was in file WeylAlgebra9.C
    ideal J = ideal(g);
    cout << "gens(J)=" << gens(J) << endl;
    cout << "TidyGens(J)=" << TidyGens(J) << endl;

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
//-x[0]*x[1] +(1),
//x[0]*x[5] +-x[2]^3 +-x[3]^3 +-x[4]^3,
//(3)*x[1]*x[2]^2*d[5] +d[2],
//(3)*x[1]*x[3]^2*d[5] +d[3],
//(3)*x[1]*x[4]^2*d[5] +d[4]]

//TidyGens(I)=[
//  (3)*v*z^2*dt +dz,
//  (3)*v*y^2*dt +dy,
//  -z^2*dy +y^2*dz,
//  (3)*v*x^2*dt +dx,
//  -z^2*dx +x^2*dz,
//  y^2*dx*dz -x^2*dy*dz,
//  (-2)*y^2*z*dx +(2)*x^2*z*dy,
//  -y^2*dx +x^2*dy,
//  u*t -x^3 -y^3 -z^3,
//  (-3)*u*v*z^2 +(3)*z^2,
//  u*dz +(-3/-1)*z^2*dt,
//  -x^3*dz -y^3*dz -z^3*dz +(3/-1)*z^2*t*dt  +(-3)*z^2,
//  (-2)*x^3*z*dy +(-2)*y^3*z*dy +(-2)*y^2*z^2*dz  +(-6)*y^2*z*t*dt +(6/-1)*y^2*z,
//  x^3*dy +y^3*dy +y^2*z*dz +(-3/-1)*y^2*t*dt +(3)*y^2,
//  (-2)*x^3*z*dx +(2/-1)*x^2*y*z*dy +(-2)*x^2*z^2*dz +(-6)*x^2*z*t*dt +(6/-1)*x^2*z,
//  x^3*dx +x^2*y*dy +x^2*z*dz +(-3/-1)*x^2*t*dt  +(3)*x^2,
//  x*dx*dz +y*dy*dz +z*dz^2 +(3)*t*dz*dt +(-4/-1)*dz,
//  (-2)*x*z*dx +(-2)*y*z*dy +(-2)*z^2*dz +(-6)*z*t*dt +(6/-1)*z,
//  -x*dx -y*dy -z*dz +(3/-1)*t*dt +(-3),
//  (3/-1)*v*x*y*dy*dt +(3/-1)*v*x*z*dz*dt  +(9/-1)*v*x*t*dt^2 +(12/-1)*v*x*dt +dx^2,
//  (-2)*v*y*z*dy*dz*dt +(-6)*v*y*t*dy*dt^2  +(6/-1)*v*z*t*dz*dt^2 +(9/-1)*v*t^2*dt^3  +(-6)*v*y*dy*dt +(6/-1)*v*z*dz*dt +(-36)*v*t*dt^2  +(16/-1)*v*dt +(1/-3)*dx^3 +(-1/-3)*dy^3 +(-1/-3)*dz^3,
//  (3/-1)*v*x*y*z*dz*dt +(9/-1)*v*x*y*t*dt^2  +(6/-1)*v*x*y*dt +y*dx^2 +x*dy^2,
//  (9/-1)*v*y*z*t*dy*dt^2 +(-27/2)*v*z*t^2*dt^3  +(3/-1)*v*y*z*dy*dt +(36/-1)*v*z*t*dt^2  +(-6)*v*z*dt +(-1/2)*z*dx^3 +(1/2)*z*dy^3 +y*dy*dz^2  +(1/2)*z*dz^3 +(3)*t*dz^2*dt +(3)*dz^2,
//  (-9)*v*y*z*t*dz*dt^2 +(-27/2)*v*y*t^2*dt^3  +(-3)*v*y*z*dz*dt +(-36)*v*y*t*dt^2 +(6/-1)*v*y*dt  +(-1/2)*y*dx^3 +(1/2)*y*dy^3 +z*dy^2*dz +(1/2)*y*dz^3  +(-3/-1)*t*dy^2*dt +(-3/-1)*dy^2,
//  (9/-1)*v*x*z*t*dz*dt^2 +(-27/2)*v*x*t^2*dt^3  +(3/-1)*v*x*z*dz*dt +(36/-1)*v*x*t*dt^2  +(6/-1)*v*x*dt +(-1/2)*y*dx^2*dy +(1/-2)*x*dy^3  +(1/2)*z*dx^2*dz +(1/2)*x*dz^3 +(3/2)*t*dx^2*dt  +(1/2)*dx^2,
//  (9/-1)*v*x*y*z*t*dt^2 +y*z*dx^2 +x*z*dy^2  +x*y*dz^2,
//  (-6)*v*y*t^2*dy*dt^3 +(6/-1)*v*t^3*dt^4  +(12/-1)*v*y*t*dy*dt^2 +(36/-1)*v*t^2*dt^3  +(4/-3)*v*y*dy*dt +(-112/3)*v*t*dt^2 +(8/-3)*v*dt  +(2/-27)*y*dx^3*dy +(-2/-27)*y*dy^4 +(4/27)*z*dx^3*dz  +(-2/9)*y*dy*dz^3 +(-4/27)*z*dz^4 +(2/-9)*t*dx^3*dt  +(-2/-3)*t*dy^3*dt +(2/-3)*t*dz^3*dt +(-4/-9)*dy^3 +(-8/9)*dz^3,
//  (-9/2)*v*z*t^2*dz*dt^3 +(-9/2)*v*t^3*dt^4  +(-9)*v*z*t*dz*dt^2 +(27/-1)*v*t^2*dt^3 -v*z*dz*dt  +(-28)*v*t*dt^2 +(-2)*v*dt +(-1/-9)*y*dx^3*dy +(1/-9)*y*dy^4  +(-1/18)*z*dx^3*dz +(1/-6)*z*dy^3*dz +(1/18)*z*dz^4  +(-1/6)*t*dx^3*dt +(-1/2)*t*dy^3*dt +(1/2)*t*dz^3*dt  +(2/-3)*dy^3 +(-1/-3)*dz^3,
//  (9/-2)*v*y*z*t^2*dt^3 +(6/-1)*v*y*z*t*dt^2  +(1/-6)*y*z*dx^3 +(-1/-6)*y*z*dy^3 +(1/3)*y^2*dy*dz^2  +(-1/-6)*y*z*dz^3 +z*t*dy^2*dt +y*t*dz^2*dt  +(1/3)*z*dy^2 +y*dz^2,
//  (-9/2)*v*x*z*t^2*dt^3 +(6/-1)*v*x*z*t*dt^2  +(-1/6)*y*z*dx^2*dy +(1/-6)*x*z*dy^3  +(-1/6)*x*y*dy*dz^2 +(1/2)*z*t*dx^2*dt  +(-1/-2)*x*t*dz^2*dt +(1/-6)*z*dx^2 +(1/-6)*x*dz^2,
//  (9/-2)*v*x*y*t^2*dt^3 +(6/-1)*v*x*y*t*dt^2  +(1/-6)*y*z*dx^2*dz +(1/-6)*x*z*dy^2*dz  +(-1/6)*x*y*dz^3 +(-1/-2)*y*t*dx^2*dt  +(-1/-2)*x*t*dy^2*dt +(1/-6)*y*dx^2 +(-1/6)*x*dy^2,
//  (3/-1)*v*z*t^3*dt^4 +(-12)*v*z*t^2*dt^3  +(-20/3)*v*z*t*dt^2 +(2/27)*y*z*dx^3*dy  +(-2/27)*y*z*dy^4 +(4/-27)*y^2*dy^2*dz^2  +(-2/27)*y*z*dy*dz^3 +(1/-9)*z*t*dx^3*dt  +(1/-3)*z*t*dy^3*dt +(-2/9)*y*t*dy*dz^2*dt  +(-1/-9)*z*t*dz^3*dt +(2/3)*t^2*dz^2*dt^2 +(2/27)*z*dx^3  +(2/-9)*z*dy^3 +(-20/27)*y*dy*dz^2 +(2/-27)*z*dz^3  +(-8/-9)*t*dz^2*dt +(4/-9)*dz^2,
//  (9/-2)*v*y*t^3*dt^4 +(-18)*v*y*t^2*dt^3  +(-10)*v*y*t*dt^2 +(-1/-9)*y*z*dx^3*dz  +(-1/9)*y*z*dy^3*dz +(-2/9)*y^2*dy*dz^3  +(1/-9)*y*z*dz^4 +(1/-6)*y*t*dx^3*dt  +(-1/-6)*y*t*dy^3*dt +(1/-3)*z*t*dy^2*dz*dt  +(1/-2)*y*t*dz^3*dt +t^2*dy^2*dt^2 +(-1/-9)*y*dx^3  +(-1/9)*y*dy^3 +(2/-9)*z*dy^2*dz +(-7/9)*y*dz^3  +(4/3)*t*dy^2*dt +(2/-9)*dy^2,
//  (-9/2)*v*x*t^3*dt^4 +(-18)*v*x*t^2*dt^3  +(10/-1)*v*x*t*dt^2 +(-1/-9)*y*z*dx^2*dy*dz  +(-1/-9)*x*z*dy^3*dz +(-1/-9)*x*y*dy*dz^3  +(1/-6)*y*t*dx^2*dy*dt +(1/-6)*x*t*dy^3*dt  +(1/-6)*z*t*dx^2*dz*dt +(-1/6)*x*t*dz^3*dt  +(-1/-2)*t^2*dx^2*dt^2 +(1/9)*y*dx^2*dy +(-1/-9)*x*dy^3  +(-1/-9)*z*dx^2*dz +(-1/-9)*x*dz^3 +(1/3)*t*dx^2*dt  +(-1/-9)*dx^2,
//  (-9/2)*v*t^4*dt^5 +(-36)*v*t^3*dt^4 +(64/-1)*v*t^2*dt^3  +(-20)*v*t*dt^2 +(-1/9)*y*z*dx^3*dy*dz  +(1/9)*y*z*dy^4*dz +(-2/-9)*y^2*dy^2*dz^3  +(1/9)*y*z*dy*dz^4 +(-1/-9)*y*t*dx^3*dy*dt  +(1/-9)*y*t*dy^4*dt +(-1/-9)*z*t*dx^3*dz*dt  +(1/3)*z*t*dy^3*dz*dt +(1/3)*y*t*dy*dz^3*dt  +(1/-9)*z*t*dz^4*dt +(-1/6)*t^2*dx^3*dt^2  +(-1/2)*t^2*dy^3*dt^2 +(-1/2)*t^2*dz^3*dt^2  +(-1/9)*y*dx^3*dy +(1/9)*y*dy^4 +(-1/9)*z*dx^3*dz  +(-1/-3)*z*dy^3*dz +(11/9)*y*dy*dz^3 +(-1/-9)*z*dz^4  +(-2/3)*t*dy^3*dt +(2/-3)*t*dz^3*dt +(-1/9)*dx^3  +(-1/-3)*dy^3 +(7/9)*dz^3,
//  (-2)*u*z*dy +(6/-1)*y^2*z*dt,
//  -u*dy +(-3)*y^2*dt,
//  (-2)*u*z*dx +(6/-1)*x^2*z*dt,
//  -u*dx +(-3)*x^2*dt,
//  (-6)*u*v*z +(6)*z,
//  -u*v +(-1/-1) ,
//  -v*x^3 -v*y^3 -v*z^3 +t]

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingWeyl1.C,v 1.12 2019/03/27 14:55:37 bigatti Exp $
// $Log: ex-RingWeyl1.C,v $
// Revision 1.12  2019/03/27 14:55:37  bigatti
// -- new syntax for assigning variables
//
// Revision 1.11  2015/07/27 11:50:50  bigatti
// -- now using "symbols(string)" for comma separated symbols
//
// Revision 1.10  2015/07/01 16:31:36  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.9  2015/06/29 15:47:57  bigatti
// -- code in namespace CoCoA
//
// Revision 1.8  2012/02/08 17:48:04  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.7  2011/03/08 18:02:07  bigatti
// -- changed size_t into long
//
// Revision 1.6  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.5  2010/10/01 16:08:07  bigatti
// -- removed function "derivation" and using RingElem(Wa, symbol("dx"))
//
// Revision 1.4  2010/05/14 11:45:07  bigatti
// -- fixed minor bug
//
// Revision 1.3  2010/05/14 09:45:29  bigatti
// -- improved syntax/style
//
// Revision 1.2  2008/11/19 09:19:48  bigatti
// -- added easy example:  ideal(x, dx)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/02/26 17:18:22  bigatti
// -- getting ready for unique ring Z: using NewZZmod(N), NewRingQQ()
//
// Revision 1.6  2007/02/12 16:15:37  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.5  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.4  2006/08/30 15:22:12  cocoa
// -- added example from WeylAlgebra9.C
//
// Revision 1.3  2006/08/17 10:03:24  cocoa
// -- creation of RingWeyl with given list of symbols
// -- updated header
//
