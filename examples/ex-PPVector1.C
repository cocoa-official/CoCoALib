// Copyright (c)  2010  John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

#include <algorithm>
using std::min;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example program shows the use of PPVector.            \n"
  "PPVector is still work in progress, so syntax might change \n";

const string LongDescription =
  "PPVector is a vector of PPWithMask.                       \n"
  "It's been designed to represent the list of generators of \n"
  "monomial ideals (e.g. to facilitate interreduction).      \n"
  "USE WITH CARE!  The interface may still change slightly.";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    const int N = 3;
    const PPMonoid PPM1 = NewPPMonoidEv(SymbolRange("x", 0, (N*N)-1), lex);
    const DivMaskRule DMR1 = NewDivMaskEvenPowers();

    const std::vector<PPMonoidElem> x = indets(PPM1);

    PPVector PPV(PPM1, DMR1);
    for (long i=0 ; i<N ; ++i )
      for (long j=0 ; j<N ; ++j )
        PushBack(PPV, x[N*i+j]*x[N*j+i]);

    for (long i=0 ; i<len(PPV) ; ++i )  cout << PPV[i] << endl;
    cout << "len(PPV) = " << len(PPV) << endl;

    interreduce(PPV);
    cout << "After interreduce(PPV): len(PPV) = " << len(PPV) << endl;
    cout << "--------------------------------------------------" << endl;

    cout << "-- working directly with PPVector -- " << endl;
    PPVector f(PPM1, DMR1);
    PPVector g(PPM1, DMR1);
    PushBack(f, x[1]*x[2]);  PushBack(f, x[3]*x[2]);
    PushBack(g, x[1]*x[5]);  PushBack(g, x[6]*x[7]);
    cout << "PP(f[0]), PP(f[1]) are: " << PP(f[0]) << ", " <<  PP(f[1]) << endl;
    cout << "PP(g[0]), PP(g[1]) are: " << PP(g[0]) << ", " <<  PP(g[1]) << endl;
  
    PPVector lcms_fg(PPM1, DMR1);
    lcms(lcms_fg, f, g);
    cout << "lcms(lcms_fg, f, g);  -->\n" << lcms_fg << endl;
    cout << "--------------------------------------------------" << endl;

    cout << "-- conversion PPVector <--> monomial ideal -- " << endl;
    const ring P = NewPolyRing(RingQQ(), symbols(PPM1));
    std::vector<RingElem> v;
    convert(v, P, lcms_fg);
    cout << "convert(v, P, lcms_fg);  -->\n" << v << endl;

    PPVector ppv(PPM1, DMR1);
    convert(ppv, v);
    cout << "convert(ppv, v);  -->\n" << ppv << endl;

    cout << "\nSame as PPVector ppv1(PPM1, DMR1, v);" << endl;
    PPVector ppv1(PPM1, DMR1, v);
    cout << "ppv1;  -->\n" << ppv1 << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-PPVector1.C,v 1.8 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-PPVector1.C,v $
// Revision 1.8  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.7  2019/10/29 16:01:17  bigatti
// -- added use of constructor
//
// Revision 1.6  2018/09/28 15:54:03  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.5  2015/06/29 15:41:17  bigatti
// *** empty log message ***
//
// Revision 1.4  2015/06/29 12:45:12  bigatti
// -- code in namespace CoCoA
//
// Revision 1.3  2014/07/14 11:47:01  abbott
// Summary: Minor tidying
// Author: JAA
//
// Revision 1.2  2014/07/03 06:44:45  bigatti
// -- improved
//
