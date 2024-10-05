// Copyright (c) 2013-2015  John Abbott,  Anna M Bigatti
// Orig author: 2013-2015  Mario Albert
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
#include <bitset>
using namespace std;
using namespace CoCoA;
using namespace CoCoA::Involutive;

//----------------------------------------------------------------------
const string ShortDescription =
  "In this file we explain how to using some functions which are using the Janet Basis.  \n";

const string LongDescription =
  " PrintNonMultVar, PrintMultVar, IsPommaretBasis, StandardRepresentation(f) \n";
//----------------------------------------------------------------------

void output(vector<RingElem> vec)
{
  for(vector<RingElem>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
    {
      cout << *iter << endl;
    }
}


void program()
{
  GlobalManager CoCoAFoundations;

  cout << ShortDescription << endl;
  cout << boolalpha; // so that bools print out as true/false
  ring Q = RingQQ();

  cout <<"/////////////////////////////////////Cyclic 4//////////////////////////////////////////" << endl;
  SparsePolyRing polyRingCyc4 = NewPolyRing(Q, SymbolRange("x",0,3));
  vector<RingElem> x = indets(polyRingCyc4);
  vector<RingElem>cyclic4;
  cyclic4.push_back(x[0]+x[1]+x[2]+x[3]);
  cyclic4.push_back(x[0]*x[1] + x[2]*x[1] + x[2]*x[3] + x[0]*x[3]);
  cyclic4.push_back(x[1]*x[2]*x[0] + x[2]*x[3]*x[1] +  x[0]*x[3]*x[1] + x[0]*x[2]*x[3]);
  cyclic4.push_back(x[0]*x[1]*x[2]*x[3] - 1);



  cout << "---------computing the Janet Base with a list of polynomials:------------------------" << endl;
  Involutive::JBMill mill(Involutive::JBMill::Builder().setInput(cyclic4));
  cout << endl;
  cout << "Janet Basis is delta-regular (there exists a finite Pommaret Basis): " << mill.IamPommaretBasis() << endl;
  cout << "Janet Basis is homogenous (Ideal is homogenous)" << mill.IamHomogenous() << endl;
  cout << "Janet Basis is monomial (Ideal is monomial)" << mill.IamMonomialIdeal() << endl;
  cout << endl;
  cout << "-----nonmultiplicative variables of the Janet-Base of cylic 4 (easy way)" << endl;
  mill.myPrintNonMultVar();

  cout << endl << endl ;

  cout << "-----nonmultiplicative variables of the Janet-base of cyclic 4 (complicate way)" << endl;
  map<PPMonoidElem, vector<bool> > nonMultVars = mill.myNonMultVars();
  for (map<PPMonoidElem, vector<bool> >::iterator i = nonMultVars.begin(); i != nonMultVars.end(); ++i)
  {
    cout << i->first << endl;
    long counter = 0;
    for (std::vector<bool>::iterator iter = i->second.begin(); iter != i->second.end(); ++iter)
    {
      if(*iter)
      {
        cout << counter << " ";
      }
      ++counter;
    }
    cout << endl;
  }

  cout << "-----multiplicative variables of the Janet-base of cyclic 4 (easy way)" << endl;
  mill.myPrintMultVar();

  cout << endl << endl ;

  cout << "-----multiplicative variables of the Janet-base of cyclic 4 (complicate way)" << endl;
  map<PPMonoidElem,vector<bool> > MultVars = mill.myMultVars();
  for (map<PPMonoidElem,vector<bool> >::iterator i = MultVars.begin(); i != MultVars.end(); ++i)
  {
    cout << i->first << endl;
    long counter = 0;
    for (std::vector<bool>::iterator iter = i->second.begin(); iter != i->second.end(); ++iter)
    {
      if(*iter)
      {
        cout << counter << " ";
      }
      ++counter;
    }
    cout << endl;
  }

  // cout << endl << endl;
  // cout << "-----checking if the Janet-base of cyclic 4 is also a Pommaret basis" << endl;

  // cout << "Janet-base of cyclic 4 is also a Pommaret basis: " << JBIsPommaretBasis(result) << endl;

  cout << endl << endl;
  RingElem f = 324*x[1]*x[1]*x[0]*x[0] + 234*power(x[0],4) - 123;

  cout << "involutive standard representation of the element f (easy way)" << endl;
  cout << "f = " << f << endl;
  mill.myOutputStandardRepresentation(f);

  cout << "involutive standard representation of the element f (complicate way)" << endl;
  pair<map<PPMonoidElem, RingElem>, RingElem> representation = mill.myStandardRepresentation(f);
  map<PPMonoidElem, RingElem> sRep = representation.first;
  RingElem res = representation.second;
  cout << "Involutive Standard Representation" << endl;
  cout << res << " =" << endl;
  cout << endl;
  for(map<PPMonoidElem, RingElem>::iterator iter(sRep.begin()); iter != sRep.end(); ++iter)
    {
      cout << iter->second <<  endl;
    }
  cout << "rest = " << res << endl;
  cout << "----------------------------------" << endl;


  cout << "the hilbert polynomial of cyclic 4 (with variable x[0])" << endl;

  ring QQ = RingQQ();
  PolyRing hilbPolyRing = NewPolyRing(QQ, symbols("t"));

  cout << "hilbert Polynomial = " << mill.myHilbertPol(indet(hilbPolyRing, 0)) << endl << endl;

  cout << "the hilbert function of cyclic 4 (with s = 10)" << endl;
  std::cout << "hilbfunc(10)=" << mill.myHilbertFunc(BigInt(10)) << endl;; //s must be of type ZZ

  cout << "hilbert function version 2 with functional expression" << endl;
  mill.myHilbertFunc();

  cout << endl << "the rational function of the hilbert series of cyclic 4" << endl;
  FractionField QR = NewFractionField(polyRingCyc4);
  RingHom embeddedHom(EmbeddingHom(QR));
  RingElem qx = embeddedHom(x[0]);

  cout << mill.myHilbertSeries(qx) << endl << endl;

  cout << "generators of the first syzygy" << endl;
  FGModule firstSyz = mill.mySyzygy();

  vector<ModuleElem> gensFirstSyz = gens(firstSyz);
  for(vector<ModuleElem>::iterator iter = gensFirstSyz.begin(); iter != gensFirstSyz.end(); ++iter)
  {
    cout << *iter << endl;
  }
  cout << endl;

  cout << "dimension of the P/I = " << mill.myDim() << endl << endl;

  cout << "Complementary Decomposition" << endl;

  // if mill.IamMonomial() is true one can use myComplementaryDecomposition as well
  std::vector< std::pair<PPMonoidElem, std::vector<bool> > > cD = mill.myComplementaryDecompositionPolynomial();
  for (std::vector< std::pair<PPMonoidElem, std::vector<bool> > >::iterator i = cD.begin(); i < cD.end(); ++i)
  {
    cout << "mon = " << i->first << endl;
    long counter = 0;
    cout << "vars = ";
    for (std::vector<bool>::iterator iter = i->second.begin(); iter != i->second.end(); ++iter)
    {
      if(*iter == true)
      {
        cout << counter << " ";
      }
      ++counter;
    }
    cout << endl;
  }
  cout << endl;
}

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
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
  return 1;
}
