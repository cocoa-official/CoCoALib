// Copyright (c) 2005 John Abbott,  Anna Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example program illustrating an implementation of Normal Remainder \n"
  "wrt a list of polynomials.                                         \n"
  "If the list is a Groebner Basis, NR returns the Normal Form.       \n";

const string LongDescription =
  "This is just an example!  If you want to compute Normal Forms \n"
  "you should use the library function \"NF\".                   \n";
//----------------------------------------------------------------------


namespace CoCoA
{

  long FindReducer(ConstRefRingElem f, const vector<RingElem>& v)
  {
    if (IsZero(f)) return -1;

    for (long j=0; j<len(v); ++j)  // 0 to len-1
      if (IsDivisible(LPP(f), LPP(v[j])))
        return j;
    return -1;
  }
  

  RingElem NRLM(ConstRefRingElem f, const vector<RingElem>& v)
  {
    const SparsePolyRing P = owner(f);
    RingElem m(P);
    RingElem r =f;
  
    long j = FindReducer(r, v);
    while (j != -1)
    {
      //   m = LM(r)/LM(v[i]);  
      P->myDivLM(raw(m), raw(r), raw(v[j])); // no checks
      r -= m * v[j];
      if (IsZero(r)) return zero(P);
      j = FindReducer(r, v);
    }
    return r;
  }


  RingElem NormalRemainder(ConstRefRingElem f, const vector<RingElem>& v)
  {
    if (IsZero(f)) return f;
    const SparsePolyRing P = owner(f);
    RingElem ansNR(P);
    RingElem tmpNR =f;
  
    tmpNR = NRLM(f, v);
    while (!IsZero(tmpNR))
    {
      P->myMoveLMToBack(raw(ansNR), raw(tmpNR));
      tmpNR = NRLM(tmpNR, v);
    }
    return ansNR;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    ring P = NewPolyRing(RingQQ(), symbols("x,y,z")); // QQ[x,y,z]

    vector<RingElem> g = RingElems(P, "x^3 -z^3,  x^2 -y");
    ideal I = ideal(g);
    cout << "gens(I)   = " << gens(I) << endl;
    vector<RingElem> GB = GBasis(I);  // it is the Groebner Basis of I
    cout << "GBasis(I) = " << GB << endl;
    cout << "When I is an ideal of polynomials \"TidyGens(I)\" returns its Groebner Basis." << endl << endl;

    RingElem f = RingElem(P, "x^12 + y^6 + z^6");
    cout << "-- f = " << f << endl;
    cout << "NormalRemainder(f, gens(I)) = " <<  NormalRemainder(f, gens(I)) << endl;
    cout << "NormalRemainder(f, GB)      = " <<  NormalRemainder(f, GB) << endl;
    cout << "NF(f, I)                    = " <<  NF(f, I) << endl;
    cout << endl;

    const vector<RingElem>& x = indets(P);
    RingElem h = g[0] + x[0]*g[1];
    cout << "-- h = " << h << endl;
    cout << "NormalRemainder(h, gens(I)) = " <<  NormalRemainder(h, gens(I)) << endl;
    cout << "NormalRemainder(h, GB)      = " <<  NormalRemainder(h, GB) << endl;
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
