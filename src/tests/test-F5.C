//   Copyright (c)  2007,2010,2024  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/BuildInfo.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpF5.H"
#include "CoCoA/symbol.H"
#include "CoCoA/verbose.H"
#include "CoCoA/VectorOps.H"

#include <algorithm> // for sort

using std::cerr;
using std::cout;
using std::endl;
using std::sort;
using std::vector;

namespace CoCoA
{

  bool LessThan(const RingElem& f, const RingElem& g)  { return LPP(f)<LPP(g); }

  
  bool EqSet(vector<RingElem>& G1, vector<RingElem>& G2)
  {
    if (len(G1) != len(G2) ) { cout << " --lengths differ"; return false; }
    sort(G1.begin(), G1.end(), LessThan);
    sort(G2.begin(), G2.end(), LessThan);
    for (long i=0; i<len(G1); ++i)  if (G1[i] != G2[i]) return false;
    return true;
  }

  
  int FindReducerIndex(const PPMonoidElem& pp, const vector<RingElem>& v)
  {
    const long nelems = len(v);
    for (long i=0; i < nelems; ++i)
      if (IsDivisible(pp, LPP(v[i])))  return i;
    return -1;
  }


  RingElem NRLM(ConstRefRingElem f, const vector<RingElem>& G)
  {
    if (IsZero(f)) return f;
    const SparsePolyRing P = owner(f);
    RingElem LMfDivLMGi(P), NRf(f);
    long i;
    while ((i = FindReducerIndex(LPP(NRf), G)) != -1)
    {
      P->myDivLM(raw(LMfDivLMGi), raw(NRf), raw(G[i]));
      NRf -= LMfDivLMGi * G[i];
      if (IsZero(NRf)) return NRf;
    }
    return NRf;
  }


  
  //-- program   --------------------------------------------------------------
  void program()
  {
    GlobalManager CoCoAFoundations;
    SparsePolyRing Qx = NewPolyRing(RingQQ(), symbols("a,b,c,d,e,f"));
    
    ideal I = ideal(RingElems(Qx, "b^6 -a*d^5 -2*b^2*c^2*f^2,  b^3 -b*c*d +d*e*f,  a^3 -b*c*e +c*e*f"));

    cout << "gens(I)   =   " << gens(I) << endl;
    vector<RingElem> GB = GBasis(I);

    SetVerbosityLevel(100);
    vector<RingElem> GB_m = F5_mat(I);
    vector<RingElem> GB_p = F5_poly(I);

    for (long i=0; i<len(GB_m); i++)  GB_m[i]=monic(GB_m[i]);
    GB_m = interreduced(GB_m);
    for (long i=0; i<len(GB_p); i++)  GB_p[i]=monic(GB_p[i]);
    GB_p = interreduced(GB_p);    
    
    if (EqSet(GB, GB_m))
      cout << "  --GB_m and GB are equal" << endl;
    else
    {
      cout << "LPP GB_m = " << endl;
      for (const auto& f:GB_m)   cout << LPP(f) << "  "  << endl;
      cout << endl;
      for (const auto& f:GB)
        cout << LPP(f) <<  "+(...)   \tvia GB_F5 reduces to " <<  NRLM(f, GB_m) << endl;
    }
    
    if (EqSet(GB, GB_p))
      cout << "  --GB_p and GB are equal" << endl;
    else
    {
      cout << "LPP GB_p = " << endl;
      for (const auto& f:GB_p)   cout << LPP(f) << "  "  << endl;
      cout << endl;
      for (const auto& f:GB)
        cout << LPP(f) <<  "+(...)   \tvia GB_F5 reduces to " <<  NRLM(f, GB_p) << endl;  
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
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
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
