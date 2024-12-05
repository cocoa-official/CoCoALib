//   Copyright (c) 2007 Alberto Arri

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/TmpF5.H"

#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/verbose.H"

#include <algorithm>
// using std::swap;
#include <vector>
#include <set>
#include <utility>
#include <set>
#include <map>
#include <cassert>
#include <iostream>
using namespace CoCoA;
using namespace std;

namespace CoCoA
{

  //-- cF5_t -----------------------------------
  class cF5_t
  {
  public: // public member fields
    int m;
    const PPMonoid& myPPM;
    SparsePolyRing env;
  public:
    cF5_t(const vector<RingElem> &I);

    
    struct module_term_t
    {
      PPMonoidElem term;
      int index;
      module_term_t(const PPMonoid &PPM, int i): term(PPM), index(i) {};
      module_term_t(ConstRefPPMonoidElem PPMel, int i): term(PPMel), index(i) {};
      bool operator < (const module_term_t & rhs) const //POS + TO
      {
        if (index != rhs.index) return index > rhs.index;
        return term < rhs.term;
      }
      module_term_t operator*(ConstRefPPMonoidElem rhs) const { return module_term_t(term * rhs, index); }
      bool operator==(const module_term_t &rhs) const { return term == rhs.term && index == rhs.index; }
    };
    

    //-- labeled_RingElem_t ----------------------------------
    class labeled_RingElem_t: public RingElem
    {
    public:
      labeled_RingElem_t(RingElem& e, int i): RingElem(e), label(PPM(owner(e)), i) {};
      labeled_RingElem_t(RingElem& e, module_term_t &mt): RingElem(e), label(mt) {};
      module_term_t label;
      RingElem & poly(){ return *this; };
    };


    //-- critpair_t --------------------------------------------
    struct critpair_t
    {
      PPMonoidElem t, u1, u2;// u,v
      int i1, i2; // k,l
      int myDeg(){ return StdDeg(t); }
      critpair_t(const PPMonoid& PPM): t(PPM), u1(PPM), u2(PPM) {};
    };


    //-- spols_order_t ------------------------------------------
    struct spols_order_t
    {
      const cF5_t &f;
      spols_order_t(const cF5_t &_f): f(_f) {};
      bool operator()(int i, int j) const
        {
          return f.r[i].label < f.r[j].label;
        }
    };

    //-- top_red_rv_t ---------------------------------------------
    struct top_red_rv_t
    {
      top_red_rv_t(): k(-1), t1(-1), t2(-1) {}
      int k, t1, t2;
    }; // K = {k}; T'=  {t1, t2} ( if k.t1.t2 != -1)


    //-- public member fields -------------------------------------------------
    vector<labeled_RingElem_t> r; //this just contains all the poly we will ever encounter
    vector<vector<pair<PPMonoidElem, int> > > Rules;
    vector<RingElem> F;
    vector<set<int> > G;
    int myNumRed0;
    void myCoreComputation(int step);//find a better name
  private:
    void spols(set<int, spols_order_t >& F, vector<critpair_t *>& B);
    critpair_t *critpair(int k, int l, int step);
    bool is_top_reducible(RefPPMonoidElem t, set<int>& Gb);
    
    void myAddRule(int j) {Rules[r[j].label.index].push_back(make_pair(r[j].label.term, j));}
    void myMonic(RingElem& p) {env->myDivByCoeff(raw(p), raw(LC(p)));}
    
    int myRewrite(ConstRefPPMonoidElem u, int k)
    {
      ConstRefPPMonoidElem v = r[k].label.term;
      int l = r[k].label.index;
      if (Rules[l].size()==0) return k;
      for (int j = Rules[l].size()-1; j >= 0; --j)
        if (IsDivisible( u*v ,  Rules[l][j].first))
          return Rules[l][j].second;
      return k;
    }

    
    bool is_rewritable(ConstRefPPMonoidElem u, int k)
    {
      int rv = myRewrite(u,k);
      return rv != k;
    }


    void myReduction(set<int, spols_order_t >& T, set<int>& Gp, vector<int>& D);

    template<typename T>
    int myFindReducer(int k, const T& Gp);

    top_red_rv_t top_reduction(int k, set<int> &Gp, vector<int>& Gpp);

    RingElem myNF(ConstRefRingElem p, set<int>& Gp);

  };    //-- cF5_t -----------------------------------


  cF5_t::cF5_t(const vector<RingElem> &I):
      myPPM(PPM(owner(I[0]))),
      env(owner(I[0])),
      F(I),
      myNumRed0(0)
       //this will make a copy of TidyGens(I) might be inefficient
  {
    m = F.size();
    r.reserve(m);
    Rules.resize(m);
    int i;
    for(i=0; i<m; ++i)
    {
      //    cout << F[i] <<  "----";
      myMonic(F[i]);
      //    cout << F[i] << endl;
      r.push_back(labeled_RingElem_t(F[i], i));
    }
    G.resize(m);
    G[m-1].insert(m-1);
  }

  
  RingElem cF5_t::myNF(ConstRefRingElem p, set<int>& Gp)
  {
    VerboseLog VERBOSE("cF5_t::myNF");
    RingElem res(p);
    VERBOSE(100) << "myNF " << p << " -> " ;
    bool red = true;
    while(red)
    {
      red = false;
      //for (set<int>::iterator it = Gp.begin(); it != Gp.end(); ++it)
      for (int j: Gp)
      {
        if (res == 0 )
        {
          VERBOSE(100) << endl;
          return res;
        }
        myMonic(res);
        PPMonoidElem oldHT = LPP(res);
        if (IsDivisible(oldHT, LPP(r[j])))
        {
          //	res.minus_eq(res.HC()/r[j].HC(),  oldHT / r[j].HT(), r[j]);
          res -= r[j] * monomial(env, LC(res)/LC(r[j]), oldHT / LPP(r[j]));
          //	  assert(res.HT()<oldHT);
          red = true;
          break;
        }
      }//end for
    }
    VERBOSE(100) << res << endl;
    return res;
  }


  cF5_t::top_red_rv_t cF5_t::top_reduction(int k, set<int>& Gp, vector<int>& Gpp)
  {
    top_red_rv_t rv;
    if (r[k] == 0)
    {
      cout << "Reduction to zero!!!!" << endl;
      ++myNumRed0;
      return rv;
    }

    RingElem p = r[k];
    int j1 = myFindReducer(k, Gp);
    int j2 = myFindReducer(k, Gpp);
    int j=0;
    if (j1 == -1 && j2 == -1)
    {
      myMonic(p);
      r[k].poly() = p;
      rv.k = k;
      return rv;
    }
    if (j1 != -1) j = j1;
    if (j2 != -1) j = j2;
    RingElem &q = r[j];
    PPMonoidElem u = LPP(p) / LPP(q);
    p -= monomial(env, LC(p)/LC(q), u) * q;
    if (r[j].label * u < r[k].label)
    {
      r[k].poly() = p;
      rv.t1 = k;
      return rv;
    }
    else
    {
      module_term_t tmp = r[j].label * u;
      r.push_back(labeled_RingElem_t(p, tmp));
      myAddRule(r.size()-1);
      rv.t1 = k;
      rv.t2 = r.size()-1;
      return rv;
    }
    return rv;
  }

  
  template<typename T>
  int cF5_t::myFindReducer(int k, const T& Gp)
  {
    PPMonoidElem  t = LPP(r[k]);
    for (typename T::const_iterator it = Gp.begin(); it != Gp.end(); ++it)
    {
      const int &j = *it;
      PPMonoidElem tp = LPP(r[j]);
      PPMonoidElem &vj = r[j].label.term;
      int kj = r[j].label.index;
      if ( IsDivisible(t, tp ))
      {
        PPMonoidElem u = t / tp;
        if ((r[j].label * u == r[k].label) ||
            is_rewritable(u,j) ||
            is_top_reducible(u*vj, G[kj]))
          continue;
        else return j;
      }
    }
    return -1;
  }


  void cF5_t::myReduction(set<int, spols_order_t>& T, set<int>& Gp, vector<int>& D)
  {
    //D = done; T = todo
    D.clear();
    while (!T.empty())
    {
      int k = *T.begin();
      T.erase(T.begin());
      RingElem h = myNF(r[k], Gp);
      r[k].poly() = h;
      top_red_rv_t trr = top_reduction(k, Gp, D);
      if (trr.k != -1) D.push_back(trr.k);
      if (trr.t1 != -1) T.insert(trr.t1);
      if (trr.t2 != -1) T.insert(trr.t2);
    }
  }


  //F is used as a return param
  void cF5_t::spols(set<int, spols_order_t>& F, vector<critpair_t *>& B)
  {
    VerboseLog VERBOSE("cF5_t::spols");
    for(vector<critpair_t *>::iterator it = B.begin(); it != B.end(); ++it)
    {
      critpair_t& cp = *(*it);
      ConstRefRingElem c1 = LC(r[cp.i1]);
      ConstRefRingElem c2 = LC(r[cp.i2]);
      RingElem s(r[cp.i1] * monomial(env, 1, cp.u1));
      //    poly_t s (cp.u1, r[cp.i1]); //s = cp.u1 * r[cp.i1]
      assert(c1 == 1 && c2 ==1);
      if (!is_rewritable(cp.u1, cp.i1) &&
          !is_rewritable(cp.u2, cp.i2))
      {
        //      s .minus_eq(c1/c2, cp.u2, r[cp.i2]);
        s -=  monomial(env, c1/c2, cp.u2) * r[cp.i2];
        if (s == 0) { cout << "**two identical polynomials found**" << endl; continue; }
        VERBOSE(100) << "S-poly( " << cp.i1 << ", " << cp.i2 << ") = " << s << endl;
        //	label_t tmp=r[cp.i2].label * cp.u1; //this what steger writes
        module_term_t tmp=r[cp.i1].label * cp.u1; //this is what Faugere writes and does work
        r.push_back(labeled_RingElem_t(s, tmp));
        myAddRule(r.size()-1);
        F.insert(r.size()-1); //F is automatically sorted in the correct way
      }
      else 
        VERBOSE(100) << "S-poly( " << cp.i1 << ", " << cp.i2
                     << ") rejected by myRewrite" << endl;
    }
  }


  bool cF5_t::is_top_reducible(RefPPMonoidElem t, set<int>& Gb)
  {
    set<int>::iterator it;
    for (it = Gb.begin(); it != Gb.end(); ++it)
      if (IsDivisible(t, LPP(r[*it]))) return true;
    return false;
  }


  cF5_t::critpair_t *cF5_t::critpair(int k, int l, int /*step*/)
  {
    VerboseLog VERBOSE("cF5_t::critpair");
    cF5_t::critpair_t *cp = new cF5_t::critpair_t(myPPM);
    VERBOSE(100) << "Considering crit pair " << k << " " << l << flush;
    cp->t = lcm(LPP(r[k]), LPP(r[l]));
    cp->u1 = cp->t / LPP(r[k]);
    cp->u2 = cp->t / LPP(r[l]);
    //
    if (r[k].label < r[l].label ) //swap the two poly to have the pair normalized
    {
      VERBOSE(100) << "  (Swapping indices)";
      std::swap(k, l);
      swap(cp->u1, cp->u2);
    }

    RefPPMonoidElem t1 = r[k].label.term;
    int k1 = r[k].label.index;
    PPMonoidElem t2 = r[l].label.term;
    int k2 = r[l].label.index;
    assert(k1<=k2);
    if ( (k1 != m-1 && is_top_reducible(cp->u1 * t1, G[k1+1])) ||
         (k2 != m-1 && is_top_reducible(cp->u2 * t2, G[k2+1])))
    {
      delete cp;
      VERBOSE(100) << "\t rejected" << endl;
      return nullptr;
    }
    cp->i1 = k;
    cp->i2 = l;
    VERBOSE(100) << "\t ok" << endl;
    return cp;
  }


  void cF5_t::myCoreComputation(int step)
  {
    VerboseLog VERBOSE("cF5_t::myCoreComputation");
    VERBOSE(100) << "Entering step " << step << endl;
    set<int> Gp = G[step+1];  //Gp = G' da steger
    Gp.insert(step);
    multimap<int, cF5_t::critpair_t *> P; //deg -> poly
    for(set<int>::iterator it = G[step+1].begin(); it != G[step+1].end(); ++it)
    {
      cF5_t::critpair_t *tmp = critpair(step, *it, step);
      if (tmp) P.insert(make_pair(tmp->myDeg(), tmp));
    }
    while (!P.empty())
    {
      int d = P.begin()->first;
      vector<critpair_t *> Pd;
      spols_order_t spo (*this);
      set<int, spols_order_t> Sd (spo);
      Pd.reserve(P.count(d));
      VERBOSE(100) << endl << "Processing " << P.count(d) << " poly of deg = " << d << endl;
      for(multimap<int, critpair_t *>::iterator it = P.begin(); it != P.end() && it->first == d; ++it)
        Pd.push_back(it->second); // Pd := { x \in P | deg(x) = d }
      P.erase(d); // P := P\Pd
      spols(Sd, Pd);
      if (VerbosityLevel()>=100)
      {
        cout << "myCoreComputation: Got spols:" << Sd.size() << endl;
        for (set<int, spols_order_t>::iterator iz = Sd.begin(); iz != Sd.end(); ++iz)
          cout << r[*iz] << endl;
        cout << endl ;
      }
      vector<int> Rd;
      Rd.clear();
      myReduction(Sd, Gp, Rd);
      for (size_t k = 0; k < Rd.size(); k++)
      {
        for (set<int>::iterator it = Gp.begin(); it != Gp.end(); ++it)
        {
          critpair_t *tmp = critpair(Rd[k], *it, step);
          VERBOSE(100) << "Trying new pair " << r[Rd[k]] << " | " << r[*it] << " \t -> ";
          if (VerbosityLevel()>=100)
          {
            if (tmp) cout << "ok" << endl;
            else cout << "rejected" << endl;
          }
          if (tmp) P.insert(make_pair(tmp->myDeg(), tmp));
        }
        Gp.insert(Rd[k]);
        VERBOSE(100) << "\t Adding new poly " << k << " = " << r[Rd[k]] << endl;
      }
      for (size_t k=0; k < Pd.size(); ++k)
        delete Pd[k];
      Pd.clear();
    }
    G[step] = Gp;
    VERBOSE(100) << "End of step " << step << " #G_" << step << " = " << Gp.size() << endl;
  }


  vector<RingElem> F5_poly(const vector<RingElem>& GensI)
  {
    VerboseLog VERBOSE("F5");
    VERBOSE(100) << "starting -- vector<RingElem> F5(const vector<RingElem>& GensI)" << endl;
    cF5_t F5i(GensI);
    int m = F5i.m;

    for(int i=m-2; i >= 0; --i)
    {
      F5i.myCoreComputation(i);
      //if 1 \in G_i -> exit
    }
    VERBOSE(100) << "F5 terminated; number of zero reductions = " << F5i.myNumRed0 << endl;
    vector<RingElem> GB;
    GB.reserve(F5i.G[0].size());
    //    for (set<int>::iterator it = F5i.G[0].begin(); it != F5i.G[0].end(); ++it)
    for (auto i:F5i.G[0])
      GB.push_back(F5i.r[i]);
    return GB;
  }


} // end of namespace CoCoA
