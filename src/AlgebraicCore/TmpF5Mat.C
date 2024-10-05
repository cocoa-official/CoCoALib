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

#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/PPMonoidOv.H"
#include "CoCoA/symbol.H"
#include "CoCoA/verbose.H"
#include "CoCoA/TmpF5.H"

#include <cmath>

#include <algorithm>
#include <utility>
#include <numeric>

#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>

#include <cassert>
#include <ctime>
#include <iostream>

#include <iterator>
#include <sstream>

#include <iomanip>


using namespace std;

namespace CoCoA{

  namespace F5ns{

    class PPdiv_t {
      typedef PPMonoidElem ppme;
      //    int nind;  // clang++ says this is unused
      vector<ppme> pp;
      vector<map<int, vector<int> > > indval;

    public:  
      PPdiv_t(int n):/*nind(n),*/indval(n){}
  
      void insert(ppme p){
        vector<long> x;
        exponents(x,p);
        pp.push_back(p);
        int cidx = pp.size() - 1;
        for(unsigned int i = 0; i < x.size(); i++){
          indval[i][x[i]].push_back(cidx);
        }
      }
  
    private:
  
      bool prj2grid(vector<long> &x)const
      {
        for (unsigned int i = 0; i < x.size(); i++){
          map<int, vector<int> >::const_iterator it = indval[i].lower_bound(x[i]);
          //if (it == indval[i].end()) return false;
          if (it != indval[i].end() && it->first == x[i]) continue;
          if (it == indval[i].begin()) return false;
          --it;
          x[i] = it->first;
        }
        return true;
      }
  
    public:
      long grid_size()const
      {
        long res = 1;
        for(unsigned int i = 0; i<indval.size(); i++) res *= indval[i].size();
        return res;
      }

      
      bool div(const ppme &p, int lc)const
      {
        vector<long> x;
        exponents(x,p);
        int o = x[lc];
        if (!prj2grid(x)) return false;    
        if (x[lc]!=o) return false;
        //    vector<int> &v = indval[lc][x[lc]];
        //    cout << "lc = " << lc << ", x[lc] = " << x[lc] << endl;
        map<int, vector<int> >::const_iterator it = indval[lc].find(x[lc]);
        assert( it != indval[lc].end() );
        const vector<int> &v = it->second;
        for (unsigned int j = 0; j < v.size(); j++)
        {
          if (IsDivisible(p, pp[v[j]])) return true;
        }
        return false;
      }
      
    };  // class PPdiv_t

  
    class PolyRingElemCmp
    {
    private:
      const SparsePolyRing& myP;
      const PPMonoid& myPPM;
    public:
      PolyRingElemCmp(const SparsePolyRing& P):myP(P),myPPM(PPM(P)) {}
      bool operator()(const RingElem &f, const RingElem &g) const {
        return myPPM->myCmp(raw(myP->myLPP(raw(f))), raw(myP->myLPP(raw(g))))<0;
      }
    };  // class PolyRingElemCmp


    F5opt_t opt;
  
    struct deg_cmp_t
    {
      bool operator()(const RingElem &x, const RingElem &y) const
      {return StdDeg(x)<StdDeg(y);}
    };

    PPMonoid *mtPPM = nullptr;

    PPMonoidElem ppe2ppm(const PPMonoidElem &pm, const PPMonoid &ppm = *mtPPM)
    {
      vector<long> x; 
      exponents(x,pm);
      return PPMonoidElem(ppm,x);
    }

    class module_term_t
    {
    public:
      PPMonoidElem term;

      int index;
      module_term_t(int i):term(*mtPPM), index(i){};
      module_term_t(ConstRefPPMonoidElem PPMel, int i):term(*mtPPM), index(i){
        term = ppe2ppm(PPMel);
      };
      bool operator < (const module_term_t & rhs) const { //REVPOS + TO
        if (index!=rhs.index) return index<rhs.index;
        return term<rhs.term;
      }  
      module_term_t operator * (const PPMonoidElem &rhs) const { return module_term_t( term * ppe2ppm(rhs), index); }
      bool operator  == (const module_term_t &rhs) const { return term == rhs.term && index == rhs.index; }
      friend ostream &operator <<( ostream &os, const module_term_t &mt);
    };  // class module_term_t

    
    class labeled_RingElem_t: public RingElem, public module_term_t
    {
    public:
      labeled_RingElem_t(ConstRefRingElem e, int i):RingElem(e), module_term_t(i) {};
      labeled_RingElem_t(ConstRefRingElem e, const module_term_t &mt):RingElem(e), module_term_t(mt) {};
   
      using module_term_t::operator <;
    }; //  class labeled_RingElem_t
    

    struct row_t
    { //inherit order from labels ie module_terms
      PPMonoidElem LT; //this is the LT of *re
      mutable bool touched;
      
    private:
      mutable RingElem *re; 
      void compute_re()const{}
    public:
      
      row_t(const row_t &r):LT(r.LT)
      {
        if (r.re) re = new RingElem(*r.re); else re = nullptr;
        touched=false;
      }

      row_t(const RingElem *_base, const RefPPMonoidElem& _coeff):
        LT(owner(_coeff)),re(nullptr)
      {
        SparsePolyRing env = owner(*_base);
        re = new RingElem (monomial(env,one(CoeffRing(env)),_coeff) * *_base);
        LT = LPP(*_base) * _coeff;
        touched=false;
      }
    
      ~row_t()      { if (re) delete re; }
    
      row_t &operator =(const row_t &r)
      {
        LT = r.LT;
        touched=r.touched;
        if (re) delete re;
        if (r.re) re = new RingElem(*r.re); else re = nullptr;
        return *this;
      }

      
      row_t &operator *=(const PPMonoidElem &ppme)
      {
        LT *= ppme;
        assert(re);
        SparsePolyRingPtr(owner(*re))->myMulByPP(raw(*re),raw(ppme));
        touched = false;
        return *this;
      }

      const RingElem &elem()const{ compute_re(); return *re; }
      RingElem &elem(){ compute_re(); return *re; }

      bool operator <(const row_t &rhs)const{ return LT<rhs.LT; }
      bool is_evaluated()const{ return true; }

      bool reduce(const row_t &r){ //returns false when result is Zero
        assert(LT==r.LT);
        touched = true;
        r.touched = true;
        SparsePolyRingPtr(owner(elem()))->myReductionStep(raw(*re),raw(r.elem()));
        if (!IsZero(elem())) {
          LT = LPP(elem());
          return true;
        } else return false;
      }
      friend ostream &operator <<( ostream &os, const row_t &r);
    };  // class row_t

    
    ostream &operator <<( ostream &os, const module_term_t &mt)
    { return os << "(" << mt.term << "," << mt.index << ")"; }

    ostream &operator <<( ostream &os, const row_t &r)
    { return os << r.elem(); }

    ostream &operator <<( ostream &os, const labeled_RingElem_t &r)
    { return os << (const module_term_t&)r << " - " << (const RingElem&) r; }

    struct cpair_t
    {
      PPMonoidElem lcm,u1,u2;
      set<labeled_RingElem_t>::iterator r1,r2;
      cpair_t(set<labeled_RingElem_t>::iterator &_r1, set<labeled_RingElem_t>::iterator &_r2):
        lcm( CoCoA::lcm(LPP(*_r1),LPP(*_r2))),
        u1(PPM(owner(*_r1))),
        u2(PPM(owner(*_r1))),
        r1(_r1),r2(_r2)
      {
        u1 = lcm / LPP(*r1);
        u2 = lcm / LPP(*r2);
      }

      module_term_t mt(int i)const
      { // returns the mt/label of ui*ri
        assert(i==1 || i==2);
        if (i==1) return static_cast<const module_term_t&>(*r1)*u1;
        return static_cast<const module_term_t&>(*r2)*u2;
      }

      module_term_t label()const { return max(mt(1),mt(2)); }

      long deg()const{ return StdDeg(lcm); }
    };

    ostream &operator << (ostream &os, const cpair_t &p){
      os << p.label() << "  " << p.lcm;
      return os;
    }

    typedef map<module_term_t,row_t> matrix_t;
    typedef matrix_t::value_type mmvt_t; //mmvt = map matrix value type = pair<module_term_t,row_t>
    typedef matrix_t::iterator mit_t;
  
    ostream &operator <<( ostream &os, const matrix_t &m){
      for (matrix_t::const_iterator it = m.begin(); it != m.end(); ++it)
      { os << it->first << ";\t " << it->second << endl; }
      return os;
    }

    class GB_t
    {
      unsigned long maxdeg;
    public:
      int nind;
      GB_t(int n/*, int m*/):maxdeg(0), nind(n), GB(),lastNRi(-1){};
      set<labeled_RingElem_t> GB;
      typedef set<labeled_RingElem_t>::iterator GBit_t;
      vector<set<PPMonoidElem> > lt_I_i; //lt_set[i] contains the leading terms of I_i

      vector<set<PPMonoidElem> > lt_Syz_i; //lt_set[i] contains the leading terms of Syz_i \ LT( I_i )

      vector<PPdiv_t> PPdiv;

      vector<vector<cpair_t> > pairs;

      void LPPinsert(const PPMonoidElem &PPM, unsigned int idx)
      {
        if (idx >= lt_I_i.size())
        { 
          lt_I_i.resize(idx+1);
          lt_Syz_i.resize(idx+1); //just to keep them the same len
        }
        if (idx >= PPdiv.size())
        {
          PPdiv.resize(idx+1,PPdiv_t(nind));
          if (idx!=0) PPdiv[idx] = PPdiv[idx-1];
        }
        lt_I_i[idx].insert(PPM);
        PPdiv[idx].insert(PPM);
      }

      
      void Syz_LPPinsert(const PPMonoidElem &PPM, unsigned int idx)
      {
        if (idx >= lt_Syz_i.size())
        { 
          lt_Syz_i.resize(idx+1); 
          lt_I_i.resize(idx+1); 
        }
        if (idx >= PPdiv.size())
        {
          PPdiv.resize(idx+1,PPdiv_t(nind));
          if (idx!=0) PPdiv[idx] = PPdiv[idx-1];
        }
        lt_Syz_i[idx].insert(PPM);
        PPdiv[idx].insert(PPM);
      }

      
      void myNR(RingElem& e, int i = 5000)
      {
        if ( lastNRi != i || (i==5000 && vecGB.size() != GB.size()) )
        {
          vecGB.clear();
          vecGB.reserve(GB.size());
          for (GBit_t it = GB.begin(); it != GB.end(); ++it)
            if (it->index <= i) vecGB.push_back(*it);
        }
        lastNRi = i;

        e = CoCoA::NR(e, vecGB);
        //      cout << " -> " << e << endl;
      }
    
    private:
      int lastNRi;
      vector<RingElem> vecGB;

      bool add_pair(GBit_t i1,GBit_t i2)
      {
        cpair_t pair(i1,i2);
        unsigned int deg = pair.deg();
        //      cout << "add_pair" << pair.mt(1) << "; " << pair.mt(2) << endl;
        if (pair.mt(1) == pair.mt(2) || is_syz_reducible(pair.mt(1)) || is_syz_reducible(pair.mt(2))) return false; //F5-pair criterion
        //      cout << "ok" << endl;
        if (deg >= pairs.size()) pairs.resize(deg+1);
        pairs[deg].push_back(pair);
        //      cout << endl << "Pair in deg = " << deg << endl;
        return true;
      }
      
      void pair_update(GBit_t &in)
      {
        unsigned int mp = max_pair();
        PPMonoidElem LPin = LPP(*in);
        for (GBit_t it = GB.begin(); it!= GB.end(); ++it)
          if (it!=in)
          {
            unsigned int ndeg = StdDeg(lcm(LPP(*it),LPin));
            if ( ndeg > mp ) 
              add_pair(it,in);
          }
      }
      
    public:
      unsigned int min_pair() const
      {
        for (unsigned int i = 0; i<pairs.size();i++) if (pairs[i].size()!=0) return i;
        return 0; //ie pairs is empty
      }
    
      unsigned int max_pair() const
      {
        if (pairs.size()==0) return 0;
        for (unsigned int i = pairs.size()-1; i>0; i--) if (pairs[i].size()!=0) return i;
        return 0;
      }

      GBit_t insert(ConstRefRingElem re, unsigned int pos)
      {
        GBit_t lit = GB.insert(labeled_RingElem_t(re,pos)).first;
        if (opt.GBLT2SYZLT) LPPinsert(ppe2ppm(LPP(re)),pos);
        if (static_cast<unsigned long>(StdDeg(re))>maxdeg) maxdeg = StdDeg(re);
        pair_update(lit);
        return lit;
      }
    
      GBit_t insert(const labeled_RingElem_t &lre)
      {
        pair<GBit_t, bool> in = GB.insert(lre);
        if (!in.second){//this is tricky.
          GB.erase(lre);
          in = GB.insert(lre);
          assert(in.second);
        }
        if (opt.GBLT2SYZLT) LPPinsert(ppe2ppm(LPP(lre)),lre.index);
        if (static_cast<unsigned long>(StdDeg(lre))>maxdeg) maxdeg = StdDeg(lre);
        pair_update(in.first);
        return in.first;
      }
    
      long max_deg() const { return maxdeg; }

      bool is_syz_reducible(const PPMonoidElem &PPE, int index)const
      {
        if (index<0) return false;

        assert(static_cast<unsigned int>(index) <= lt_I_i.size());
        if (static_cast<unsigned int>(index) < lt_Syz_i.size())
          for (set<PPMonoidElem>::const_iterator it = lt_Syz_i[index].begin(); it!= lt_Syz_i[index].end(); ++it) 
            if (IsDivisible(PPE,*it))  return true; 
      
        if (static_cast<unsigned int>(index) < lt_I_i.size())
          for (set<PPMonoidElem>::const_iterator it = lt_I_i[index].begin(); it!= lt_I_i[index].end(); ++it) 
            if (IsDivisible(PPE,*it)) return true;
      
        return false;
      }

      
      bool is_syz_reducible_gen(const module_term_t& mt, int lc)const{
        if (mt.index<1) return false;
        return PPdiv[mt.index-1].div(mt.term,lc);
      }

      
      bool is_syz_reducible(const module_term_t& mt)const{ return mt.index==0?false:is_syz_reducible(mt.term,mt.index-1);}

      
      const labeled_RingElem_t *is_reducible(const PPMonoidElem &PPE, int idx = 5000)const
      {
        for (GBit_t it = GB.begin(); it!= GB.end(); ++it)
          if (it->index <= idx && IsDivisible(PPE, LPP(*it))) return &*it; 
        return nullptr;
      }

      
      void LT_print(){
        for(unsigned int i = 0; i<lt_I_i.size(); i++){
          cout << "Set #" << i << " = { ";
          for (set<PPMonoidElem>::const_iterator it = lt_I_i[i].begin(); it!=lt_I_i[i].end(); ++it){
            cout << *it << ", ";
          }
          cout << " }" << endl;
        }
      }
    };  // class GB_t

    
    class matrF5_t
    {
    private:
      const PPMonoid& myPPM;
      SparsePolyRing env;
    public:
      unsigned int curr_deg,curr_poly,red2zero,n_indets,max_gen_deg;
      const vector<PPMonoidElem> &PPMindets;
      vector<PPMonoidElem> mtPPMindets;
      vector<RingElem> gens;
      GB_t GB;
      
      matrF5_t(const vector<RingElem> &I);
      void do_it();
    
    private:
      void do_poly_step();
      map<int, map<int, matrix_t> > matrix_hist;
      //    PPmap<RingElem *> fast_reductor;
      bool do_deg_step();

      mmvt_t * matrix_insert_row(matrix_t &mtr, module_term_t mt, row_t row)const;

      bool gen_macaulay_rows(matrix_t &mtr, const set<labeled_RingElem_t>::iterator &set_it, RefPPMonoidElem& ppme, 
                             unsigned int d, unsigned int ri = 0);
      void generate_macaulay( matrix_t &matrix, int mini = -1, int max1 = 5000 );

      void generate_macaulay_xx( matrix_t &matrix );

      struct gauss_stat_t{
        int r2z,sums;
        gauss_stat_t():r2z(0),sums(0){};
      };

      gauss_stat_t gauss( matrix_t &m);
      int pNF(RingElem &p, int idx);
      map<PPMonoidElem,int> red_frq;
    };  // class matrF5_t

    
    matrF5_t::matrF5_t(const vector<RingElem> &I): myPPM(PPM(owner(I[0]))),
                                                   env(owner(I[0])),
                                                   PPMindets(indets(PPM(owner(I[0])))),
                                                   gens(I),
                                                   GB(NumIndets(PPM(owner(I[0])))/*,I.size()*/)
    {
      red2zero = curr_deg = curr_poly = 0;
      n_indets = NumIndets(myPPM);

      mtPPM = new PPMonoid(NewPPMonoidOv(symbols(myPPM), StdDegRevLex)); // :(
      mtPPMindets = indets(*mtPPM);
      sort(gens.begin(),gens.end(),deg_cmp_t());
      for(unsigned int i = 0; i<gens.size(); i++)  assert(IsHomog(gens[i]));
      max_gen_deg = 0;
      for(unsigned int i=0;i<gens.size();i++) max_gen_deg = max<unsigned int>(max_gen_deg,StdDeg(gens[i]));
      GB.lt_I_i.reserve(gens.size());
      if (opt.incremental) { GB.insert(gens[0],0); /*fast_reductor[LPP(gens[0])] = &gens[0];*/ }
    }

    void matrF5_t::do_it()
    {
      //      VerboseLog VERBOSE("matrF5_t::do_it");
      if (VerbosityLevel()>=100)
        cout << "poly\t" << "deg\t"  << "Rows\t" << "CoeffLen\t"
             << "rsteps\t" << "Ctime\t" << "Gtime\t" << "R2Z\t" << "Sums\t" 
             << "Touched\t" << "New\t" << endl;
      if (opt.incremental)
      {
        curr_poly = 1;
        while(curr_poly < gens.size()) do_poly_step();
      }
      else
      {//non inc.
        GB.lt_I_i.resize(gens.size());
        GB.lt_Syz_i.resize(gens.size());
        //      GB.PPdiv.resize(gens.size(),PPdiv_t(GB.nind));
        for (unsigned int i=0; i<gens.size();i++)
        {
          RingElem pnf(env); //pnf = poly normal form
          pnf = gens[i];
          GB.myNR(pnf);
          GB.insert(pnf,i);
          gens[i] = pnf;
          //	fast_reductor[LPP(gens[i])] = &gens[i];
        }
        curr_poly = gens.size();
        curr_deg = 100;
        for (unsigned int i=0; i<gens.size();i++) curr_deg = min<int>(curr_deg,StdDeg(gens[i]));
        //      GB.pairs.resize(max_gen_deg+1);
        matrix_t &m = matrix_hist[curr_poly][curr_deg];

        generate_macaulay(m,0,curr_poly);
        curr_deg++;

        while(curr_deg <= max(GB.min_pair(),max_gen_deg))
        { 
          do_deg_step();
          if (VerbosityLevel()>=100) cout << "GBmin = " << GB.min_pair() << endl;
          if (GB.min_pair()==0) break;
        }
      }
    }

  
    void matrF5_t::do_poly_step(){
      RingElem pnf(env); //pnf = poly normal form
      do
      {
        pnf = gens[curr_poly];
        GB.myNR(pnf);
        if (IsZero(pnf))
        { 
          gens.erase(gens.begin()+curr_poly); 
          if (curr_poly == gens.size()) return;
        }
      }
      while(IsZero(pnf));
      GB.insert(pnf,curr_poly);
      curr_deg = StdDeg(pnf);
      matrix_t &M = matrix_hist[curr_poly][curr_deg];
      generate_macaulay(M, opt.skip_rows ? curr_poly:0, curr_poly);
      curr_deg++;
      while (GB.min_pair()!=0)  if (!do_deg_step()) break;
      matrix_hist[curr_poly].clear();
      curr_poly++;
    }

  
    struct matrixrow_sorter /*: public binary_function<bool,mmvt_t *,mmvt_t *>*/{ 
      const SparsePolyRing &env;
      const PPMonoid& myPPM;
      matrixrow_sorter(const SparsePolyRing &e):env(e),myPPM(PPM(e)){}
      bool operator () (const mmvt_t *x,const mmvt_t *y){
        int cmp = myPPM->myCmp(raw(x->second.LT),raw(y->second.LT));
        //      if (x->second.LT != y->second.LT) return x->second.LT < y->second.LT;
        if (cmp!=0) return cmp<0;
        return !(x->first < y->first);
      }
    };

  
    int matrF5_t::pNF(RingElem &p, int idx){
      const labeled_RingElem_t *rrle = nullptr;
      int loop_cnt = 0;
      PPMonoidElem LPPp = LPP(p);

      for(;;)
      {
        rrle = GB.is_reducible(LPPp, idx);
        if (!rrle) return loop_cnt;
        env->myReductionStep(raw(p), raw(*rrle));
        if (IsZero(p)) return loop_cnt;
        LPPp = LPP(p);
        loop_cnt++;
      }

      //     while ( (!IsZero(p)) && (rrle = GB.is_reducible(LPPp, idx) )) { 
      //       bool didit = false;
      //       //red_frq[LPP(p)]++;
      //       if (opt.prev_red){
      // 	int idl = rrle->index;
      // 	for(; idl>=0; idl--) {
      // 	  map<int, matrix_t>::iterator it = matrix_hist[idl].find(curr_deg);
      // 	  if (it != matrix_hist[idl].end() ) {
      // 	    for(mit_t mit = it->second.begin(); mit!=it->second.end(); ++mit) if (!IsZero(mit->second.elem())){
      // 	      if (mit->second.LT == LPPp) { 
      // 		env->myReductionStep(raw(p),raw(mit->second.elem())); 
      // 		if (IsZero(p)) return loop_cnt;
      // 		LPPp = LPP(p);
      // 		didit = true; 
      // 		break; 
      // 	      }
      // 	    }
      // 	    if (didit) break;
      // 	  }
      // 	}
      //       }
      //       if (!didit) env->myReductionStep(raw(p),raw(*rrle));
      //       if (IsZero(p)) return loop_cnt;
      //       LPPp = LPP(p);
      //       loop_cnt++; 
      //     }
      return loop_cnt;
    }

    matrF5_t::gauss_stat_t matrF5_t::gauss( matrix_t &m)
    {
      VerboseLog VERBOSE("gauss::do_it");
      gauss_stat_t gs;
      if (m.size() < 1) { VERBOSE(100) << "gauss: empty matrix"; return gs; }
      matrixrow_sorter mrs(env);
      priority_queue< mmvt_t *, vector<mmvt_t*>, matrixrow_sorter> pq(mrs);
      int realred = 0, rsteps = 0;
      long ml  = 0;
      if (VerbosityLevel()>=100) cout << "\t&" << ml;
      for (mit_t it = m.begin(); it!=m.end(); ++it)
      { 
        if (opt.skip_rows)
        {
          int loop_cnt=0;
          if (opt.use_NR)
            GB.myNR(it->second.elem(),curr_poly - 1);
          else
            loop_cnt=pNF(it->second.elem(),curr_poly - 1);
          if (loop_cnt!=0 || opt.use_NR)
          {
            it->second.touched = true;
            realred++;
            rsteps+=loop_cnt;
            if (IsZero(it->second.elem())) { gs.r2z++; continue; }
            it->second.LT = LPP(it->second.elem());
          }
        }
        pq.push(&*it); 
      }
      //    cout << "\t" << realred;
      if (pq.size() <= 1) {  VERBOSE(100) << "one liner"; return gs; }
      do
      {
        mmvt_t *i = pq.top();
        pq.pop();
        PPMonoidElem cLT = i->second.LT;
        //???SET BUT NOT USED (JAA,2013-03-13)      bool did_loop = false;
        while ( ( !pq.empty()) && pq.top()->second.LT == cLT ){
          //???SET BUT NOT USED (JAA,2013-03-13)        did_loop = true;
          mmvt_t *j = pq.top();
          //cout << endl << "P1 = " << i->second.elem() << ";" << endl << "P2 = " << j->second.elem() << endl;
          pq.pop(); //	cout << i->first << "   " << j->first << endl;
          assert(i->first < j->first); //*i has smaller label
          gs.sums++;
          if (j->second.reduce(i->second)){ //*i won't be changed
            assert(j->second.LT < i->second.LT);
            if (opt.skip_rows && !opt.use_NR) rsteps += pNF(j->second.elem(),curr_poly - 1);
            if (!IsZero(j->second.elem())) {
              j->second.LT = LPP(j->second.elem());
              assert(j->second.elem()!=0);
              pq.push(j);
            }
            else
            {  gs.r2z++; }
          }
          else
          {  gs.r2z++; } //we had a reduction to zero :(
        }
      }
      while(! pq.empty());
      if (VerbosityLevel()>=100) cout << "\t&" << rsteps + realred; 
      return gs;
    }

  
    mmvt_t * matrF5_t::matrix_insert_row(matrix_t &mtr, module_term_t mt, row_t row)const
    {
      mit_t it = mtr.find(mt);
    
      PPMonoidElem rwc(myPPM);
      if (it == mtr.end())
      { //new row
        it = mtr.insert(make_pair(mt,row)).first;
        return &*it;
      }
      else return nullptr;
      return &*it;
    }

  
    bool matrF5_t::gen_macaulay_rows(matrix_t &mtr, const set<labeled_RingElem_t>::iterator &set_it, 
                                     RefPPMonoidElem& ppme, 
                                     unsigned int d, unsigned int ri)
    {  
      if (GB.is_syz_reducible(ppme,set_it->index - 1))
      { 
        //      cout << "F5 criteria killed " << ppme << " idx " << set_it->index << endl;
        return false;
      }
      if (d == 0)
      { //try adding row
        module_term_t mt(ppme, set_it->index);
        row_t row(&*set_it,ppe2ppm(ppme/set_it->term,myPPM));
        return matrix_insert_row(mtr,mt,row);
      }
      unsigned int i=d;
      if (ri != n_indets - 1) {
        if (gen_macaulay_rows(mtr,set_it,ppme,d,ri+1)) 
          for (i = 1; i<=d; i++) {
            ppme *= mtPPMindets[ri];
            if (!gen_macaulay_rows(mtr,set_it,ppme,d-i,ri+1)) break;
          }
      } else {
        ppme *= power(mtPPMindets[ri],d);
        gen_macaulay_rows(mtr,set_it,ppme,0,ri+1);
        i=d;
      }
      ppme /= power(mtPPMindets[ri], i>=d ? d : i);
      return true;
    }

    void matrF5_t::generate_macaulay( matrix_t &matrix, int mini, int maxi){ 
      if (GB.pairs.size()>=curr_deg) GB.pairs[curr_deg].clear();
      for(GB_t::GBit_t it = GB.GB.begin(); it != GB.GB.end(); ++it) if (mini <= it->index && it->index <= maxi){
          int delta_deg = curr_deg - StdDeg(*it);
          PPMonoidElem wt = it->term;
          if (delta_deg >= 0) gen_macaulay_rows(matrix,it,wt,delta_deg);
        }
    }

    void matrF5_t::generate_macaulay_xx( matrix_t &matrix)
    { 
      VerboseLog VERBOSE("matrF5_t::generate_macaulay_xx");
      if (GB.pairs.size()>=curr_deg) GB.pairs[curr_deg].clear();
      matrix_t &oldm = matrix_hist[curr_poly][curr_deg-1];
      matrix.clear();
      if (!opt.incremental)
      {
        for(unsigned int i = 0; i < gens.size(); ++i )
          if ((!IsZero(gens[i])) && static_cast<unsigned long>(StdDeg(gens[i])) == curr_deg) {
            if (!IsZero(gens[i])){
              //	  cout << "New gen " << i << endl;
              GB_t::GBit_t it = GB.insert(gens[i],i);
              matrix.insert(mmvt_t(*it,row_t(&*it,PPMonoidElem(myPPM))));
            }
            else
              VERBOSE(100) << "Gen #" << i << " is useless" << endl;
          }
      }
    
      RingElem TheOne(one(env));
      PPMonoidElem TheId(one(myPPM));
      row_t TheOtherOne(&TheOne,TheId);
    
      for (mit_t it = oldm.begin(); it!=oldm.end(); ++it) if (!IsZero(it->second.elem())) {
          for (unsigned int i=0; i<n_indets; i++){
            module_term_t mt = it->first;
            assert(static_cast<unsigned long>(StdDeg(it->second.LT)) + 1 == curr_deg);
            mt.term *= mtPPMindets[i];
            if (!GB.is_syz_reducible_gen(mt,i)){
              mit_t ip = matrix.find(mt);
              if (ip==matrix.end()) { //new row/label
                ip = matrix.insert(mmvt_t(mt,TheOtherOne )).first;
                //ip = matrix.insert(mmvt_t(mt,it->second )).first;
                ip->second = it->second;
                ip->second *= PPMindets[i];
              }else {
                //if (1 + StdDeg(it->second.coeff) <= StdDeg(ip->second.coeff)){
                if (PPMindets[i]*it->second.LT < ip->second.LT){
                  ip->second = it->second;
                  ip->second *= PPMindets[i];
                  //	    cout << "** delta = " << 1 + StdDeg(it->second.coeff) << "; Old lt = " << ip->second.LT 
                  //	 << "; NewLT = " << nrw.LT << endl;
                } //else
              }
              // 	} else {
              // 	  cout << mt << i << "  killed" << endl;
            }//if us syz etc..
          }//end for indet
        }//end for it    
      oldm.clear();
    }

    
    bool matrF5_t::do_deg_step()
    {
      //      VerboseLog VERBOSE("matrF5_t::do_deg_step");
      matrix_t &matrix = matrix_hist[curr_poly][curr_deg];
      if (curr_deg==0) return false;
      if (VerbosityLevel()>=100) cout << curr_poly << "\t&" << setw(2) << curr_deg; //<< "\t" << GB.pairs[curr_deg].size();

      //    copy(GB.pairs[curr_deg].begin(),GB.pairs[curr_deg].end(),ostream_iterator<cpair_t>(cout, "\n"));
      clock_t timeb4 = clock(),timeafter;
      double gtime = 0;
      generate_macaulay_xx(matrix);
      if (VerbosityLevel()>=100)
      {
        timeafter = clock();
        gtime = static_cast<double>(timeafter - timeb4) / CLOCKS_PER_SEC;
        cout << "\t&" << setw(5) << matrix.size() << flush;
        timeb4 = clock();
      }
      gauss_stat_t gs = gauss(matrix);
      if (VerbosityLevel()>=100)
      {
        timeafter = clock();
        cout << "\t&" << gtime << "\t& " << static_cast<double>(timeafter - timeb4) / CLOCKS_PER_SEC ;
      }
      int newel=0,touched=0;
      if (VerbosityLevel()>=100) cout << " \t&" << gs.r2z << " \t&" << setw(5) << gs.sums << flush;
      //        cout << endl << "After Gauss:" << endl << matrix;

      for( mit_t it = matrix.begin(); it != matrix.end(); ++it) 
        if (it->second.is_evaluated()) {
          if (it->second.touched) touched++;
          RingElem &rre = it->second.elem();
          if (rre == 0) { 
            GB.Syz_LPPinsert(it->first.term,it->first.index-1); // *non Faugere extension*
            continue;
          }
          if (!GB.is_reducible(LPP(rre))) {
            //	  cout << "New poly :" << it->first << " : " << LPP(rre) << endl;
            newel++;
            GB.insert(labeled_RingElem_t(rre,it->first));
          }
          //	fast_reductor[LPP(rre)] = &it->second.elem();
        }
      if (VerbosityLevel()>=100) cout << "\t& " << setw(5) << touched << "\t& " << newel << "\\\\" << endl;
      curr_deg++;
      //        if (touched==0) return false; //?? why ??
      return true;
    }


    // void interreduce(vector<RingElem> &v)
    // {
    //   SparsePolyRing env = owner(v[0]);
    //   for(unsigned int i = 0; i < v.size(); i++)
    //   {
    //     vector<RingElem> rl;
    //     for(unsigned int j = 0; j < v.size(); j++)
    //       if (i!=j && !IsZero(v[j])) rl.push_back(v[j]);
    //     v[i] = NR(v[i],rl);
    //   }
    //   vector<RingElem>::iterator nend = remove(v.begin(),v.end(),zero(env));
    //   v.erase(nend,v.end());
    //   for(unsigned int i = 0; i < v.size(); i++) v[i] /= monomial(env,LC(v[i]),one(PPM(env)));
    // }

    

  }//end of namespace F5ns

  
  vector<RingElem> F5_mat(const vector<RingElem>& GensI, const F5ns::F5opt_t &/*F5opt*/ )
  {
    VerboseLog VERBOSE("F5");
    VERBOSE(100) << "starting" << endl;
    vector<RingElem> GB;
    
    F5ns::matrF5_t F5i(GensI);

    F5i.do_it();
    GB.clear();
    GB.reserve(F5i.GB.GB.size());
    VERBOSE(100) << "GB size = " << F5i.GB.GB.size() << endl;
    //    copy(F5i.GB.GB.begin(),F5i.GB.GB.end(),back_inserter(GB));
    GB.insert(GB.end(), F5i.GB.GB.begin(), F5i.GB.GB.end());
    
    // //check answer
    // if (F5ns::opt.checkGB)
    // {
    //   cout << "Checking GB .... " << flush;
    //   F5ns::interreduce(GB);
    //   cout << (F5ns::GBtest(GB) ? "ok" : "not ok") << endl;
    // }
    return GB;
  }

}

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpF5Mat.C,v 1.23 2024/07/02 12:39:58 bigatti Exp $
// $Log: TmpF5Mat.C,v $
// Revision 1.23  2024/07/02 12:39:58  bigatti
// Summary: Making F5 more readable with CoCoALib conventions
//
// Revision 1.22  2023/03/27 13:59:46  abbott
// Summary: Commented out unused var (to silence compiler warning)
//
// Revision 1.21  2023/03/13 19:51:56  abbott
// Summary: Commented out needless binary_function<...>
//
// Revision 1.20  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.19  2022/02/08 20:21:37  abbott
// Summary: Commented out unused param in GB_t (to avoid clang warning)
//
// Revision 1.18  2022/02/05 20:01:33  abbott
// Summary: Minor changes to avoid warnings from clang (redmine 1528)
//
// Revision 1.17  2019/03/19 11:07:08  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.16  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.15  2018/03/15 14:39:20  bigatti
// -- new file SparsePolyOps-ideal.H
//
// Revision 1.14  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.13  2014/07/07 14:12:17  abbott
// Summary: Corrected silly typo
// Author: JAA
//
// Revision 1.12  2014/07/07 12:51:09  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.11  2014/04/30 16:14:12  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.10  2014/01/16 16:13:57  abbott
// Replaced  RingElem(env,0)  by  zero(env)
//
// Revision 1.9  2013/03/15 17:49:57  abbott
// Commented out unused variable.
//
// Revision 1.8  2012/10/16 10:15:43  abbott
// Replaced  RefRingElem  by  RingElem&  (in an old style cast?!?)
//
// Revision 1.7  2012/01/26 16:47:49  bigatti
// -- changed back_inserter into insert
//
// Revision 1.6  2012/01/26 16:37:10  abbott
// Removed an apparently totally pointless CPP symbol definition (NDEBUG).
//
// Revision 1.5  2008/06/30 17:15:41  abbott
// Used const_iterator instead of iterator (where appropriate).
//
// Revision 1.4  2008/05/30 12:50:31  abbott
// Comment out some unused formal parameters (to avoid compiler warnings)
//
// Revision 1.3  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.2  2007/12/11 14:39:25  bigatti
// --
//
// Revision 1.1  2007/11/20 10:01:26  bigatti
// -- change: TmpF5.C --> TmpF5Mat.C  (by Alberto Arri)
// -- updated and improved test-F5.C
//
