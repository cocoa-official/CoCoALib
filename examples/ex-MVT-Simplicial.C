// Copyright (c) 2019  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

// Authors: Anna M. Bigatti, Eduardo Saenz-De-Cabezon

#include "CoCoA/library.H"
#include <fstream>
#include <string>
using std::string;

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program contains the code described in ISSAC 2019 paper: \n"
  "  Bigatti, Heras, Saenz-De-Cabezon\n"
  "  \"Monomial resolutions for efficient computation of simplicial homology\"\n"
  "  https://doi.org/10.1145/3326229.3326266";


const string LongDescription = 
  "The paper examples are available at:\n"
  "  http://www.dima.unige.it/~bigatti/cocoa/ex-MVT-Simplicial-tests.zip\n"
  "This program reads from stdin.  For your convenience run it as:\n"
  "  ./ex-MVT-Simplicial < ex-MVT-Simplicial-tests/<...>\n";


//----------------------------------------------------------------------

namespace CoCoA
{

  //----------------------------------------------------------------------------
  bool DEBUG = false;

  typedef std::vector<DynamicBitset> PPVectorSqFr;
  typedef std::set<int> intset;
  //----------------------------------------------------------------------------
  // Taylor  
  //  void InterreduceSort(PPVectorSqFr& g);
  bool DBLessThan(const DynamicBitset& f, const DynamicBitset& g);
  bool divides00(const PPVectorSqFr& g, DynamicBitset b);
  DynamicBitset listLCM(list<int> s, const PPVectorSqFr& g, const DynamicBitset& l);
  int binomial(int n, int k);
  //----------------------------------------------------------------------------
  // MVT-Simplicial

  class MVTSquareFree
  {
  public:
    explicit MVTSquareFree(const ideal& I);
    void myComputeBetti(std::vector<int>& bettis);
    std::ostream& myOutput(std::ostream& out) const;

  private:
    explicit MVTSquareFree(const PPVectorSqFr& I, int n);
    void myComputeRec(PPVectorSqFr& I, const DynamicBitset& AccPivot, const int dim, const BigInt& pos);
    void myInsert(const int&, const BigInt&);
    void myInsertGens(const PPVectorSqFr& I, const int dim, const BigInt& pos);
    typedef std::list<BigInt> positions;
    
  private: // data members
    std::vector<positions> myDimVec;
    bool IamEmpty;
  };


  //*******************************************************//
  void InterreduceSort(PPVectorSqFr& g);
  void tilde(PPVectorSqFr& g, DynamicBitset& AccPivot, PPVectorSqFr& I);
  bool hasFullSupport(const PPVectorSqFr& g, const DynamicBitset& AccPivot);

  

  //********************************************************//
  // FUNCTIONS and MEMBER FUNCTIONS for MVTSquareFree class //
  //********************************************************//
  
  std::ostream& operator<<(std::ostream& out, const MVTSquareFree& M)
  { return M.myOutput(out); }

  //-- public member functions -------------------------------------

  //-- public constructor
  MVTSquareFree::MVTSquareFree(const ideal& I):
    myDimVec(len(gens(I)))
  {
    IamEmpty = true;
    PPMonoid PPMI = PPM(RingOf(I));
    const std::vector<RingElem>& gensI = gens(I);
    //    if(deg(product(gensI))<NumIndets(RingOf(I))) CoCoA_THROW_ERROR("the ideal has not full support","MVTSquareFree Simplicial"); 
    //      else
    PPVectorSqFr J;
    J.reserve(len(gensI));
    for (long i=0; i<len(gensI); ++i)
	  {
	    // ANNA: should be a CoCoA_ASSERT
	    if (!IsMonomial(gensI[i])) CoCoA_THROW_ERROR("is not monomial","MVTSquareFree");
	    if (!IsOne(LC(gensI[i]))) CoCoA_THROW_ERROR("coeff is not 1","MVTSquareFree");
	    //if (deg(LPP(gensI[i]))>2) CoCoA_THROW_ERROR("deg>2","MVTSquareFree");
	    if (!IsSqFree(LPP(gensI[i]))) CoCoA_THROW_ERROR("not squarefree","MVTSquareFree");
	    J.push_back(DynamicBitset(LPP(gensI[i])));
	  }
    InterreduceSort(J);
    //cout<<"The ideal: "<<J<<endl;
    DynamicBitset AccumulatingPivot(NumIndets(PPMI));
    if(!hasFullSupport(J, AccumulatingPivot))
      CoCoA_THROW_ERROR("the ideal has not full support","MVTSquareFree Simplicial"); 
    BigInt pos(1);
    if (len(gens(I))==1) myInsertGens(J, 0, pos);
    myComputeRec(J, AccumulatingPivot, 0, pos);
  }


  std::ostream& MVTSquareFree::myOutput(std::ostream& out) const
  {
    if (IamEmpty) out << " [empty]"; else out << endl;
    for (long i=0; i<len(myDimVec); ++i)
    {
      const positions& m=myDimVec[i];
      if (!m.empty())
      {
        out << "tree in dim = " << i << endl;
        for (positions::const_iterator it=m.begin();it!=m.end();++it)
        {
          //          out<<"  Multidegree "<< it->first << " (deg="<<count(it->first)<<"):";
          out<<"  pos "<< *it << endl;
        }
      }
    }
    return out;
  }


  void MVTSquareFree::myComputeBetti(vector<int>& bettis)
  {
    for (long i=0; i<len(myDimVec); ++i)
      bettis[i] = (myDimVec[i]).size();
  }
  

  //-- private member functions -------------------------------------

  void MVTSquareFree::myInsert(const int& dim, const BigInt& pos)
  {
    IamEmpty = false;
    //    (myDimVec[dim]).insert(pos);
    (myDimVec[dim]).push_back(pos);
  }    


  void MVTSquareFree::myInsertGens(const PPVectorSqFr& I, const int dim, const BigInt& pos)
  {
    for (PPVectorSqFr::const_iterator it=I.begin(); it!=I.end(); ++it)
      myInsert(dim, pos);
  }
  

  void MVTSquareFree::myComputeRec(PPVectorSqFr& I, const DynamicBitset& AccPivot, const int dim, const BigInt& pos_orig)
  {
    BigInt pos = pos_orig;
    //if (dim!=0) myInsertGens(I, dim, pos);
    //cout<<"pos: "<<pos<<" dim: "<<dim<<" ideal: "<<I<<endl;
    while (len(I)>2 && hasFullSupport(I,AccPivot))
    {
      //cout << " -- MVT     pos = " << pos   << "  I = " << I << " has fullsupport= "<<hasFullSupport(I,AccPivot)<<endl;
      PPVectorSqFr TildeI;
      DynamicBitset TildeAccPivot = AccPivot;
      tilde(TildeI, TildeAccPivot, I); // modifies I
      // vertical branch
      myComputeRec(TildeI, TildeAccPivot, 1+dim, 2*pos);
      // horizontal branch
      pos = 2*pos+1;
    }  // while
    if (len(I)==2 && hasFullSupport(I,AccPivot)) myInsert(1+dim, 2*pos);
    if (len(I)==1 && hasFullSupport(I,AccPivot)) myInsert(dim, pos);

  }

  //********************************************************//
  //********************************************************//

  //  namespace // anonymous for functions local to this file/compilation unit.
  //  {
   
  bool DBLessThan(const DynamicBitset& f, const DynamicBitset& g)
  { return f<g; }


  bool divides00(const PPVectorSqFr& g, DynamicBitset b)
  {
    for (PPVectorSqFr::const_iterator it=g.begin(); it!=g.end(); ++it)
      if (IsSubset(*it,b)) return true;
    return false;
  }

  bool hasFullSupport(const PPVectorSqFr& g, const DynamicBitset& AccPivot)
  {
    DynamicBitset s(len(g[0]));
    s=AccPivot;
    for (PPVectorSqFr::const_iterator it=g.begin(); it!=g.end(); ++it)
      s=s|*it;
    return(s.IamAll1s());
  }

  void InterreduceSort(PPVectorSqFr& g)
  {
    PPVectorSqFr newg;
    newg.reserve(len(g));
    sort(g.begin(), g.end(), DBLessThan);
    for (PPVectorSqFr::const_iterator it=g.begin(); it!=g.end(); ++it)
      if (!divides00(newg, *it)) newg.push_back(*it);
    swap(newg,g);
  }


  void tilde1(PPVectorSqFr& g, const PPVectorSqFr& I, const DynamicBitset& pivot)
  {
    long count=0;
    for (long i=0; i<len(I); ++i)
      if (IsDisjoint(I[i], pivot))
        g.push_back(I[i]);
      else
      {
        g.push_back(I[i]-pivot);
        ++count;
      }
    if (count!=0)
      InterreduceSort(g);
  }
  

  void tilde2(PPVectorSqFr& g, const PPVectorSqFr& I, const DynamicBitset& pivot)
  {
    PPVectorSqFr Coprime_g;
    PPVectorSqFr Coprime_g_reduced;
    for (long i=0; i<len(I); ++i)
      if (IsDisjoint(I[i], pivot))
        Coprime_g.push_back(I[i]);
      else
        g.push_back(I[i]-pivot);
    InterreduceSort(g);
    for (long i=0; i<len(Coprime_g); ++i)
      if (!divides00(g, Coprime_g[i]))
        Coprime_g_reduced.push_back(Coprime_g[i]);    
    for (long i=0; i<len(Coprime_g_reduced); ++i)
      g.push_back(Coprime_g_reduced[i]);
    sort(g.begin(), g.end(), DBLessThan);
  }


  void tilde(PPVectorSqFr& g, DynamicBitset& AccPivot, PPVectorSqFr& I)
  {
    if (!g.empty()) CoCoA_THROW_ERROR("g: not empty", "tilde");
    g.reserve(len(I));
    DynamicBitset pivot = I.back();
    I.pop_back();
    if (DEBUG) cout << "tilde: pivot = " << pivot << endl;
    AccPivot |= pivot;
    if (len(I) > 20)
      tilde2(g, I, pivot);
    else
      tilde1(g, I, pivot);
    if (DEBUG) cout << "tilde: g = " << g << endl;
  }
  


  //----------------------------------------------------------------------------
  // Taylor  

  //   void InterreduceSort(PPVectorSqFr& g)
  //   {
  //     PPVectorSqFr newg;
  //     newg.reserve(len(g));
  //     sort(g.begin(), g.end(), DBLessThan);
  //     for (PPVectorSqFr::const_iterator it=g.begin(); it!=g.end(); ++it)
  //       if (!divides00(newg, *it)) newg.push_back(*it);
  //     swap(newg,g);
  //   }

  //   bool DBLessThan(const DynamicBitset& f, const DynamicBitset& g)
  //   { return f<g; }


  //   bool divides00(const PPVectorSqFr& g, DynamicBitset b)
  //   {
  //     for (PPVectorSqFr::const_iterator it=g.begin(); it!=g.end(); ++it)
  //       if (IsSubset(*it,b)) return true;
  //     return false;
  //   }

  DynamicBitset listLCM(list<int> s, const PPVectorSqFr& g, DynamicBitset& l)// the lcm of the generators whose indices are in the list s
  {
    DynamicBitset newl(len(*(g.begin())));
    newl=l;
    for (list<int>::const_iterator it=s.begin(); it!=s.end(); ++it)
    {
      newl=newl|g[*it];
    }
    return(newl); 
  }

  int binomial(int n, int k)
  {
    if(k >  n)
      return 0;
    if(k == 0)
      return 1;
    if(k > n/2)
      return binomial(n,n-k);
    return (n * binomial(n-1,k-1)/k);
  }

  void print( list<int> l)//only for printing lists of integers
  {
    for(list<int>::iterator it=l.begin(); it!=l.end() ; ++it)
      cout << " " << *it;
    cout<<endl;
  }

  void subset(vector<int> arr, int size, int left, int index, list<int> &l, vector<list<int> > &bl)//to compute k-subsets of a range 1..r
  {
    if(left==0)
    {
      // print(l);
      bl.push_back(l);
      return;
    }
    for(int i=index; i<size;i++)
    {
      l.push_back(arr[i]);
      subset(arr,size,left-1,i+1,l,bl);
      l.pop_back();
    }
  }  


  int orderSubset(list<int>& s, int n, int k)
  {
    int r=binomial(n,k);
    int m=0;
    list<int>::iterator it=s.end();
    for(m=0;m<k;m++)
    {
      it--;
      r=r-binomial(n-(*it)-1,m+1);
    }
    return r;
  }


  int orderSubsetExceptOne(list<int>& s, int n, int k, int e)
  {
    int r=binomial(n,k-1);
    int m=0;
    list<int>::iterator it=s.end();
    for(int mm=0;mm<k;mm++)
    {
      it--;
      if(mm==k-e-1)
      {}
      else {r=r-binomial(n-(*it)-1,m+1);m++;}
    }
    return r;
  }


  int oneStrand(vector<vector<list<int> > >& st, const PPVectorSqFr& J, PPMonoid PPMI, long lo, long hi)
  // computes the skeleton of the Taylor 1-strand of sqfree ideal J
  // it cuts computation if at some point every k-set LCM is 1
  {
    vector<int> range;
    int r=J.size();
    for(int i=0; i<r ; ++i) range.push_back(i);
    int strandsize=0;
    list<int> lt;
    vector<list<int> > biglt;
    bool complete=0;
    int i=lo;
    while (i<=hi+1 && !complete)
    {
      subset(range,r,i,0,lt,biglt);
      DynamicBitset l(len(J[0]));
      DynamicBitset k(len(J[0]));    
      l=DynamicBitset(NumIndets(PPMI));
      vector<list<int> > thisList;
      for(std::vector<list<int> >::iterator it=biglt.begin();it!=biglt.end(); ++it)
      {  
        k=listLCM(*(it),J,l);
        if (k.IamAll1s())
        {
          thisList.push_back(*(it));
          if(thisList.size()==binomial(r,i)) //once every k-set LCM is 1, all the next ones will be 1 too and we don't need to compute them
          {
            //cout<<"Complete at i= "<<i<<endl;
            complete=1;
          };
          //print(*it);
          //cout<<k<<endl;
          strandsize++;
        }
      }
      biglt.clear();
      st.push_back(thisList);
      i++;
    }
    return strandsize;
  } 


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << LongDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    double t0; 

    cout << endl;
    cout << "Reading from stdin... (type C-c or 0 to quit)" << endl; // 2019-04

    long N,r;
    string CommentLine;
    while (cin.peek()=='-')  // skip initial comment lines, if any
      getline(cin, CommentLine);

    cin >> N;
    if (N==0) return;

    // We read the monomials and construct an ideal in a polynomial ring
    SparsePolyRing P = NewPolyRing(RingQQ(), SymbolRange("x",0,N-1)); 
    RingElem m = product(indets(P));

    std::vector<RingElem> f;
    t0 = CpuTime();
    f = RingElems(P, cin);
    r = len(f);
    cout << "Finished reading: "
         << r << " monomials in " << N << " variables"
         << " (sec " << CpuTime()-t0 << ")" << endl << endl;
    
    for(int i=0; i<len(f); ++i)      f[i] = m/f[i];
    ideal I=ideal(f);
    //cout<<I<<endl;

    // We transform the ideal in the polynomial ring to a vector of DynamicBitsets in a PPMOnoid
    PPMonoid PPMI = PPM(RingOf(I));
    const std::vector<RingElem>& gensI = gens(I);
    PPVectorSqFr J;
    J.reserve(len(gensI));
    for (long i=0; i<len(gensI); ++i)
    {
      if (!IsMonomial(gensI[i])) CoCoA_THROW_ERROR("is not monomial","MVTSquareFree");
      if (!IsOne(LC(gensI[i]))) CoCoA_THROW_ERROR("coeff is not 1","MVTSquareFree");
      if (!IsSqFree(LPP(gensI[i]))) CoCoA_THROW_ERROR("not squarefree","MVTSquareFree");
      J.push_back(DynamicBitset(LPP(gensI[i])));
    }
    InterreduceSort(J);
    //cout << "PPVectorSqFr: converted & interreduced" << endl; // 2019-04
    //cout<<J<<endl;


    //--- MVT-Simplicial -----------
    // MVTSquareFree
    t0 = CpuTime();
    MVTSquareFree MVTSqFr(I);
    std::vector<int> bettis2(NumIndets(P),0);
    cout << "Construction of MV Tree SqFr (sec " << CpuTime()-t0 << ")" << endl;

    t0 = CpuTime();
    MVTSqFr.myComputeBetti(bettis2);
    cout << "Reading the Betti numbers SqFr (sec " << CpuTime()-t0 <<")" <<endl;
    OutputRange(cout, bettis2.begin(), bettis2.end());  cout << endl;

    //-------------------------------
    long MaxDeg=0;
    for (long i=0; i<len(gens(I)); ++i)
      MaxDeg = max(MaxDeg, deg(gens(I)[i])); // max dim of facets
    cout << "MaxDeg = " << MaxDeg << ";  " << flush;
    //-------------------------------
    long lo=0, hi=0;
    for (long i=0; i<len(bettis2); ++i)
      if (bettis2[i] != 0)
      {
        lo=i-1;
        cout << "lo = " << lo << ";  " << flush;
        break;
      }
    for (long i=len(bettis2)-1; i>=0; --i)
      if (bettis2[i] != 0)
      {
        hi=i+1;
        cout << "hi = " << hi << endl;
        break;
      }
    //-------------------------------
    if (hi == lo+2)
      cout<<"Betti["<<lo+2<<"] = "<<bettis2[lo+1]<<endl;
    else
    { 
      //-------------------------------
      long HI = min(hi, MaxDeg);
      //-------------------------------

      //We compute the skeleton of the Taylor 1-Strand
      vector<vector<list<int> > > st;

      t0=CpuTime();
      int stSize=oneStrand(st,J,PPMI,lo,HI);
      cout << "size of the skeleton of the Taylor 1-strand: " << stSize;
      cout << " (sec " << CpuTime()-t0 << ")" << endl;
    
      t0=CpuTime();
      cout << "actual Taylor 1-strand with its differentials " << flush;
      ring QQ=RingQQ();
      vector<matrix> diff;
      vector<int> diffRanks;
      diffRanks.push_back(0);
      unsigned int rows,columns;
      for(unsigned int i=1;i<st.size();i++)
      {
        rows=st[i-1].size();
        columns=st[i].size();
        if(rows>0 && columns>0)
        {
          //we compute the index-subset of each generator to use it as an identifier
          vector<int> ti,ti1;
          //for (int j=0;j<st[i].size();j++)//WARNING: may be this is useless and we should also compute ti1 (and thusavoid the double computation)
          //	     ti.push_back(orderSubset(st[i][j],r,i));
          for (int j=0;j<st[i-1].size();j++) //WARNING: by now we have already computed this,  we should avoid this double computation
            ti1.push_back(orderSubset(st[i-1][j],r,i+lo-1));
          //cout<<endl<<"$$$"<<ti<<endl<<"$$$"<<ti1<<endl;
          //we create zero matrices to store the differentials
          matrix D(NewDenseMat(QQ,rows,columns));
          diff.push_back(D);
          //cout<<D;
          //Now we fill the matrices
          for (unsigned int e=0;e<st[i].size();e++)
          {
            list<int> p;
            p=st[i][e];
            list<int>::iterator it=p.begin();
            //print(p);
            int j=0;
            int sign=1;
            for(it=p.begin();it!=p.end();it++)
            {
              int entryColumn=e;
              int entryRow=j;
              int image=orderSubsetExceptOne(p,r,i+lo,j);
              int imageIndex=0;
              bool found=0;
              while(imageIndex<ti1.size() && found==0)
              {
                if(ti1[imageIndex]==image)
                {
                  found=1;
                }
                imageIndex++;
              }
              if(found)
              {
                entryRow=imageIndex-1;
                SetEntry(D,entryRow,entryColumn,sign);
              }
              sign=sign*(-1);
              //cout<<(*it)<<"-->"<<j<<"-->"<<orderSubsetExceptOne(p,r,i,j)<<endl;
              j++;
            }
          }
          //  cout<<D;
          int thisRank=rk(D);
          diffRanks.push_back(thisRank);
          //cout<<i<<" RANK= "<<rk(D)<<endl;
        }
        else
          if (rows==0) diffRanks.push_back(0);
      }
      diffRanks.push_back(binomial(r-1,st.size()-1));
      cout << " (sec " << CpuTime()-t0 << ")" << endl;
    
      //Finally, we compute the Betti numbers:
      vector<int> bettis;
      for(int i=1; i<st.size()-1; ++i)
        if (st.size()==0)
          bettis.push_back(0);
        else
          bettis.push_back(st[i].size()-diffRanks[i]-diffRanks[i+1]);

      for (int i=0; i<st.size(); ++i)
      {
        cout << "i: " << i+lo
             << " \tmodule rank: " << st[i].size()
             << " \tdiffRank: " << diffRanks[i] << endl;
      }
    
      for (int i=1; i<bettis.size(); ++i)
        cout << "Betti[" << i+lo+1 << "] = " << bettis[i] << endl;
    }
    

    //--------------
    //PRUEBAS
    //---------------
    //   cout<<endl<<"***PRUEBAS***"<<endl;

    // for(unsigned int i=0;i<st.size();i++)
    //{cout<<i<<": "<<st[i].size()<<endl;
    //for(std::vector<list<int> >::iterator it=st[i].begin();it!=st[i].end();it++)
    //  {cout<<i<<": "; print(*it);cout<<orderSubset(*it,r,i)<<endl;}
    //}

    
    //list<int> p;
    //p=st[4][3];
    //list<int>::iterator it=p.begin();
    //print(p);
    //advance(it,2);
    //p.erase(it);
    //print(p);
    //int j=0;
    //for(it=p.begin();it!=p.end();it++)
    //{cout<<(*it)<<"-->"<<j<<"-->"<<orderSubsetExceptOne(p,r,4,j)<<endl;
    // j++;
    // }
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
  // catch (const CoCoA::InterruptReceived& intr)
  // {
  //   cerr << endl
  //        << "------------------------------" << endl
  //        << ">>>  CoCoALib interrupted  <<<" << endl
  //        << "------------------------------" << endl
  //        << "-->>  " << intr << "  <<--" << endl;
  //   return 2;
  // }
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

