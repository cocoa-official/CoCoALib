// Copyright (c) 2006  John Abbott, Anna Bigatti
// Orig author: 2006  Eduardo Saenz-de-Cabezon
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <fstream>
// using std::ifstream; using std::ofstream;
#include <cstdlib>
// using exit

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example of use of the Mayer-Vietoris trees.  \n";

const string LongDescription =
  "Example of use of the Mayer-Vietoris trees.  \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {

    // ANNA: moved here from MVT_LEX_u.C, global
    GlobalManager CoCoAFoundations;

    ifstream fich;
    int N;
    int r;
    fich.open("ex-MVT.in");
    if (fich.fail())
    {
      cerr << "unable to open file ex-MVT.in for reading" << endl;
      exit(1);
    }
    fich>>N;
    fich>>r;
    cout<<"Number of variables of the ring: "<<N<<endl;
    cout<<"Number of generators of the ideal in the file: "<<r<<endl;
  
    cout<<endl;


    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x",0,N-1), lex);
    DivMaskRule DMR = NewDivMaskEvenPowers();
    const vector<PPMonoidElem>& x = indets(PPM);

    PPVector f(PPM, DMR);
  
    for (int i=0; i<r; ++i)
    {
      PPMonoidElem otro=one(PPM);
      int pow;				
      for (long j=0; j<N; ++j)
      {
        fich>>pow;
        otro=otro*power(x[j],pow);
      };
      PushBack(f, otro);
    };
  
  
    double start, finish;

    start = CpuTime();
    PPVector n1betti(PPM,DMR);
    MayerVietorisTreeN1(n1betti, f);
    finish = CpuTime();
    cout << "Construction of N-1 MV Tree: Time in total (seconds): "  << finish-start;
    cout<<endl<<len(n1betti)<<endl;
    cout<<endl;


    start = CpuTime();
    MultidegreeMap output_list, undecided_list;
    MayerVietorisTree(output_list, f);
    finish = CpuTime();
    cout << "Construction of MV Tree: Time in total (seconds): "  << finish - start;
    cout<<endl;
  
  
  
  
    vector<int> bettis(N,0);
    bettis[0]=r;
    start = CpuTime();
    Bettis(bettis, output_list);
    finish = CpuTime();
    cout << "Reading the Betti numbers: Time in total (seconds): "  << finish - start;
    cout<<endl;
    cout<<bettis<<"Size of resolution: "<<ResSize(bettis)<<endl;

    //put(output_list);
    //cout<<endl;


    vector<int> bettis2(N,0);
    bettis2[0]=r;
    start = CpuTime();
    ReduceMVTree(output_list, undecided_list);
    finish = CpuTime();
    cout << "Reducing the tree: Time in total (seconds): "  << finish - start;

    Bettis(bettis2, output_list);

    cout<<endl;

    //put(output_list);
    //cout<<endl;
    cout<<bettis2<<"Size of resolution: "<<ResSize(bettis2)<<endl;


    vector<int> bettis3(N,0);
    bettis3[0]=0;
    start = CpuTime();

    Bettis(bettis3, undecided_list);

    cout<<endl;
    cout<<bettis3<<"Size of undecided: "<<ResSize(bettis3)<<endl;
//put(output_list);
//cout<<endl;
    int reg,pd;

    cout<<endl;
    cout<<"Betti diagram of the DECIDED ones "<<endl;
    BettiPseudoDiagram bettidiagram;
    GradedBettis(output_list,bettidiagram,reg,pd);

//PrintBettiDiagram(bettidiagram);

    cout<<endl<<"C-M regularity: "<<reg<<endl;
    cout<<endl<<"Projective dimension: "<<pd<<endl;

    bettidiagram.clear();
    cout<<endl;
    cout<<"Betti diagram of the UNDECIDED ones "<<endl;
    GradedBettis(undecided_list,bettidiagram,reg,pd);

//PrintBettiDiagram(bettidiagram);
    cout<<endl<<"C-M regularity: "<<reg<<endl;
    cout<<endl<<"Projective dimension: "<<pd<<endl;

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
