// Copyright (c) 2015 John Abbott, Anna Bigatti, Anders Nedergaard Jensen
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows what can be computed in CoCoALib using GFan:  \n"
  "a library for computations in ..., \n"
  "....";

const string LongDescription = 
  "This program shows what can be computed in CoCoALib using GFan:  \n"
  "a library for computations in ..., \n"
  "....";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

#ifndef CoCoA_WITH_GFAN
    cout << "GFan library is not available to CoCoALib." << endl;
#else // GFan is available
    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    matrix CM(NewDenseMat(RingQQ(),2,3));
    SetEntry(CM,0,0, 300);  //     M[0][0] = 300; 
    SetEntry(CM,0,1, 2);    //     M[0][1] = 2;   
    SetEntry(CM,0,2, -2);   //     M[0][2] = -2;  
    SetEntry(CM,1,0, 30);   //     M[1][0] = 30;  
    SetEntry(CM,1,1, 5);    //     M[1][1] = 5;   
    SetEntry(CM,1,2, 1);    //     M[1][2] = 1;

    matrix CM2(NewDenseMat(RingQQ(),1,3));
    SetEntry(CM,0,0, 4);
    SetEntry(CM,0,1, 1);
    SetEntry(CM,0,2, 5);

    GFan::cone CC(CM, CM2);
    cout << "The GFan cone is " << CC << endl;

    cout << "equations(CC) = " << equations(CC) << endl;
    cout << "inequalities(CC) = " << inequalities(CC) << endl;
    cout << NumRows(inequalities(CC)) << "x"
         << NumCols(inequalities(CC)) << endl;
    
    cout << RelativeInteriorPoint(CC) << endl;


#endif // CoCoA_WITH_GFAN
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
