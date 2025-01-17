// Copyright (c) 2015  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This simple example shows how to use the SmallFp export conventions\n"
  "\"SymmResidues\" and \"NonNegResidues\", and the effect they have. \n";

const string LongDescription =
  "This simple example shows how to use the SmallFp export conventions\n"
  "\"SymmResidues\" and \"NonNegResidues\", and the effect they have. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void PrintFpElems(const SmallFpImpl& ModP)
  {
    cout << "Exported values from " << ModP << endl;
    SmallFpImpl::value a = zero(SmallFp);
    do
    {
      cout << ModP.myExport(a) << "  ";
      a = ModP.myAdd(a,one(SmallFp));
    } while (!IsZero(a));
    cout << endl << endl;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    SmallFpImpl FF5(5); // uses default setting stored in GlobalManager
    SmallFpImpl FF5symm(5, GlobalSettings::ResidueRepr::symmetric);
    SmallFpImpl FF5nonneg(5, GlobalSettings::ResidueRepr::NonNegative);
    PrintFpElems(FF5);
    PrintFpElems(FF5symm);
    PrintFpElems(FF5nonneg);
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
