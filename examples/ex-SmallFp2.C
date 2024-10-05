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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-SmallFp2.C,v 1.3 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-SmallFp2.C,v $
// Revision 1.3  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.2  2021/03/04 21:03:45  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.1  2015/11/04 10:12:02  abbott
// Summary: New example
//
//
