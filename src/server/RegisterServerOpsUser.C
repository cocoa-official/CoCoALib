//   Copyright (c)  2010  Anna Bigatti,  John Abbott

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

#include "CoCoA/library.H"
#include "CoCoA5io.H"
#include "ServerOp.H"

// #include <iostream> --- already included in library.H
using std::ostream;
using std::istream;
using std::endl;
// #include <list>   --- already included in library.H
using std::list;
// #include <memory> --- already included in library.H
using std::unique_ptr;
// #include <string> --- already included in library.H
using std::string;
// #include <vector> --- already included in library.H
using std::vector;

namespace CoCoA
{

  // write the name your sublibrary like this
  const ServerOpBase::LibraryInfo& CoCoALib_user()
  {
    static const ServerOpBase::LibraryInfo UniqueValue("CoCoALib", BuildInfo::version(), "user");
    return UniqueValue;
  }


  // your function should be a class like this
  class EmptyOperation: public ServerOpBase
  {
  public:
    EmptyOperation(): ServerOpBase(CoCoALib_user()) {};
    ~EmptyOperation() {};
    void myOutputSelf(std::ostream& out) const { out << "EmptyOperation"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutEmpty = 1234; /* your function call */; }
    void myWriteResult(std::ostream& out) const { out << ourVarName4 << " := " << myOutEmpty << ";" << endl; }
    void myClear() { myInEmpty = 0;  myOutEmpty = 0; }
  private:
    int myInEmpty;
    int myOutEmpty;
  };

  // write your function as a class like this
  void EmptyOperation::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and myIn
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    /* ..... */
  }



class MonomialIsIrreducible: public ServerOpBase
  {
  public:
    MonomialIsIrreducible(): ServerOpBase(CoCoALib_user()) {};
    ~MonomialIsIrreducible() {};
    void myOutputSelf(std::ostream& out) const {out << "Monomial_IsIrreducible";}
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutReg = IsIrreducible(*myInPPV); }
    void myWriteResult(std::ostream& out) const{out << ourVarName4 <<" := "<< myOutReg << ";"<<endl;}
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    bool myOutReg;
  };

  void MonomialIsIrreducible::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and PPList
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x", 0, NumIndets(P)-1), lex);
    vector<PPMonoidElem> V;
    ReadPPs(in, V, PPM, GetTag);
    DivMaskRule DMR = NewDivMaskEvenPowers();
    vector<long> exps(NumIndets(PPM));
    myInPPV.reset(new PPVector(PPM, DMR));
    //myOutPPsPtr.reset(new int);

    for(size_t i=0; i<V.size(); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }



class MonomialIsPrimary: public ServerOpBase
  {
  public:
    MonomialIsPrimary(): ServerOpBase(CoCoALib_user()) {};
    ~MonomialIsPrimary() {};
    void myOutputSelf(std::ostream& out) const {out << "Monomial_Is_Primary";}
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutReg = IsPrimary(*myInPPV); }
    void myWriteResult(std::ostream& out) const{out << ourVarName4 <<" := "<< myOutReg << ";"<<endl;}
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    bool myOutReg;
  };

  void MonomialIsPrimary::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and PPList
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x", 0, NumIndets(P)-1), lex);
    vector<PPMonoidElem> V;
    ReadPPs(in, V, PPM, GetTag);
    DivMaskRule DMR = NewDivMaskEvenPowers();
    vector<long> exps(NumIndets(PPM));
    myInPPV.reset(new PPVector(PPM, DMR));
    //myOutPPsPtr.reset(new int);

    for(size_t i=0; i<V.size(); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }


class MonomialIsPrime: public ServerOpBase
  {
  public:
    MonomialIsPrime(): ServerOpBase(CoCoALib_user()) {};
    ~MonomialIsPrime() {};
    void myOutputSelf(std::ostream& out) const {out << "Monomial_IsPrime";}
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutReg = IsPrime(*myInPPV); }
    void myWriteResult(std::ostream& out) const{out << ourVarName4 <<" := "<< myOutReg << ";"<<endl;}
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    bool myOutReg;
  };

  void MonomialIsPrime::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and PPList
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x", 0, NumIndets(P)-1), lex);
    vector<PPMonoidElem> V;
    ReadPPs(in, V, PPM, GetTag);
    DivMaskRule DMR = NewDivMaskEvenPowers();
    vector<long> exps(NumIndets(PPM));
    myInPPV.reset(new PPVector(PPM, DMR));
    //myOutPPsPtr.reset(new int);

    for(size_t i=0; i<V.size(); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }




//   namespace
//   {
//     class IdealOperationRegistrationClass
//     {
//     public:
//       IdealOperationRegistrationClass()
//       {
//       }
//     } FakeVariable1;
//   }



  namespace CoCoAServerOperationsFromUser
  {
    bool RegisterOps()
    {
      RegisterOp("Empty_Operation_Name",   ServerOp(new EmptyOperation()));
      RegisterOp("Monomial_Is_Irreducible",   ServerOp(new MonomialIsIrreducible()));
      RegisterOp("Monomial_Is_Primary",   ServerOp(new MonomialIsPrimary()));
      RegisterOp("Monomial_Is_Prime",   ServerOp(new MonomialIsPrime()));
      return true;
    }


    bool RegisterOpsOnce()
    {
      static bool EvalOnce = RegisterOps();
      return EvalOnce;
    }
  }


}
