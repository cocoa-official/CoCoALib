//   Copyright (c)  2007-2009  John Abbott and Anna Bigatti

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
#include "CoCoA5io.H" // #include "CoCoA4io.H"
#include "GlobalIO.H"
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

  // sublibrary of CoCoALib for groebner related operations
  // by M.Caboara
  const ServerOpBase::LibraryInfo& CoCoALib_groebner()
  {
    static const ServerOpBase::LibraryInfo UniqueValue("CoCoALib", BuildInfo::version(), "groebner");
    return UniqueValue;
  }

  // sublibrary of CoCoALib for approx points
  // by J Abbott, M-L Torrente
  const ServerOpBase::LibraryInfo& CoCoALib_approx()
  {
    static const ServerOpBase::LibraryInfo UniqueValue("CoCoALib", BuildInfo::version(), "approx");
    return UniqueValue;
  }

  // sublibrary of CoCoALib for monomial (squarefree) ideals
  // by M.Caboara, E.Saenz-de-Cabezon
  const ServerOpBase::LibraryInfo& CoCoALib_combinatorics()
  {
    static const ServerOpBase::LibraryInfo UniqueValue("CoCoALib", BuildInfo::version(), "combinatorics");
    return UniqueValue;
  }


  // ---- Verbosity Level as optional last argument

//   int TryReadingVerbosityLevel(istream& in)
//   {
//     SkipTag(GlobalInput(), "<verbosity_level>");
//     return ReadVerbosityLevel(in, TagWasRead);
//   }


//   // ---- TestSocket ----
//   class TestSocket: public ServerOpBase
//   {
//   public:
//     TestSocket(): ServerOpBase(CoCoALib_groebner()) {};
//     ~TestSocket() {};
//     void myOutputSelf(std::ostream& out) const { out << "TestSocket"; }
//     void myReadArgs(std::istream& in, int NumArgs);
//     void myCompute()  { /**/ }
//     void myWriteResult(std::ostream& out) const {out << ourVarName4 << " := True"; }
//     void myClear() {};
//   private:
//     vector<RingElem> myInPL;  // ANNA: this is totally useless, I'll clean it up..
//   };


//   void TestSocket::myReadArgs(std::istream& in, int NumArgs)
//   {
//     CoCoA_ASSERT(NumArgs==2);
//     const SparsePolyRing P(ReadPolyRing(in, GetTag));
//     ReadPolyList(in, myInPL, P, GetTag);
//   }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ModuleGBasis: public ServerOpBase
  {
  public:
    ModuleGBasis(): ServerOpBase(CoCoALib_groebner()) {};
    ~ModuleGBasis() {};
    void myOutputSelf(std::ostream& out) const { out << "ModuleGBasis"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute()  {VectorList NoUse; ComputeGBasis(myResVL, NoUse, myInVL);}
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInVL.clear(); myResVL.clear(); }
  private:
    VectorList myInVL, myResVL;
  };


  void ModuleGBasis::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    ReadVectorList(in, myInVL, FM, GetTag);
  }


  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ModuleLT: public ServerOpBase
  {
  public:
    ModuleLT(): ServerOpBase(CoCoALib_groebner()) {};
    ~ModuleLT() {};
    void myOutputSelf(std::ostream& out) const { out << "ModuleLT"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeLT(myResVL, myInVL); }
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInVL.clear(); myResVL.clear(); }
  private:
    VectorList myInVL, myResVL;
  };


  void ModuleLT::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    ReadVectorList(in, myInVL, FM, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ModuleSyzygy: public ServerOpBase
  {
  public:
    ModuleSyzygy(): ServerOpBase(CoCoALib_groebner()) {};
    ~ModuleSyzygy() {};
    void myOutputSelf(std::ostream& out) const { out << "ModuleSyzygy"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeSyz(myResVL, NewFreeModuleForSyz(myInVL), myInVL); }
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInVL.clear(); myResVL.clear(); }
  private:
    VectorList myInVL, myResVL;
  };


  void ModuleSyzygy::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    ReadVectorList(in, myInVL, FM, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ModuleIntersection: public ServerOpBase
  {
  public:
    ModuleIntersection(): ServerOpBase(CoCoALib_groebner()) {};
    ~ModuleIntersection() {};
    void myOutputSelf(std::ostream& out) const { out << "ModuleIntersection"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeIntersection(myResVL, myInVL1, myInVL2); }
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInVL1.clear(); myInVL2.clear(); myResVL.clear(); }
  private:
    VectorList myInVL1, myInVL2, myResVL;
  };


  void ModuleIntersection::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    ReadVectorList(in, myInVL1, FM, GetTag);
    ReadVectorList(in, myInVL2, FM, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ColonModMod: public ServerOpBase
  {
  public:
    ColonModMod(): ServerOpBase(CoCoALib_groebner()) {};
    ~ColonModMod() {};
    void myOutputSelf(std::ostream& out) const { out << "ColonModMod"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeColon(myResPL, myInVL1, myInVL2); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInVL1.clear(); myInVL2.clear(); myResPL.clear(); }
  private:
    VectorList myInVL1, myInVL2;
    PolyList myResPL;
  };


  void ColonModMod::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    ReadVectorList(in, myInVL1, FM, GetTag);
    ReadVectorList(in, myInVL2, FM, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ModuleSaturation: public ServerOpBase
  {
  public:
    ModuleSaturation(): ServerOpBase(CoCoALib_groebner()) {};
    ~ModuleSaturation() {};
    void myOutputSelf(std::ostream& out) const { out << "ModuleSaturation"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeSaturation(myResVL, myInVL, myInPL); }
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInVL.clear(); myInPL.clear(); myResVL.clear(); }
  private:
    VectorList myInVL, myResVL;
    PolyList myInPL;
  };


  void ModuleSaturation::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    const SparsePolyRing P = RingOf(FM);
    ReadVectorList(in, myInVL, FM, GetTag);
    ReadPolyList(in, myInPL, P, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealSaturation: public ServerOpBase
  {
  public:
    IdealSaturation(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealSaturation() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealSaturation"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeSaturation(myResPL, myInPL1, myInPL2); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL1.clear(); myInPL2.clear(); myResPL.clear(); }
  private:
    PolyList myInPL1, myInPL2, myResPL;
  };


  void IdealSaturation::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL1, P, GetTag);
    ReadPolyList(in, myInPL2, P, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  // class ColonModId: public ServerOpBase
  // {
  // public:
  //   ColonModId(): ServerOpBase(CoCoALib_groebner()) {};
  //   ~ColonModId() {};
  //   void myOutputSelf(std::ostream& out) const { out << "ColonModId"; }
  //   void myReadArgs(std::istream& in, int NumArgs);
  //   void myCompute() { ComputeCColon(myResVL, myInVL, myInPL); }
  //   void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
  //  void myClear() { myInVL.clear(); myInPL.clear(); myResVL.clear(); }
  // private:
  //   VectorList myInVL, myResVL;
  //   PolyList myInPL;
  // };


  // void ColonModId::myReadArgs(std::istream& in, int NumArgs)
  // {
  //   const FreeModule FM(ReadFreeModule(in, GetTag));
  //   const SparsePolyRing P = RingOf(FM);
  //   ReadVectorList(myInVL, FM, GetTag);
  //   ReadPolyList(myInPL, P, GetTag);
  // }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ColonIdId: public ServerOpBase
  {
  public:
    ColonIdId(): ServerOpBase(CoCoALib_groebner()) {};
    ~ColonIdId() {};
    void myOutputSelf(std::ostream& out) const { out << "ColonIdId"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeColon(myResPL, myInPL1, myInPL2); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL1.clear(); myInPL2.clear(); myResPL.clear(); }
  private:
    PolyList myInPL1, myInPL2, myResPL;
  };


  void ColonIdId::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL1, P, GetTag);
    ReadPolyList(in, myInPL2, P, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ModuleElim: public ServerOpBase
  {
  public:
    ModuleElim(): ServerOpBase(CoCoALib_groebner()) {};
    ~ModuleElim() {};
    void myOutputSelf(std::ostream& out) const { out << "ModuleElim"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeElim(myResVL, myInVL, *myInElimPPPtr); }
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInVL.clear(); myInElimPPPtr.reset(0); myResVL.clear(); }
  private:
    VectorList myInVL, myResVL;
    unique_ptr<PPMonoidElem> myInElimPPPtr;
  };


  void ModuleElim::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    const SparsePolyRing P = RingOf(FM);
    ReadVectorList(in, myInVL, FM, GetTag);
    PolyList PL;
    myInElimPPPtr.reset(new PPMonoidElem(PPM(P)));
    ReadPolyList(in, PL, P, GetTag);
    for (PolyList::const_iterator it=PL.begin(); it != PL.end(); ++it)
      (*myInElimPPPtr) *= LPP(*it);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class ModuleHomog: public ServerOpBase
  {
  public:
    ModuleHomog(): ServerOpBase(CoCoALib_groebner()) {};
    ~ModuleHomog() {};
    void myOutputSelf(std::ostream& out) const { out << "ModuleHomog"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeHomogenization(myResVL, myInVL, myInHomIndets); }
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInVL.clear(); myInHomIndets.clear(); myResVL.clear(); }
  private:
    VectorList myInVL, myResVL;
    PolyList myInHomIndets;
  };


  void ModuleHomog::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const FreeModule FM(ReadFreeModule(in, GetTag));
    const SparsePolyRing P = RingOf(FM);
    ReadVectorList(in, myInVL, FM, GetTag);
    ReadPolyList(in, myInHomIndets, P, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealGBasis: public ServerOpBase
  {
  public:
    IdealGBasis(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealGBasis() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealGBasis"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { PolyList tmp; ComputeGBasis(myResPL, tmp, myInPL); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL.clear(); myResPL.clear(); }
  private:
    PolyList myInPL, myResPL;
  };


  void IdealGBasis::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs==2);
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL, P, GetTag);
  }

 // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealSATGBasis: public ServerOpBase
  {
  public:
    IdealSATGBasis(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealSATGBasis() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealSATGBasis"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeGBasisSelfSatCore(myResPL, myInPL); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL.clear(); myResPL.clear(); }
  private:
    PolyList myInPL, myResPL;
  };

  void IdealSATGBasis::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs==2);
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL, P, GetTag);
  }

// ---- CoCoA/GOperations.H  by  M.Caboara ----
//   class IdealSATMixGBasis: public ServerOpBase
//   {
//   public:
//     IdealSATMixGBasis(): ServerOpBase(CoCoALib_groebner()) {};
//     ~IdealSATMixGBasis() {};
//     void myOutputSelf(std::ostream& out) const { out << "IdealSATMixGBasis"; }
//     void myReadArgs(std::istream& in, int NumArgs);
//     void myCompute() { ComputeSATMixGBasis(myResPL, myInPL); }
//     void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
//     void myClear() { myInPL.clear(); myResPL.clear(); }
//   private:
//     PolyList myInPL, myResPL;
//   };

//   void IdealSATMixGBasis::myReadArgs(std::istream& in, int NumArgs)
//   {
//     //    CoCoA_ASSERT(NumArgs==2 || NumArgs==3);
//     CoCoA_ASSERT(NumArgs==2);
//     const SparsePolyRing P(ReadPolyRing(in, GetTag));
//     ReadPolyList(in, myInPL, P, GetTag);
//   }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealLT: public ServerOpBase
  {
  public:
    IdealLT(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealLT() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealLT"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeLT(myResPL, myInPL); }
    void myWriteResult(std::ostream& out) const;
    void myClear() { myInPL.clear(); myResPL.clear(); }
  private:
    PolyList myInPL, myResPL;
    //    ideal myOutIdeal;
  };


  void IdealLT::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs==2);
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL, P, GetTag);
  }


  void IdealLT::myWriteResult(std::ostream& out) const
  {
    // this function will be nicer when ComputeLT "returns" an ideal
    if (myResPL.empty())
      out << ourVarName4 << " := ideal()";
    else
    {
      out << ourVarName4 << " := ";
      WriteIdeal(out, ideal(myResPL));
    }
  }
  

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealSyzygy: public ServerOpBase
  {
  public:
    IdealSyzygy(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealSyzygy() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealSyzygy"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeSyz(myResVL, NewFreeModuleForSyz(myInPL), myInPL); }
    void myWriteResult(std::ostream& out) const { WriteVectorListInVar(out, ourVarName4, myResVL); }
    void myClear() { myInPL.clear(); myResVL.clear(); }
  private:
    PolyList myInPL;
    VectorList myResVL;
  };


  void IdealSyzygy::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs==2);
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL, P, GetTag);
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealIntersection: public ServerOpBase
  {
  public:
    IdealIntersection(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealIntersection() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealIntersection"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeIntersection(myResPL, myInPL1, myInPL2); }
    //    void myCompute() { myOutIdeal := intersect(myInI1, myInI2); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL1.clear(); myInPL2.clear(); myResPL.clear(); }
  private:
    PolyList myInPL1, myInPL2, myResPL;
    //    ideal myInI1, myInI2, myOutI;
  };


  void IdealIntersection::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL1, P, GetTag);
    ReadPolyList(in, myInPL2, P, GetTag);
  }


//   void IdealIntersection::myWriteResult(std::ostream& out) const
//   {
//     out << ourVarName4 << " := ";
//     WriteIdeal(out, myOutIdeal);
//   }
  


  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealElim: public ServerOpBase
  {
  public:
    IdealElim(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealElim() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealElim"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeElim(myResPL, myInPL, *myInElimPPPtr); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL.clear(); myResPL.clear(); myInElimPPPtr.reset(0); }
  private:
    PolyList myInPL, myResPL;
    unique_ptr<PPMonoidElem> myInElimPPPtr;
  };


  void IdealElim::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL, P, GetTag);
    PolyList PL;
    PPMonoidElem t(PPM(P));
    ReadPolyList(in, PL, P, GetTag);
    for (PolyList::const_iterator it=PL.begin(); it != PL.end(); ++it)
    {
      if ( !IsIndet(*it) )
        CoCoA_THROW_ERROR("Expected indet", "IdealElim::myReadArgs");
      t *= LPP(*it);
    }
    myInElimPPPtr.reset(new PPMonoidElem(t));
  }

  // ---- CoCoA/GOperations.H  by  M.Caboara ----
  class IdealHomog: public ServerOpBase
  {
  public:
    IdealHomog(): ServerOpBase(CoCoALib_groebner()) {};
    ~IdealHomog() {};
    void myOutputSelf(std::ostream& out) const { out << "IdealHomog"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { ComputeHomogenization(myResPL, myInPL, myInHomIndets); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL.clear(); myInHomIndets.clear(); myResPL.clear(); }
  private:
    PolyList myInPL, myInHomIndets, myResPL;
  };


  void IdealHomog::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL, P, GetTag);
    ReadPolyList(in, myInHomIndets, P, GetTag);
  }

  // ---- CoCoA/TmpF5.H  by  A.Arri ----
  class F5GBasis: public ServerOpBase
  {
  public:
    F5GBasis(): ServerOpBase(CoCoALib_groebner()) {};
    ~F5GBasis() {};
    void myOutputSelf(std::ostream& out) const { out << "F5GBasis"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myResPL = F5_poly(myInPL); /*F5(myResPL, myInPL);*/ }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, myResPL); }
    void myClear() { myInPL.clear(); myResPL.clear(); }
  private:
    PolyList myInPL, myResPL;
  };


  void F5GBasis::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    ReadPolyList(in, myInPL, P, GetTag);
  }

  // ---- CoCoA/TmpIsTree.H  by  M.Caboara ----
  class IsTreeNoOpt: public ServerOpBase
  {
  public:
    IsTreeNoOpt(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~IsTreeNoOpt() {};
    void myOutputSelf(std::ostream& out) const { out << "IsTreeNoOpt"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutCycle = myInComplex.myIsTreeNoOpt(); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, FacetList2PolyList(*myRingPtr,myOutCycle)); }
    void myClear() { myInComplex.myClear(); myOutCycle.clear(); }
  private:
    FacetComplex myInComplex;
    list<facet> myOutCycle;
    unique_ptr<SparsePolyRing> myRingPtr;
  };


  void IsTreeNoOpt::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    PolyList PL;
    myRingPtr.reset(new SparsePolyRing(ReadPolyRing(in, GetTag)));
    ReadPolyList(in, PL, *myRingPtr, GetTag);
    myInComplex = FacetComplex(*myRingPtr, PL);
  }

  // ---- CoCoA/TmpIsTree.H  by  M.Caboara ----
  class IsTreeOpt: public ServerOpBase
  {
  public:
    IsTreeOpt(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~IsTreeOpt() {};
    void myOutputSelf(std::ostream& out) const { out << "IsTreeOpt"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutCycle = myInComplex.myIsTreeOpt(); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, FacetList2PolyList(*myRingPtr,myOutCycle)); }
    void myClear() { myInComplex.myClear(); myOutCycle.clear(); }
  private:
    FacetComplex myInComplex;
    list<facet> myOutCycle;
    unique_ptr<SparsePolyRing> myRingPtr;
  };


  void IsTreeOpt::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    PolyList PL;
    myRingPtr.reset(new SparsePolyRing(ReadPolyRing(in, GetTag)));
    ReadPolyList(in, PL, *myRingPtr, GetTag);
    myInComplex = FacetComplex(*myRingPtr, PL);
  }

  // ---- CoCoA/TmpIsTree.H  by  M.Caboara ----
  class IsTreeCBNoOpt: public ServerOpBase
  {
  public:
    IsTreeCBNoOpt(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~IsTreeCBNoOpt() {};
    void myOutputSelf(std::ostream& out) const { out << "IsTreeCBNoOpt"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutCycle = myInComplex.myIsTreeCBNoOpt(); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, FacetList2PolyList(*myRingPtr,myOutCycle)); }
    void myClear() { myInComplex.myClear(); myOutCycle.clear(); }
  private:
    FacetComplex myInComplex;
    list<facet> myOutCycle;
    unique_ptr<SparsePolyRing> myRingPtr;
  };


  void IsTreeCBNoOpt::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    PolyList PL;
    myRingPtr.reset(new SparsePolyRing(ReadPolyRing(in, GetTag)));
    ReadPolyList(in, PL, *myRingPtr, GetTag);
    myInComplex = FacetComplex(*myRingPtr, PL);
  }

  // ---- CoCoA/TmpIsTree.H  by  M.Caboara ----
  class IsTreeCBOpt: public ServerOpBase
  {
  public:
    IsTreeCBOpt(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~IsTreeCBOpt() {};
    void myOutputSelf(std::ostream& out) const { out << "IsTreeCBOpt"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutCycle = myInComplex.myIsTreeCBOpt(); }
    void myWriteResult(std::ostream& out) const { WritePolyListInVar(out, ourVarName4, FacetList2PolyList(*myRingPtr,myOutCycle)); }
    void myClear() { myInComplex.myClear(); myOutCycle.clear(); }
  private:
    FacetComplex myInComplex;
    list<facet> myOutCycle;
    unique_ptr<SparsePolyRing> myRingPtr;
  };


  void IsTreeCBOpt::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs<10);  // ??? fix with allowed values
    PolyList PL;
    myRingPtr.reset(new SparsePolyRing(ReadPolyRing(in, GetTag)));
    ReadPolyList(in, PL, *myRingPtr, GetTag);
    myInComplex = FacetComplex(*myRingPtr, PL);
  }


  //----------------------------------------------------
  // ---- CoCoA/ApproxPts.H  by  L.Torrente and J.Abbott ----
  class StableBorder: public ServerOpBase
  {
  public:
    StableBorder(): ServerOpBase(CoCoALib_approx()) {};
    ~StableBorder() {};
    void myOutputSelf(std::ostream& out) const { out << "StableBorder"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute()  { ApproxPts::SOITwinFloat(myOutSOI, myOutBBasis, myOutAlmostVanishing, *myPolyRingPtr, myInPts, myInTolerance, *myInGammaPtr); }
    void myWriteResult(std::ostream& out) const;
    void myClear() { myOutSOI.clear(); myOutBBasis.clear(); myOutAlmostVanishing.clear(); myInPts.clear(); myInTolerance.clear(); myInGammaPtr.reset(); }
  private:
    unique_ptr<ring> myPolyRingPtr;
    vector<ApproxPts::PointR> myInPts;
    vector<RingElem> myInTolerance;
    unique_ptr<RingElem> myInGammaPtr; /// BUG BUG BUG  cannot (yet?) use a plain RingElem (causes "race cond" type problems with RingQQ)
    vector<PPMonoidElem> myOutSOI;
    vector<RingElem> myOutBBasis;
    vector<RingElem> myOutAlmostVanishing;
  };


  void StableBorder::myWriteResult(std::ostream& out) const
  {
    out << ourVarName4 << " := Record[\n"
      " SOI := " << myOutSOI << ",\n"
      " AlmostVanishing := " << myOutAlmostVanishing << ",\n";
    if (myOutBBasis.empty())
      out << " StableBBasisFound := FALSE";
    else
      out << " StableBBasisFound := TRUE,\n"
             " BBasis := " << myOutBBasis;

    out << "\n];" << endl;
  }


  void StableBorder::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 3+1);

    myPolyRingPtr.reset(new SparsePolyRing(ReadPolyRing(in, GetTag)));
    matrix  PtsMat = ReadRationalMatrix(in, GetTag);
    matrix  TolMat  = ReadRationalMatrix(in, GetTag);
    matrix gamma = ReadRationalMatrix(in, GetTag);
    if (NumRows(TolMat) != 1) CoCoA_THROW_ERROR("NumRows(TolMat) should be 1","StableBorder");
    if (NumCols(TolMat) != NumCols(PtsMat))
      CoCoA_THROW_ERROR("PtsMat and TolMat should have same NumCols","StableBorder");
    if (NumRows(gamma) != 1 || NumCols(gamma) != 1)
      CoCoA_THROW_ERROR("gamma should be 1x1 matrix","StableBorder");

    vector<ApproxPts::PointR> pts(NumRows(PtsMat), ApproxPts::PointR(NumCols(PtsMat), zero(RingQQ())));
    for (long i=0; i < NumRows(PtsMat); ++i)
      for (long j=0; j < NumCols(PtsMat); ++j)
        pts[i][j] = PtsMat(i,j);

    vector<RingElem> tolerance(NumCols(PtsMat), zero(RingQQ()));
    for (long j=0; j < NumCols(PtsMat); ++j)
      tolerance[j] = TolMat(0,j);

    myInGammaPtr.reset(new RingElem(gamma(0,0)));
    swap(myInPts, pts);
    swap(myInTolerance, tolerance);
  }



 // ---- CoCoA/ApproxPts.H  by  L.Torrente and J.Abbott ----
  class NumBMBorder: public ServerOpBase
  {
  public:
    NumBMBorder(): ServerOpBase(CoCoALib_approx()) {};
    ~NumBMBorder() {};
    void myOutputSelf(std::ostream& out) const { out << "NumBMBorder"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute()  { ApproxPts::NBMTwinFloat(myOutQB, myOutBBasis, myOutAlmostVanishing, myInPts, myInTolerance); }
    void myWriteResult(std::ostream& out) const;
    void myClear() { myOutQB.clear(); myOutBBasis.clear(); myOutAlmostVanishing.clear(); myInPts.clear(); myInTolerance.clear(); }
  private:
    vector<ApproxPts::PointR> myInPts;
    vector<RingElem> myInTolerance;
    vector<PPMonoidElem> myOutQB;
    vector<RingElem> myOutBBasis;
    vector<RingElem> myOutAlmostVanishing;
  };


  void NumBMBorder::myWriteResult(std::ostream& out) const
  {
    out << ourVarName4 << " := Record[\n"
      " QB := " << myOutQB << ",\n"
      " AlmostVanishing := " << myOutAlmostVanishing << ",\n";
    if (myOutBBasis.empty())
      out << " StableBBasisFound := FALSE";
    else
      out << " StableBBasisFound := TRUE,\n"
             " BBasis := " << myOutBBasis;

    out << "\n];" << endl;
  }


  void NumBMBorder::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2+1);

    /*const SparsePolyRing NoUse = */ ReadPolyRing(in, GetTag);
    matrix  PtsMat = ReadRationalMatrix(in, GetTag);
    matrix  TolMat  = ReadRationalMatrix(in, GetTag);
    if (NumRows(TolMat) != 1) CoCoA_THROW_ERROR("NumRows(TolMat) should be 1","StableBorder");
    if (NumCols(TolMat) != NumCols(PtsMat))
      CoCoA_THROW_ERROR("PtsMat and TolMat should have same NumCols","StableBorder");

    vector<ApproxPts::PointR> pts(NumRows(PtsMat), ApproxPts::PointR(NumCols(PtsMat), zero(RingQQ())));
    for (long i=0; i < NumRows(PtsMat); ++i)
      for (long j=0; j < NumCols(PtsMat); ++j)
        pts[i][j] = PtsMat(i,j);

    vector<RingElem> tolerance(NumCols(PtsMat), zero(RingQQ()));
    for (long j=0; j < NumCols(PtsMat); ++j)
      tolerance[j] = TolMat(0,j);

    swap(myInPts, pts);
    swap(myInTolerance, tolerance);
  }


  // ---- CoCoA/ApproxPts.H  by  L.Torrente and J.Abbott ----
  class PreprocessBase: public ServerOpBase
  {
  protected:
    PreprocessBase(): ServerOpBase(CoCoALib_approx()) {};
    virtual ~PreprocessBase() {};
    virtual void myOutputSelf(std::ostream& out) const { out << myAlgName(); }
    virtual void myReadArgs(std::istream& in, int NumArgs);
    virtual void myWriteResult(std::ostream& out) const;
    virtual void myClear() { myInPts.clear(); myInTolerance.clear(); myOutPts.clear(); myOutWeights.clear(); }
    virtual const char* myAlgName() const = 0;
  protected:
    vector<ApproxPts::PointR> myInPts, myOutPts;
    vector<RingElem> myInTolerance;
    vector<long> myOutWeights;
  };


  void PreprocessBase::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 3);
    /*const SparsePolyRing NoUse = */ ReadPolyRing(in, GetTag);
    matrix  PtsMat = ReadRationalMatrix(in, GetTag);
    matrix  TolMat  = ReadRationalMatrix(in, GetTag);
    if (NumRows(TolMat) != 1) CoCoA_THROW_ERROR("NumRows(TolMat) should be 1", myAlgName());
    if (NumCols(TolMat) != NumCols(PtsMat))
      CoCoA_THROW_ERROR("PtsMat and TolMat should have same NumCols", myAlgName());

    vector<ApproxPts::PointR> pts(NumRows(PtsMat), ApproxPts::PointR(NumCols(PtsMat), zero(RingQQ())));
    for (long i=0; i < NumRows(PtsMat); ++i)
      for (long j=0; j < NumCols(PtsMat); ++j)
        pts[i][j] = PtsMat(i,j);

    vector<RingElem> tolerance(NumCols(PtsMat), zero(RingQQ()));
    for (long j=0; j < NumCols(PtsMat); ++j)
      tolerance[j] = TolMat(0,j);

    swap(myInTolerance, tolerance);
    swap(myInPts, pts);
  }


  void PreprocessBase::myWriteResult(std::ostream& out) const
  {
    out << ourVarName4 << " := Record[\n"
      " Points := " << myOutPts << ",\n"
      " Weights := " << myOutWeights << "\n];" << endl;
  }


  // ---- CoCoA/ApproxPts.H  by  L.Torrente and J.Abbott ----
  class PreprocessAggr: public PreprocessBase
  {
  public:
    void myCompute()  { PreprocessPtsAggr(myOutPts, myOutWeights, myInPts, myInTolerance); }
    const char* myAlgName() const { return "PreprocessAggr"; }
  };

  // ---- CoCoA/ApproxPts.H  by  L.Torrente and J.Abbott ----
  class PreprocessSubdiv: public PreprocessBase
  {
  public:
    void myCompute()  { PreprocessPtsSubdiv(myOutPts, myOutWeights, myInPts, myInTolerance); }
    const char* myAlgName() const { return "PreprocessSubdiv"; }
  };

  // ---- CoCoA/ApproxPts.H  by  L.Torrente and J.Abbott ----
  class PreprocessGrid: public PreprocessBase
  {
  public:
    void myCompute()  { PreprocessPtsGrid(myOutPts, myOutWeights, myInPts, myInTolerance); }
    const char* myAlgName() const { return "PreprocessGrid"; }
  };



// ---- CoCoA/TmpMayerVietorisTree.H  by  E. Saenz-de-Cabezon ----

class MVT: public ServerOpBase
  {
  public:
    MVT(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVT() {};
    void myOutputSelf(std::ostream& out) const { out << "MV_Tree"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { MayerVietorisTree(*myOutMdMpPtr, *myInPPsPtr); }
    void myWriteResult(std::ostream& out) const;
    void myClear() { myInPPsPtr.reset(); myOutMdMpPtr.reset(); }
  private:
    unique_ptr<PPVector> myInPPsPtr;
    unique_ptr<MultidegreeMap>  myOutMdMpPtr;
  };


  void MVT::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and PPList
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    //    PPOrdering ord = NewLexOrdering(NumIndets(P));
    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x", 0, NumIndets(P)-1), lex);
    vector<PPMonoidElem> V;
    ReadPPs(in, V, PPM, GetTag);
    DivMaskRule DMR = NewDivMaskEvenPowers();
    vector<long> exps(NumIndets(PPM));
    myInPPsPtr.reset(new PPVector(PPM, DMR));
    myOutMdMpPtr.reset(new MultidegreeMap);

    for (long i=0 ; i<len(V); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPsPtr->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    convert(*myInPPsPtr, V);  // ANNA: fix this!
    
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }


  void MVT::myWriteResult(std::ostream& out) const
  { MultidegreeMap::const_iterator pos;
    out << ourVarName4 << " := [];";
    for (pos= (*myOutMdMpPtr).begin(); pos!=(*myOutMdMpPtr).end();++pos)
	{
	out<< "Append(" << ourVarName4<< ", ["<<pos->first;
	out<<",[";
	ListOfDims::const_iterator inner_pos;
 	for (inner_pos= (pos->second).begin(); inner_pos!=(pos->second).end();++inner_pos)
  		{
 		out<<"["<<inner_pos->first;
		list<position_t>::const_iterator my_pos;
		out<<",[";
		for (my_pos= (inner_pos->second).begin(); my_pos!=(inner_pos->second).end();++my_pos)
                {
		  ++my_pos;
		  if(my_pos==(inner_pos->second).end())
			{--my_pos;  out<<*my_pos;}
		  else
			{--my_pos;  out<<*my_pos<<",";}
                }
		++inner_pos;
		if(inner_pos==(pos->second).end())
			{--inner_pos;out<<"]]";}
		  else
			{--inner_pos;  out<<"]],";}
  		}
	out<<"]";
	out<< "]);" <<endl;
	}
  }

// void PrintMultidegreeMap(const MultidegreeMap& myMap)
// {
// MultidegreeMap::const_iterator pos;
// 
// }

  class MVTN1: public ServerOpBase
  {
  public:
    MVTN1(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVTN1() {};
    void myOutputSelf(std::ostream& out) const { out << "MVT_N-1"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { MayerVietorisTreeN1(*myOutPPsPtr, *myInPPsPtr); }
    void myWriteResult(std::ostream& out) const;
    void myClear() { myInPPsPtr.reset(); myOutPPsPtr.reset(); }
  private:
    unique_ptr<PPVector> myInPPsPtr, myOutPPsPtr;
  };


  void MVTN1::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and PPList
    const SparsePolyRing P(ReadPolyRing(in, GetTag));
    PPMonoid PPM = NewPPMonoidEv(SymbolRange("x", 0, NumIndets(P)-1), lex);
    vector<PPMonoidElem> V;
    ReadPPs(in, V, PPM, GetTag);
    DivMaskRule DMR = NewDivMaskEvenPowers();
    vector<long> exps(NumIndets(PPM));
    myInPPsPtr.reset(new PPVector(PPM, DMR));
    myOutPPsPtr.reset(new PPVector(PPM, DMR));

    for (long i=0; i < len(V); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPsPtr->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }


  void MVTN1::myWriteResult(std::ostream& out) const
  {
    out << ourVarName4 << " := [];";
    for (long i=0; i < len(*myOutPPsPtr); ++i)
      out<< "Append(" << ourVarName4<< ", "<< PP((*myOutPPsPtr)[i]) << ");" <<endl;
  }


  class MVTReg: public ServerOpBase
  {
  public:
    MVTReg(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVTReg() {};
    void myOutputSelf(std::ostream& out) const { out << "MVT_Reg"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutReg = MVTRegularity(*myInPPV); }
    void myWriteResult(std::ostream& out) const { out << ourVarName4 << " := " << myOutReg << ";" << endl; }
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    int myOutReg;
  };

  void MVTReg::myReadArgs(std::istream& in, int NumArgs)
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

    for (long i=0; i<len(V); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }

class MVTRegUB: public ServerOpBase
  {
  public:
    MVTRegUB(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVTRegUB() {};
    void myOutputSelf(std::ostream& out) const { out << "MVT_Reg_Upper_Bound"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutReg = MVTRegularityUpperBound(*myInPPV); }
    void myWriteResult(std::ostream& out) const { out << ourVarName4 << " := " << myOutReg << ";" << endl; }
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    int myOutReg;
  };

  void MVTRegUB::myReadArgs(std::istream& in, int NumArgs)
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

    for (long i=0; i<len(V); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }

class MVTRegLB: public ServerOpBase
  {
  public:
    MVTRegLB(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVTRegLB() {};
    void myOutputSelf(std::ostream& out) const { out << "MVT_Reg_Lower_Bound"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutReg = MVTRegularityLowerBound(*myInPPV); }
    void myWriteResult(std::ostream& out) const { out << ourVarName4 << " := " << myOutReg << ";" << endl; }
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    int myOutReg;
  };

  void MVTRegLB::myReadArgs(std::istream& in, int NumArgs)
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

    for (long i=0; i<len(V); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }


class MVTProjDimension: public ServerOpBase
  {
  public:
    MVTProjDimension(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVTProjDimension() {};
    void myOutputSelf(std::ostream& out) const { out << "MVTProjDimension"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutProjDim = MVTProjDim(*myInPPV); }
    void myWriteResult(std::ostream& out) const { out << ourVarName4 << " := " << myOutProjDim << ";" << endl; }
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    int myOutProjDim;
  };

  void MVTProjDimension::myReadArgs(std::istream& in, int NumArgs)
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

    for (long i=0; i<len(V); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }



class MVTProjDimUB: public ServerOpBase
  {
  public:
    MVTProjDimUB(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVTProjDimUB() {};
    void myOutputSelf(std::ostream& out) const { out << "MVT_ProjDim_Upper_Bound"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutProjDim = MVTProjDimUpperBound(*myInPPV); }
    void myWriteResult(std::ostream& out) const { out << ourVarName4 << " := " << myOutProjDim << ";" << endl; }
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    int myOutProjDim;
  };

  void MVTProjDimUB::myReadArgs(std::istream& in, int NumArgs)
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

    for (long i=0 ; i<len(V); ++i)
    {
      exponents(exps, V[i]);
      // anna: should find a way to avoid making copies the pps
      myInPPV->myPushBack(PPMonoidElem(PPM,exps));
    }
    //    GlobalOutput() << "Print " << myInPPsPtr->mySize() << ";";
  }


class MVTProjDimLB: public ServerOpBase
  {
  public:
    MVTProjDimLB(): ServerOpBase(CoCoALib_combinatorics()) {};
    ~MVTProjDimLB() {};
    void myOutputSelf(std::ostream& out) const { out << "MVT_ProjDim_Lower_Bound"; }
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute() { myOutProjDim = MVTProjDimLowerBound(*myInPPV); }
    void myWriteResult(std::ostream& out) const { out << ourVarName4 << " := " << myOutProjDim << ";" << endl; }
    void myClear() { myInPPV.reset(); }
  private:
    unique_ptr<PPVector> myInPPV;
    int myOutProjDim;
  };

  void MVTProjDimLB::myReadArgs(std::istream& in, int NumArgs)
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

    for (long i=0 ; i<len(V); ++i)
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



  namespace CoCoAServerOperationsFromCoCoALib
  {
    bool RegisterOps()
    {
      //      RegisterOp("TestSocket",           ServerOp(new TestSocket()));
      // ideal
      RegisterOp("F5",                   ServerOp(new F5GBasis()));
      RegisterOp("ideal_LT",             ServerOp(new IdealLT()));
      RegisterOp("ideal_colon_ideal",    ServerOp(new ColonIdId()));
      RegisterOp("ideal_elim",           ServerOp(new IdealElim()));
      RegisterOp("ideal_groebner",       ServerOp(new IdealGBasis()));
      RegisterOp("ideal_sat_groebner",   ServerOp(new IdealSATGBasis()));
      //RegisterOp("ideal_satmix_groebner", ServerOp(new IdealSATMixGBasis()));
      RegisterOp("ideal_homogenization", ServerOp(new IdealHomog()));
      RegisterOp("ideal_intersection",   ServerOp(new IdealIntersection()));
      RegisterOp("ideal_saturation",     ServerOp(new IdealSaturation()));
      RegisterOp("ideal_syzygy",         ServerOp(new IdealSyzygy()));
      // module
      RegisterOp("module_LT",            ServerOp(new ModuleLT()));
      RegisterOp("module_colon_module",  ServerOp(new ColonModMod()));
      RegisterOp("module_elim",          ServerOp(new ModuleElim()));
      RegisterOp("module_groebner",      ServerOp(new ModuleGBasis()));
      RegisterOp("module_homogenization", ServerOp(new ModuleHomog()));
      RegisterOp("module_intersection",  ServerOp(new ModuleIntersection()));
      RegisterOp("module_saturation",    ServerOp(new ModuleSaturation()));
      RegisterOp("module_syzygy",        ServerOp(new ModuleSyzygy()));
      // approx points
      RegisterOp("StableBorder",         ServerOp(new StableBorder()));
      RegisterOp("NumBMBorder",          ServerOp(new NumBMBorder()));
      RegisterOp("PreprocessAggr",       ServerOp(new PreprocessAggr()));
      RegisterOp("PreprocessGrid",       ServerOp(new PreprocessGrid()));
      RegisterOp("PreprocessSubdiv",     ServerOp(new PreprocessSubdiv()));
      // IsTree
      RegisterOp("IsTree_NoOpt",         ServerOp(new IsTreeNoOpt()));
      RegisterOp("IsTree_Opt",           ServerOp(new IsTreeOpt()));
      RegisterOp("IsTree_CBOpt",         ServerOp(new IsTreeCBOpt()));
      RegisterOp("IsTree_CBNoOpt",       ServerOp(new IsTreeCBNoOpt()));
      // Mayer Vietoris Trees
      RegisterOp("Mayer_Vietoris_Tree",  ServerOp(new MVT()));
      RegisterOp("MVT_N_minus_one",      ServerOp(new MVTN1()));
      RegisterOp("MVT_Regularity",       ServerOp(new MVTReg()));
      RegisterOp("MVT_Regularity_Upper_Bound", ServerOp(new MVTRegUB()));
      RegisterOp("MVT_Regularity_Lower_Bound", ServerOp(new MVTRegLB()));
      RegisterOp("MVT_ProjDim",          ServerOp(new MVTProjDimension()));
      RegisterOp("MVT_ProjDim_Upper_Bound",    ServerOp(new MVTProjDimUB()));
      RegisterOp("MVT_ProjDim_Lower_Bound",    ServerOp(new MVTProjDimLB()));
      return true;
    }


    bool RegisterOpsOnce()
    {
      static bool EvalOnce = RegisterOps();
      return EvalOnce;
    }
  }


}

// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/RegisterServerOps.C,v 1.18 2024/07/14 10:29:47 abbott Exp $
// $Log: RegisterServerOps.C,v $
// Revision 1.18  2024/07/14 10:29:47  abbott
// Summary: Changed call to F5() into call to F5_poly (hope this is correct)
//
// Revision 1.17  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.16  2020/06/17 15:49:31  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.15  2019/03/20 16:20:54  bigatti
// -- updates for the new server
//
// Revision 1.14  2019/03/19 11:05:02  abbott
// Summary: In StableBorder replaced data mem myPolyRing by myPolyRingPtr (o/w tries to create ring before GlobalManager exists)
//
// Revision 1.13  2019/03/04 11:37:14  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.12  2017/11/28 13:27:30  bigatti
// -- fixed IdealSATGBasis
// -- commented out IdealSATMixGBasis
//
// Revision 1.11  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.10  2017/04/18 16:23:20  bigatti
// -- removed use of stats level (to be replaced by VerbosityLevel)
//
// Revision 1.9  2017/04/07 20:06:41  abbott
// Summary: Hacked to make it compile (added member myPolyRing); no idea if it works!!!
//
// Revision 1.8  2015/12/04 15:22:24  bigatti
// -- renamed ComputeSSaturation into ComputeSaturation
//
// Revision 1.7  2015/05/20 13:27:08  bigatti
// -- ComputeCColon --> ComputeColon
//
// Revision 1.6  2014/07/30 15:07:13  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.5  2014/07/08 09:18:04  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.4  2014/05/15 12:31:34  abbott
// Summary: Now using new files server/GlobalIO.HC (previously in CoCoA/io.H)
// Author: JAA
//
// Revision 1.3  2014/03/26 16:53:38  bigatti
// -- added unused arg in GBasis
//
// Revision 1.2  2013/06/12 08:55:48  bigatti
// -- added unused arg (in ComputeGBasis)
//
// Revision 1.1  2013/05/27 12:57:39  abbott
// Moved all server-related code into src/server/
//
// Revision 1.46  2013/03/27 18:24:33  abbott
// Added approx point preprocessing to C5; also changed names of the fns, and updated doc.
//
// Revision 1.45  2013/02/21 17:15:06  bigatti
// -- changes syntax for ComputeSyz
//
// Revision 1.44  2012/10/02 10:36:12  abbott
// Revised interface to BuildInfo information strings.
// Several consequential changes.
//
// Revision 1.43  2012/02/08 17:09:42  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.42  2011/03/11 11:05:24  bigatti
// -- changed size_t --> long
// -- changed size --> len
//
// Revision 1.41  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.40  2010/09/15 21:24:04  abbott
// Corrected assertion for number of args to StableBorder.
//
// Revision 1.39  2010/02/03 18:05:15  bigatti
// -- minor fix
//
// Revision 1.38  2010/02/03 17:04:53  bigatti
// -- first use of "lex" in definition of PPMonoidEv
// -- added commented out code for TestSocket (wip)
//
// Revision 1.37  2009/11/03 17:50:16  bigatti
// -- added MVTProjDim functions by Eduardo Saenz-de-Cabezon
// -- code cleaning
//
// Revision 1.36  2009/11/03 14:57:46  abbott
// Inserted some space to improve code readability.
// Replaced several endls with '\n' to avoid needless flushing,
// and removed some << operators between string literals.
//
// Revision 1.35  2009/10/28 16:09:09  bigatti
// -- modified class name MVTProjDim into MVTProjDimension to avoid name clashing
//
// Revision 1.34  2009/10/28 14:46:37  bigatti
// -- new MVT functions by Eduardo Saenz-de-Cabezon
//
// Revision 1.33  2009/10/27 13:06:12  bigatti
// -- just a "\n"
//
// Revision 1.32  2009/10/26 17:21:27  bigatti
// -- trying to avoid sending small strings to sockets...
//
// Revision 1.31  2009/09/24 15:23:01  abbott
// Added some missing using commands.
//
// Revision 1.30  2009/09/24 14:46:06  abbott
// Removed some unnecessary "std::" prefixes.
// Changed two "unsigned int" into "size_t".
//
// Revision 1.29  2009/07/30 15:46:19  bigatti
// -- IdealLT::myWriteResult using new ideal ctor
//
// Revision 1.28  2009/07/24 14:48:03  abbott
// Corrected an assertion about number of args (for NBM).
//
// Revision 1.27  2009/07/06 12:32:33  abbott
// Removed some useless consts on return type for myAlgName
// (should keep the compiler quiet in debugging mode).
//
// Revision 1.26  2009/04/23 15:20:48  bigatti
// -- changed "=" into ":=" for Records
// -- fixed LT5(Ideal(x))  (WriteIdeal looks at item [0], was [1])
// -- fixed number of arg for approx functions
//
// Revision 1.25  2009/01/26 15:57:48  bigatti
// -- added "const" to libraries
// -- first use of WriteIdeal
//
// Revision 1.24  2009/01/08 13:11:08  bigatti
// -- fixed MVTN1 bug (myOutPPsPtr is now initialized in myReadArgs)
//
// Revision 1.23  2008/11/24 17:11:12  abbott
// Final tidying of code for preprocessing and SOI/NBM.
//
// Revision 1.22  2008/11/23 18:58:32  abbott
// Major overhaul to preprocessing and SOI/NBM code.
// Split SOI/NBM off into a separate file.
// Preprocessing is now "rational" (but internally guided by C++ doubles).
// SOI/NBM now each have 3 similar interfaces: one purely rational, one for
// input which is represented as doubles, and one which converts the input
// to RingTwinFloat values and produces a result which is over some RingTwinFloat
// (the precision is increased automatically until an answer is obtained).
//
// Revision 1.21  2008/11/21 21:17:35  abbott
// Added 3 new "rational" preprocessing algms.
//
// Revision 1.20  2008/11/20 09:58:53  abbott
// Cleaned up code for Preprocessing: introduced another base class.
// Change the way it sends results to CoCoA4 -- sends both points and weights.
//
// Revision 1.19  2008/11/13 12:13:17  bigatti
// -- NumArgs is now used in every myReadArgs, but usually just a dumb check
//
// Revision 1.18  2008/09/22 16:42:48  bigatti
// -- small fix
//
// Revision 1.17  2008/09/22 16:07:02  bigatti
// -- tested (and fixed) communication with cocoa-4 (number of arguments
//    passed to myReadArgs   )
//
// Revision 1.16  2008/09/19 12:20:22  bigatti
// -- fix call to SatGB
//
// Revision 1.15  2008/09/19 11:41:23  bigatti
// -- first fix for new mechanism for passing verbosity
//
// Revision 1.14  2008/09/19 11:34:15  bigatti
// -- new mechanism for passing verbosity level (or StatLevel)
//    [only partially tested]
//
// Revision 1.13  2008/09/12 13:28:43  bigatti
// -- new: NBM implementation
//
// Revision 1.12  2008/07/04 09:11:04  bigatti
// -- new PPVector class
//
// Revision 1.11  2008/06/04 18:27:37  abbott
// Modified the server interface for "SOI": it now accepts a 3rd arg (gamma).
//
// Revision 1.10  2008/05/30 14:20:43  abbott
// SOI now returns also the "almost vanishing" polynomials.
//
// Revision 1.9  2008/05/29 15:46:29  bigatti
// -- added Approximate Border Basis (by Abbott,Torrente)
//
// Revision 1.8  2008/05/28 16:21:05  bigatti
// -- using the new function ReadPPs fom MVTN1
//
// Revision 1.7  2008/05/27 16:22:04  bigatti
// -- added MayerVietorisTreeN1
//
// Revision 1.6  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.5  2007/11/09 10:45:52  bigatti
// -- [caboara] preparation for self-saturating algorithm
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/09/25 16:32:30  abbott
// Several minor changes to silence gcc-4.3:
//    more #includes,
//    and fixed a template problemm in RegisterServerOps.C
//
// Revision 1.2  2007/09/25 16:28:31  bigatti
// -- ServerOp includes infos about the library it is defined in
// -- CoCoAServer no longer prints its own version
// -- CoCoAServer prints all offered operations
// -- CoCoAServer can change stat_level (Max's verbosity) at runtime
//
// Revision 1.1  2007/04/27 14:54:22  bigatti
// -- content of CoCoAServer.C split into dedicated files
// -- new registration mechanism (through include "RegisterServerOps.H")
//
