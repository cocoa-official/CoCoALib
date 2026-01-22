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
    void myCompute() { myResPL = ComputeColon(myInVL1, myInVL2); }
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
    void myCompute() { myResPL = ComputeSaturation(myInPL1, myInPL2); }
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
    void myCompute() { myResPL = ComputeColon(myInPL1, myInPL2); }
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
    void myCompute() { myResPL = ComputeIntersection(myInPL1, myInPL2); }
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
    void myCompute() { myResPL = ComputeElim(myInPL, *myInElimPPPtr); }
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
