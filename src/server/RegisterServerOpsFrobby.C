//   Copyright (c)  2009  John Abbott,  Anna M Bigatti
//   Orig authors: (2009) Anna Bigatti and Bjarke Hammersholt Roune

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

#ifdef CoCoA_WITH_FROBBY
#include "CoCoA5io.H"
#include "ServerOp.H"

// #include <iostream>
using std::endl;
using std::clog;
#include <sstream>
using std::stringstream;
// #include <memory>
using std::unique_ptr;
// #include <string>
using std::string;
// #include <vector>
using std::vector;

namespace CoCoA
{
  // sublibrary of CoCoALib for integration with Frobby
  // by Bjarke Hammersholt Roune.
  const ServerOpBase::LibraryInfo& CoCoALib_frobby()
  {
    static const ServerOpBase::LibraryInfo
      UniqueValue("CoCoALib", BuildInfo::version(), "frobby");
    return UniqueValue;
  }

  // ---- CoCoA/ExternalLibs-Frobby.H by Bjarke Hammersholt Roune

  // Common base class for Frobby operations. Contains some code that is
  // useful for Frobby operations.
  class FrobbyOpBase : public ServerOpBase {
  public:
	FrobbyOpBase(): ServerOpBase(CoCoALib_frobby()) {
	}

  protected:
	unique_ptr<ideal> myReadMonomialIdeal
	  (std::istream& in, const SparsePolyRing& ring);
	unique_ptr<PPMonoidElem> myReadMonomial
	  (std::istream& in, const SparsePolyRing& ring);
  };

  unique_ptr<ideal> FrobbyOpBase::myReadMonomialIdeal
  (std::istream& in, const SparsePolyRing& ring) {
    PolyList polyList;
    ReadPolyList(in, polyList, ring, GetTag);
    return unique_ptr<ideal>(new ideal(ring, polyList));
  }

  unique_ptr<PPMonoidElem> FrobbyOpBase::myReadMonomial
  (std::istream& in, const SparsePolyRing& ring) {
    PolyList polyList;
    ReadPolyList(in, polyList, ring, GetTag);

    if (polyList.size() != 1 ||	!IsMonomial(polyList.front()))
      CoCoA_THROW_ERROR("Expected a single monomial", "FrobbyBaseOp::myReadMonomial");

    return unique_ptr<PPMonoidElem>(new PPMonoidElem(LPP(polyList.front())));
  }

#ifdef CoCoA_WITH_FROBBY

  class AlexanderDualOp : public FrobbyOpBase
  {
  public:
    void myOutputSelf(std::ostream& out) const;
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute();
    void myWriteResult(std::ostream& out) const;
    void myClear();

  private:
    unique_ptr<ideal> myIdeal;
    unique_ptr<PPMonoidElem> myPP;
    unique_ptr<ideal> myDual;
  };

  void AlexanderDualOp::myOutputSelf(std::ostream& out) const {
	out << "AlexanderDualFrobby";
  }

  void AlexanderDualOp::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 3);  // ring, ideal, pp
    const SparsePolyRing ring(ReadPolyRing(in, GetTag));
	myIdeal = myReadMonomialIdeal(in, ring);
	myPP = myReadMonomial(in, ring);
  }

  void AlexanderDualOp::myCompute() {
	myDual.reset(new ideal(FrbAlexanderDual(*myIdeal, *myPP)));
  }

  void AlexanderDualOp::myWriteResult(std::ostream& out) const
  {
    WritePolyListInVar(out, ourVarName4, gens(*myDual));
  }

  void AlexanderDualOp::myClear() {
	myIdeal.reset();
	myPP.reset();
	myDual.reset();
  }

  class IrreducibleDecomOp : public FrobbyOpBase
  {
  public:
    void myOutputSelf(std::ostream& out) const;
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute();
    void myWriteResult(std::ostream& out) const;
    void myClear();

  private:
    unique_ptr<ideal> myIdeal;
    unique_ptr<vector<ideal> > myDecom;
  };

  void IrreducibleDecomOp::myOutputSelf(std::ostream& out) const {
	out << "IrreducibleDecom";
  }

  void IrreducibleDecomOp::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and ideal
    const SparsePolyRing ring(ReadPolyRing(in, GetTag));
	myIdeal = myReadMonomialIdeal(in, ring);
  }

  void IrreducibleDecomOp::myCompute() {
	myDecom.reset(new vector<ideal>());
	FrbIrreducibleDecomposition(*myDecom, *myIdeal);
  }

  void IrreducibleDecomOp::myWriteResult(std::ostream& out) const
  {
	out << ourVarName4 << " := [];\n";
	for (size_t component = 0; component < myDecom->size(); ++component) {
	  WritePolyListInVar(out, "Tmp", gens((*myDecom)[component]));
	  out << "Append(" << ourVarName4 << ", Tmp);\n";
	}
  }

  void IrreducibleDecomOp::myClear() {
	myIdeal.reset();
	myDecom.reset();
  }

  class MaximalStandardMonomialsOp : public FrobbyOpBase
  {
  public:
    void myOutputSelf(std::ostream& out) const;
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute();
    void myWriteResult(std::ostream& out) const;
    void myClear();

  private:
    unique_ptr<ideal> myIdeal;
    unique_ptr<ideal> myMsms;
  };

  void MaximalStandardMonomialsOp::myOutputSelf(std::ostream& out) const {
	out << "MaximalStandardMonomials";
  }

  void MaximalStandardMonomialsOp::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and ideal
    const SparsePolyRing ring(ReadPolyRing(in, GetTag));
	myIdeal = myReadMonomialIdeal(in, ring);
  }

  void MaximalStandardMonomialsOp::myCompute() {
	myMsms.reset(new ideal(FrbMaximalStandardMonomials(*myIdeal)));
  }

  void MaximalStandardMonomialsOp::myWriteResult(std::ostream& out) const
  {
    WritePolyListInVar(out, ourVarName4, gens(*myMsms));
  }

  void MaximalStandardMonomialsOp::myClear() {
	myIdeal.reset();
	myMsms.reset();
  }

  class PrimaryDecomOp : public FrobbyOpBase
  {
  public:
    void myOutputSelf(std::ostream& out) const;
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute();
    void myWriteResult(std::ostream& out) const;
    void myClear();

  private:
    unique_ptr<ideal> myIdeal;
    unique_ptr<vector<ideal> > myDecom;
  };

  void PrimaryDecomOp::myOutputSelf(std::ostream& out) const {
	out << "PrimaryDecompositionFrobby";
  }

  void PrimaryDecomOp::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and ideal
    const SparsePolyRing ring(ReadPolyRing(in, GetTag));
	myIdeal = myReadMonomialIdeal(in, ring);
  }

  void PrimaryDecomOp::myCompute() {
	myDecom.reset(new vector<ideal>());
	FrbPrimaryDecomposition(*myDecom, *myIdeal);
  }

  void PrimaryDecomOp::myWriteResult(std::ostream& out) const
  {
	out << ourVarName4 << " := [];\n";
	for (size_t component = 0; component < myDecom->size(); ++component) {
	  WritePolyListInVar(out, "Tmp", gens((*myDecom)[component]));
	  out << "Append(" << ourVarName4 << ", Tmp);\n";
	}
  }

  void PrimaryDecomOp::myClear() {
	myIdeal.reset();
	myDecom.reset();
  }

  class AssociatedPrimesOp : public FrobbyOpBase
  {
  public:
    void myOutputSelf(std::ostream& out) const;
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute();
    void myWriteResult(std::ostream& out) const;
    void myClear();

  private:
    unique_ptr<ideal> myIdeal;
    unique_ptr<vector<ideal> > myPrimes;
  };

  void AssociatedPrimesOp::myOutputSelf(std::ostream& out) const {
	out << "AssociatedPrimes";
  }

  void AssociatedPrimesOp::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and ideal
    const SparsePolyRing ring(ReadPolyRing(in, GetTag));
	myIdeal = myReadMonomialIdeal(in, ring);
  }

  void AssociatedPrimesOp::myCompute() {
	myPrimes.reset(new vector<ideal>());
	FrbAssociatedPrimes(*myPrimes, *myIdeal);
  }

  void AssociatedPrimesOp::myWriteResult(std::ostream& out) const
  {
	out << ourVarName4 << " := [];\n";
	for (size_t component = 0; component < myPrimes->size(); ++component) {
	  WritePolyListInVar(out, "Tmp", gens((*myPrimes)[component]));
	  out << "Append(" << ourVarName4 << ", Tmp);\n";
	}
  }

  void AssociatedPrimesOp::myClear() {
	myIdeal.reset();
	myPrimes.reset();
  }

  class DimensionOp : public FrobbyOpBase
  {
  public:
	DimensionOp();

    void myOutputSelf(std::ostream& out) const;
    void myReadArgs(std::istream& in, int NumArgs);
    void myCompute();
    void myWriteResult(std::ostream& out) const;
    void myClear();

  private:
    unique_ptr<ideal> myIdeal;
	long myDimension;
  };

  DimensionOp::DimensionOp():
	myDimension(0) {
  }

  void DimensionOp::myOutputSelf(std::ostream& out) const {
	out << "Dimension";
  }

  void DimensionOp::myReadArgs(std::istream& in, int NumArgs)
  {
    CoCoA_ASSERT(NumArgs == 2);  // ring and ideal
    const SparsePolyRing ring(ReadPolyRing(in, GetTag));
	myIdeal = myReadMonomialIdeal(in, ring);
  }

  void DimensionOp::myCompute() {
	myDimension = FrbDimension(*myIdeal);
  }

  void DimensionOp::myWriteResult(std::ostream& out) const
  {
	out << ourVarName4 << " := " << myDimension << ";\n";
  }

  void DimensionOp::myClear() {
	myIdeal.reset();
	myDimension = 0;
  }
#else // if no Frobby
  class NoFrobbyErrorOp : public FrobbyOpBase
  {
  public:
    void myOutputSelf(std::ostream& out) const {myError();}
    void myReadArgs(std::istream& in, int NumArgs) {myError();}
    void myCompute() {myError();}
    void myWriteResult(std::ostream& out) const {myError();}
    void myClear() {myError();}

  private:
	void myError() const {
	  CoCoA_THROW_ERROR("Frobby not present. Build CoCoA with Frobby to enable "
				  "this function.", "NoFRobbyErrorOp");
	}
  };

  typedef NoFrobbyErrorOp AlexanderDualOp;
  typedef NoFrobbyErrorOp DimensionOp;
  typedef NoFrobbyErrorOp IrreducibleDecomOp;
  typedef NoFrobbyErrorOp PrimaryDecomOp;
  typedef NoFrobbyErrorOp MaximalStandardMonomialsOp;
  typedef NoFrobbyErrorOp AssociatedPrimesOp;

#endif

//----------------------------------------------------------------------

  namespace CoCoAServerOperationsFromFrobby
  {
    bool RegisterOps()
    {
      // integration with Frobby
      RegisterOp("AlexanderDual_Frobby", ServerOp(new AlexanderDualOp()));
      RegisterOp("Dimension_Frobby", ServerOp(new DimensionOp()));
      RegisterOp("IrreducibleDecom_Frobby", ServerOp(new IrreducibleDecomOp()));
      RegisterOp("PrimaryDecom_Frobby", ServerOp(new PrimaryDecomOp()));
      RegisterOp("MaximalStandardMonomials_Frobby",
                 ServerOp(new MaximalStandardMonomialsOp()));
      RegisterOp("AssociatedPrimes_Frobby", ServerOp(new AssociatedPrimesOp()));

      return true;
    }


    bool RegisterOpsOnce()
    {
      static bool EvalOnce = RegisterOps();
      return EvalOnce;
    }
  }
}
#endif
