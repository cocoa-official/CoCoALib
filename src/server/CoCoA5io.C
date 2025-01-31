//   Copyright (c)  2019  Anna Bigatti,  John Abbott

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

#include "CoCoA5io.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpGReductor.H"
#include "CoCoA/degree.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/module.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"
#include "GlobalIO.H"

// #include<algorithm>
// #include<cstdlib>
// #include<fstream>
// #include<functional>
// #include<typeinfo>

//#include <vector>
using std::vector;
#include <iostream>
using std::endl;
using std::flush;
//#include <string>
using std::string;
#include<sstream>
using std::ostringstream;

namespace CoCoA
{
  string ReadOperationString(std::istream& in)
  {
    string s;
    SkipTag(in, "<operation>");
    in >> s;
    InputFailCheck(in, "ReadOperationString");
    SkipTag(in, "</operation>");
    return s;
  }


  RingElem ReadPoly(std::istream& in, const SparsePolyRing& P)
  {
    ring R = CoeffRing(P);
    long NumParams=0;
    const long NumInds = NumIndets(P);
    if (IsFractionField(R) && IsPolyRing(BaseRing(R)))
    {
      NumParams = NumIndets(BaseRing(R));
    }
    vector<long> v(NumInds);
    vector<long> ParV(NumParams);

    long NumSummands;
    in >> NumSummands;
    InputFailCheck(in, "ReadPoly -- NumSummands");
    BigInt IntCoeff; // to avoid ctor-dtor inside loop
    RingElem c(R), ans(P), mon(P);
    for (long NS=0; NS<NumSummands; ++NS)
    {
      in >> IntCoeff;
      for (long i=0 ; i<NumParams ; ++i)  in >> ParV[i];
      for (long i=0 ; i<NumInds ; ++i)  in >> v[i];
      InputFailCheck(in, "ReadPoly -- summand");
      if (NumParams==0)
        c = IntCoeff;
      else
      {
        const SparsePolyRing BR = BaseRing(R);
        c = CanonicalHom(BR,R)(monomial(BR, IntCoeff, ParV));
      }
      mon = monomial(P, c, v);
      P->myAddClear(raw(ans), raw(mon)); // ANNA: use geobucket
    }
    return ans;
  }


  PPMonoidElem ReadPP(std::istream& in, const PPMonoid& PPM)
  {
    long NumInds = NumIndets(PPM);
    vector<long> v(NumInds);
    for (long i=0; i<NumInds; ++i)  in >> v[i];
    InputFailCheck(in, "ReadPP");
    return PPMonoidElem(PPM, v);
  }


  int ReadVerbosityLevel(std::istream& in, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<verbosity_level>");
    int vl;
    in >> vl;
    InputFailCheck(in, "ReadVerbosityLevel");
    SkipTag(in, "</verbosity_level>");
    return vl;
  }


  matrix ReadRationalMatrix(std::istream& in, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<rational_matrix>");
    ring QQ = RingQQ();
    long rows, cols;
    in >> rows;
    in >> cols;
    InputFailCheck(in, "ReadRationalMatrix -- rows/cols");
    matrix M(NewDenseMat(QQ, rows, cols));
    BigInt N, D;
    for (long i = 0; i < rows; ++i)
      for (long j = 0; j < cols; ++j)
      {
        in >> N;
        in >> D;
        InputFailCheck(in, "ReadRationalMatrix -- entry");
        SetEntry(M, i, j, RingElem(QQ, N)/D);
      }
    SkipTag(in, "</rational_matrix>");
    return M;
  }


  matrix ReadIntegerMatrix(std::istream& in, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<integer_matrix>");
    ring ZZ = RingZZ();
    long rows, cols;
    in >> rows;
    in >> cols;
    InputFailCheck(in, "ReadIntegerMatrix -- rows/cols");
    matrix M(NewDenseMat(ZZ, rows, cols));
    BigInt N;
    for (long i = 0; i < rows; ++i)
      for (long j = 0; j < cols; ++j)
      {
        in >> N;
        InputFailCheck(in, "ReadIntegerMatrix -- entry");
        SetEntry(M, i, j, N);
      }
    SkipTag(in, "</integer_matrix>");
    return M;
  }


  // Expects in to be positioned on the number of polys in the PolyList
  void ReadPolyList(std::istream& in, PolyList& PL, const SparsePolyRing& P, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<polynomial_list>");
    long NumPolys;
    PL.clear();
    in >> NumPolys;
    InputFailCheck(in, "ReadPolyList -- NumPolys");
    PL.reserve(NumPolys);
    for (long i=0 ; i<NumPolys ; ++i)  PL.push_back(ReadPoly(in, P));
    SkipTag(in, "</polynomial_list>");
  }


  // Expects in to be positioned on the number of PPs
  void ReadPPs(std::istream& in, std::vector<PPMonoidElem>& PPs, const PPMonoid& PPM, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<powerproduct_list>");
    long NumPP;
    PPs.clear();
    in >> NumPP;
    InputFailCheck(in, "ReadPPs -- NumPP");
    PPs.reserve(NumPP);
    for (long i=0 ; i<NumPP ; ++i)  PPs.push_back(ReadPP(in, PPM));
    SkipTag(in, "</powerproduct_list>");
  }


  void ReadVectorList(std::istream& in, VectorList& VL, const FreeModule& FM, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<vector_list>");
    const SparsePolyRing P = RingOf(FM);
    ModuleElem TmpVector(FM);
    const vector<ModuleElem>& e = gens(FM);
    const long NumComponents = NumCompts(FM);

    VL.erase(VL.begin(), VL.end());
    //PL.reserve(NumPolys);
    long NumVectors;
    in >> NumVectors;
    InputFailCheck(in, "ReadVectorList -- NumVectors");
    for (long i=0 ; i<NumVectors ; ++i)
    {
      TmpVector = zero(FM);
      for (long j=0 ; j<NumComponents ; ++j)
      {
        long ComponentIndex;
        in >> ComponentIndex;  // ANNA: does this really make sense??
        if (ComponentIndex != j+1)
          CoCoA_THROW_ERROR("wrong component index","ReadVectorList");
        TmpVector += e[j]*ReadPoly(in, P);
      }
      VL.push_back(TmpVector);
    }
    SkipTag(in, "</vector_list>");
  }


  degree ReadDegree(std::istream& in, long GradingDim)
  {
    degree ans(GradingDim);
    BigInt tmp;
    for (long i=0; i<GradingDim; ++i)
    {
      in >> tmp;
      InputFailCheck(in, "ReadDegree");
      ans.mySetComponent(i, tmp);  // ans[i] is read only.
    }
    return ans;
  }


  FreeModule ReadFreeModule(std::istream& in, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<graded_free_module>");
    bool IsPosTO=false;
    vector<degree> InputShifts;
    vector<size_t> ComponentOrdering;

    const SparsePolyRing P(ReadPolyRing(in, GetTag));

    SkipTag(in, "number_module_components");
    long NumComponents;
    in >> NumComponents;
    InputFailCheck(in, "ReadFreeModule -- NumComponents");

    SkipTag(in, "<ordering_grading>");
    string tag;
    in >> tag;
    InputFailCheck(in, "ReadFreeModule -- tag");
    while (tag!="</ordering_grading>")
    {
      switch (tag[1])
      {
      case 'm':
        AssertTag(tag, "<module_shifts>");
        for (long i=0 ; i<NumComponents ; ++i)
          InputShifts.push_back(ReadDegree(in, GradingDim(P)));
        SkipTag(in, "</module_shifts>");
        break;
      case 'o':
        AssertTag(tag, "<ordering_type>");
        in >> tag;
        InputFailCheck(in, "ReadFreeModule -- tag2");
        while (tag!="</ordering_type>")
        {
          if      (tag=="wdeg_pos_to" || tag=="pos_to") IsPosTO = true;
          else if (tag=="wdeg_to_pos" || tag=="to_pos") IsPosTO = false;
          else if (tag=="<module_components_ordering>")
          {
            for (long i=0 ; i<NumComponents ; ++i)
            {
              long c;
              in >> c;
              ComponentOrdering.push_back(c);
            }
            InputFailCheck(in, "ReadFreeModule -- ComponentOrdering");
            SkipTag(in, "</module_components_ordering>");
          }
          else ThrowInputError(tag);
          in >> tag;
          InputFailCheck(in, "ReadFreeModule -- tag3");
        }
        break;
      }
      in >> tag;
      InputFailCheck(in, "ReadFreeModule -- tag4");
    }
    SkipTag(in, "</graded_free_module>");

    if (InputShifts.empty()) // GradimDim is 0
    {
      CoCoA_ASSERT(GradingDim(P)==0);
      for (long i=0; i != NumComponents; ++i) InputShifts.push_back(degree(0));
    }

    if (IsPosTO)
      return NewFreeModule(P, InputShifts, WDegPosnOrd);
    else
      return NewFreeModule(P, InputShifts, OrdPosn);
  }


  SparsePolyRing ReadPolyRing(std::istream& in, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<polynomial_ring>");
    long NumInds = 0;
    long NumPars = 0;
    BigInt charact(32003);
    long PrecBits = 0;

    string PolyRingName;
    string tag;
    in >> tag;
    InputFailCheck(in, "ReadPolyRing -- tag");
    if (tag=="name")
    {
      in >> PolyRingName;
      in >> tag;
      InputFailCheck(in, "ReadPolyRing -- PolyRingName");
    }
    AssertTag(tag, "<coefficient_ring>");
    in >> tag;
    InputFailCheck(in, "ReadPolyRing -- tag2");
    while (tag!="</coefficient_ring>")
    {
      switch (tag[0])
      {
      case 'f': AssertTag(tag,"float_precision"); in >> PrecBits;break;
      case 'c': AssertTag(tag, "characteristic"); in >> charact; break;
      default:  ThrowInputError(tag);
      }
      in >> tag;
      InputFailCheck(in, "ReadPolyRing -- tag3");
    }
    // Read Indets
    SkipTag(in, "<indeterminates>");
    in >> tag;
    InputFailCheck(in, "ReadPolyRing -- tag4");
    while (tag!="</indeterminates>")
    {
      switch (tag[7])
      {
      case 'i': AssertTag(tag, "number_indets"); in >> NumInds; break;
      case 'p': AssertTag(tag, "number_params"); in >> NumPars; break;
      default:  ThrowInputError(tag);
      }
      in >> tag;
      InputFailCheck(in, "ReadPolyRing -- tag5");
    }
    const PPOrdering O(ReadPolyRingOrdering(in, NumInds, GetTag));
    SkipTag(in, "</polynomial_ring>");

    SparsePolyRing P = NewPolyRingServer(charact, PrecBits, NumPars, O);
    WritePolyRing(GlobalLogput(), P);
    return P;
  }


  void WriteMatrix(std::ostream& out, const matrix &M)
  {
    if (NumCols(M)==0) { out << "Mat([[]])" << endl;  return; }
    ostringstream buf;
    buf << "Mat([";
    for (long i=0 ; i < NumRows(M); ++i)
    {
      if (i > 0) buf << ",";
      buf  << "\n  [" << M(i,0);
      for (long j=1 ; j < NumCols(M); ++j)  buf << ", " << M(i,j);
      buf << "]";
    }
    buf << "\n])" << flush;
    out << buf.str();
  }


//   void WriteMatrixInVar(std::ostream& out, const std::string& VarName, const matrix &M)
//   {
//     out << VarName << " := ";
//     WriteMatrix(out, M);
//     out << ";" << endl;
//   }


  void WritePolyRing(std::ostream& out, const SparsePolyRing& P)
  {
    const BigInt charact = characteristic(P);
    long NumPars = 0;
    long prec = 0;
    PPOrdering PPO = ordering(PPM(P));

    if (IsRingTwinFloat(CoeffRing(P)))
      prec = PrecisionBits(CoeffRing(P));
    if (IsQuotientRing(CoeffRing(P)) &&
        IsSparsePolyRing(BaseRing(CoeffRing(P))))
    { // this is just enough for the rings CoCoA-4 can send
      const SparsePolyRing R = BaseRing(CoeffRing(P));
      if (IsRingTwinFloat(CoeffRing(R)))
        prec = PrecisionBits(CoeffRing(R));
      NumPars = NumIndets(R);
    }
    out << "Use P::=";
    if      (charact > 0) out << "ZZ/(" << charact << ")";
    else if (prec == 0)   out << "QQ";
    else                  out << "RR(" << prec << ")";
    if (NumPars > 0)
      out << "(a[0.." << NumPars-1 << "])";
    out << "[x[0.." << NumIndets(P)-1 << "]]";
    if (IsLex(PPO))          out << ",Lex;\n";
    else if (IsStdDegLex(PPO))    out << ",DegLex;\n";
    else if (IsStdDegRevLex(PPO)) out << ";\n";
    else
    {
      out << ",Mat([..]);\n";
      out << "-- GrDim = " << GradingDim(PPO) << endl;
    }
    out << flush;
  }


//   void WriteMonomials(std::ostream& out, ConstRefRingElem f)
//   {
//     // prints the polynomial writing [term,..,term] (called by WritePolyList)
//     if (IsZero(f)) // trivial special case
//     { out << "[0]" << endl; return; }

//     ostringstream buf;
//     SparsePolyIter i=BeginIter(f);
//     buf << "[ (" << coeff(i) << ")*" << PP(i);
//     for (++i; !IsEnded(i); ++i)
//       buf << ", (" << coeff(i) << ")*" << PP(i);
//     buf << " ]" << flush;
//     out << buf.str() << flush;
//   }


  void WritePoly(std::ostream& out, ConstRefRingElem f)
  { out << f; }


  void WritePolyList(std::ostream& out, const PolyList& PL)
  {
    out << "<poly_list>" << endl;
    for (PolyList::const_iterator it=PL.begin(); it!=PL.end(); ++it)
    { 
      WritePoly(out, *it);
      out << ";\n";
    }
    out << "</poly_list>" << endl;
  }


  void WriteIdeal(std::ostream& out, const ideal& I)
  {
    out << "<ideal>" << endl;
    WritePolyList(out, gens(I));
    out << "</ideal>" << endl;
  }


  void WritePolyListInVar(std::ostream& out,
                          const std::string& VarName,
                          const PolyList& PL)
  {
    //    out << VarName << " := ";
    WritePolyList(out, PL);
    //    out << ";" << endl;
  }


  void WriteVectorListInVar(std::ostream& out,
                            const std::string& VarName,
                            const VectorList& VL)
  {
    out << VarName << " := [];" << endl;
    for (VectorList::const_iterator it=VL.begin(); it!=VL.end(); ++it)
      out << "Append(" << VarName << ", Vector(" << *it << "));" << endl;
  }


  void PrintTimeToLog(double T)
  {
    GlobalLogput() << "[log]              TotalTime="
                   << T << endl;
  }


  void PrintTimeToCoCoA5(double T)
  {
    GlobalOutput() << "--CoCoAServer: computing Cpu Time = "
                   << T << endl;
  }


  void PrintVersionToCoCoA5()
  {
    GlobalOutput() << "--CoCoAServer: CoCoALib version "
                   << BuildInfo::version() << endl;
  }

  void EndOfTransmissionToCoCoA5()
  {
    // The CoCoA4 parser interprets the character '\001' as meaning end of file;
    // we use this character to mark the end of a transmission.
    GlobalOutput() << char(1) << flush;
  }

  //----------------------------------------------------------------------

  PPOrdering ReadPolyRingOrdering(std::istream& in, long NumInds, SkipTagType ST)
  {
    if (ST == GetTag) SkipTag(in, "<ordering_grading>");
    string tag;
    SkipTag(in, "ordering_type");
    in >> tag;
    InputFailCheck(in, "ReadPolyRingOrdering -- tag");
    switch (tag[9])
    {
    case 'p': AssertTag(tag, "lexicographic");
      SkipTag(in, "</ordering_grading>");
      return lex(NumInds);
    case 'v': AssertTag(tag, "graded_reverse_lexicographic");
      SkipTag(in, "</ordering_grading>");
      return StdDegRevLex(NumInds);
    case 'x': AssertTag(tag, "graded_lexicographic");
      SkipTag(in, "</ordering_grading>");
      return StdDegLex(NumInds);
    case 'd': AssertTag(tag, "matrix_ordering");
      {
        long GrDim;
        SkipTag(in, "grading_dim");
        in >> GrDim;
        InputFailCheck(in, "ReadPolyRingOrdering -- GrDim");
        matrix M = ReadIntegerMatrix(in, GetTag);
        if (NumRows(M) != NumInds || NumCols(M) != NumInds)
          CoCoA_THROW_ERROR(ERR::BadMatrixSize, "ReadPolyRingOrdering");
        SkipTag(in, "</ordering_grading>");
        return NewMatrixOrdering(M, GrDim);
      }
    default:  ThrowInputError(tag);
    }
    return StdDegRevLex(NumInds); // to stop compilation warning
  }


  ring NewCoeffRing(const BigInt& charact, long FloatPrecision, long NumParams)
  {
    if (NumParams == 0)
    {
      if (IsZero(charact) && FloatPrecision==0)  return RingQQ();
      if (IsZero(charact))  return NewRingTwinFloat(FloatPrecision);
      return NewZZmod(charact);
    }

    //-------- Parameter Ring --------//
    const vector<symbol> ParamNames(SymbolRange("a", 0, NumParams-1));
    if (IsZero(charact) && FloatPrecision!=0)
      return NewPolyRing_DMPI(NewRingTwinFloat(FloatPrecision), NumParams);
    if (IsZero(charact))
      return NewFractionField(NewPolyRing_DMPI(RingQQ(), ParamNames));
    ring Zmod = NewZZmod(charact);
    if (IsRingFp(Zmod))
      return NewFractionField(NewPolyRing_DMPII(Zmod, ParamNames));
    else
      return NewFractionField(NewPolyRing_DMPI(Zmod, ParamNames));
  }


  SparsePolyRing NewPolyRingServer(const BigInt& charact, long FloatPrecision, long NumParams, const PPOrdering& O)
  {
    //-------- coefficient ring --------//
    ring R = NewCoeffRing(charact, FloatPrecision, NumParams);
    const vector<symbol> IndetNames(SymbolRange("x", 0, NumIndets(O)-1));

    //-------- PolyRing --------//
    if (IsRingFp(R))
      return NewPolyRing_DMPII(R, IndetNames, O);
    /*else*/
    return NewPolyRing_DMPI(R, IndetNames, O);
  }


  //----------------------------------------------------------------------

  void ThrowInputError(const std::string& unknown_tag)
  {
    CoCoA_THROW_ERROR("UNKNOWN TAG: " + unknown_tag, "ThrowInputError");
  }


  void AssertTag(const std::string& input_tag, const std::string& expected_tag)
  {
    if (input_tag == expected_tag) return;
    CoCoA_THROW_ERROR("Expected tag: " + expected_tag +
                "; Read tag: " + input_tag,
                "CoCoAServer");
  }


  void SkipTag(std::istream& in, const std::string& expected_tag)
  {
    string tag;
    in >> tag;
    InputFailCheck(in, "SkipTag "+expected_tag);
    AssertTag(tag, expected_tag);
  }

}
