//   Copyright (c)  2014-2017 John Abbott and Anna M. Bigatti
//   Author:  2014-2017 John Abbott and Anna M. Bigatti

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


// Source code for input of polynomials and polynomial expressions

#include "CoCoA/RingElemInput.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/geobucket.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"
#include "CoCoA/time.H"

#include <algorithm>
using std::sort;
#include <istream>
using std::istream;
#include <iostream> // just for debugging
#include <map>
using std::map;
#include <sstream>
using std::istringstream;
using std::ostringstream;
#include <string>
using std::string;
#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    // returns 'a' for alpha, '0' for digit, '\0' for EOF
    // std expr chars are returned as themselves +-*/^()
    // other chars are returned as '?'
    char WhatsNext(std::istream& in)
    { // recognized characters
      CoCoA_ASSERT(!in.bad());
      in >> std::ws; // skip white spaces
      //    std::cout << " WhatsNext: testing " << char(in.peek()) << std::endl;
      if (in.eof())  return '\0';
      const char ch = in.peek();
      if (isdigit(ch))  return '0'; // 0 for digit
      if (isalpha(ch))  return 'a'; // a for alpha
      switch (ch)
      {
      case '+':
      case '-':
      case '*':
      case '/':
      case '^':
      case '(':
      case ')':
        return ch;
      case ';':
      case ',':
        return ch;
      }
      // if (ch == '+') return '+';
      // if (ch == '-') return '-';
      // if (ch == '*') return '*';
      // if (ch == '/') return '/';
      // if (ch == '^') return '^';
      // if (ch == '(') return '(';
      // if (ch == ')') return ')';
      // if (ch == ';') return ';';
      // if (ch == ',') return ',';
   //   CoCoA_THROW_ERROR("Illegal char \'"+string(1,ch)+'\'', "WhatsNext");
      return '?'; // never reached -- just to keep compiler quiet
    }


    void CoCoA_THROW_ERROR_NextChar(std::istream& in, const string& FuncName)
    {
      if (in.eof())
        CoCoA_THROW_ERROR("Unexpected EOF", FuncName);
//???    in.clear();
      const char ch = in.peek();
      const int ByteCode = static_cast<int>(static_cast<unsigned char>(ch));
      // 
      if (ByteCode >= 32 && ByteCode < 127)
        CoCoA_THROW_ERROR("Unexpected \'"+string(1,ch)+'\'', FuncName);
      std::ostringstream ByteCodeDecimal; ByteCodeDecimal << ByteCode;
      CoCoA_THROW_ERROR("Unexpected char with byte code "+ByteCodeDecimal.str(), FuncName);
    }


    RingElem ReadExpr(const ring& P, std::istream& in, const map<symbol, long>& SymTable); // fwd decl

    RingElem ReadAtom(const ring& P, std::istream& in, const map<symbol, long>& SymTable)
    {
      CoCoA_ASSERT(in.good());
      CoCoA_ASSERT(in.flags() & std::ios::dec); // fails if std:ios:oct or std::ios::hex
////      const int base = (in.flags() & std::ios::oct) ? 8 : (in.flags() & std::ios::hex) ? 16 : 10;
////      RingElem tmp(P);

      const char ch = WhatsNext(in);
      if (ch == '0')
      {
        BigInt N;
        in >> N;
        // No need to check status of in, since there was at least 1 digit
        if (in.peek() != '.')  { return RingElem(P,N); }
        in.ignore(); // skip over "."
        // Next lines copied from BigRat.C
        const string AfterDot = ScanUnsignedIntegerLiteral(in);
        const long NumPlaces = len(AfterDot);
        if (NumPlaces == 0)  { return RingElem(P,N); }
        istringstream FracDigits(AfterDot);
///        if (base == 8) FracDigits >> std::oct;
///        if (base == 16) FracDigits >> std::hex;
        BigRat FracPart;
        FracDigits >> FracPart;
        FracPart /= power(10, NumPlaces); // ASSUMES DECIMAL!!!
        return  RingElem(P, N+FracPart);
      }

      if (ch == 'a') //???? (symbol::IsValidFirstChar(ch))
      {
        symbol s("dummy");
        in >> s;  if (!in)  CoCoA_THROW_ERROR_NextChar(in, "ReadFactor: malformed symbol"); // can happen if input is "a[1"
        if (SymTable.find(s) == SymTable.end())
        {
          std::ostringstream os;
          os << "symbol \"" << s << "\" not in ring";
          CoCoA_THROW_ERROR(std::string(os.str()), "ReadFactor: symbol");
        }
        //        return indet(P, SymTable.find(s)->second);
        return RingElem(P, s);
      }
      if (ch != '(')  CoCoA_THROW_ERROR_NextChar(in, "ReadAtom");

//      if (ch == '(')
//      {
///      case '(':
        in.ignore(); // '('
        RingElem SubExpr = ReadExpr(P, in, SymTable); // looks wrong, why not readexpr(blah,symboltable)??????????????????????//??????????
        if ( in.peek() != ')' )
          CoCoA_THROW_ERROR_NextChar(in, "ReadFactor: expected close bracket for (...)");
        in.ignore(); // ')'
        return SubExpr;

/////      default:
    }


    RingElem ReadFactor(const ring& P, std::istream& in, const map<symbol, long>& SymTable)
    {
      CoCoA_ASSERT(in.good());
      CoCoA_ASSERT(in.flags() & std::ios::dec); // fails if std:ios:oct or std::ios::hex
      const RingElem base = ReadAtom(P, in, SymTable);

      if ( WhatsNext(in) != '^' )  return base;
      in.ignore(); // '^'
      // Just read a '^', so expect either a non-negative integer literal or '(' integer-literal ')'
      in >> std::ws;
      const char AfterUpArrow = in.peek();
      if (!in.good())
        CoCoA_THROW_ERROR_NextChar(in, "ReadFactor: missing exponent");
      if (AfterUpArrow != '(' && !isdigit(AfterUpArrow)) // ASSUMES DECIMAL!!!!
        CoCoA_THROW_ERROR_NextChar(in, "ReadFactor: exponent must be non-neg integer or in brackets");
      bool IsNegative = false;
      if (AfterUpArrow == '(')
      {
        in.ignore(); // skip over the '('
        in >> std::ws;
        const char SignOrDigit = in.peek();
        if (!in.good() || (SignOrDigit != '+' && SignOrDigit != '-' && !isdigit(SignOrDigit)))
          CoCoA_THROW_ERROR_NextChar(in, "ReadFactor: malformed exponent");
        if (!isdigit(SignOrDigit))  { IsNegative = (SignOrDigit == '-'); in.ignore(); } // skip '+' or '-'
      }
      BigInt exp;
      in >> exp;
      if (!in.good())
        CoCoA_THROW_ERROR_NextChar(in, "ReadFactor: malformed exponent");
      if (IsNegative) negate(exp);
      if (AfterUpArrow == '(')
      {
        // Now check that there was the matching ')'
        // We move past the next char only if it is ')', so that the error message makes sense
        in >> std::ws;
        const char CloseBracket = in.peek();
        if (in.eof() || CloseBracket != ')')
          CoCoA_THROW_ERROR_NextChar(in, "ReadFactor: missing close bracket in exponent");
        in.ignore();
      }
      return power(base, exp);
    }


    RingElem ReadProduct(const ring& P, std::istream& in, const map<symbol, long>& SymTable)
    {
      CoCoA_ASSERT(in.good());
      CoCoA_ASSERT(in.flags() & std::ios::dec); // fails if std:ios:oct or std::ios::hex
      char w = '*';
      // char PrevOp = ' ';
      RingElem resProd(one(P));

      while (true)
      {
        if ( w == '*' )  resProd *= ReadFactor(P, in, SymTable);
        else             resProd /= ReadFactor(P, in, SymTable); // ( w == '/' )
        w = WhatsNext(in);
        if ( w != '*' &&  w != '/' )  break;
        in.ignore();
        // uncomment the following to prevent a/b*c, a/b/c
        // if ( PrevOp == '/' ) 
        //   CoCoA_THROW_ERROR_NextChar(in, "ReadProduct: ambiguous after \'/\' - evaluating left-to-right");
// PrevOp = w;
      }
      return resProd;
    }


    //-------- input from istream --------------------

    bool CmpLPP(const RingElem& f1, const RingElem& f2)
    {
      return LPP(f1) < LPP(f2);
    }

    // Faster version if input is a large sum of terms -- uses geobuckets.
    RingElem ReadExprInSparsePolyRing(const ring& P, std::istream& in, const std::map<symbol, long>& SymTable)
    {
      CoCoA_ASSERT(in.good());
      char w = WhatsNext(in);
      if (w == '\0')
        CoCoA_THROW_ERROR_NextChar(in, "ReadExprInSparsePolyRing");
      geobucket gbk(P);
      // Put single monomials in a vector; later I'll sort them and sum them.
      // This works well if the order of the monomials is different from that of P.
      // JAA 2016-05-18 Still not happy with this solution; a std::list may be better;
      // perhaps use MoveLM instead of += in loop when summing vector entries???
      vector<RingElem> TermList;
      
      while (true)
      {
        if ( w == '-' || w == '+' )  in.ignore(); // '+/-'
        RingElem summand = ReadProduct(P,in, SymTable);
        if ( w == '-' )
          summand = -summand; // inefficient?
        if (IsMonomial(summand))
          TermList.push_back(summand); // SLUG, makes wasteful copy; use MoveLM???
        else
          gbk.myAddClear(summand, NumTerms(summand));
        w = WhatsNext(in);
/////        if ( (w == ',') || (w == ';') || (w == ')') || (w == '\0') )  break;
        if ( (w != '-') && (w != '+') )  break; /////CoCoA_THROW_ERROR_NextChar(in, "ReadExprInSparsePolyRing");
      }
      RingElem ans(P);
      if (!TermList.empty())
      {
        sort(TermList.begin(), TermList.end(), CmpLPP);
        const long nterms = len(TermList);
        for (long i=0; i < nterms; ++i)
          ans += TermList[i];
        if (IsZero(gbk))
          return ans;
        gbk.myAddClear(ans, NumTerms(ans));
      }

      // Get answer out of geobucket -- awkward syntax :-(
      AddClear(ans, gbk);
      return ans;
    }


    RingElem ReadExpr(const ring& P, std::istream& in, const map<symbol, long>& SymTable)
    {
      CoCoA_ASSERT(in.good());
      if (IsSparsePolyRing(P))
        return ReadExprInSparsePolyRing(P, in, SymTable); // usefully faster if input is large sum of terms
      char w = WhatsNext(in);
      if (w == '\0')
        CoCoA_THROW_ERROR_NextChar(in, "ReadExpr");
      RingElem f(P);

      while (true)
      {
        if ( w == '-' || w == '+' )  in.ignore(); // '+/-'
        if ( w == '-' )
          f -= ReadProduct(P,in, SymTable);
        else  f += ReadProduct(P,in, SymTable);
        w = WhatsNext(in);
////        if ( (w == ',') || (w == ';') || (w == ')') || (w == '\0') )  break;
        if ( (w != '-') && (w != '+') )  break; ///???CoCoA_THROW_ERROR_NextChar(in, "ReadExpr");
      }
      return f;
    }


  } // namespace anonymous -------------------------------------------------


  RingElem ReadExpr(const ring& P, std::istream& in)
  {
    if (!in.good())
      CoCoA_THROW_ERROR("istream is not good", "ReadExpr");
    if (!IsDecimal(in))
      CoCoA_THROW_ERROR("istream is not in \"decimal\" mode", "ReadExpr");

    const vector<symbol> syms = symbols(P);    // makes a copy of vec
    map<symbol, long> SymTable;
    for (int i=0; i < len(syms); ++i)
      SymTable[syms[i]] = i;
    
    const RingElem ans = ReadExpr(P, in, SymTable);
//    if (WhatsNext(in) != '\0') CoCoA_THROW_ERROR("Extra chars after reading ringelem expression", "ReadExpr");
    return ans;
  }


  RingElem ReadExprSemicolon(const ring& P, std::istream& in)
  {
    if (!in.good())
      CoCoA_THROW_ERROR("istream is not good", "ReadExprSemicolon");
    if (!IsDecimal(in))
      CoCoA_THROW_ERROR("istream is not in \"decimal\" mode", "ReadExprSemicolon");

    const RingElem ans = ReadExpr(P, in);
    if ( WhatsNext(in) != ';' )
      CoCoA_THROW_ERROR_NextChar(in, "ReadExprSemicolon");
    in.ignore(); // ';'
    return ans;
  }


  namespace // anonymous
  {
    
    std::vector<RingElem> ReadCommaSeparatedListOfRingElems(const ring& P, std::istream& in)
    {
      CoCoA_ASSERT(in.good());
      // Skip initial whitespace; if we hit EOF, return empty vector
      in >> std::ws;
      std::vector<RingElem> v;
      if (in.eof()) return v;
    
      const vector<symbol> syms = symbols(P);    // makes a copy of vec
      map<symbol, long> SymTable;
      for (int i=0; i < len(syms); ++i)
        SymTable[syms[i]] = i;

      while (true)
      {
        v.push_back(ReadExpr(P, in, SymTable));
        in >> std::ws;
        const char ch = in.peek();
        if (in.eof() || ch != ',')  return v;
        in.ignore(); // ignore the comma
      }
    }

  } // end of namespace anonymous


  std::vector<RingElem> RingElems(const ring& P, std::istream& in)
  {
    if (!in.good())
      CoCoA_THROW_ERROR("istream is not good",  "RingElems");
    if (!IsDecimal(in))
      CoCoA_THROW_ERROR("istream is not in \"decimal\" mode", "RingElems");

    /*const*/ vector<RingElem> v = ReadCommaSeparatedListOfRingElems(P, in);
    return v;
  }
  

  std::vector<RingElem> RingElemVec(const ring& P, std::istream& in)
  {
    static const char* const FnName = "RingElemVec";
    if (!in.good())
      CoCoA_THROW_ERROR("istream is not good", FnName);
    if (!IsDecimal(in))
      CoCoA_THROW_ERROR("istream is not in \"decimal\" mode", "RingElemVec");

    char ch;
    in >> ch;
    if (in.eof() || ch != '[')
      CoCoA_THROW_ERROR("Expected initial '['", FnName);
    in >> std::ws;
    if (in.eof())  CoCoA_THROW_ERROR_NextChar(in, FnName);
    if (in.peek() == ']')  { in.ignore(); return vector<RingElem>(0); }
    // now we expect the list to be non-empty (or else junk input)
    /*const*/ vector<RingElem> v = ReadCommaSeparatedListOfRingElems(P, in);
    if (in.eof())
      CoCoA_THROW_ERROR_NextChar(in, FnName); // missing close ]
    ch = in.peek();
    if (ch != ']')
      CoCoA_THROW_ERROR_NextChar(in, FnName); // expected close ], but got sthg else
    in.ignore();
    return v;
  }


  //-------- input from string --------------------

  RingElem ReadExprSemicolon(const ring& P, const std::string& s)
  {
    istringstream is(s);
    /*const*/ RingElem ans =  ReadExprSemicolon(P, is);
    if (!is.eof())
      CoCoA_THROW_ERROR("Extra chars after RingElem expression", "ReadExprSemicolon(string)");
    return ans;
  }

  namespace {    // anonymous

    std::string snippet(std::istringstream& is)
    {
      std::string s;
      if (is.eof()) return s; // should never happen
      s = "..";
      char c;
      is >> std::noskipws; // no skip white spaces
      for (long l=0; l<10; ++l)
      {
        is >> c;
        if (is.eof()) return s;
        s.append(1,c);
      }
      return s+"..";
    }
    
  }
  
  
  RingElem ReadExpr(const ring& P, const std::string& s)
  {
    istringstream is(s);
    /*const*/ RingElem ans = ReadExpr(P, is);
    if (!is.eof())
    {
      // There are extra chars after the end of valid expr: try to give helpful err mesg.
      const char ch = is.peek();
      // Special handling for "double exponents" (redmine 1579)
      if (ch == '^')
        CoCoA_THROW_ERROR("Double exponent not allowed; use brackets like (a^b)^c", "ReadExpr(string)");
      std::ostringstream os;
      os << "Unexpected char \'" << ch << "\' in string from \"" << snippet(is) << "\"";
      if (ch == '(' || isalnum(ch)) os << " -- maybe forgotten '*'?";
      CoCoA_THROW_ERROR(os.str(), "ReadExpr(string)");
    }
    return ans;
  }


  std::vector<RingElem> RingElems(const ring& P, const std::string& s)
  {
    istringstream is(s);
    /*const*/ vector<RingElem> v = RingElems(P, is);
    if (!is.eof())
      CoCoA_THROW_ERROR("Extra chars after last RingElem expression", "RingElems(string)");
////    if (WhatsNext(is) != '\0') CoCoA_THROW_ERROR_NextChar(is, "RingElems");
    return v;
  }
  
  std::vector<RingElem> RingElemVec(const ring& P, const std::string& s)
  {
    istringstream is(s);
    /*const*/ vector<RingElem> v = RingElemVec(P, is);
    if (!is.eof())  is >> std::ws; // flush any trailing whitespace
    if (!is.eof())
      CoCoA_THROW_ERROR("Extra chars after list inside []", "RingElemVec");
////    if (WhatsNext(is) != '\0') CoCoA_THROW_ERROR_NextChar(is, "RingElemVec");
    return v;
  }


} // end of namespace CoCoA
