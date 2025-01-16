//   Copyright (c)  2005,2007,2008,2012  John Abbott and Anna M. Bigatti

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

#include "CoCoA/symbol.H"

#include "CoCoA/utils.H" // must come before VectorOps!!
#include "CoCoA/VectorOps.H"
#include "CoCoA/error.H"

#include <algorithm>
using std::min;
using std::sort;
#include <atomic>
#include <cstddef>
using std::size_t;
#include <ios>
#include <iostream>
using std::ostream;
using std::istream;
#include <limits>
using std::numeric_limits;
#include <sstream>
using std::istringstream; // for symbols(string)
//#include <string>
using std::string;
//#include <vector>
using std::vector;


namespace CoCoA
{

  namespace // anonymous, for file local defns
  {
    std::atomic<long> NewSymbolCounter;

    // Printable head for anonymous symbols
    static const string AnonHead = "#";

    // symbol::myInput ASSUMES this fn accepts a (weak) superset of IsValidFirstChar
    inline bool IsValidAfterFirstChar(char ch) noexcept
    {
      return (ch == '_') || isalnum(ch);
    }

  } // end of anonymous namespace 


  symbol NewSymbol()
  {
    return symbol(/*anonymous,*/ NewSymbolCounter++); // post-incr so first symbol has index 0
  }


  // Not fully exception-safe: loop could throw, but NewSymbolCounter was updated
  std::vector<symbol> NewSymbols(const long NumSyms)
  {
    if (NumSyms <= 0)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "NewSymbols");
    vector<symbol> ans; ans.reserve(NumSyms);
    // Next 2 lines should be threadsafe (I hope)
    const long LAST = (/*atomic*/NewSymbolCounter += NumSyms);  // BUG no overflow check!
    const long FIRST = LAST - NumSyms;
    for (long j = FIRST; j < LAST; ++j)
      ans.push_back(symbol(/*anonymous,*/ j));
    return ans;
  }


  // head must be non empty & begin with a letter and be comprised of letters or underscore.
  bool symbol::IsValidHead(const std::string& head) noexcept
  {
    if (head.empty() || !IsValidFirstChar(head[0])) return false;
    const long n = len(head);
    for (long i=1; i < n; ++i)
    {
      if (!IsValidAfterFirstChar(head[i]))
        return false;
    }
    return true;
  }


  symbol::symbol(long subscript):
      myHead(/*anonymous*/),
      mySubscripts(1)
  {
    mySubscripts[0] = subscript;
  }



  symbol::symbol(const std::string& head):
      myHead(head),
      mySubscripts(0)
  {
    if (!IsValidHead(myHead))
      CoCoA_THROW_ERROR(ERR::BadSymbolHead, "symbol(head)");
  }


  symbol::symbol(const std::string& head, long subscript):
      myHead(head),
      mySubscripts(1)
  {
    if (!IsValidHead(myHead))
      CoCoA_THROW_ERROR(ERR::BadSymbolHead, "symbol(head, subscript)");
    mySubscripts[0] = subscript;
  }



  symbol::symbol(const std::string& head, long subscript1, long subscript2):
      myHead(head),
      mySubscripts(2)
  {
    if (!IsValidHead(myHead))
      CoCoA_THROW_ERROR(ERR::BadSymbolHead, "symbol(head, subscript)");
    mySubscripts[0] = subscript1;
    mySubscripts[1] = subscript2;
  }



  symbol::symbol(const std::string& head, const std::vector<long>& subscripts):
      myHead(head),
      mySubscripts(subscripts)
  {
    if (!IsValidHead(myHead))
      CoCoA_THROW_ERROR(ERR::BadSymbolHead, "symbol(head, subscripts)");
  }


  void symbol::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out));

    if (myHead.empty()) out << AnonHead; else out << myHead;
    if (mySubscripts.empty()) return;
    out << '[' << mySubscripts[0];
    const long n = len(mySubscripts);
    for (long i=1; i < n; ++i)
    {
      out << "," << mySubscripts[i];
    }
    out << ']';
  }

//    void myOutputSelf(OpenMath::ostream& out) const;????


  bool symbol::IsValidFirstChar(char ch) noexcept
  {
    return isalpha(ch);
  }



  // This is semi-exception clean: value of object is modified only
  // if the entire read is successful.
  // This version does not handle white space as I'd like -- see symbol.txt.
  // Normal input fns merely set the failbit of the stream for parse errors.
  void symbol::myInput(istream& in)
  {
    static const char* const FnName = "symbol::myInput";
    if (!in.good()) CoCoA_THROW_ERROR("istream is not good", FnName);
    if (!IsDecimal(in)) CoCoA_THROW_ERROR("istream is not in \"decimal\" mode", FnName);
    in >> std::ws;
    if (!in.good() || !IsValidFirstChar(in.peek())) { CoCoA_THROW_ERROR("No symbol head found", FnName); }
    string head;
    while (IsValidAfterFirstChar(in.peek())) // ASSUMES that IsValidFirstChar ==> IsValidAfterFirstChar
    {
      char ch;
      in >> ch;
      head.push_back(ch);
    }
    if (in.eof() || in.peek() != '[')
    {
      // symbol without indexes/subscripts
      in.clear();
      swap(myHead, head);
      mySubscripts.clear();
      return;
    }

    in.ignore(); // ignore the initial open square bracket
    vector<long> subscripts;
    while (true)
    {
      long subscript;
      in >> subscript >> std::ws;
      if (!in) { CoCoA_THROW_ERROR("Invalid symbol index", FnName); }
      subscripts.push_back(subscript);
      if (in.peek() == ']') break;
      if (!in || in.peek() != ',') { CoCoA_THROW_ERROR("Symbol index list not closed, or contains unexpected char", FnName); } 
      in.ignore(); // ignore the comma
    }
    in.ignore(); // ignore the final ']'

    swap(myHead, head);
    swap(mySubscripts, subscripts);
  }


  int symbol::myCmp(const symbol& sym2) const noexcept
  {
    const int HeadCmp = LexCmp3(myHead.begin(), myHead.end(),
                                sym2.myHead.begin(), sym2.myHead.end());
    if (HeadCmp != 0) return HeadCmp;
    return LexCmp3(mySubscripts.begin(), mySubscripts.end(),
                   sym2.mySubscripts.begin(), sym2.mySubscripts.end());
  }


  namespace // anon namespace for file local fn
  {
    void FillRange(vector<symbol>& ans,
                   const string& head,
                   const vector<long>& subscripts1,
                   const vector<long>& subscripts2,
                   vector<long>& workspace,
                   long i)
    {
//      CoCoA_ASSERT(symbol::IsValidHead(head));
//      CoCoA_ASSERT(len(subscripts1) == len(subscripts2));
      if (i == len(workspace))
      {
        ans.push_back(symbol(head, workspace));
        return;
      }

//      CoCoA_ASSERT(i < len(subscripts1));
//      CoCoA_ASSERT(subscripts1[i] <= subscripts2[i]);
      for (long j=subscripts1[i]; /*j <= subscripts2[i]*/; ++j)
      {
        workspace[i] = j;
        FillRange(ans, head, subscripts1, subscripts2, workspace, i+1);
        if (j == subscripts2[i]) break; // do it this way in case subscripts2[i] == MaxLong
      }
    }
  }


  // This does a sort of cartesian product; not tricky, just long and tedious.
  // !!! ASSUMES that numeric_limits<size_t>::max() >= numeric_limits<long>::max() !!!
  // JAA: I don't see how to get rid of size_t from this fn, because max_size() produces size_t.
  std::vector<symbol> SymbolRange(const symbol& sym1, const symbol& sym2)
  {
    const long n = NumSubscripts(sym1);
    if (head(sym1) != head(sym2) || NumSubscripts(sym2) != n)
      CoCoA_THROW_ERROR(ERR::BadSymbolRange, "SymbolRange(sym1,sym2)");

    // Boring trivial case; maybe it should be an error?
    if (sym1 == sym2)
      return vector<symbol>(1, sym1);

    const long MaxLong = numeric_limits<long>::max();
    vector<symbol> ans;
    const size_t MAX = min(size_t(MaxLong), ans.max_size());
    // Conduct sanity check on subscript ranges, and compute
    // number of elements in final range; complain if too big.
    size_t RangeSize = 1;
    for (long i=0; i < n; ++i)
    {
      if (subscript(sym1,i) > subscript(sym2,i))
        CoCoA_THROW_ERROR(ERR::BadSymbolRange, "SymbolRange(sym1,sym2)");
      const size_t dim = ULongDiff(subscript(sym2,i), subscript(sym1,i));
      // Next "if" checks whether dim+1 would overflow:
      if (dim >= static_cast<size_t>(MaxLong))
        CoCoA_THROW_ERROR(ERR::ArgTooBig, "SymbolRange(sym1,sym2)");
      // Next "if" checks whether product Range*(dim+1) would overflow.
      if (MAX/(dim+1) < RangeSize)
        CoCoA_THROW_ERROR(ERR::ArgTooBig, "SymbolRange(sym1,sym2)");
      RangeSize *= dim+1; // cannot overflow thanks to checks above.
    }
    ans.reserve(RangeSize);
    vector<long> subscripts1(n);
    vector<long> subscripts2(n);
    vector<long> workspace(n);
    for (long i=0; i < n; ++i)
    {
      subscripts1[i] = subscript(sym1,i);
      subscripts2[i] = subscript(sym2,i);
    }
    FillRange(ans, head(sym1), subscripts1, subscripts2, workspace, 0);
    return ans;
  }

  std::vector<symbol> SymbolRange(const std::string& head, long lo, long hi)
  {
    const char* const FnName = "SymbolRange(hd,ind1,ind2)";
    if (!symbol::IsValidHead(head)) CoCoA_THROW_ERROR(ERR::BadSymbolHead, FnName);
    return SymbolRange(symbol(head,lo), symbol(head,hi));
  }


  long subscript(const symbol& sym, long n)
  {
    if (n < 0 || n >= NumSubscripts(sym))
      CoCoA_THROW_ERROR(ERR::BadSymbolSubscript, "subscript(symbol,n)");
    return sym.mySubscripts[n];
  }


  bool AreDistinct(const std::vector<symbol>& syms) noexcept
  {
    vector<symbol> s = syms;
    sort(s.begin(), s.end());
    const long NumSyms = len(s);
    for (long i=0; i < NumSyms-1; ++i)
      if (s[i] == s[i+1]) return false;
    return true;
  }


  bool AreArityConsistent(const std::vector<symbol>& syms) noexcept
  {
    vector<symbol> s = syms;
    sort(s.begin(), s.end());
    const long NumSyms = len(s);
    for (long i=0; i < NumSyms-1; ++i)
      if (head(s[i]) == head(s[i+1]) && NumSubscripts(s[i]) != NumSubscripts(s[i+1]))
        return false;
    return true;
  }


  std::ostream& operator<<(std::ostream& out, const symbol& sym)
  {
    if (!out) return out;  // short-cut for bad ostreams
    sym.myOutputSelf(out);
    return out;
  }


  std::istream& operator>>(std::istream& in, symbol& sym)
  {
    sym.myInput(in);
    return in;
  }


  namespace // anonymous, for file local fns related to symbols(...)
  {

    void CoCoA_THROW_ERROR_NextChar(std::istream& in, const string& FuncName)
    {
      in.clear();
      CoCoA_THROW_ERROR("Unexpected \'"+ string(1,char(in.peek())) +"\'", FuncName);
    }


    // void ReadSymbolRange(std::istream& in, std::vector<symbol>& symbs)
    // {
    //   CoCoA_THROW_ERROR(ERR::NYI, "ReadSymbolRange");
    // }
    
  } // end of namespace anonymous
  

  std::vector<symbol> symbols(const std::string& str)
  {
    static const char* const FnName = "symbols";
    vector<symbol> symbs;
    if (str.empty()) return symbs;
    istringstream in(str);
    char ch;
    while (true)
    {
      in >> std::ws; ch = in.peek();
      if (!in.good() || !symbol::IsValidFirstChar(ch))
        CoCoA_THROW_ERROR_NextChar(in, FnName);
      symbol s("dummy");
      in >> s; // throws if there is a problem
      symbs.push_back(s);

      in >> std::ws; ch = in.peek();
      if (in.eof()) break;
      if (ch != ',') CoCoA_THROW_ERROR("Expected comma in list of symbols", FnName);
      in.ignore(); // skip the comma
    }
    return symbs;
  }


  std::vector<symbol> symbols(const std::vector<std::string>& s)
  {
    vector<symbol> IndetNames;
    for (long i=0; i < len(s); ++i)
      IndetNames.push_back(symbol(s[i]));
    return IndetNames;
  }


} // end of namespace CoCoA
