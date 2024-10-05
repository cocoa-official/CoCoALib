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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/symbol.C,v 1.50 2024/07/02 15:37:50 bigatti Exp $
// $Log: symbol.C,v $
// Revision 1.50  2024/07/02 15:37:50  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.49  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.48  2021/02/17 17:50:36  abbott
// Summary: Now asserts ostream is in decimal mode (redmine 1547)
//
// Revision 1.47  2021/02/10 19:42:56  abbott
// Summary: Added noexcept (redmine 1572); change index type to long from MachineInte (redmine 925)
//
// Revision 1.46  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.45  2020/12/04 21:38:28  abbott
// Summary: Major revision to symbol::myInput and to symbols (redmine 1529, 1523)
//
// Revision 1.44  2020/06/17 19:01:09  abbott
// Summary: Renamed two functions to better align with CoCoA_THROW_ERROR
//
// Revision 1.43  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.42  2020/02/11 16:56:43  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.41  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.40  2019/03/27 13:51:07  bigatti
// (abbott) removed ThreadsafeCounter --> using atomic
//
// Revision 1.39  2018/06/13 13:11:42  abbott
// Summary: Symbols now accepts empty string, and returns empty list
//
// Revision 1.38  2018/06/13 12:34:15  abbott
// Summary: Added arg check to symbols (empty string)
//
// Revision 1.37  2018/05/18 12:25:54  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.36  2018/05/17 15:56:32  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.35  2018/04/18 14:29:04  abbott
// Summary: Commented out unused fn "ReadSymbolRange"
//
// Revision 1.34  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.33  2016/02/01 12:42:55  abbott
// Summary: Removed cruft - the old defns of symbols.
//
// Revision 1.32  2015/07/28 11:55:20  bigatti
// -- symbols(string) a bit more robust
//
// Revision 1.31  2015/07/27 11:47:50  bigatti
// -- function "symbols(string)" now reads comma separated symbols
// -- function "symbols(string,string,..)" now disabled
//
// Revision 1.30  2014/07/31 14:45:19  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.29  2014/07/14 15:10:41  abbott
// Summary: Added include of UtilsTemplate.H (after moving LexCmp3 there)
// Author: JAA
//
// Revision 1.28  2014/06/14 19:27:56  abbott
// Summary: Modified to allow digits in symbol heads
// Author: JAA
//
// Revision 1.27  2014/04/03 15:33:25  abbott
// Summary: Now uses ULongDiff
// Author: JAA
//
// Revision 1.26  2014/03/21 17:12:26  abbott
// Summary: Corrected so it now reads negative indices too
// Author: JAA
//
// Revision 1.25  2014/03/21 15:45:27  abbott
// Summary: Corrected handling of EOF right after symbol head
// Author: JAA
//
// Revision 1.24  2014/03/21 15:15:55  abbott
// Summary: symbol::myInput now sets failbit if the head cannot be read
// Author: JAA
//
// Revision 1.23  2014/01/29 11:22:20  abbott
// Improved input for symbols.
//
// Revision 1.22  2012/07/04 12:26:50  abbott
// Removed obsolete include directives (for BOOST stuff).
//
// Revision 1.21  2012/05/29 14:56:15  abbott
// Separated ThreadsafeCounter from symbol; also employed it in ring.C.
//
// Revision 1.20  2012/05/29 07:48:13  abbott
// Simplified private ctor for anonymous symbols.
// Some general cleaning.
//
// Revision 1.19  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.18  2012/05/24 14:46:20  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.17  2012/05/10 14:43:26  abbott
// Added new fns NewSymbol & NewSymbols.
// Also added ThreadsafeCounter - will be separated off later.
//
// Revision 1.16  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.15  2011/09/06 15:22:12  abbott
// Added "const" to a local variable.
//
// Revision 1.14  2011/03/14 10:30:08  abbott
// Changed size_t into long (in fn interfaces).
//
// Revision 1.13  2011/02/16 14:58:10  bigatti
// -- added    symbols(vector<string> s)
//
// Revision 1.12  2010/11/26 15:23:55  bigatti
// -- removed space in subscripts: c[1, 2] is now printed c[1,2]
//
// Revision 1.11  2008/12/17 10:49:22  abbott
// Change loop control variable from int into long (as it should always have been).
//
// Revision 1.10  2008/12/12 11:29:47  abbott
// Fixed a bug in SymbolRange.  Added example and test for symbols.
//
// Revision 1.9  2008/12/11 20:12:32  abbott
// Now using MachineInt type in ctors to specify the values of symbol subscripts
// (but we check that the values fit into a long).  Corrected two bugs in
// SymbolRange (an undetected overflow and an off-by-1 error).
//
// Revision 1.8  2008/12/11 10:54:01  abbott
// Now subscripts are [long] values rather than [int].  The extra range will probably
// never be useful, but I prefer to offer it all the same.
//
// Revision 1.7  2008/11/21 21:16:34  abbott
// Added ctor for symbols with 2 subscripts.
//
// Revision 1.6  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.5  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/09/25 16:32:30  abbott
// Several minor changes to silence gcc-4.3:
//    more #includes,
//    and fixed a template problemm in RegisterServerOps.C
//
// Revision 1.3  2007/05/31 14:50:09  bigatti
// -- added AreDistinct and AreArityConsistent
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.5  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.4  2006/10/27 19:09:45  cocoa
// Replaced some member functions of CoCoA::symbol by friend functions.
// Removed some include dependency on symbol.H
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
