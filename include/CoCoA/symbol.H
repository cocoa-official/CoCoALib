#ifndef CoCoA_symbol_H
#define CoCoA_symbol_H

//   Copyright (c)  2005,2008,2012  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/utils.H"


#include <vector>
//using std::vector;
#include <string>
//using std::string;
#include <iosfwd>
//using std::ostream;


namespace CoCoA
{

///  class MachineInt; // fwd decl -- defined in MachineInt.H

  class symbol
  {
  private:
    explicit symbol(/*anonymous*/ long subscript); // called only by NewSymbol & NewSymbols
    friend symbol NewSymbol();
    friend std::vector<symbol> NewSymbols(const long NumSyms);
  public:
    explicit symbol(const std::string& head);                             // symbolic name with 0 subscripts
    symbol(const std::string& head, long subscript);                      // symbolic name with 1 subscript
    symbol(const std::string& head, long subscript1, long subscript2);    // symbolic name with 2 subscripts
    symbol(const std::string& head, const std::vector<long>& subscripts); // symbolic name with many subscripts
    // Default copy ctor and assignment are fine.
    // Default dtor is fine.

    void myOutputSelf(std::ostream& out) const;
    void myInput(std::istream& in);

    int myCmp(const symbol& sym2) const noexcept;
    static bool IsValidHead(const std::string& head) noexcept;
    static bool IsValidFirstChar(char ch) noexcept;
  private: // data members
    std::string myHead;
    std::vector<long> mySubscripts;
    friend const std::string& head(const symbol& sym) noexcept;
    friend long NumSubscripts(const symbol& sym) noexcept;
    friend long subscript(const symbol& sym, long n);
  };


  symbol NewSymbol();  ///< pseudo ctor for anonymous symbols
  std::vector<symbol> NewSymbols(const long NumSyms); // several anon symbols
  const std::string& head(const symbol& sym) noexcept;     // inline
  long NumSubscripts(const symbol& sym) noexcept;          // inline
  long subscript(const symbol& sym, long n);
  bool AreDistinct(const std::vector<symbol>& syms) noexcept;
  bool AreArityConsistent(const std::vector<symbol>& syms) noexcept;
  std::ostream& operator<<(std::ostream& out, const symbol& sym);
  std::istream& operator>>(std::istream& in, symbol& sym);
  int cmp(const symbol& sym1, const symbol& sym2) noexcept;         // inline
  bool operator==(const symbol& sym1, const symbol& sym2) noexcept; // inline
  bool operator!=(const symbol& sym1, const symbol& sym2) noexcept; // inline
  bool operator<(const symbol& sym1, const symbol& sym2) noexcept;  // inline
  bool operator<=(const symbol& sym1, const symbol& sym2) noexcept; // inline
  bool operator>(const symbol& sym1, const symbol& sym2) noexcept;  // inline
  bool operator>=(const symbol& sym1, const symbol& sym2) noexcept; // inline

  std::vector<symbol> SymbolRange(const symbol& sym1, const symbol& sym2);
  std::vector<symbol> SymbolRange(const std::string& head, long lo, long hi); // <==> SymbolRange(symbol(head,lo),symbol(head,hi));
  std::vector<symbol> symbols(const std::string& str);
  std::vector<symbol> symbols(const std::vector<std::string>& s);

  //---------------------------------------------------------------------------
  // Inline implementations

  inline const std::string& head(const symbol& sym) noexcept
  {
    return sym.myHead;
  }


  inline long NumSubscripts(const symbol& sym) noexcept
  {
    return len(sym.mySubscripts);
  }


  inline int cmp(const symbol& sym1, const symbol& sym2) noexcept
  {
    return sym1.myCmp(sym2);
  }

  inline bool operator==(const symbol& sym1, const symbol& sym2) noexcept
  {
    return sym1.myCmp(sym2) == 0;
  }

  inline bool operator!=(const symbol& sym1, const symbol& sym2) noexcept
  {
    return sym1.myCmp(sym2) != 0;
  }

  inline bool operator<(const symbol& sym1, const symbol& sym2) noexcept
  {
    return sym1.myCmp(sym2) < 0;
  }

  inline bool operator<=(const symbol& sym1, const symbol& sym2) noexcept
  {
    return sym1.myCmp(sym2) <= 0;
  }

  inline bool operator>(const symbol& sym1, const symbol& sym2) noexcept
  {
    return sym1.myCmp(sym2) > 0;
  }

  inline bool operator>=(const symbol& sym1, const symbol& sym2) noexcept
  {
    return sym1.myCmp(sym2) >= 0;
  }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/symbol.H,v 1.21 2022/02/18 14:11:52 abbott Exp $
// $Log: symbol.H,v $
// Revision 1.21  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.20  2021/02/10 19:34:25  abbott
// Summary: Added noexcept (redmine 1572), and changed indices to long from MachineInt (redmine 925)
//
// Revision 1.19  2021/01/07 14:48:34  abbott
// Summary: Corrected copyright
//
// Revision 1.18  2020/12/04 21:35:09  abbott
// Summary: Added IsValidFirstChar
//
// Revision 1.17  2016/02/01 12:42:55  abbott
// Summary: Removed cruft - the old defns of symbols.
//
// Revision 1.16  2015/07/27 11:47:50  bigatti
// -- function "symbols(string)" now reads comma separated symbols
// -- function "symbols(string,string,..)" now disabled
//
// Revision 1.15  2012/05/29 14:56:15  abbott
// Separated ThreadsafeCounter from symbol; also employed it in ring.C.
//
// Revision 1.14  2012/05/29 07:48:13  abbott
// Simplified private ctor for anonymous symbols.
// Some general cleaning.
//
// Revision 1.13  2012/05/24 14:46:20  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.12  2012/05/10 14:43:27  abbott
// Added new fns NewSymbol & NewSymbols.
// Also added ThreadsafeCounter - will be separated off later.
//
// Revision 1.11  2011/11/09 13:51:29  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.10  2011/03/14 10:30:08  abbott
// Changed size_t into long (in fn interfaces).
//
// Revision 1.9  2011/02/16 14:58:10  bigatti
// -- added    symbols(vector<string> s)
//
// Revision 1.8  2008/12/11 20:12:32  abbott
// Now using MachineInt type in ctors to specify the values of symbol subscripts
// (but we check that the values fit into a long).  Corrected two bugs in
// SymbolRange (an undetected overflow and an off-by-1 error).
//
// Revision 1.7  2008/12/11 10:54:01  abbott
// Now subscripts are [long] values rather than [int].  The extra range will probably
// never be useful, but I prefer to offer it all the same.
//
// Revision 1.6  2008/11/21 21:16:14  abbott
// Added ctor for symbols with 2 subscripts.
//
// Revision 1.5  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.4  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
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
// Revision 1.4  2006/11/27 14:20:09  cocoa
// -- removed #include "CoCoA/config.H"
//
// Revision 1.3  2006/10/27 19:09:45  cocoa
// Replaced some member functions of CoCoA::symbol by friend functions.
// Removed some include dependency on symbol.H
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
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
// Revision 1.1  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//

#endif
