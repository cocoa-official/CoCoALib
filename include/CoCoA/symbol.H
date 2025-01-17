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

#endif
