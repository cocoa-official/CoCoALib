#ifndef CoCoA_OpenMathXML_H
#define CoCoA_OpenMathXML_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/OpenMath.H"
#include "CoCoA/BigInt.H"

#include <iostream>
// using std::ostream;
#include <stack>
// using std::stack;
#include <string>
// using std::string


namespace CoCoA
{

  class OpenMathOutputXML: public OpenMathOutputBase
  {
  public:
    OpenMathOutputXML(std::ostream& out);
    ~OpenMathOutputXML();
    void mySend(const MachineInt& n) override;
    void mySend(const BigInt& N) override;
    void mySend(const OpenMathSymbol& s) override;
    void mySendApplyStart() override;
    void mySendApplyEnd() override;
    void mySendObjectStart() override;
    void mySendObjectEnd() override;
    //void OutputAttribute();//??????
  private:
    std::ostream& myOut;  // UGLY, RISKY!!???!!
    long myLevel;  // nesting level
    std::stack<OpenMathTag> myTagStack;  // only really needed for debugging
    const char* myIndent();
  };


  class OpenMathInputXML: public OpenMathInputBase
  {
  public:
    OpenMathInputXML(std::istream& in); // just sets ref count to zero.
    virtual ~OpenMathInputXML();

    void advance() override;
    OpenMathTag myCurrTag() override;
    long NumDescendants() const override;
    bool myRecv(long & n) override; // SHORT INTEGER
//????    virtual bool myRecv(unsigned long & n); // UNSIGNED SHORT INTEGER
    bool myRecv(BigInt& N) override;    // LONG INTEGER
    bool myRecv(OpenMathSymbol& s) override;
  private: // data members
    enum class ReadStatus { NotYetRead, AlreadyRead };
    static constexpr ReadStatus NotYetRead = ReadStatus::NotYetRead;
    static constexpr ReadStatus AlreadyRead = ReadStatus::AlreadyRead;
    ReadStatus myReadStatus;
    OpenMathTag myCurrTagType;
    BigInt myIntValue;
    OpenMathSymbol mySymbol;
    std::istream& myIn;  // UGLY, RISKY!!???!!
    long myLevel;  // nesting level
    std::stack<OpenMathTag> myTagStack;  // only really needed for debugging
  private:
    char ReadChar();
    char SkipWSReadChar();
    void ReadNextNode();
    bool SkipMatch(const std::string& expected);
    bool ReadDecimalString(std::string& DecimalDigits);
    bool ReadStringInQuotes(std::string& QuotedString);
  };


} // end of namespace CoCoA

#endif
