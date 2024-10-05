#ifndef CoCoA5_LineProviders_H
#define CoCoA5_LineProviders_H

//   Copyright (c) 2009-2017  John Abbott,  Anna M. Bigatti
//   Original author: 2009 Giovanni Lagorio (lagorio@disi.unige.it)
//
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


#include <boost/intrusive_ptr.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <exception>
#include <vector>
#include <string>

namespace CoCoA
{

  namespace ParserNS
  {

    class ParserStatus //: boost::noncopyable
    {
      friend class Parser;
      bool inTheMiddleOfAStatement;
      bool inRecoveryMode;
      std::vector<std::string> openBlocks;
    public:
      ParserStatus() :
          inTheMiddleOfAStatement(false),
          inRecoveryMode(false)
	{}
      ParserStatus(const ParserStatus&) = delete;             // noncopyable
      ParserStatus& operator=(const ParserStatus&) = delete;  // noncopyable
      bool isInTheMiddleOfAStatement() const { return inTheMiddleOfAStatement; }
      bool isInRecoveryMode() const { return inRecoveryMode; }
      std::string promptHeader() const { return (openBlocks.empty() ? "" : openBlocks.back()); }
    };

    struct AskingForNewInteractiveInputDuringRecoveryException : std::exception {};

  } // namespace ParserNS

  namespace LexerNS
  {

    class LexerStatus //: boost::noncopyable
    {
      friend class Lexer;
      bool inStringLiteral;
      bool inMultiLineComment;
      bool inCoCoA4CompatibilityMode;
    public:
      LexerStatus() :
          inStringLiteral(false),
          inMultiLineComment(false),
          inCoCoA4CompatibilityMode(false)
	{}
      LexerStatus(const LexerStatus&) = delete;            // noncopyable
      LexerStatus& operator=(const LexerStatus&) = delete; // noncopyable
      bool isInStringLiteral() const { return inStringLiteral; }
      bool isInMultiLineComment() const { return inMultiLineComment; }
      bool isInCoCoA4CompatibilityMode() const { return inCoCoA4CompatibilityMode; }
    };

    class ReferenceCountedObject //: boost::noncopyable
    {
      mutable int refCount;
      friend void intrusive_ptr_add_ref(const ReferenceCountedObject *);
      friend void intrusive_ptr_release(const ReferenceCountedObject *);
    public:
      int getRefCount() const { return this->refCount; }
      ReferenceCountedObject() : refCount(0) {}
      ReferenceCountedObject(const ReferenceCountedObject&) = delete;             // noncopyable
      ReferenceCountedObject& operator=(const ReferenceCountedObject&) = delete;  // noncopyable
#ifdef PRINT_REFCOUNTS
      inline std::string refCountAsString() const
        {
        std::ostringstream os;
        os << " /* refCount=" << this->refCount << " */";
        return os.str();
      }
#endif // #ifdef PRINT_REFCOUNTS
      inline void dumpRefCountAsString(std::ostream &
#ifdef PRINT_REFCOUNTS
                                       out
#endif // #ifdef PRINT_REFCOUNTS
                                       ) const
        {
#ifdef PRINT_REFCOUNTS
        out << this->refCountAsString();
#endif // #ifdef PRINT_REFCOUNTS
      }
      virtual ~ReferenceCountedObject() = 0;
    };

    inline void intrusive_ptr_add_ref(const ReferenceCountedObject *o)
    {
      ++o->refCount;
      assert(o->refCount<INT_MAX); // *not* <= (otherwise, the next increment would bring it to INT_MIN, which is indeed <= INT_MAX!)
    }

    inline void intrusive_ptr_release(const ReferenceCountedObject *o)
    {
      assert(o->refCount>=0);
      if (--o->refCount==0)
        delete o;
    }

    class LineProvider : public ReferenceCountedObject
    {
      int lineNumber;
    protected:
      std::string filename; // used only in FileLineProvider & FileRegionLineProvider, o/w left empty.
      virtual bool doReadNextLine(const LexerStatus &, const ParserNS::ParserStatus &, std::string &chars) = 0;
    public:
      LineProvider() : lineNumber(0) {}
      inline int currentLineNumber() const { return lineNumber; }
      // Knowing the current state of the Lexer and the Parser, the ReadlineProvider (that is, the one
      // that fetches lines from the console) can provide meaningful prompts (indicating, for instance,
      // the nesting level of the current statement) and warning (for instance, newline inside a string literal).
      const std::string& myFileName() const { return filename; }
      bool IamReadingFromFile() const { return !myFileName().empty(); }
      inline bool readNextLine(const LexerStatus &ls, const ParserNS::ParserStatus &ps, std::string &chars)
      {
        const bool b = doReadNextLine(ls, ps, chars);
        if (b)
          ++lineNumber;
        return b;
      }
      virtual ~LineProvider() {}
    };


    class InteractiveLineProvider : public LineProvider
    {
    public:
      static std::string promptSuffix;
    protected:
      static std::string prompt(const LexerStatus &, const ParserNS::ParserStatus &);
    };

    class GetlineLineProvider : public InteractiveLineProvider
    {
    protected:
      bool doReadNextLine(const LexerStatus &, const ParserNS::ParserStatus &, std::string &chars);
    };

#ifdef CoCoA_WITH_READLINE
    class ReadlineLineProvider : public InteractiveLineProvider
    {
    protected:
      bool doReadNextLine(const LexerStatus &, const ParserNS::ParserStatus &, std::string &chars);
    };
#endif

    class FileLineProvider : public LineProvider {
      boost::iostreams::filtering_istream input;
    protected:
      bool doReadNextLine(const LexerStatus &, const ParserNS::ParserStatus &, std::string &chars);
    public:
#ifdef C5IDE
      std::string wholeFile;
#endif // #ifdef C5IDE
      FileLineProvider(const std::string &filename);
    };

    class FileRegionLineProvider : public LineProvider {
      boost::iostreams::filtering_istream input;
      const long myFromLine;
      const long myFromChar;
      const long myToLine;
      const long myToChar;
      long myCurrLine;
    protected:
      bool doReadNextLine(const LexerStatus &, const ParserNS::ParserStatus &, std::string &chars);
    public:
#ifdef C5IDE
      std::string wholeFile;
#endif // #ifdef C5IDE
      FileRegionLineProvider(const std::string &filename, long FromLine, long FromChar, long ToLine, long ToChar);
    };

    class BaseException : public std::exception
    {
    public:
      const std::string reason;
      explicit BaseException(const std::string &reason) : reason(reason.empty()?">>>> no description given <<<<":reason)
      {}
      const char *what() const throw () { return reason.c_str(); }
      ~BaseException() throw () {}
    };

    struct IOException : public BaseException
    {
      explicit IOException(const std::string &reason) :
          BaseException(reason)
	{}
    };

    struct FailedAssertionException : public BaseException
    {
      explicit FailedAssertionException(const std::string &reason) :
          BaseException(reason)
	{}
    };

  } // namespace LexerNS
} // namespace CoCoA


// Close include guard:
#endif
