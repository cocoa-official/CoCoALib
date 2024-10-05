#ifndef CoCoA5_Lexer_H
#define CoCoA5_Lexer_H

//   Copyright (c) 2009-2020   John Abbott,  Anna M. Bigatti
//   Original author:  2009 Giovanni Lagorio (lagorio@disi.unige.it)
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


#include <string>
#include <cassert>
#include <stack>

#include "LineProviders.H"

namespace CoCoA {

namespace InterpreterNS {
	class OSTREAM;
	class CppOSTREAM;
	class RuntimeException;
	class InterruptException;
} // namespace InterpreterNS {

namespace LexerNS {

class CharPointer;
class Lexer;

class Line : public ReferenceCountedObject {
private:
	friend class Lexer;
	friend class CharPointer;
	boost::intrusive_ptr<const Line> nextLine;
	Line() :
		number(0),
		chars(1, '\0')
	{}
public:
	const int number;
	const std::string chars;
	const boost::intrusive_ptr<const LineProvider> provider;
	Line(const int number, const std::string &chars, const boost::intrusive_ptr<const LineProvider> provider) :
		number(number),
		chars(chars),
		provider(provider)
	{
		assert(number>=0);
	}
	const boost::intrusive_ptr<const Line> getNextLineInBuffer() const { return nextLine; }
	//~Line() { std::cout << "Releasing line: " << chars; }
};

class CharPointer {
	friend class Lexer;
	boost::intrusive_ptr<const Line> line;
	size_t positionInLine;
	static boost::intrusive_ptr<const Line> emptyLine;
	CharPointer() :
		line(emptyLine),
		positionInLine(0)
	{}
public:
	static const CharPointer Null;
	CharPointer(boost::intrusive_ptr<const Line> line, const size_t positionInLine) :
		line(line),
		positionInLine(positionInLine)
	{
		assert(line);
		assert(positionInLine<line->chars.length());
	}
	inline boost::intrusive_ptr<const Line> getLine() const { return line; }
	inline size_t getPositionInLine() const { return positionInLine; }
	inline char operator *() const {
		return line->chars[positionInLine];
	}
	std::string toPrettyString() const;
	std::string stringTo(const CharPointer &end) const;
	bool operator==(const CharPointer &other) const {
		return positionInLine==other.positionInLine && line==other.line;
	}
	bool operator!=(const CharPointer &other) const {
		return positionInLine!=other.positionInLine || line!=other.line;
	}
};

enum TokenType {
	TT_ALIAS,
	TT_AND,
	TT_ASSIGN, // :=
	TT_BLOCK,
	TT_BREAK,
	TT_CART_PROD, // ><
	TT_CATCH,
	TT_CIAO,
	TT_CLEAR,
	//TT_CLOSE_COMPAT, // }$
	TT_CLOSED_ROUND, // )
	TT_CLOSED_SQUARE, // ]
	TT_COLON,
	TT_COLONCOLON, // ::
	TT_COMMA,
	TT_COMPATIBILITY_MARKER,
	TT_CONTINUE,
	TT_DEFINE,
	TT_DEGLEX,
	TT_DEGREVLEX,
	TT_DELETE,
	TT_DESCRIBE,
	TT_DESTROY,
	TT_DO,
	TT_DOT,
	TT_DOTDOT,
	TT_ELLIPSIS,
	TT_ELIF,
	TT_ELSE,
	TT_END,
	TT_ENDALIAS,
	TT_ENDBLOCK,
	TT_ENDCATCH,
	TT_ENDDEFINE,
	TT_ENDFOR,
	TT_ENDFOREACH,
	TT_ENDIF,
	TT_ENDLAMBDA,
	TT_ENDPACKAGE,
	TT_ENDREPEAT,
	TT_ENDTRY,
	TT_ENDUSING,
	TT_ENDWHILE,
	TT_EOF,
	TT_EQUAL,
	TT_EXIT,
	TT_EXPORT,
	TT_FALSE,
	TT_FLOAT_LITERAL,
	TT_FOR,
	TT_FOREACH,
	TT_GE, // >=
	TT_GT, // >
	TT_HELP, // ? .... <newline>
	TT_IDENTIFIER,
	TT_IF,
	TT_IMPORTBYREF,
	TT_IMPORTBYVALUE,
	TT_IN, // In
	TT_INT_LITERAL,
	TT_ISDEFINED,
	TT_ISIN, // IsIn
	TT_LAMBDA,
	TT_LE, // <=
	TT_LEX,
	TT_LOAD,
	TT_LT, // <
	TT_MINUS,
	TT_MOD, // %
	TT_NOTEQUAL, // <>
	TT_ON,
	//TT_OPEN_COMPAT, // ${
	TT_OPEN_ROUND, // (
	TT_OPEN_SQUARE, // [
	TT_OPT,
	TT_OR,
	//TT_ORD,
	TT_PACKAGE,
	TT_PACKAGENAME,
	TT_PIPE, // |
	TT_PLUS,
	TT_POSTO,
	TT_POWER, // ^
	TT_PRINT,
	TT_PRINTLN,
	TT_PROTECT,
	TT_QUIT,
	TT_RECORD,
	TT_REF,
	TT_REPEAT,
	TT_RETURN,
	TT_RING_ASSIGN, // ::=
	TT_SEMICOLON,
	TT_SET,
	TT_SKIP,
	TT_SLASH,
	TT_SOURCE,
	TT_SOURCE_AS_LSHIFT, // <<
	TT_SOURCEREGION,
	TT_STAR, // *
	TT_STEP,
	TT_STRING_LITERAL,
	TT_THEN,
	TT_TIME,
	TT_TO,
	TT_TOPLEVEL,
	TT_TOPOS,
	TT_TRUE,
	TT_TRY,
	TT_UNPROTECT,
	TT_UNSET,
	TT_UNTIL,
	TT_UPONERROR,
	TT_USE,
	TT_USING,
	TT_VAR,
  //	TT_WEIGHTS,  // disabled 2015-03-09
	TT_WHILE,
	TT_XEL
};

class Token {
	TokenType type;
	CharPointer begin, end;
	Token() :
		type(TT_EOF),
		begin(CharPointer::Null),
		end(CharPointer::Null)
	{}
public:
	static Token EndOfFile;
	Token(const CharPointer &begin, const CharPointer &end, TokenType type) :
		type(type),
		begin(begin),
		end(end)
	{
		assert(type!=TT_EOF);
		assert(*begin);
		assert(*end);
	}
	inline TokenType getType() const { return type; }
	inline const CharPointer &getBegin() const { return begin; }
	inline const CharPointer &getEnd() const { return end; }
        std::string LexemeClass() const;
	std::string lexeme() const;
};

enum WarningSeverity {
	WS_NONE=0,
	WS_LOW=1,
	WS_NORMAL=2,
	WS_PEDANTIC=3
};

class ErrorReporter : public ReferenceCountedObject {
	int errorCount, warningCount;
protected:
	static const std::string ContextPrefix, WherePrefix, ErrorPrefix, WarningPrefix, CalledbyPrefix;
	WarningSeverity warningLevel;
	virtual void implReportWarning(const std::string &msg)=0;
	virtual void implReportWarning(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to)=0;
	virtual void implReportInterrupt(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to)=0;
	virtual void implReportError(const std::string &msg)=0;
	virtual void implReportError(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to)=0;
public:
	const boost::intrusive_ptr<InterpreterNS::OSTREAM> outputStream;
	ErrorReporter(WarningSeverity warningLevel, const boost::intrusive_ptr<InterpreterNS::OSTREAM> outputStream);
	virtual void printContext();
	virtual void printWhere();
	virtual void printError();
	virtual void printWarning();
	virtual void printCalledBy();
	virtual void printBold(const std::string &s);
	void reportWarning(const std::string &msg, WarningSeverity severity = WS_NORMAL) {
		if (this->warningLevel>=severity) {
			++this->warningCount;
			this->implReportWarning(msg);
		}
	}
	void reportWarning(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to,  WarningSeverity severity = WS_NORMAL) {
		if (this->warningLevel>=severity) {
			++this->warningCount;
			this->implReportWarning(msg, from, to);
		}
	}
	void reportWarning(const std::string &msg, const LexerNS::Token &tok) {
		this->reportWarning(msg, tok.getBegin(), tok.getEnd());
	}
	void reportInterrupt(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) {
		++this->errorCount;
		this->implReportInterrupt(from, to);
	}
        void reportInterrupt(const InterpreterNS::InterruptException &intr);
	void reportError(const std::string &msg) {
		++this->errorCount;
		this->implReportError(msg);
	}
	void reportError(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) {
		++this->errorCount;
		this->implReportError(msg, from, to);
	}
	void reportError(const std::string &msg, const LexerNS::Token &tok) {
		this->reportError(msg, tok.getBegin(), tok.getEnd());
	}
	virtual void reportError(const InterpreterNS::RuntimeException &exception);
	int getErrorCount() const { return this->errorCount; }
	int getWarningCount() const { return this->warningCount; }
	void resetErrorCounts() { this->errorCount = this->warningCount = 0; }
	int getWarningLevel() const { return this->warningLevel; }
	virtual bool reportLineNumberWhenMeaningful(const CharPointer &fromPos, const CharPointer &toPos, bool printColumns, bool includeHeader);
	virtual ~ErrorReporter() /*throw ()*/;
};

class DefaultErrorReporter : public ErrorReporter {
	void outputUnderlinedChars(std::ostream &out, const CharPointer &from, const CharPointer &to);
protected:
	void implReportWarning(const std::string &msg);
	void implReportWarning(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void implReportInterrupt(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void implReportError(const std::string &msg);
	void implReportError(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
public:
	DefaultErrorReporter(WarningSeverity warningLevel, const boost::intrusive_ptr<InterpreterNS::CppOSTREAM> outputStream);
	~DefaultErrorReporter() /*throw ()*/;
};

class Lexer : public ReferenceCountedObject {
	boost::intrusive_ptr<ErrorReporter> errorReporter;
	boost::intrusive_ptr<LineProvider> lineProvider;
	boost::intrusive_ptr<Line> previousLine, currentLine;
	size_t currentPositionInLine;
	std::stack<CharPointer> ungetCPs;
	std::stack<Token> ungetTokens;
	LexerStatus status;
	inline void ungetCP(CharPointer cp) { this->ungetCPs.push(cp); }
	void getCP(const ParserNS::ParserStatus &ps, CharPointer &result);
	inline CharPointer getCP(const ParserNS::ParserStatus &ps);
	inline CharPointer peekCP(const ParserNS::ParserStatus &ps);
	Token getIdentifierOrKeyword(const CharPointer &beginCP, const ParserNS::ParserStatus &ps);
	Token getPackageName(const CharPointer &beginCP, const ParserNS::ParserStatus &ps);
	Token getNumberLiteral(const CharPointer &beginCP, const ParserNS::ParserStatus &ps);
	Token getStringLiteral(const CharPointer &beginCP, const ParserNS::ParserStatus &ps);
	Token getHelpStatement(const CharPointer &beginCP, const ParserNS::ParserStatus &ps);
	void consumeSingleLineComment(const ParserNS::ParserStatus &ps);
	void checkErrorCounter();
	CharPointer lastReportedErrorPosition;
	bool warnAboutFeaturesThatMayBeRemoved;
	bool interactiveInput;
public:
	explicit Lexer(boost::intrusive_ptr<ErrorReporter> errorReporter, boost::intrusive_ptr<LineProvider> lineProvider, bool warnAboutFeaturesThatMayBeRemoved, bool interactiveInput) :
		errorReporter(errorReporter),
		lineProvider(lineProvider),
		currentPositionInLine(0),
		lastReportedErrorPosition(CharPointer::Null),
		warnAboutFeaturesThatMayBeRemoved(warnAboutFeaturesThatMayBeRemoved),
		interactiveInput(interactiveInput)
	{
		assert(this->errorReporter);
		assert(this->lineProvider);
	}
	static bool isKeyword(const std::string &);
	static bool isValidIdentifier(const std::string &);
	bool isInteractive() { return interactiveInput; }
	inline const LexerStatus & getStatus() const { return status; }
	inline void ungetToken(const Token &t) {
		ungetTokens.push(t);
	}
	inline void startingTopLevelStatement();
	Token getToken(const ParserNS::ParserStatus &ps);
	void resetCompatibilityMode();
	void reportError(const std::string &msg);
	void reportError(const std::string &msg, const CharPointer &from, const CharPointer &to, bool canThrow = true);
	void reportWarning(const std::string &msg, WarningSeverity severity);
	void reportWarning(const std::string &msg, const CharPointer &from, const CharPointer &to, WarningSeverity severity);
	void reportWarningMayBeRemoved(const std::string &prefixMsg, const CharPointer &from, const CharPointer &to);
};

class ExceptionWithSourcePosition : public BaseException {
public:
	const CharPointer from, to;
	const bool needsRecovery;
	ExceptionWithSourcePosition(const std::string &reason, const CharPointer &from, const CharPointer &to, bool needsRecovery) :
		BaseException(reason),
		from(from),
		to(to),
		needsRecovery(needsRecovery)
	{}
	ExceptionWithSourcePosition(const std::string &reason, const Token &token, bool needsRecovery) :
		BaseException(reason),
		from(token.getBegin()),
		to(token.getEnd()),
		needsRecovery(needsRecovery)
	{}
	~ExceptionWithSourcePosition() throw () {}
};

class LexerException : public ExceptionWithSourcePosition {
public:
	LexerException(const std::string &reason, const CharPointer &from, const CharPointer &to, bool needsRecovery) :
		ExceptionWithSourcePosition(reason, from, to, needsRecovery)
	{}
	LexerException(const std::string &reason, const Token &token, bool needsRecovery) :
		ExceptionWithSourcePosition(reason, token.getBegin(), token.getEnd(), needsRecovery)
	{}
	~LexerException() throw () {}
};

class TooManyErrorsException : public BaseException {
	static std::string humanReadable(LexerNS::TokenType tt);
public:
	static const int MAX_ERRORS = 15;
	TooManyErrorsException() :
		BaseException("Too many errors, bailing out")
	{}
	~TooManyErrorsException() throw () {}
};

inline void Lexer::startingTopLevelStatement() {
	this->previousLine = 0;
}


inline CharPointer Lexer::getCP(const ParserNS::ParserStatus &ps) {
	CharPointer result(CharPointer::Null);
	this->getCP(ps, result);
	return result;
}

inline CharPointer Lexer::peekCP(const ParserNS::ParserStatus &ps) {
	CharPointer cp = this->getCP(ps);
	this->ungetCP(cp);
	return cp;
}

} // namespace LexerNS
} // namespace CoCoA


// Close include guard:
#endif
