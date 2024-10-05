//   Copyright (c)  2009  John Abbott,  Anna M Bigatti
//   Orig author: 2009 Giovanni Lagorio (lagorio@disi.unige.it)
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

#include <iostream>
#include <boost/algorithm/string/case_conv.hpp>
// case_conv:  several calls to to_lower_copy
#include <boost/scope_exit.hpp>
//#include <cstring>
//#include <cctype>
#include <sstream>
//#include "CoCoA/library.H" // already included by AST->Parser->Interpreter
#include "C5Utils.H"
#include "Interpreter.H"


namespace CoCoA {
namespace LexerNS {

using namespace std;
using namespace boost;
using namespace CoCoA;

void Lexer::reportWarningMayBeRemoved(const string &prefixMsg, const CharPointer &from, const CharPointer &to)
{
	if (warnAboutFeaturesThatMayBeRemoved)
		this->reportWarning(prefixMsg+" will be probably removed in the final version of CoCoA 5", from, to, WS_LOW);
}

// used only by ctor for UnexpectedTokenException (see Parser.H:158)
  std::string Token::LexemeClass() const
  {
    if (type==TT_EOF) return string(""); // empty string
    if (type == TT_STRING_LITERAL) return string("string literal");
    if (type == TT_IDENTIFIER) return string("identifier");
    if (type == TT_INT_LITERAL || type == TT_FLOAT_LITERAL)
      return string("number");
    const std::string tok = lexeme(); // wasteful copy
    if (tok.empty()) return string("???NOTHING???"); // should never do this
    const char FirstChar = tok[0];
    if (isalpha(FirstChar)) return string("keyword");
    return string("operator/symbol");
  }

string Token::lexeme() const {
	return type==TT_EOF ? begin.toPrettyString() : begin.stringTo(end);
}

void Lexer::resetCompatibilityMode() {
	status.inCoCoA4CompatibilityMode=false;
}

void Lexer::consumeSingleLineComment(const ParserNS::ParserStatus &ps) {
	for(;;) {
		char c = *this->getCP(ps);
		if (c=='\n' || !c)
			return;
	}
}

Token Lexer::getHelpStatement(const CharPointer &begin, const ParserNS::ParserStatus &ps) {
	CharPointer end=begin;
	for(;;) {
		CharPointer cp = getCP(ps);
		const char c = *cp;
		if (c=='\n' || !c)
			return Token(begin, end, TT_HELP);
		end = cp;
	}
}


void Lexer::getCP(const ParserNS::ParserStatus &ps, CharPointer &result) {
	if (!this->ungetCPs.empty()) {
		result = this->ungetCPs.top();
		this->ungetCPs.pop();
		return;
	}
fetch:
	while (!this->currentLine || this->currentPositionInLine>=this->currentLine->chars.length()) {
                if (this->status.isInMultiLineComment())
                  this->reportWarning("Multiline comments are OBSOLESCENT; please comment each line separately", WS_PEDANTIC);
		std::string chars;
		if (!this->lineProvider->readNextLine(this->status, ps, chars)) {
			result = CharPointer::Null;
			return;
		}
		chars.append("\n");
		this->currentLine = new Line(lineProvider->currentLineNumber(), chars, lineProvider);
		if (this->previousLine)
			this->previousLine->nextLine = currentLine;
		this->previousLine = currentLine;
		this->currentPositionInLine = 0;
	}
	if (!this->previousLine) {
		const char currChar = this->currentLine->chars[this->currentPositionInLine];
		if (!IsBlank(currChar)) {
			if ( (currChar=='-' || currChar=='/') && this->currentPositionInLine<(this->currentLine->chars.length()-1) && currChar==this->currentLine->chars[this->currentPositionInLine+1]) {
				this->currentLine = nullptr;
				goto fetch;
			}
			this->previousLine = this->currentLine;
		}
	}
	result.line = this->currentLine;
	result.positionInLine = this->currentPositionInLine++;
}


Token Lexer::getToken(const ParserNS::ParserStatus &ps) {
	if (!this->ungetTokens.empty()) {
		const Token t = this->ungetTokens.top();
		this->ungetTokens.pop();
		return t;
	}
	for(;;) {
		CharPointer begin = getCP(ps), end=begin;
		const char c = *begin;
		switch (c) {
		case ' ':
		case '\n':
		case '\r':
		case '\t':
		case '\f':
			continue;
		case '\0':
			return Token::EndOfFile;
		case ';':
			if (this->status.inCoCoA4CompatibilityMode) {
				this->status.inCoCoA4CompatibilityMode = false;
				this->reportError("Semicolons are not allowed inside backward-compatibility markers (exiting backward-compatibility mode)", begin, begin);
			}
			return Token(begin, begin, TT_SEMICOLON);
		case '|':
			return Token(begin, begin, TT_PIPE);
		case '^':
			return Token(begin, begin, TT_POWER);
		case ',':
			return Token(begin, begin, TT_COMMA);
		case '%':
			return Token(begin, begin, TT_MOD);
		case '+':
			return Token(begin, begin, TT_PLUS);
		case '(':
			return Token(begin, begin, TT_OPEN_ROUND);
		case ')':
			return Token(begin, begin, TT_CLOSED_ROUND);
		case '[':
			return Token(begin, begin, TT_OPEN_SQUARE);
		case ']':
			return Token(begin, begin, TT_CLOSED_SQUARE);
		case '*':
			end = this->getCP(ps);
			if (*end=='*') {
				if (*this->peekCP(ps)=='*') {
					end = this->getCP(ps);
					assert(*end=='*');
					this->status.inCoCoA4CompatibilityMode = !this->status.inCoCoA4CompatibilityMode;
					return Token(begin, end, TT_COMPATIBILITY_MARKER);
				}
			}
			this->ungetCP(end);
			return Token(begin, begin, TT_STAR);
		case '=':
			return Token(begin, begin, TT_EQUAL);
		case '<':
			end = this->getCP(ps);
			switch (*end) {
			case '=':
				return Token(begin, end, TT_LE);
			case '>':
				return Token(begin, end, TT_NOTEQUAL);
			case '<':
				return Token(begin, end, TT_SOURCE_AS_LSHIFT);
			}
			this->ungetCP(end);
			return Token(begin, begin, TT_LT);
		case '>':
			end = this->getCP(ps);
			switch (*end) {
			case '=':
				return Token(begin, end, TT_GE);
			case '<':
				return Token(begin, end, TT_CART_PROD);
			}
			this->ungetCP(end);
			return Token(begin, begin, TT_GT);
		case '$':
			/*end = this->getCP(ps);
			switch (*end) {
			case '{':
				{
					if (this->status.isInCoCoA4CompatibilityMode())
						throw LexerException("Compatibility parentheses cannot be nested", begin, end, false);
					this->status.inCoCoA4CompatibilityMode = true;
					return Token(begin, end, TT_OPEN_COMPAT);
				}
			case '/':
				goto closeCompatibilityParenthesis;
			}
			this->ungetCP(end);*/
			return this->getPackageName(begin, ps);
		/*case '}':
			end = this->getCP(ps);
			if (*end=='$') {
				closeCompatibilityParenthesis:
				if (!this->status.isInCoCoA4CompatibilityMode())
					throw LexerException("Trying to close an unopened compatibility parenthesis", begin, end, false);
				this->status.inCoCoA4CompatibilityMode = false;
				return Token(begin, end, TT_CLOSE_COMPAT);
			}
			this->ungetCP(end);
			break;*/
		case '.':
			end = this->getCP(ps);
			if (*end=='.') {
				CharPointer cp = this->getCP(ps);
				if (*cp=='.')
					return Token(begin, cp, TT_ELLIPSIS);
				this->ungetCP(cp);
				return Token(begin, end, TT_DOTDOT);
			}
			this->ungetCP(end);
			return Token(begin, begin, TT_DOT);
		case '-':
			if (*this->peekCP(ps)=='-') {
				this->consumeSingleLineComment(ps);
				continue;
			}
			return Token(begin, begin, TT_MINUS);
		case '?':
			return this->getHelpStatement(begin, ps);
		case '/':
			switch(*this->peekCP(ps)) {
			case '/':
				this->consumeSingleLineComment(ps);
				continue;
			case '*':
				assert(!this->status.inMultiLineComment);
				this->status.inMultiLineComment=true;
				{
					BOOST_SCOPE_EXIT( (&status) ) {
						status.inMultiLineComment=false;
					} BOOST_SCOPE_EXIT_END
					this->getCP(ps); // consumes the '*'

					CharPointer lastReadCP = begin;
					CharPointer cp(CharPointer::Null);
                                        char prev = ' '; // any value except '/' and '*'
					for(;;) {
						getCP(ps, cp);
						const char ch = *cp;
						if (!ch)
                                                  throw LexerException("Unclosed multi-line comment /*...*/", begin, lastReadCP, false);
                                                if (prev == '/' && ch == '*')
                                                  this->reportWarning("Possible nested comment!", lastReadCP, cp, WS_LOW);
						lastReadCP = cp;
                                                if (prev == '*' && ch == '/')
							break;
                                                prev = ch;
					}
				}
				continue;
			/*case '$':
				{
					this->getCP(ps);
					if (this->status.isInCoCoA4CompatibilityMode())
						throw LexerException("Compatibility parenthesis already open", begin, end, false);
					this->status.inCoCoA4CompatibilityMode = true;
					return Token(begin, end, TT_OPEN_COMPAT);
				}
			*/
			}
			return Token(begin, begin, TT_SLASH);
		case '\'':
			this->reportWarning("OBSOLESCENT string literal; please use double quotes as delimiters", WS_LOW);
                        // fall through
		case '\"':
			return this->getStringLiteral(begin, ps);
		case ':':
			end = this->getCP(ps);
			switch(*end) {
			case '=':
				return Token(begin, end, TT_ASSIGN);
			case ':':
				if (*this->peekCP(ps)=='=')
					return Token(begin, getCP(ps), TT_RING_ASSIGN);
				return Token(begin, end, TT_COLONCOLON);
			}
			this->ungetCP(end);
			return Token(begin, begin, TT_COLON);
		}
		if (isalpha(c) || c=='_')
			return this->getIdentifierOrKeyword(begin, ps);
		if (isdigit(c))
			return this->getNumberLiteral(begin, ps);
		//if (c=='{' || c=='}')
		//	throw LexerException("Unknown symbol "+begin.toPrettyString()+", did you forget a $? Backward compatibility parentheses are ${ and }$", begin, begin, true);
		throw LexerException("Unknown symbol "+begin.toPrettyString(), begin, begin, true);
	}
}

namespace {
	struct Keyword {
		const string realId;
		const string lowerId;
		const TokenType tokenType;
		Keyword(const string & realId, TokenType tokenType) :
			realId(realId),
			lowerId(to_lower_copy(realId)),
			tokenType(tokenType)
		{}
	} keywordTable [] = {
			Keyword("Alias", TT_ALIAS),
			Keyword("And", TT_AND),
			Keyword("Block", TT_BLOCK),
			Keyword("Break", TT_BREAK),
			Keyword("Catch", TT_CATCH),
			Keyword("Ciao", TT_CIAO),
			Keyword("Clear", TT_CLEAR),
			Keyword("Continue", TT_CONTINUE),
			Keyword("Define", TT_DEFINE),
			Keyword("DegLex", TT_DEGLEX),
			Keyword("DegRevLex", TT_DEGREVLEX),
			Keyword("Delete", TT_DELETE),
			Keyword("Describe", TT_DESCRIBE),
			Keyword("Destroy", TT_DESTROY),
			Keyword("Do", TT_DO),
			Keyword("Elif", TT_ELIF),
			Keyword("Else", TT_ELSE),
			Keyword("End", TT_END),
			Keyword("EndAlias", TT_ENDALIAS),
			Keyword("EndBlock", TT_ENDBLOCK),
			Keyword("EndCatch", TT_ENDCATCH),
			Keyword("EndDefine", TT_ENDDEFINE),
			Keyword("EndForeach", TT_ENDFOREACH),
			Keyword("EndFor", TT_ENDFOR),
			Keyword("EndFunc", TT_ENDLAMBDA),
			Keyword("EndIf", TT_ENDIF),
			Keyword("EndPackage", TT_ENDPACKAGE),
			Keyword("EndRepeat", TT_ENDREPEAT),
			Keyword("EndTry", TT_ENDTRY),
			Keyword("EndUsing", TT_ENDUSING),
			Keyword("EndWhile", TT_ENDWHILE),
			Keyword("Exit", TT_EXIT),
			Keyword("Export", TT_EXPORT),
			Keyword("False", TT_FALSE),
			Keyword("Foreach", TT_FOREACH),
			Keyword("For", TT_FOR),
			Keyword("Func", TT_LAMBDA),
			Keyword("If", TT_IF),
			Keyword("ImportByRef", TT_IMPORTBYREF),
			Keyword("ImportByValue", TT_IMPORTBYVALUE),
			Keyword("In", TT_IN),
			Keyword("IsDefined", TT_ISDEFINED),
			Keyword("IsIn", TT_ISIN),
			Keyword("Lex", TT_LEX),
			Keyword("Load", TT_LOAD),
			Keyword("On", TT_ON),
			Keyword("Opt", TT_OPT),
			Keyword("Or", TT_OR),
			Keyword("Package", TT_PACKAGE),
			Keyword("PosTo", TT_POSTO),
			Keyword("PrintLn", TT_PRINTLN),
			Keyword("Print", TT_PRINT),
			Keyword("Protect", TT_PROTECT),
			Keyword("Quit", TT_QUIT),
			Keyword("Record", TT_RECORD),
			Keyword("Ref", TT_REF),
			Keyword("Repeat", TT_REPEAT),
			Keyword("Return", TT_RETURN),
			Keyword("Set", TT_SET),
			Keyword("Skip", TT_SKIP),
			Keyword("Source", TT_SOURCE),
			Keyword("SourceRegion", TT_SOURCEREGION),
			Keyword("Step", TT_STEP),
			Keyword("Then", TT_THEN),
			Keyword("Time", TT_TIME),
			Keyword("To", TT_TO),
			Keyword("TopLevel", TT_TOPLEVEL),
			Keyword("ToPos", TT_TOPOS),
			Keyword("True", TT_TRUE),
			Keyword("Try", TT_TRY),
			Keyword("Unset", TT_UNSET),
			Keyword("UponError", TT_UPONERROR),
			Keyword("Until", TT_UNTIL),
			Keyword("Unprotect", TT_UNPROTECT),
			Keyword("Use", TT_USE),
			Keyword("Using", TT_USING),
			Keyword("Var", TT_VAR),
      //			Keyword("Weights", TT_WEIGHTS),  // disabled 2015-03-09
			Keyword("While", TT_WHILE),
			Keyword("Xel", TT_XEL)
	};
}

Token Lexer::getNumberLiteral(const CharPointer &begin, const ParserNS::ParserStatus &ps) {
	CharPointer end = begin;
	for(;;) {
		CharPointer cp = this->getCP(ps);
		const char c = *cp;
		if (c=='.') {
			if (isdigit(*this->peekCP(ps))) {
				for(;;) {
					cp = this->getCP(ps);
					if (!isdigit(*cp)) {
						this->ungetCP(cp);
						break;
					}
					end = cp;
				}
				return Token(begin, end, TT_FLOAT_LITERAL);
			}
		}
		if (!isdigit(c)) {
			this->ungetCP(cp);
			break;
		}
		end = cp;
	}
	return Token(begin, end, TT_INT_LITERAL);
}

Token Lexer::getStringLiteral(const CharPointer &begin, const ParserNS::ParserStatus &ps)
{
  this->status.inStringLiteral = true;
  BOOST_SCOPE_EXIT( (&status) )
  {
    status.inStringLiteral = false;
  } BOOST_SCOPE_EXIT_END
  const char ClosingQuote = *begin;
  CharPointer end = begin;
  CharPointer previousCharP = begin;
  bool lastWasEscape=false;
  for(;;)
  {
    CharPointer cp = this->getCP(ps);
    const char ch = *cp;
    if (cp == CharPointer::Null) // was (!ch)   JAA changed 2020-03-03 to allow NULs in string literals
      unclosed: throw LexerException("Unclosed string literal", begin, end, false);
    previousCharP = end;
    end = cp;
    if (ch==ClosingQuote && !lastWasEscape)
      break;
    if (lastWasEscape)
    {
      switch (ch)
      {
      case 'n': // newline
      case 'r': // carriage return (without newline)
      case 't': // tab
      case 'a': // alert (i.e. beep)
      case '\'':
      case '\"':
      case '\\':
        break;
      case 'x': // hex character code
      {
        CharPointer d1 = this->getCP(ps);
        char ch = *d1;
        if (!ch)
          goto unclosed;
        if (!isxdigit(ch))
          throw LexerException("After \\x you must insert two hex-digits", end, d1, true);
        CharPointer d2 = this->getCP(ps);
        ch = *d2;
        if (!ch)
          goto unclosed;
        if (!isxdigit(ch))
          throw LexerException("After \\x you must insert two hex-digits", end, d2, true);
        end = d2;
        break;
      }
      default:
        this->reportError("Unrecognized escape sequence", previousCharP, end);
      }
    }
    lastWasEscape = !lastWasEscape && ch=='\\';
    if (ch=='\n')
      this->reportError("Unclosed string literal?  Or use \\n to put a newline in a string", end, end);
  }
  return Token(begin, end, TT_STRING_LITERAL);
}


Token Lexer::getPackageName(const CharPointer &begin, const ParserNS::ParserStatus &ps) {
	CharPointer end = begin;
	CharPointer cp = begin;
	char c;
	for(;;) {
		cp = this->getCP(ps);
		c = *cp;
		if (!(isalnum(c) || c=='_' || c=='/')) {
			this->ungetCP(cp);
			break;
		}
		end = cp;
	}
	return Token(begin, end, TT_PACKAGENAME);
}

Token Lexer::getIdentifierOrKeyword(const CharPointer &begin, const ParserNS::ParserStatus &ps) {
	CharPointer end = begin;
	if (this->status.isInCoCoA4CompatibilityMode() && islower(*begin)) {
		CharPointer cp = this->peekCP(ps);
		if (isupper(*cp))
			throw LexerException("This form of implicit multiplication is not allowed (anymore)", begin, cp, false);
		return Token(begin, begin, TT_IDENTIFIER);
	}
	for(;;) {
		CharPointer cp = this->getCP(ps);
		char c = *cp;
		if (!(isalnum(c) || c=='_')) {
			this->ungetCP(cp);
			break;
		}
		end = cp;
	}
	const size_t beginPos = begin.getPositionInLine();
	const size_t endPos = end.getPositionInLine();
	assert(begin.getLine()==end.getLine());
	const string chars = begin.getLine()->chars.substr(beginPos, endPos-beginPos+1);
	if (chars=="RECORD")
		return Token(begin, end, TT_IDENTIFIER);
	if (this->errorReporter->getWarningLevel()!=WS_PEDANTIC) {
		// upper cases "TRUE" and "FALSE" are for backward compatibility (they should be "True" and "False")
		if (chars=="TRUE")
			return Token(begin, end, TT_TRUE);
		if (chars=="FALSE")
			return Token(begin, end, TT_FALSE);
	}
	const string lowerChars = to_lower_copy(chars);
	TokenType type = TT_IDENTIFIER;
	for (const Keyword &k: keywordTable)
        {
          if (chars==k.realId || chars==k.lowerId) {
            type = k.tokenType;
            break;
          }
          if (lowerChars==k.lowerId) {
            string msg("\""+chars+"\" recognized as the keyword \""+k.realId+"\"; use \""+k.realId+"\" ");
            msg += status.isInCoCoA4CompatibilityMode() ?
              "to avoid this warning" :
              "or \""+k.lowerId+"\" to avoid this warning";
            //					  "to avoid this warning (note: outside CoCoA4 compatibility mode, \""+k.lowerId+"\" is allowed too)" :
            //					"or \""+k.lowerId+"\" to avoid this warning (note: in CoCoA4 compatibility mode, only the former casing is allowed)";
            //					"or \""+k.lowerId+"\" to avoid this warning (note: in CoCoA4 compatibility mode, only \""+k.realId+"\" is allowed)";
            this->reportWarning(msg, begin, end, WS_LOW);
            type = k.tokenType;
            break;
          }
	}
	if (type==TT_LAMBDA && this->status.inCoCoA4CompatibilityMode) {
		this->status.inCoCoA4CompatibilityMode = false;
		this->reportError("anonymous functions are not allowed inside backward-compatibility markers (exiting backward-compatibility mode)", begin, end);
	}
	return Token(begin, end, type);
}

bool Lexer::isKeyword(const string &s) {
	const string lowerChars = to_lower_copy(s);
	for (const Keyword &k: keywordTable)
        {
          if (lowerChars==k.lowerId)
            return true;
	}
	return false;
}

bool Lexer::isValidIdentifier(const string &s) {
	const int length = s.length();
	if (length==0)
		return false;
	if (!isalpha(s[0]))
		return false;
	for(int a=1; a<length; ++a) {
		const char c = s[a];
		if (!(isalnum(c) || c=='_'))
			return false;
	}
        return !isKeyword(s);
}

namespace {
	constexpr char HICHAR = '^';
	constexpr size_t MAX_PREFIX = 40;
	constexpr size_t MAX_SUFFIX = 30;

	inline void outputChar(ostream &out, char c, size_t n) {
		for(size_t i=0; i<n; ++i)
			out << c;
	}

	void outputUnderlinedLine(ostream &out, intrusive_ptr<const Line> line, size_t from, size_t to) {
		string prefix, suffix;
		const string &chars = line->chars;
        size_t printingStart=0, printingEnd=chars.length()==0 ? 0 : chars.length()-1;
		if (from>=MAX_PREFIX) {
			printingStart = from-MAX_PREFIX;
			prefix = "... ";
		}
		if ( (printingEnd-to)>=MAX_SUFFIX ) {
			printingEnd = to+MAX_SUFFIX;
			suffix = " ...";
		}
		assert(printingStart<=printingEnd);
		assert(printingEnd<chars.length());
		out << "--> " << prefix; // "--> " for emacs #567 issue
		for (size_t i=printingStart; i<=printingEnd; ++i)
                {
                  const char c = chars[i];
                  if (c == '\n') break;
                  if (c == '\t') { out << ' '; continue; }   // replace TAB by SPACE to avoid alignment problems
                  if (c < 32 || c > 126) { out << '#'; continue; }  // replace unprintable chars by #
                  out << c;  // c is a "normal, printable" char, so print it
		}
		out << suffix << '\n' << "--> "; // "--> " for emacs #567 issue
		outputChar(out, ' ', prefix.length()+from-printingStart);
		outputChar(out, HICHAR, to-from+1);
		out << '\n';
	}
}

void DefaultErrorReporter::outputUnderlinedChars(ostream &out, const CharPointer &from, const CharPointer &to) {
	if (!*from) {
//JAA 2013-05-06		out << "<End of input>" << endl;
		return;
	}
	intrusive_ptr<const Line> fromLine = from.getLine();
	intrusive_ptr<const Line> toLine = to.getLine();
	if (fromLine==toLine)
		return outputUnderlinedLine(out, fromLine, from.getPositionInLine(), to.getPositionInLine());
	if (!*to) {
//JAA 2013-05-06  ever get here???
		out << "... <End of input>" << endl;
		return;
	}
	outputUnderlinedLine(out, fromLine, from.getPositionInLine(), fromLine->chars.length()-2); // the last char of every line is \n and we don't want to underline it
	int skippedLines=0;
	intrusive_ptr<const Line> currLine = fromLine;
	while ( (currLine = currLine->getNextLineInBuffer())!= toLine )
		++skippedLines;
	if (skippedLines) {
		out << "--> ... (";
		if (skippedLines==1)
			out << "another source line";
		else
			out << "other " << skippedLines << " source lines";
		out << ") ...\n";
	}
	outputUnderlinedLine(out, toLine, 0, to.getPositionInLine());
}

intrusive_ptr<const Line> CharPointer::emptyLine(new Line());

namespace {
	string stripped(const string &filename) {
		const string::size_type pos = filename.rfind('/');
		return pos==string::npos ? filename : filename.substr(pos+1, string::npos);
	}
}

bool ErrorReporter::reportLineNumberWhenMeaningful(const CharPointer &fromPos, const CharPointer &toPos, bool printColumns, bool includeHeader) {
  if (fromPos == CharPointer::Null) return false; //JAA 2013-04-29 WORKAROUND HACK!!!!
	intrusive_ptr<const Line> l = fromPos.getLine();
        if (l->provider->IamReadingFromFile())
        {
		if (includeHeader)
			this->printWhere();
		const int l1 = l->number;
		const int l2 = toPos.getLine()->number;
		if (l1==l2) {
			this->outputStream->print(" at line ")->print(lexical_cast<string>(l1));
			if (printColumns)
				this->outputStream->print(" (column ")->print(lexical_cast<string>(fromPos.getPositionInLine()+1))->print(")");
		} else {
			this->outputStream->print(" at lines ")->print(lexical_cast<string>(l1));
			if (printColumns)
				this->outputStream->print(" (column ")->print(lexical_cast<string>(fromPos.getPositionInLine()+1))->print(")");
			this->outputStream->print(" .. ")->print(lexical_cast<string>(l2));
			if (printColumns)
				this->outputStream->print(" (column ")->print(lexical_cast<string>(toPos.getPositionInLine()+1))->print(")");
		}
		this->outputStream->print(" of ")->print(stripped(l->provider->myFileName()));
		return true;
	}
	return false;
}

ErrorReporter::ErrorReporter(WarningSeverity warningLevel, const intrusive_ptr<InterpreterNS::OSTREAM> outputStream) :
	errorCount(0),
	warningCount(0),
	warningLevel(warningLevel),
	outputStream(outputStream)
{}

  ErrorReporter::~ErrorReporter() /*throw ()*/ {}

DefaultErrorReporter::DefaultErrorReporter(WarningSeverity warningLevel, const intrusive_ptr<InterpreterNS::CppOSTREAM> outputStream) :
	ErrorReporter(warningLevel, outputStream)
{}

DefaultErrorReporter::~DefaultErrorReporter() {}

void Lexer::checkErrorCounter() {
	if (this->errorReporter->getErrorCount() >= TooManyErrorsException::MAX_ERRORS)
		throw TooManyErrorsException();
}

void Lexer::reportError(const string &msg) {
	this->errorReporter->reportError(msg);
	this->checkErrorCounter();
}

void Lexer::reportError(const string &msg, const CharPointer &from, const CharPointer &to, bool canThrow) {
	if (from!=CharPointer::Null && from==lastReportedErrorPosition) // may happen because of parsing error recovery
		return;
	if (canThrow && this->isInteractive())
		throw LexerException(msg, from, to, false);
	this->lastReportedErrorPosition = from;
	this->errorReporter->reportError(msg, from, to);
	this->checkErrorCounter();
}

void Lexer::reportWarning(const string &msg, WarningSeverity severity) {
	this->errorReporter->reportWarning(msg, severity);
}

void Lexer::reportWarning(const string &msg, const CharPointer &from, const CharPointer &to, WarningSeverity severity) {
	this->errorReporter->reportWarning(msg, from, to, severity);
}

string CharPointer::toPrettyString() const {
	const unsigned char c = **this;
	if (c>' ' && c <='~') // visible chars
		return '\"'+std::string(1, c)+'\"';
	switch(c) {
	case '\0': return "<end of input>";
	case '\t': return "<tab>";
	case '\n': return "<newline>";
	case '\r': return "<carriage return>";
	case ' ' : return "<space>";
	default  : ostringstream os;
			   os << "<unprintable char (bytecode=" << static_cast<int>(c) << ")>";
			   return os.str();
	}
}

string CharPointer::stringTo(const CharPointer &end) const {
	assert(line);
	assert(end.line);
	std::string result;
	boost::intrusive_ptr<const Line> l = line;
	size_t pos = positionInLine;
	for(; l!=end.line; l=l->getNextLineInBuffer()) {
		assert(l);
		assert(pos<l->chars.length());
		result.append(l->chars.substr(pos));
		pos=0;
	}
	assert(pos<=end.positionInLine);
	result.append(l->chars.substr(pos, end.positionInLine-pos+1));
	return result;
}

const CharPointer CharPointer::Null;

Token Token::EndOfFile;

const string
  ErrorReporter::ContextPrefix("--> CONTEXT: "), // "--> " for issue #148
  ErrorReporter::WherePrefix("\n--> WHERE:"),
  ErrorReporter::ErrorPrefix("--> ERROR: "), // "--> " for issue #567
  ErrorReporter::WarningPrefix("--> WARNING: "),
  ErrorReporter::CalledbyPrefix("CALLED BY: ");

void ErrorReporter::printWarning() {
	this->outputStream->print(WarningPrefix);
}

void ErrorReporter::printBold(const string &s) {
	this->outputStream->print(s);
}

void ErrorReporter::printContext() {
	this->outputStream->print(ContextPrefix);
}

void ErrorReporter::printError() {
	this->outputStream->print(ErrorPrefix);
}

void ErrorReporter::printCalledBy() {
	this->outputStream->print(CalledbyPrefix);
}

void ErrorReporter::printWhere() {
	this->outputStream->print(WherePrefix);
}

void DefaultErrorReporter::implReportWarning(const string &msg) {
	this->printWarning();
	this->outputStream->print(msg)->newline();
}

void DefaultErrorReporter::implReportWarning(const string &msg, const CharPointer &from, const CharPointer &to) {
	this->printWarning();
	this->outputStream->print(msg);
	this->reportLineNumberWhenMeaningful(from, to, true, true);
	this->outputStream->newline();
	this->outputUnderlinedChars(intrusive_ptr_cast<InterpreterNS::CppOSTREAM>(this->outputStream)->out, from, to);
}

  void DefaultErrorReporter::implReportInterrupt(const CharPointer &from, const CharPointer &to) {
    this->outputStream->print("\n-----------------------------\n>>>  CoCoA-5 interrupted  <<<\n-----------------------------")->newline();
	this->reportLineNumberWhenMeaningful(from, to, true, true);
	this->outputStream->newline();
	this->outputUnderlinedChars(intrusive_ptr_cast<InterpreterNS::CppOSTREAM>(this->outputStream)->out, from, to);
}

void DefaultErrorReporter::implReportError(const string &msg) {
	this->printError();
	this->outputStream->print(msg)->newline();
}

void DefaultErrorReporter::implReportError(const string &msg, const CharPointer &from, const CharPointer &to) {
	this->printError();
	this->outputStream->print(msg);
	this->reportLineNumberWhenMeaningful(from, to, true, true);
	this->outputStream->newline();
	this->outputUnderlinedChars(intrusive_ptr_cast<InterpreterNS::CppOSTREAM>(this->outputStream)->out, from, to);
}

} // namespace LexerNS
} // namespace CoCoA
