//   Copyright (c)  2009  John Abbott,  Anna M Bigatti
//   Origj author: 2009 Giovanni Lagorio (lagorio@disi.unige.it)
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
////#include <algorithm>
#include <exception>
#include <boost/scope_exit.hpp>
#include <boost/algorithm/string/case_conv.hpp>
// case_conv:  just 1 call to to_lower_copy

#include "Parser.H"

namespace CoCoA {
namespace ParserNS {

using namespace std;
using namespace boost;
using namespace CoCoA::AST;
using namespace CoCoA::LexerNS;

namespace {
	WarningSeverity WSEndInsteadOfErrBlabla = WS_NORMAL;
}

void Parser::reportError(const string &msg, const CharPointer &from, const CharPointer &to) {
	if (this->lexer->isInteractive())
		throw ExceptionWithSourcePosition(msg, from, to, false);
	this->lexer->reportError(msg, from, to);
}

void Parser::reportWarning(const string &msg, const CharPointer &from, const CharPointer &to, WarningSeverity severity) {
	this->lexer->reportWarning(msg, from, to, severity);
}

void Parser::reportWarningMayBeRemoved(const string &prefixMsg, const CharPointer &from, const CharPointer &to) {
	this->lexer->reportWarningMayBeRemoved(prefixMsg, from, to);
}

Token Parser::expecting(const TokenType expectedTT) {
	const Token token = this->getToken();
	const TokenType type = token.getType();
	if (type!=expectedTT) {
		this->ungetToken(token);
		throw UnexpectedTokenException(expectedTT, token);
	}
	return token;
}

std::string UnexpectedTokenException::humanReadable(const TokenType type) {
	switch (type) {
	case TT_ALIAS: return "\"Alias\"";
	case TT_AND: return "\"And\"";
	case TT_ASSIGN: return "an assignment operator (\":=\")";
	case TT_BLOCK: return "\"Block\"";
	case TT_BREAK: return "\"Break\"";
	case TT_CART_PROD: return "a \"><\" operator";
	case TT_CATCH: return "\"Catch\"";
	case TT_CIAO: return "\"Ciao\"";
	case TT_CLEAR: return "\"Clear\"";
	//case TT_CLOSE_COMPAT: return "the backward-compatibility close marker";
	case TT_CLOSED_ROUND: return "\")\"";
	case TT_CLOSED_SQUARE: return "\"]\"";
	case TT_COLONCOLON: return "a scope operator (\"::\")";
	case TT_COLON: return "colon";
	case TT_COMMA: return "a comma";
	case TT_COMPATIBILITY_MARKER: return "a compatibility marker";
	case TT_CONTINUE: return "\"Continue\"";
	case TT_DEFINE: return "\"Define\"";
	case TT_DEGLEX: return "\"DegLex\"";
	case TT_DEGREVLEX: return "\"DegRevLex\"";
	case TT_DELETE: return "\"Delete\"";
	case TT_DESCRIBE: return "\"Describe\"";
	case TT_DESTROY: return "\"Destroy\"";
	case TT_DO: return "\"Do\"";
	case TT_DOTDOT: return "a \"..\" operator";
	case TT_DOT: return "a dot";
	case TT_ELLIPSIS: return "an ellipsis (\"...\")";
	case TT_ELIF: return "\"Elif\"";
	case TT_ELSE: return "\"Else\"";
	case TT_END: return "\"End\"";
	case TT_ENDALIAS: return "\"EndAlias\"";
	case TT_ENDBLOCK: return "\"EndBlock\"";
	case TT_ENDCATCH: return "\"EndCatch\"";
	case TT_ENDDEFINE: return "\"EndDefine\"";
	case TT_ENDFOREACH: return "\"EndForeach\"";
	case TT_ENDFOR: return "\"EndFor\"";
	case TT_ENDIF: return "\"EndIf\"";
	case TT_ENDLAMBDA: return "\"EndFunc\"";
	case TT_ENDPACKAGE: return "\"EndPackage\"";
	case TT_ENDREPEAT: return "\"EndRepeat\"";
	case TT_ENDTRY: return "\"EndTry\"";
	case TT_ENDUSING: return "\"EndUsing\"";
	case TT_ENDWHILE: return "\"EndWhile\"";
	case TT_EOF: return "<EOF>";
	case TT_EQUAL: return "an \"=\" operator";
	case TT_EXIT: return "\"Exit\"";
	case TT_EXPORT: return "\"Export\"";
	case TT_FALSE: return "\"False\"";
	case TT_FLOAT_LITERAL: return "a floating-point literal";
	case TT_FOREACH: return "\"Foreach\"";
	case TT_FOR: return "\"For\"";
	case TT_GE: return "a \">=\" operator";
	case TT_GT: return "a \">\" operator";
	case TT_HELP: return "a help statement";
	case TT_IDENTIFIER: return "an identifier";
	case TT_IMPORTBYREF: return "\"ImportByRef\"";
	case TT_IMPORTBYVALUE: return "\"ImportByValue\"";
	case TT_ISDEFINED: return "\"IsDefined\"";
	case TT_IF: return "\"If\"";
	case TT_IN: return "\"In\"";
	case TT_INT_LITERAL: return "an integer literal";
	case TT_ISIN: return "\"IsIn\"";
	case TT_LAMBDA: return "\"Func\"";
	case TT_LE: return "a \"<=\" operator";
	case TT_LEX: return "\"Lex\"";
	case TT_LOAD: return "\"Load\"";
	case TT_LT: return "a \"<\" operator";
	case TT_MINUS: return "a \"-\" operator";
	case TT_MOD: return "a \"%\" operator";
	case TT_NOTEQUAL: return "a \"<>\" operator";
	case TT_ON: return "\"On\"";
	//case TT_OPEN_COMPAT: return "the backward-compatibility open marker";
	case TT_OPEN_ROUND: return "\"(\"";
	case TT_OPEN_SQUARE: return "\"[\"";
	case TT_OPT: return "\"Opt\"";
	case TT_OR: return "\"Or\"";
	case TT_PACKAGE: return "\"Package\"";
	case TT_PACKAGENAME: return "a package name (\"$PkgName\")";
	case TT_PIPE: return "a pipe (\"|\")";
	case TT_PLUS: return "a \"+\" operator";
	case TT_POSTO: return "\"PosTo\"";
	case TT_POWER: return "a \"^\" operator";
	case TT_PRINTLN: return "\"PrintLn\"";
	case TT_PRINT: return "\"Print\"";
	case TT_PROTECT: return "\"Protect\"";
	case TT_QUIT: return "\"Quit\"";
	case TT_RECORD: return "\"Record\"";
	case TT_REF: return "\"Ref\"";
	case TT_REPEAT: return "\"Repeat\"";
	case TT_RETURN: return "\"Return\"";
	case TT_RING_ASSIGN: return "a ring-assignment operator (\"::=\")";
	case TT_SEMICOLON: return "a semicolon";
	case TT_SET: return "\"Set\"";
	case TT_SKIP: return "\"Skip\"";
	case TT_SLASH: return "a \"/\" operator";
	case TT_SOURCE: return "\"Source\"";
	case TT_SOURCE_AS_LSHIFT: return "an obsolete source operator (\"<<\")";
	case TT_SOURCEREGION: return "\"SourceRegion\"";
	case TT_STEP: return "\"Step\"";
	case TT_STAR: return "a \"*\" operator";
	case TT_STRING_LITERAL: return "a string literal";
	case TT_THEN: return "\"Then\"";
	case TT_TIME: return "\"Time\"";
	case TT_TO: return "\"To\"";
	case TT_TOPLEVEL: return "\"TopLevel\"";
	case TT_TOPOS: return "\"ToPos\"";
	case TT_TRUE: return "\"True\"";
	case TT_TRY: return "\"Try\"";
	case TT_UNPROTECT: return "\"Unprotect\"";
	case TT_UNSET: return "\"Unset\"";
	case TT_UNTIL: return "\"Until\"";
	case TT_UPONERROR: return "\"UponError\"";
	case TT_USE: return "\"Use\"";
	case TT_USING: return "\"Using\"";
	case TT_VAR: return "\"Var\"";
    //	case TT_WEIGHTS: return "\"Weights\"";   // disabled 2015-03-09
	case TT_WHILE: return "\"While\"";
	case TT_XEL: return "\"Xel\"";
	}
	assert(false);
	throw FailedAssertionException("Unknown token type");
}

void Parser::reportSkippedInput(const CharPointer &from, const CharPointer &to) {
	if (from!=to)
		this->reportWarning("Skipping some input symbols while trying to continue the parsing after the error", from, to, WS_PEDANTIC);
	this->lastRecoveredPosition = to;
}

namespace {
	// the semicolon is rather special: it might start a statement, but only the empty one!
	// So, having the semicolon here would do more harm than good
	bool mayStartOrFollowAStatement(const TokenType tt) {
		switch (tt) {
		case TT_ALIAS:
		case TT_BLOCK:
		case TT_BREAK:
		case TT_CATCH:
		case TT_CIAO:
		case TT_CLEAR:
		case TT_DEFINE:
		case TT_DELETE:
		case TT_DESCRIBE:
		case TT_DESTROY:
		case TT_ELIF:
		case TT_ELSE:
		case TT_END:
		case TT_ENDALIAS:
		case TT_ENDBLOCK:
		case TT_ENDCATCH:
		case TT_ENDDEFINE:
		case TT_ENDFOR:
		case TT_ENDFOREACH:
		case TT_ENDIF:
		case TT_ENDPACKAGE:
		case TT_ENDREPEAT:
		case TT_ENDTRY:
		case TT_ENDUSING:
		case TT_ENDWHILE:
		case TT_EXIT:
		case TT_FOR:
		case TT_FOREACH:
		case TT_HELP:
		case TT_IF:
		case TT_LAMBDA:
		case TT_LOAD:
		case TT_PACKAGE:
		case TT_PRINT:
		case TT_PRINTLN:
		case TT_QUIT:
		case TT_REPEAT:
		case TT_RETURN:
		case TT_SET:
		case TT_SKIP:
		case TT_SOURCE:
		case TT_SOURCE_AS_LSHIFT:
		case TT_SOURCEREGION:
		case TT_TRY:
		case TT_UNSET:
		case TT_UPONERROR:
		case TT_USE:
		case TT_USING:
		case TT_WHILE:
			return true;
		default:
			return false;
		}
	}
}

namespace {
	string aSemicolonShouldBeInsertedHere("A semicolon should be inserted here");
}

Token Parser::getTokenWithVirtualSemicolon(const CharPointer &position) {
	const Token token = this->getToken();
	switch (token.getType()) {
	case TT_ELIF:
	case TT_ELSE:
	case TT_END:
	case TT_ENDALIAS:
	case TT_ENDBLOCK:
	case TT_ENDCATCH:
	case TT_ENDDEFINE:
	case TT_ENDFOR:
	case TT_ENDFOREACH:
	case TT_ENDIF:
	case TT_ENDLAMBDA:
	case TT_ENDPACKAGE:
	case TT_ENDREPEAT:
	case TT_ENDTRY:
	case TT_ENDUSING:
	case TT_ENDWHILE:
	case TT_IN:
	case TT_UPONERROR:
		this->ungetToken(token);
		this->reportWarning(aSemicolonShouldBeInsertedHere, position, token.getBegin(), WS_PEDANTIC);
		return Token(position, position, TT_SEMICOLON);
	default:
		return token;
	}
}

Token Parser::expectingSemicolon(const CharPointer &position, const char * const errorMsg) {
	const Token token = this->getToken();
	switch (token.getType()) {
	case TT_SEMICOLON:
		return token;
	case TT_ELIF:
	case TT_ELSE:
	case TT_END:
	case TT_ENDALIAS:
	case TT_ENDBLOCK:
	case TT_ENDCATCH:
	case TT_ENDDEFINE:
	case TT_ENDFOR:
	case TT_ENDFOREACH:
	case TT_ENDIF:
	case TT_ENDLAMBDA:
	case TT_ENDPACKAGE:
	case TT_ENDREPEAT:
	case TT_ENDTRY:
	case TT_ENDUSING:
	case TT_ENDWHILE:
	case TT_IN:
	case TT_UPONERROR:
		this->ungetToken(token);
		this->reportWarning(aSemicolonShouldBeInsertedHere, position, token.getBegin(), WS_PEDANTIC);
		return Token(position, position, TT_SEMICOLON);
	default:
		this->ungetToken(token);
		if (errorMsg)
			throw UnexpectedTokenException(errorMsg, token);
		throw UnexpectedTokenException(TT_SEMICOLON, token);
	}
}

void Parser::tryToRecover(const CharPointer &beginPosition) {
	this->lexer->resetCompatibilityMode();
	this->status.inRecoveryMode = true;
	BOOST_SCOPE_EXIT( (&status) ) {
		status.inRecoveryMode = false;
	} BOOST_SCOPE_EXIT_END
	CharPointer endPosition = beginPosition;
	for(;;) {
		try	{
			const Token token = this->getToken();
			switch (TokenType tt = token.getType()) {
			case TT_SEMICOLON:
				endPosition = token.getEnd();
				// fall through
			case TT_EOF:
				this->reportSkippedInput(beginPosition, endPosition);
				return;
			default:
				if (mayStartOrFollowAStatement(tt) && endPosition!=lastRecoveredPosition) {
					this->ungetToken(token);
					this->reportSkippedInput(beginPosition, endPosition);
					return;
				}
			}
			endPosition = token.getEnd();
		} catch (const AskingForNewInteractiveInputDuringRecoveryException &) {
			// we've consumed all (interactive) input after an error, we're good to go
			this->reportSkippedInput(beginPosition, endPosition);
			return;
		} catch (const LexerException &pe) {
			// we ignore all lexer errors (for instance, unbalanced compatibility-parentheses) until we're synchronized
		}
	}
}

intrusive_ptr<UseStatement> Parser::parseUseStatement(const Token &tokUse) {
	intrusive_ptr<Identifier> identifier;
	intrusive_ptr<RingDefinition> ringDefinition;
	const Token tokIdenfier = this->expecting(TT_IDENTIFIER);
	const Token lookahead = this->getToken();
	const TokenType laTT = lookahead.getType();
	if (laTT==TT_RING_ASSIGN || laTT==TT_ASSIGN)
        {
          if (laTT == TT_ASSIGN)
            this->reportWarning("Write ::= instead of := to avoid this warning", lookahead.getBegin(), lookahead.getEnd(), WS_LOW);
		identifier = new Identifier(tokIdenfier);
        }
	else {
		this->ungetToken(lookahead);
		this->ungetToken(tokIdenfier);
	}
	ringDefinition = this->parseRingDefinition();
	const Token tokSemicolon = this->expectingSemicolon(ringDefinition->getEnd());
	return new UseStatement(tokUse.getBegin(), identifier, ringDefinition, tokSemicolon.getEnd());
}

void Parser::parseSetStatement(const Token &tokSet) {
	intrusive_ptr<Identifier> identifier = new Identifier(this->expecting(TT_IDENTIFIER));
	intrusive_ptr<Expression> exp;
	const Token t = this->getToken();
	if (t.getType()==TT_ASSIGN)
		exp = this->parseExpression();
	else
		this->ungetToken(t);
	const Token tokSemicolon = this->expectingSemicolon(exp ? exp->getEnd() : identifier->getEnd());
	this->reportError("The statement Set has been removed", tokSet.getBegin(), tokSemicolon.getEnd());
}

void Parser::parseUnsetStatement(const Token &tokUnset) {
	intrusive_ptr<Identifier> identifier = new Identifier(this->expecting(TT_IDENTIFIER));
	const Token tokSemicolon = this->expectingSemicolon(identifier->getEnd());
	this->reportError("The statement Unset has been removed", tokUnset.getBegin(), tokSemicolon.getEnd());
}

intrusive_ptr<PrintStatement> Parser::parsePrintStatement(const Token &tokPrintOrPrintln) {
	vector<intrusive_ptr<Expression> > exps;
	CharPointer position = tokPrintOrPrintln.getEnd();
	Token t = this->getTokenWithVirtualSemicolon(position);
	TokenType tt = t.getType();
	if (tt!=TT_ON && tt!=TT_SEMICOLON) {
		this->ungetToken(t);
		for(;;) {
			intrusive_ptr<Expression> e = this->parseExpression();
			position = e->getEnd();
			exps.push_back(e);
			t = this->getTokenWithVirtualSemicolon(position);
			tt = t.getType();
			if (tt!=TT_COMMA)
				break;
		}
	}
	if (tt==TT_ON) {
		const intrusive_ptr<Expression> onExp = this->parseExpression();
		const Token tokSemicolon = this->expectingSemicolon(onExp->getEnd());
		return new PrintStatement(tokPrintOrPrintln, exps, onExp, tokSemicolon.getEnd());
	}
	if (tt!=TT_SEMICOLON)
		throw UnexpectedTokenException("Expecting an \"On\" or a semicolon to end the Print(Ln) statement", t);
	return new PrintStatement(tokPrintOrPrintln, exps, intrusive_ptr<Expression>(), t.getEnd());
}

intrusive_ptr<ForeachStatement> Parser::parseForeachStatement(const Token &tokForeach, const string &label) {
	const Token tokId = expecting(TT_IDENTIFIER);
	this->expecting(TT_IN);
	const intrusive_ptr<Expression> inExp = this->parseExpression();
	const Token tokDo = this->expecting(TT_DO);
	this->enteringBlock("foreach", true, label);
	this->currentScopeData().iterationVars.push_back(tokId.lexeme());
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.currentScopeData().iterationVars.pop_back();
		This.exitingBlock(true);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> statements = this->parseStatements(tokDo);
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndForeach", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDFOREACH) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDFOREACH, t);
	}
	return new ForeachStatement(tokForeach.getBegin(), tokId, inExp, statements, t.getEnd(), label);
}

intrusive_ptr<TryStatement> Parser::parseTryStatement(const Token &tokTry) {
	Parser &This = *this;
	this->enteringBlock("try", false);
	BOOST_SCOPE_EXIT( (&This) ) {
			This.exitingBlock(false);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> tryStatements = this->parseStatements(tokTry);
	this->exitingBlock(false);
	this->enteringBlock("UponError", false);
	expecting(TT_UPONERROR);
	const Token tokId = expecting(TT_IDENTIFIER);
	this->currentScopeData().iterationVars.push_back(tokId.lexeme());
	BOOST_SCOPE_EXIT( (&This) ) {
		This.currentScopeData().iterationVars.pop_back();
	} BOOST_SCOPE_EXIT_END
	const Token tokDo = this->expecting(TT_DO);
	intrusive_ptr<Statements> uponErrorStatements = this->parseStatements(tokDo);
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndTry", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDTRY) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDTRY, t);
	}
	return new TryStatement(tokTry.getBegin(), tryStatements, tokId, uponErrorStatements, t.getEnd());
}

intrusive_ptr<Statements> Parser::parseStatements(const Token &tokBeforeStatements) {
	vector<intrusive_ptr<Statement> > statements;
	Token t = tokBeforeStatements; // because we cannot default-construct it
	for(;;) {
		this->status.inTheMiddleOfAStatement = false; // plainly false, but it gives a better feedback
		t = this->getToken();
		const TokenType tt = t.getType();
		this->ungetToken(t);
		switch(tt) {
		case TT_END:
		case TT_ENDFOR:
		case TT_ENDFOREACH:
		case TT_ELIF:
		case TT_ELSE:
		case TT_ENDIF:
		case TT_ENDUSING:
		case TT_UNTIL:
		case TT_ENDWHILE:
		case TT_ENDALIAS:
		case TT_IN:
		case TT_ENDCATCH:
		case TT_ENDBLOCK:
		case TT_ENDPACKAGE:
		case TT_UPONERROR:
		case TT_ENDTRY:
		case TT_ENDDEFINE:
		case TT_ENDLAMBDA:
		case TT_ENDREPEAT:
			return new Statements(tokBeforeStatements, statements);
		default:
			intrusive_ptr<Statement> stmt = this->tryParseFunBodyStatement();
			if (stmt)
				statements.push_back(stmt);
		}
	}
}

intrusive_ptr<BlockStatement> Parser::parseBlockStatement(const Token &tokBlock) {
	this->enteringBlock("block", false);
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.exitingBlock(false);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> statements = this->parseStatements(tokBlock);
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndBlock", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDBLOCK) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDBLOCK, t);
	}
	return new BlockStatement(tokBlock.getBegin(), statements, t.getEnd());
}


intrusive_ptr<PackageStatement> Parser::parsePackageDefinition(const Token &tokPackage) {
	const Token tokPackageName = expecting(TT_PACKAGENAME);
	const string pkgName(tokPackageName.lexeme());
	if (pkgName=="$")
		this->reportError("Invalid package name", tokPackageName.getBegin(), tokPackageName.getEnd());
	if (pkgName==TopLevelPackageName)
		this->reportError("Reserved package name", tokPackageName.getBegin(), tokPackageName.getEnd());
	this->scopeDataStack.push(ScopeData(ScopeData::ST_PACKAGE));
	this->enteringBlock("package", false);
	Parser &This = *this;
	this->currentPackage = pkgName;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.scopeDataStack.pop();
		This.exitingBlock(false);
		This.currentPackage.clear();
	} BOOST_SCOPE_EXIT_END
	std::vector<boost::intrusive_ptr<Statement> > statements;
	set<string> aliases;
	set<string> definedNames;
	set<string> exports;
	vector<Token> exportTokens;
	bool okAliases = true;
	bool okExports = true;
        bool skipExports = false;
	for(;;) {
		const Token t = this->getToken();
		const TokenType tt = t.getType();
		switch (tt) {
		case TT_PACKAGE:
			throw UnexpectedTokenException("Package declarations cannot be nested", t);
		case TT_EXPORT:
                {
                  if (skipExports)
                    this->reportError("Cannot export a name after declaring `export skip'", t.getBegin(), t.getEnd());
                  if (!okExports)
                    this->reportError("All exports must be declared at the beginning of the package", t.getBegin(), t.getEnd());
                  const Token CheckForSkip = this->getToken();
                  const TokenType type = CheckForSkip.getType();
                  if (type == TT_SKIP)
                  {
                    if (!exportTokens.empty())
                      this->reportError("Cannot use `export skip' after exporting some names", t.getBegin(), t.getEnd());
                    skipExports = true;  // flag in case there is another "export" command
                    this->expectingSemicolon(CheckForSkip.getEnd());
                    break;
                  }
                  this->ungetToken(CheckForSkip); // oops, wasn't "skip", put it back
                  /*const*/ Token tokId = this->expecting(TT_IDENTIFIER);
                  const string identifier(tokId.lexeme());
                  if (!exports.insert(identifier).second)
                    this->reportError("Duplicated export", tokId.getBegin(), tokId.getEnd());
                  else
                    exportTokens.push_back(tokId);
                  this->expectingSemicolon(CheckForSkip.getEnd());
                  break;
		}
		case TT_ALIAS: {
			okExports = false;
			intrusive_ptr<AliasStatement> alias = this->parseAliasStatement(t);
			if (!okAliases)
				this->reportError("All aliases must be declared at the beginning of the package, after the export declarations (if any)", alias->getBegin(), alias->getEnd());
			else
				statements.push_back(alias);
			for (const Binding &b: alias->bindings)
				if (!aliases.insert(b.identifier->identifier).second)
					this->reportError("Duplicated alias", b.identifier->getBegin(), b.identifier->getEnd());
			break;
		}
		case TT_SEMICOLON:
			continue;
		case TT_DEFINE: {
			okAliases = okExports = false;
			intrusive_ptr<DefineStatement> define = this->parseDefineStatement(t);
			if (aliases.find(define->name)!=aliases.end())
				this->reportError("The function name conflicts with an alias", define->expId->getBegin(), define->expId->getEnd());
			if (!definedNames.insert(define->name).second)
				this->reportError("The function name has been already used in this package", define->expId->getBegin(), define->expId->getEnd());
			statements.push_back(define);
			break;
		}
		case TT_IDENTIFIER: {
			okAliases = okExports = false;
			this->expecting(TT_ASSIGN);
			intrusive_ptr<Expression> exp = this->parseExpression();
			const Token semicolon = this->expectingSemicolon(exp->getEnd());
			const intrusive_ptr<Identifier> identifier(new Identifier(t));
			if (aliases.find(identifier->identifier)!=aliases.end())
				this->reportError("The variable name conflicts with an alias", identifier->getBegin(), identifier->getEnd());
			else if (!definedNames.insert(identifier->identifier).second)
				this->reportError("The variable name has been already used in this package", identifier->getBegin(), identifier->getEnd());
			else
				statements.push_back(new AssignmentStatement(identifier, exp, semicolon.getBegin()));
			break;
		}
		case TT_END:
			this->reportWarning("Interpreting End as EndPackage", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
			// fall through
		case TT_ENDPACKAGE:
			if (skipExports==false && exportTokens.empty())
				this->reportWarning("No export declarations found inside the package declaration; no members will be exported to top-level", tokPackage.getBegin(), tokPackageName.getEnd(), WS_NORMAL);
			else {
                          for (const Token &t: exportTokens)
                          {
                            const string name(t.lexeme());
                            if (definedNames.find(name)==definedNames.end())
                              this->reportError("Cannot find a definition for the exported name \""+name+"\"", t.getBegin(), t.getEnd());
                          }
			}
			return new PackageStatement(tokPackage.getBegin(), tokPackageName, new Statements(tokPackageName, statements), t.getEnd(), exportTokens, definedNames);
                case TT_EOF:
                  	throw UnexpectedTokenException("End-of-input while looking for EndPackage", t);
		default:
			throw UnexpectedTokenException("Package members can be Alias, Define and simple assignments", t);
		}
	}
}

void Parser::parseCatchStatement(const Token &tokCatch) {
	this->enteringBlock("catch", false);
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.exitingBlock(false);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> statements = this->parseStatements(tokCatch);
	boost::intrusive_ptr<Expression> leftValue;
	Token t = this->getToken();
	TokenType tt = t.getType();
	if (tt==TT_IN) {
		leftValue = this->parseExpression();
		if (intrusive_ptr<Identifier> id = dynamic_pointer_cast<Identifier>(leftValue)) {
			if (this->currentScopeData().isIterationVar(id->identifier))
				this->reportError("This variable cannot be used to catch exceptions", id->getBegin(), id->getEnd());
			else
				this->maybeLocal(id->identifier);
		}
		if (!leftValue->isLeftValue())
			this->reportError("A left-value (an identifier, a field access or an indexed-access) is required here", leftValue->getBegin(), leftValue->getEnd());
		t = this->getToken();
		tt = t.getType();
	}
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndCatch", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDCATCH) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDCATCH, t);
	}
	this->reportError("The CoCoA-4 Catch construct is buggy and no more supported, please use the new Try construct instead", tokCatch.getBegin(), t.getEnd());
}

vector<Param> Parser::parseParamList(int &nMandatoryParameters, set<string> &paramNames) {
	static const string cannotFollowError("Mandatory parameters cannot follow optional ones");
	bool thereAreOptParameters = false;
	vector<Param> params;
	for(;;) {
		Param p;
		Token t = this->getToken();
		const TokenType tt = t.getType();
		if (tt==TT_CLOSED_ROUND) {
			this->ungetToken(t);
			return params;
		}
		if (tt==TT_VAR || tt==TT_REF) {
			if (tt==TT_VAR)
				this->reportWarning("Interpreting Var as Ref (use the latter one to avoid this warning)", t.getBegin(), t.getEnd(), WS_LOW);
			if (thereAreOptParameters) {
				this->reportError(cannotFollowError, t.getBegin(), t.getEnd());
				goto skip_pushback; // avoids to build inconsistent ASTs
			}
			p.byRef = true;
			t = this->expecting(TT_IDENTIFIER);
			p.name = t.lexeme();
			++nMandatoryParameters;
		} else if (tt==TT_OPT) {
			thereAreOptParameters = true;
			p.opt = true;
			t = this->expecting(TT_IDENTIFIER);
			p.name = t.lexeme();
		} else if (tt==TT_IDENTIFIER) {
			if (thereAreOptParameters) {
				this->reportError(cannotFollowError, t.getBegin(), t.getEnd());
				goto skip_pushback; // avoids to build inconsistent ASTs
			}
			p.name = t.lexeme();
			++nMandatoryParameters;
		} else
			throw UnexpectedTokenException("Expecting an identifier or a closed round parenthesis", t);
		if (!paramNames.insert(p.name).second)
			this->reportError("Duplicated parameter name", t.getBegin(), t.getEnd());
		params.push_back(p);
skip_pushback:
		const Token maybeComma = this->getToken();
		if (maybeComma.getType()!=TT_COMMA) {
			this->ungetToken(maybeComma);
			return params;
		}
	}
}

void Parser::findLocalsImportsAndLambdas(intrusive_ptr<FunctionDeclaration> funDecl, const string &name) {
	ScopeData &sd = this->currentScopeData();
	assert(sd.type==ScopeData::ST_DEFINE || sd.type==ScopeData::ST_LAMBDA);
	int nParams=0;
	if (funDecl->thereIsEllipsis)
		funDecl->localNames.insert(make_pair(Parser::ARGV, nParams++));
	else
          for (const Param &param: funDecl->params)
            funDecl->localNames.insert(make_pair(param.name, nParams++));
	if (name.length())
		sd.maybeImplicitImports.insert(name);
	set<string> explicitlyImported;
	for (const Import &import: funDecl->imports)
		explicitlyImported.insert(import.name);
	for (const string &s: sd.maybeLocals)
        {
          if (explicitlyImported.find(s)==explicitlyImported.end())
            funDecl->localNames.insert(make_pair(s, funDecl->localNames.size())); // note: this insert may fail, but it's ok
	}
	for (const string &ii: sd.maybeImplicitImports)
        {
          if (explicitlyImported.find(ii)!=explicitlyImported.end())
            continue;
          if (funDecl->localNames.find(ii)!=funDecl->localNames.end())
            continue;
          funDecl->imports.push_back(Import(intrusive_ptr<Identifier>(), ii, Import::IT_TOPLEVEL, true));
	}
}

intrusive_ptr<DefineStatement> Parser::parseDefineStatement(const Token &tokDefine) {
	const Token tokIdentifier = this->expecting(TT_IDENTIFIER);
	this->scopeDataStack.push(ScopeData(ScopeData::ST_DEFINE));
	this->enteringBlock("define", false);
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.exitingBlock(false);
		This.scopeDataStack.pop();
	} BOOST_SCOPE_EXIT_END
	const string fnName(tokIdentifier.lexeme());
	intrusive_ptr<FunctionDeclaration> funDecl(parseFunctionDeclaration(fnName));
	this->findLocalsImportsAndLambdas(funDecl, fnName);
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndDefine", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDDEFINE) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDDEFINE, t);
	}
	if (funDecl->thereIsEllipsis && fnName==ARGV)
		this->reportError("An fn-proc with a variable number of arguments cannot be called "+ARGV, tokIdentifier.getBegin(), tokIdentifier.getEnd());
	assert(this->scopeDataStack.size()>=1);
	assert(this->scopeDataStack.top().type==ScopeData::ST_DEFINE);
	return new DefineStatement(tokDefine.getBegin(), new Identifier(tokIdentifier), funDecl, t.getEnd());
}

intrusive_ptr<LambdaExpression> Parser::parseLambdaExpression(const LexerNS::Token &tokLambda) {
	this->scopeDataStack.push(ScopeData(ScopeData::ST_LAMBDA));
	this->enteringBlock("func", false);
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.exitingBlock(false);
		This.scopeDataStack.pop();
	} BOOST_SCOPE_EXIT_END
	Token peek = this->getToken();
	if (peek.getType()==TT_IDENTIFIER)
		this->reportError("Unexpected identifier; parameter(s) must be between brackets", tokLambda.getBegin(), peek.getEnd());
	else
		this->ungetToken(peek);
	intrusive_ptr<FunctionDeclaration> funDecl(parseFunctionDeclaration(string()));
	this->findLocalsImportsAndLambdas(funDecl, string());
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndFunc", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDLAMBDA) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDLAMBDA, t);
	}
	assert(this->scopeDataStack.size()>=1);
	assert(this->scopeDataStack.top().type==ScopeData::ST_LAMBDA);
	return new LambdaExpression(tokLambda.getBegin(), funDecl, t.getEnd());
}

intrusive_ptr<FunctionDeclaration> Parser::parseFunctionDeclaration(const string &fnName) {
	ScopeData &sd = this->currentScopeData();
	assert(sd.type==ScopeData::ST_DEFINE || sd.type==ScopeData::ST_LAMBDA);
	const Token tokOpenRound = this->expecting(TT_OPEN_ROUND);
	vector<Param> params;
	int nMandatoryParameters = 0;
	const Token maybeEllipsis = this->getToken();
	set<string> paramNames;
	if (maybeEllipsis.getType()==TT_ELLIPSIS) {
		sd.thereIsEllipsisParameter = true;
		paramNames.insert(ARGV);
	} else {
		this->ungetToken(maybeEllipsis);
		params = this->parseParamList(nMandatoryParameters, paramNames);
	}
	const Token tokClosedRound = this->expecting(TT_CLOSED_ROUND);
	vector<Import> imports;
	this->parseImports(fnName, paramNames, imports);
	intrusive_ptr<Statements> statements = this->parseStatements(tokClosedRound);
	return new FunctionDeclaration(tokOpenRound.getBegin(), sd.thereAreReturnsWithExpr, sd.thereIsEllipsisParameter, nMandatoryParameters, params, imports, statements);
}

void Parser::parseImports(const string &fnName, const set<string> &paramNames, vector<Import> &imports) {
	set<string> importedNames;
	int nImportByValue = 0;
	for(;;) {
		Token t = this->getToken();
		const TokenType tt = t.getType();
		Import::ImportType type;
		switch (tt) {
		case TT_TOPLEVEL:
			type = Import::IT_TOPLEVEL;
			break;
		case TT_IMPORTBYREF:
			type = Import::IT_BYREF;
			break;
		case TT_IMPORTBYVALUE:
			type = Import::IT_BYVALUE;
			break;
		default:
			this->ungetToken(t);
			return;
		}
		Token tokId(t);
		do {
			tokId = this->expecting(TT_IDENTIFIER);
			const string name = tokId.lexeme();
			if (name==fnName && type==Import::IT_BYVALUE)
				this->reportError("Importing by value the fn-proc name is not allowed", tokId.getBegin(), tokId.getEnd());
			if (paramNames.find(name)!=paramNames.end())
				this->reportError("The name conflicts with a parameter name", tokId.getBegin(), tokId.getEnd());
			if (!importedNames.insert(name).second)
				this->reportError("Duplicated import", tokId.getBegin(), tokId.getEnd());
			imports.push_back(Import(new Identifier(tokId), name, type, false));
			if (type==Import::IT_BYVALUE)
				imports.back().byValueIndex = nImportByValue++;
			t = this->getToken();
		} while (t.getType()==TT_COMMA);
		this->ungetToken(t);
		this->expectingSemicolon(tokId.getEnd());
	}
}

const string Parser::ARGV("ARGV");
const string Parser::TopLevelPackageName("$TopLevel");

intrusive_ptr<AliasStatement> Parser::parseAliasStatement(const Token &tokAlias) {
	const ScopeData::ScopeType scopeType = this->currentScopeType();
	CharPointer position = tokAlias.getEnd();
	vector<Binding> bindingList = this->parseBindingList(position);
///NO LONGER USED 20140616	bool thereIsIn = false;
	intrusive_ptr<Statements> statements;
	Token t = this->getToken();
	TokenType tt = t.getType();
	if (tt==TT_SEMICOLON) {
		if (scopeType!=ScopeData::ST_TOPLEVEL && scopeType!=ScopeData::ST_PACKAGE)
			this->reportError("Aliases can be declared at top- or package-level only; inside fn-procs you can use only Alias blocks", t.getBegin(), t.getEnd());
	} else {
		if (tt!=TT_IN)
			throw UnexpectedTokenException("Binding lists must be followed by \"In\" or a semicolon", t);
///NO LONGER USED 20140616		thereIsIn = true;
		this->enteringBlock("alias", false);
		Parser &This = *this;
		BOOST_SCOPE_EXIT( (&This) ) {
			This.exitingBlock(false);
		} BOOST_SCOPE_EXIT_END
		statements = this->parseStatements(t);
		t = this->getToken();
		tt = t.getType();
		if (tt==TT_END)
			this->reportWarning("Interpreting End as EndAlias", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
		else if (tt!=TT_ENDALIAS) {
			this->ungetToken(t);
			throw UnexpectedTokenException(TT_ENDALIAS, t);
		}
		this->reportError("The statement Alias...In...EndAlias, which was declared obsolescent in 4.7.5, is no longer supported", tokAlias.getBegin(), t.getEnd());
		return nullptr;
	}
	return new AliasStatement(tokAlias.getBegin(), bindingList, /* thereIsIn, statements, */ t.getEnd());
}

intrusive_ptr<WhileStatement> Parser::parseWhileStatement(const Token &tokWhile, const string &label) {
	intrusive_ptr<Expression> exp = parseExpression();
	const Token tokDo = expecting(TT_DO);
	this->enteringBlock("while", true, label);
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.exitingBlock(true);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> statements = this->parseStatements(tokDo);
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndWhile", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDWHILE) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDWHILE, t);
	}
	return new WhileStatement(tokWhile.getBegin(), exp, statements, t.getEnd(), label);
}

intrusive_ptr<RepeatUntilStatement> Parser::parseRepeatUntilStatement(const Token &tokRepeat, const string &label) {
	this->enteringBlock("repeat", true, label);
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.exitingBlock(true);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> statements = this->parseStatements(tokRepeat);
	const Token shouldBeUntilOrEndRepeat = this->getToken();
	const TokenType tt = shouldBeUntilOrEndRepeat.getType();
	intrusive_ptr<Expression> exp;
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndRepeat", shouldBeUntilOrEndRepeat.getBegin(), shouldBeUntilOrEndRepeat.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt==TT_UNTIL) {
		exp = this->parseExpression();
		const Token t = this->expectingSemicolon(exp->getEnd());
		return new RepeatUntilStatement(tokRepeat.getBegin(), statements, exp, t.getEnd(), label);
	} else if (tt!=TT_ENDREPEAT) {
		this->ungetToken(shouldBeUntilOrEndRepeat);
		throw UnexpectedTokenException("I was expecting Until or EndRepeat but I've found "+shouldBeUntilOrEndRepeat.lexeme(), shouldBeUntilOrEndRepeat);
	}
	return new RepeatUntilStatement(tokRepeat.getBegin(), statements, exp, shouldBeUntilOrEndRepeat.getEnd(), label);
}

void Parser::enteringBlock(const string &blockName, bool isLoop, const string label) {
	assert(isLoop || label.length()==0);
	if (isLoop) {
		ScopeData &sd = this->currentScopeData();
		++sd.nNestingLoops;
		sd.loopLabels.push_back(label);
	}
	vector<string> &blocks = this->status.openBlocks;
	if (blocks.empty())
		blocks.push_back(blockName);
	else {
		const string& current = blocks.back();
		blocks.push_back(current + " >> " + blockName);
	}
}

void Parser::exitingBlock(bool isLoop) {
	assert(this->status.openBlocks.size()>=1);
	if (isLoop) {
		ScopeData &sd = this->currentScopeData();
		--sd.nNestingLoops;
		assert(sd.loopLabels.size()>=1);
		sd.loopLabels.pop_back();
	}
	this->status.openBlocks.pop_back();
}

intrusive_ptr<ForStatement> Parser::parseForStatement(const Token &tokFor, const string &label) {
	const Token tokId = this->expecting(TT_IDENTIFIER);
	this->expecting(TT_ASSIGN);
	const intrusive_ptr<Expression> beginExp = this->parseExpression();
	this->expecting(TT_TO);
	const intrusive_ptr<Expression> endExp = this->parseExpression();
	intrusive_ptr<Expression> stepExp;
	const Token maybeTokStep = this->getToken();
	if (maybeTokStep.getType()==TT_STEP)
		stepExp = this->parseExpression();
	else
		this->ungetToken(maybeTokStep);
	const Token tokDo = this->expecting(TT_DO);
	this->enteringBlock("For", true, label);
	this->currentScopeData().iterationVars.push_back(tokId.lexeme());
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.currentScopeData().iterationVars.pop_back();
		This.exitingBlock(true);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> statements = this->parseStatements(tokDo);
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndFor", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDFOR) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDFOR, t);
	}
	return new ForStatement(tokFor.getBegin(), tokId, beginExp, endExp, stepExp, statements, t.getEnd(), label);
}

intrusive_ptr<UsingStatement> Parser::parseUsingStatement(const Token &tokUsing) {
	const intrusive_ptr<Identifier> identifier = new Identifier(this->expecting(TT_IDENTIFIER));
	const Token tokDo = this->expecting(TT_DO);
	this->enteringBlock("using", false);
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This) ) {
		This.exitingBlock(false);
	} BOOST_SCOPE_EXIT_END
	intrusive_ptr<Statements> statements = this->parseStatements(tokDo);
	const Token t = this->getToken();
	const TokenType tt = t.getType();
	if (tt==TT_END)
		this->reportWarning("Interpreting End as EndUsing", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
	else if (tt!=TT_ENDUSING) {
		this->ungetToken(t);
		throw UnexpectedTokenException(TT_ENDUSING, t);
	}
	return new UsingStatement(tokUsing.getBegin(), identifier, statements, t.getEnd());
}

intrusive_ptr<ReturnStatement> Parser::parseReturnStatement(const Token &tokReturn) {
	Token t = this->getTokenWithVirtualSemicolon(tokReturn.getEnd());
	intrusive_ptr<Expression> exp;
	ScopeData &scopeData = this->currentScopeData();
	const ScopeData::ScopeType scopeType = scopeData.type;
	if (scopeType!=ScopeData::ST_LAMBDA && scopeType!=ScopeData::ST_DEFINE)
		throw ExceptionWithSourcePosition("Return statements can be used only inside functions", tokReturn.getBegin(), t.getEnd(), false);
	if (t.getType()!=TT_SEMICOLON) {
		scopeData.thereAreReturnsWithExpr = true;
		this->ungetToken(t);
		exp = this->parseExpression();
		t = this->expectingSemicolon(exp ? exp->getEnd():tokReturn.getEnd());
	} else
		scopeData.thereAreReturnsWithoutExpr = true;
	if (scopeData.thereAreReturnsWithExpr && scopeData.thereAreReturnsWithoutExpr)
		this->reportError("Inside a function definition all Return statements must be either with or without an expression", tokReturn.getBegin(), t.getEnd());
	return new ReturnStatement(tokReturn.getBegin(), exp, t.getEnd());
}

intrusive_ptr<SourceStatement> Parser::parseSourceStatement(const Token &tokSource) {
	if (this->currentScopeType()!=ScopeData::ST_TOPLEVEL)
		this->reportError("Files can be Source-d at top-level only", tokSource.getBegin(), tokSource.getEnd());
	const intrusive_ptr<Expression> expFileName = this->parseExpression();
        const Token tokSemicolon = this->expectingSemicolon(expFileName->getEnd());
        return new SourceStatement(tokSource, expFileName, tokSemicolon.getEnd());
}

intrusive_ptr<SourceRegionStatement> Parser::parseSourceRegionStatement(const Token &tokSourceRegion) {
	if (this->currentScopeType()!=ScopeData::ST_TOPLEVEL)
		this->reportError("Files can be SourceRegion-d at top-level only", tokSourceRegion.getBegin(), tokSourceRegion.getEnd());
	const intrusive_ptr<Expression> expFromLine = this->parseExpression();
        const Token comma1 = this->expecting(TT_COMMA);
	const intrusive_ptr<Expression>	expFromChar = this->parseExpression();

	const Token tokTo = this->expecting(TT_TO);
        const intrusive_ptr<Expression> expToLine = this->parseExpression();
        const Token comma2 = this->expecting(TT_COMMA);
        const intrusive_ptr<Expression> expToChar = this->parseExpression();
	const Token tokIn = this->expecting(TT_IN);
        const intrusive_ptr<Expression> expFileName = this->parseExpression();
        const Token tokSemicolon = this->expectingSemicolon(expFileName->getEnd());

	return new SourceRegionStatement(tokSourceRegion, expFromLine,expFromChar, expToLine,expToChar, expFileName, tokSemicolon.getEnd());
}

void Parser::parseClearStatement(const Token &tokClear) {
	CharPointer position = tokClear.getEnd();
	const vector<intrusive_ptr<Identifier> > list = this->parseIdList(position);
	const Token tokSemicolon = this->expectingSemicolon(position);
	this->reportError("The Clear statement, which has been declared obsolescent in 4.7.5, is no more supported", tokClear.getBegin(), tokSemicolon.getEnd());
}

void Parser::parseDeleteStatement(const Token &tokDelete) {
	CharPointer position = tokDelete.getEnd();
	const vector<intrusive_ptr<Identifier> > list = this->parseIdList(position);
	const Token tokSemicolon = this->expectingSemicolon(position);
	this->reportError("The Delete statement, which has been declared obsolescent in 4.7.5, is no more supported", tokDelete.getBegin(), tokSemicolon.getEnd());
}

void Parser::parseDestroyStatement(const Token &tokDestroy) {
	CharPointer position = tokDestroy.getEnd();
	const vector<intrusive_ptr<Identifier> > list = this->parseIdList(position);
	const Token tokSemicolon = this->expectingSemicolon(position);
	this->reportError("The Destroy statement, which has been declared obsolescent in 4.7.5, is no more supported", tokDestroy.getBegin(), tokSemicolon.getEnd());
}

intrusive_ptr<DescribeStatement> Parser::parseDescribeStatement(const Token &tokDescribe) {
	intrusive_ptr<Expression> exp;
	const Token t = this->getToken();
	if (t.getType()==TT_PACKAGENAME) {
		Token maybeSemicolon = this->getTokenWithVirtualSemicolon(t.getEnd());
		this->ungetToken(maybeSemicolon);
		if (maybeSemicolon.getType()==TT_SEMICOLON)
			exp = new Identifier(t);
		else
			goto notPackageName;
	} else {
notPackageName:
		this->ungetToken(t);
		exp = this->parseExpression();
	}
	const Token tokSemicolon = this->expectingSemicolon(exp->getEnd());
	return new DescribeStatement(tokDescribe.getBegin(), exp, tokSemicolon.getEnd());
}

void Parser::checkLoopLabel(const Token &t) {
	const string label(t.lexeme());
	ScopeData &sd = this->currentScopeData();
	if (find(sd.loopLabels.begin(), sd.loopLabels.end(), label)!=sd.loopLabels.end())
		this->reportError("Duplicate loop label \""+label+"\"", t.getBegin(), t.getEnd());
}

intrusive_ptr<Statement> Parser::parseBreakOrContinueStatement(const LexerNS::Token &tokBreakOrContinue) {
	string label;
	const Token t = this->getToken();
	if (t.getType()==TT_IDENTIFIER) {
		label = t.lexeme();
		ScopeData &sd = this->currentScopeData();
		if (find(sd.loopLabels.begin(), sd.loopLabels.end(), label)==sd.loopLabels.end())
			this->reportError("There is no outer loop labelled \""+label+"\"", t.getBegin(), t.getEnd());
	} else
		this->ungetToken(t);
	const Token tokSemicolon = this->expectingSemicolon(tokBreakOrContinue.getEnd());
	assert(this->scopeDataStack.size()>=1);
	if (!this->scopeDataStack.top().nNestingLoops)
		this->reportError("`break' and `continue' statements can be used only inside loops", tokBreakOrContinue.getBegin(), tokSemicolon.getEnd());
	if (tokBreakOrContinue.getType()==TT_BREAK)
		return new BreakStatement(tokBreakOrContinue.getBegin(), tokSemicolon.getEnd(), label);
	assert(tokBreakOrContinue.getType()==TT_CONTINUE);
	return new ContinueStatement(tokBreakOrContinue.getBegin(), tokSemicolon.getEnd(), label);
}

intrusive_ptr<SkipStatement> Parser::parseSkipStatement(const Token &tokSkip) {
	const Token tokSemicolon = this->expectingSemicolon(tokSkip.getEnd());
	return new SkipStatement(tokSkip.getBegin(), tokSemicolon.getEnd());
}

intrusive_ptr<CiaoOrQuitStatement> Parser::parseCiaoOrQuitStatement(const Token &tokCiaoOrQuit) {
	const Token tokSemicolon = this->expectingSemicolon(tokCiaoOrQuit.getEnd());
	return new CiaoOrQuitStatement(tokCiaoOrQuit, tokSemicolon.getEnd());
}

intrusive_ptr<IfStatement> Parser::parseIfStatement(const Token &tokIf) {
	vector<IfBranch> branches;
	bool foundElse = false;
	bool firstTime = true;
	Parser &This = *this;
	BOOST_SCOPE_EXIT( (&This)(&firstTime) ) {
		if (!firstTime)
			This.exitingBlock(false);
	} BOOST_SCOPE_EXIT_END
	Token t = tokIf;
	TokenType tt = TT_IF; // anything != TT_END, TT_ENDIF (to enter the loop) and TT_ELSE (to read the expression) would work
	while (tt!=TT_END && tt!=TT_ENDIF) {
		IfBranch branch;
		if (tt!=TT_ELSE) {
			branch.optExp = parseExpression();
			this->expecting(TT_THEN);
			if (firstTime) {
				firstTime = false;
				this->enteringBlock("if", false);
			}
		}
		branch.statements = this->parseStatements(t);
		branches.push_back(branch);
		t = this->getToken();
		tt = t.getType();
		if (tt==TT_END) {
			this->reportWarning("Interpreting End as EndIf", t.getBegin(), t.getEnd(), WSEndInsteadOfErrBlabla);
			break;
		}
		if (tt==TT_ENDIF)
			break;
		if (tt==TT_ELIF) {
			if (foundElse) {
				branches.pop_back();
				this->reportError("The Else branch must be the last one, I was expecting Endif", t.getBegin(), t.getEnd());
			}
			continue;
		}
		if (tt==TT_ELSE) {
			 if (foundElse) {
				 branches.pop_back();
				 this->reportError("An If statement cannot have more than one Else branch, I was expecting Endif", t.getBegin(), t.getEnd());
			 }
			 foundElse = true;
		} else
			throw UnexpectedTokenException(TT_ENDIF, t);
	}
	return new IfStatement(tokIf.getBegin(), branches, t.getEnd());
}

intrusive_ptr<Statement> Parser::tryParseFunBodyStatement() {
	if (this->lexer->isInteractive())
		return this->parseFunBodyStatement();
	try {
		return this->parseFunBodyStatement();
	} catch (const UnexpectedTokenException &ete) {
		const Token &foundToken = ete.found;
		if (foundToken.getType()==TT_EOF)
			throw;
		this->reportError(ete.reason, foundToken.getBegin(), foundToken.getEnd());
		this->tryToRecover(foundToken.getBegin());
		return intrusive_ptr<Statement>();
	} catch (const ExceptionWithSourcePosition &ewsp) {
		this->reportError(ewsp.reason, ewsp.from, ewsp.to);
		if (ewsp.needsRecovery)
			this->tryToRecover(ewsp.from);
		return intrusive_ptr<Statement>();
	}
}

intrusive_ptr<ProtectStatement> Parser::parseProtectStatement(const Token &tokProtect) {
	const Token tokIdentifier = this->expecting(TT_IDENTIFIER);
	const Token maybeColon = this->getToken();
	intrusive_ptr<Expression> optExp;
	if (maybeColon.getType()!=TT_COLON)
		this->ungetToken(maybeColon);
	else
		optExp = this->parseExpression();
	const Token tokSemicolon = this->expectingSemicolon(tokIdentifier.getEnd());
	return new ProtectStatement(tokProtect.getBegin(), new Identifier(tokIdentifier), optExp, tokSemicolon.getEnd());
}

intrusive_ptr<UnprotectStatement> Parser::parseUnprotectStatement(const Token &tokUnProtect) {
	const Token tokIdentifier = this->expecting(TT_IDENTIFIER);
	const Token tokSemicolon = this->expectingSemicolon(tokIdentifier.getEnd());
	return new UnprotectStatement(tokUnProtect.getBegin(), new Identifier(tokIdentifier), tokSemicolon.getEnd());
}

intrusive_ptr<Statement> Parser::parseTopLevelStatement() {
	this->status.inTheMiddleOfAStatement = false;
	Token t = this->getToken();
	this->status.inTheMiddleOfAStatement = true;
	BOOST_SCOPE_EXIT( (&status) ) {
		status.inTheMiddleOfAStatement = false;
	} BOOST_SCOPE_EXIT_END
	switch (t.getType()) {
	case TT_USE:
		return this->parseUseStatement(t);
	case TT_PACKAGE:
		return this->parsePackageDefinition(t);
	case TT_SOURCE_AS_LSHIFT:
          // issue warning that "<<" is obsolescent ???
	case TT_SOURCE:
	case TT_LOAD:
		return this->parseSourceStatement(t);
	case TT_SOURCEREGION:
		return this->parseSourceRegionStatement(t);
	case TT_ALIAS:
		return this->parseAliasStatement(t);
	case TT_DEFINE:
		return this->parseDefineStatement(t);
	default:
		this->ungetToken(t);
		return this->parseFunBodyStatement();
	}
}

intrusive_ptr<Statement> Parser::parseFunBodyStatement() {
	const Token t = this->getToken();
	switch (t.getType()) {
	case TT_SEMICOLON:
		//this->reportWarning("Empty statement", t.getBegin(), t.getEnd(), WS_PEDANTIC);
		return new EmptyStatement(t.getBegin());
	case TT_IF:
		return this->parseIfStatement(t);
	case TT_RETURN:
		return this->parseReturnStatement(t);
	case TT_USING:
		return this->parseUsingStatement(t);
	case TT_PRINT:
	case TT_PRINTLN:
		return this->parsePrintStatement(t);
	case TT_BREAK:
	case TT_CONTINUE:
		return this->parseBreakOrContinueStatement(t);
	case TT_DESCRIBE:
		return this->parseDescribeStatement(t);
	case TT_CIAO:
	case TT_EXIT:
	case TT_QUIT:
		return this->parseCiaoOrQuitStatement(t);
	case TT_SKIP:
		return this->parseSkipStatement(t);
	case TT_BLOCK:
		return this->parseBlockStatement(t);
	case TT_HELP:
		return new HelpStatement(t);
	case TT_TRY:
		return this->parseTryStatement(t);
	case TT_PROTECT:
		return this->parseProtectStatement(t);
	case TT_UNPROTECT:
		return this->parseUnprotectStatement(t);
	case TT_IMPORTBYREF:
	case TT_IMPORTBYVALUE:
	case TT_TOPLEVEL:
		throw UnexpectedTokenException("Imports (TopLevel, ImportByRef and ImportByValue) must be placed at the beginning of function bodies", t);
	case TT_SET:
		this->parseSetStatement(t);
		return nullptr;
	case TT_UNSET:
		this->parseUnsetStatement(t);
		return nullptr;
	case TT_CLEAR:
		this->parseClearStatement(t);
		return nullptr;
	case TT_DELETE:
		this->parseDeleteStatement(t);
		return nullptr;
	case TT_DESTROY:
		this->parseDestroyStatement(t);
		return nullptr;
	case TT_CATCH:
		this->parseCatchStatement(t);
		return nullptr;
	case TT_TIME:
		{
			intrusive_ptr<Statement> stmt = this->parseFunBodyStatement();
			if (!stmt)
				return stmt;
			if (intrusive_ptr<EvalStatement> evalStmt = dynamic_pointer_cast<EvalStatement>(stmt))
				this->reportWarning("This use of Time statement is discouraged; please use an assignment to avoid this warning", evalStmt->getBegin(), evalStmt->getEnd(), WS_LOW);
			return new TimeStatement(stmt, t.getBegin());
		}
	case TT_IDENTIFIER: {
		const Token maybeColon = this->getToken();
		if (maybeColon.getType()==TT_COLON) {
			const string label(t.lexeme());
			const Token maybeLoopKeyword = this->getToken();
			switch (maybeLoopKeyword.getType()) {
			case TT_FOR:
				this->checkLoopLabel(t);
				return this->parseForStatement(maybeLoopKeyword, label);
			case TT_FOREACH:
				this->checkLoopLabel(t);
				return this->parseForeachStatement(maybeLoopKeyword, label);
			case TT_WHILE:
				this->checkLoopLabel(t);
				return this->parseWhileStatement(maybeLoopKeyword, label);
			case TT_REPEAT:
				this->checkLoopLabel(t);
				return this->parseRepeatUntilStatement(t, label);
			default:
				this->ungetToken(maybeLoopKeyword);
			}
		}
		this->ungetToken(maybeColon);
		break;
	}
	case TT_SOURCE:  // Anna 2014-07-24 allow source in top-level loops
    return this->parseSourceStatement(t);
	case TT_SOURCEREGION:  // Anna 2014-07-24 allow source in top-level loops
    return this->parseSourceRegionStatement(t);
	case TT_FOR:
		return this->parseForStatement(t, "");
	case TT_FOREACH:
		return this->parseForeachStatement(t ,"");
	case TT_WHILE:
		return this->parseWhileStatement(t, "");
	case TT_REPEAT:
		return this->parseRepeatUntilStatement(t, "");
	case TT_DEFINE:
		this->reportError("Named fn-procs can be declared only at top-level or inside a package", t.getBegin(), t.getEnd());
		this->parseDefineStatement(t);
		return nullptr;
	case TT_USE:
		this->reportError("Statement Use can only be used at top-level", t.getBegin(), t.getEnd());
		this->parseUseStatement(t);
		return nullptr;
	case TT_PACKAGE:
		this->reportError("Packages can be declared only at top-level", t.getBegin(), t.getEnd());
		this->parsePackageDefinition(t);
		return nullptr;
	case TT_ALIAS:
		this->reportError("Aliases can only be declared at top or package level", t.getBegin(), t.getEnd());
		this->parseAliasStatement(t);
		return nullptr;
	default:
		break; // apparently useless, but it shuts the compiler (blabla not handled in switch)
	}
	this->ungetToken(t);
	intrusive_ptr<Expression> exp = this->parseExpression();
	const Token t2 = this->getTokenWithVirtualSemicolon(exp->getEnd());
	const TokenType tt2 = t2.getType();
	if (intrusive_ptr<Identifier> id = dynamic_pointer_cast<Identifier>(exp)) {
		if (this->currentScopeData().isIterationVar(id->identifier))
			this->reportError("This variable cannot be assigned", id->getBegin(), id->getEnd());
		else
			this->maybeLocal(id->identifier);
	}
	if (tt2==TT_RING_ASSIGN) {
		intrusive_ptr<RingDefinition> rdef = this->parseRingDefinition();
		Token tokSemicolon = this->expectingSemicolon(rdef->getEnd(), "Expecting a semicolon after a ring assignment");
		if (!exp->isLeftValue())
			this->reportError("Invalid left-value for ring-assignment", exp->getBegin(), exp->getEnd());
		return new RingAssignStatement(exp, rdef, tokSemicolon.getEnd());
	}
	if (tt2==TT_ASSIGN) {
		intrusive_ptr<Expression> rightExp = this->parseExpression();
		Token tokSemicolon = this->expectingSemicolon(rightExp->getEnd(), "Expecting a semicolon after an assignment");
		if (!exp->isLeftValue())
			this->reportError("Invalid left-value for assignment", exp->getBegin(), exp->getEnd());
		return new AssignmentStatement(exp, rightExp, tokSemicolon.getEnd());
	}
	if (tt2!=TT_SEMICOLON)
		throw UnexpectedTokenException("Expecting an operator, or a semicolon to end the statement", t2);
	const ScopeData::ScopeType scopeType = this->currentScopeType();
	if (scopeType==ScopeData::ST_DEFINE || scopeType==ScopeData::ST_LAMBDA) {
		if (!dynamic_pointer_cast<InvocationExpression>(exp))
			this->reportError("Expression-statements inside fn-procs must invoke another fn-proc", exp->getBegin(), t2.getEnd());
	}
	return new EvalStatement(exp, t2.getEnd());
}

void Parser::maybeLocal(const string &identifier) {
	ScopeData &sd = this->currentScopeData();
	const bool insideFun = sd.type==ScopeData::ST_DEFINE || sd.type==ScopeData::ST_LAMBDA;
	if (insideFun)
		sd.maybeLocals.insert(identifier);
}

intrusive_ptr<RingDefinition> Parser::parseRingDefinition() {
	const Token tokRingIdentifier = this->expecting(TT_IDENTIFIER);
	Token lastToken = tokRingIdentifier;
	intrusive_ptr<Expression> optionalExp;
	const Token maybeSlash = this->getToken();
	if (maybeSlash.getType()!=TT_SLASH)
		this->ungetToken(maybeSlash);
	else {
		this->expecting(TT_OPEN_ROUND);
		optionalExp = this->parseExpression();
		lastToken = this->expecting(TT_CLOSED_ROUND);
	}
	vector<intrusive_ptr<IndeterminateDeclaration> > indeterminates;
	const Token maybeOpenSquare = this->getToken();
	if (maybeOpenSquare.getType()!=TT_OPEN_SQUARE)
		this->ungetToken(maybeOpenSquare);
	else {
		for(;;) {
			indeterminates.push_back(parseIndeterminateDecl(false));
			const Token t = this->getToken();
			if (t.getType()!=TT_COMMA) {
				this->ungetToken(t);
				break;
			}
		}
		lastToken = this->expecting(TT_CLOSED_SQUARE);
		bool isFirst=true;
		for(;;) {
			static const char * const ignoredOrders = "We currently support only Lex, DegLex, DegRevLex, and Xel; others are recognized but *ignored*";
			Token maybeComma = this->getToken();
			if (maybeComma.getType()!=TT_COMMA) {
				this->ungetToken(maybeComma);
				break;
			} else {
				lastToken = this->getToken();
				if (isFirst)
					isFirst = false;
				else
					this->reportWarning("This modifier will be ignored (only the first one is currently considered)", lastToken.getBegin(), lastToken.getEnd(), WS_LOW);
				switch (lastToken.getType()) {
				case TT_LEX:
				case TT_DEGREVLEX:
				case TT_DEGLEX:
				case TT_XEL:
					break;
				case TT_POSTO:
				case TT_TOPOS:
					this->reportWarning(ignoredOrders, lastToken.getBegin(), lastToken.getEnd(), WS_LOW);
					break;
// 				case TT_WEIGHTS:    // disabled 2015-03-09
// 					this->reportWarning(ignoredOrders, lastToken.getBegin(), lastToken.getEnd(), WS_LOW);
// 					this->expecting(TT_OPEN_ROUND);
// 					this->parseExpressionList();
// 					lastToken = this->expecting(TT_CLOSED_ROUND);
// 					break;
				case TT_IDENTIFIER: {
					const string maybeFormerKeyword(to_lower_copy(lastToken.lexeme()));
					if (maybeFormerKeyword=="elim" || maybeFormerKeyword=="ord") { /* horrible backward-compatibility hack (needed since Elim and Ord are no longer keywords) */
						this->reportWarning(ignoredOrders, lastToken.getBegin(), lastToken.getEnd(), WS_LOW);
						this->expecting(TT_OPEN_ROUND);
						this->parseIndeterminateDecl(true);
						Token maybeDotdot = this->getToken();
						if (maybeDotdot.getType()==TT_DOTDOT)
							this->parseIndeterminateDecl(true);
						else
							this->ungetToken(maybeDotdot);
						lastToken = this->expecting(TT_CLOSED_ROUND);
						break;
					}
					/* no break here */
				}
                                  // fall through
				default:
					throw UnexpectedTokenException("Expecting PosTo, ToPos, Lex, Xel, DegLex, DegRevLex, Ord or Elim", lastToken);
				}
			}
		}
	}
	return new RingDefinition(new Identifier(tokRingIdentifier), optionalExp, indeterminates, lastToken, lastToken.getEnd());
}

intrusive_ptr<IndeterminateDeclaration> Parser::parseIndeterminateDecl(bool isInsideElim) {
	const intrusive_ptr<Identifier> identifier = new Identifier(this->expecting(TT_IDENTIFIER));
	CharPointer endPosition = identifier->getEnd();
	vector<pair <intrusive_ptr<Expression>, intrusive_ptr<Expression> > > ranges;
	Token t = this->getToken();
	if (t.getType()!=TT_OPEN_SQUARE)
		this->ungetToken(t);
	else {
		for(;;) {
			intrusive_ptr<Expression> e1 = this->parseAdditiveExpression();
			intrusive_ptr<Expression> e2;
			if (!isInsideElim) {
				const Token maybeDotDot = this->getToken();
				if (maybeDotDot.getType()==TT_DOTDOT)
					e2 = this->parseAdditiveExpression();
				else
					this->ungetToken(maybeDotDot);
			}
			ranges.push_back(make_pair(e1, e2));
			t = this->getToken();
			if (t.getType()!=TT_COMMA) {
				this->ungetToken(t);
				break;
			}
		}
		Token closedSquare = this->expecting(TT_CLOSED_SQUARE);
		endPosition = closedSquare.getEnd();
	}
	return new IndeterminateDeclaration(identifier, ranges, endPosition);
}

vector<intrusive_ptr<Expression> > Parser::parseExpressionList() {
	vector<intrusive_ptr<Expression> > list;
	for(;;) {
		const intrusive_ptr<Expression> exp = this->parseExpression();
		list.push_back(exp);
		const Token t = this->getToken();
		if (t.getType()!=TT_COMMA) {
			this->ungetToken(t);
			break;
		}
	}
	return list;
}

vector<Argument> Parser::parseArgumentList() {
	vector<Argument> list;
	for(;;) {
		bool byRef = false;
		const Token maybeRef = this->getToken();
		if (maybeRef.getType()==TT_REF)
			byRef = true;
		else
			this->ungetToken(maybeRef);
		const intrusive_ptr<Expression> exp = this->parseExpression();
		if (byRef)
			if (intrusive_ptr<const Identifier> id = dynamic_pointer_cast<const Identifier>(exp)) {
				if (this->currentScopeData().isIterationVar(id->identifier))
					this->reportError("This variable cannot be passed by-reference", id->getBegin(), id->getEnd());
				else
					this->maybeLocal(id->identifier);
			}
		list.push_back(Argument(byRef, exp));
		const Token t = this->getToken();
		if (t.getType()!=TT_COMMA) {
			this->ungetToken(t);
			break;
		}
	}
	return list;
}

vector<intrusive_ptr<Identifier> > Parser::parseIdList(CharPointer &position) {
	vector<intrusive_ptr<Identifier> > list;
	for(;;) {
		const intrusive_ptr<Identifier> identifier = new Identifier(this->expecting(TT_IDENTIFIER));
		position = identifier->getEnd();
		list.push_back(identifier);
		const Token t = this->getToken();
		if (t.getType()!=TT_COMMA) {
			this->ungetToken(t);
			break;
		}
	}
	return list;
}

vector<Binding> Parser::parseBindingList(CharPointer &position) {
	vector<Binding> list;
	set<string> aliases;
	for(;;) {
		Binding b;
		b.identifier = new Identifier(this->expecting(TT_IDENTIFIER));
		this->expecting(TT_ASSIGN);
		Token pn = this->expecting(TT_PACKAGENAME);
		position = pn.getEnd();
		b.packageName = this->normalizePackageName(pn);
		if (b.packageName==TopLevelPackageName)
			this->reportError(TopLevelPackageName+" it's not a real package, use TopLevel declarations to import top-level names", pn.getBegin(), position);
		if (!aliases.insert(b.identifier->identifier).second)
			this->reportError("Duplicated alias", b.identifier->getBegin(), b.identifier->getEnd());
		else
			list.push_back(b);
		const Token t = this->getToken();
		if (t.getType()!=TT_COMMA) {
			this->ungetToken(t);
			break;
		}
	}
	return list;
}

intrusive_ptr<Expression> Parser::parseExpression() {
	return this->parseConditionalOrExpression();
}

intrusive_ptr<Expression> Parser::parseConditionalOrExpression() {
	intrusive_ptr<Expression> returnValue = this->parseConditionalAndExpression();
	for(;;) {
		const Token t = this->getToken();
		if (t.getType()!=TT_OR) {
			this->ungetToken(t);
			break;
		}
		intrusive_ptr<Expression> exp = this->parseExpression();
		returnValue = new OrExpression(returnValue, t.getBegin(), t.getEnd(), exp);
	}
	return returnValue;
}

intrusive_ptr<Expression> Parser::parseConditionalAndExpression() {
	intrusive_ptr<Expression> returnValue = this->parseEqualityExpression();
	for(;;) {
		const Token t = this->getToken();
		if (t.getType()!=TT_AND) {
			this->ungetToken(t);
			break;
		}
		const intrusive_ptr<Expression> exp = this->parseEqualityExpression();
		returnValue = new AndExpression(returnValue, t.getBegin(), t.getEnd(), exp);
	}
	return returnValue;
}

intrusive_ptr<Expression> Parser::parseEqualityExpression() {
	intrusive_ptr<Expression> returnValue = this->parseRelationalExpression();
	for(;;) {
		const Token t = this->getToken();
		const TokenType type = t.getType();
		if (type!=TT_EQUAL && type!=TT_NOTEQUAL && type!=TT_ISIN) {
			this->ungetToken(t);
			break;
		}
		const intrusive_ptr<Expression> exp = this->parseRelationalExpression();
		switch (type) {
		case TT_EQUAL:
			returnValue = new EqualExpression(returnValue, t.getBegin(), t.getEnd(), exp);
			break;
		case TT_NOTEQUAL:
			returnValue = new NotEqualExpression(returnValue, t.getBegin(), t.getEnd(), exp);
			break;
		case TT_ISIN:
			returnValue = new IsInExpression(returnValue, t.getBegin(), t.getEnd(), exp);
			break;
		default:
			assert(false);
		}
	}
	return returnValue;
}

intrusive_ptr<Expression> Parser::parseRelationalExpression() {
	intrusive_ptr<Expression> returnValue = this->parseCartesianProductExpression();
	for(;;) {
		const Token t = this->getToken();
		const TokenType type = t.getType();
		if (type!=TT_LE && type!=TT_LT && type!=TT_GE && type!=TT_GT) {
			this->ungetToken(t);
			break;
		}
		const intrusive_ptr<Expression> exp = this->parseCartesianProductExpression();
		switch (type) {
		case TT_LE:
			returnValue = new LessOrEqualExpression(returnValue, t.getBegin(), t.getEnd(), exp);
			break;
		case TT_LT:
			returnValue = new LessThanExpression(returnValue, t.getBegin(), t.getEnd(), exp);
			break;
		case TT_GE:
			returnValue = new GreaterOrEqualExpression(returnValue, t.getBegin(), t.getEnd(), exp);
			break;
		case TT_GT:
			returnValue = new GreaterThanExpression(returnValue, t.getBegin(), t.getEnd(), exp);
			break;
		default:
			assert(false);
		}
	}
	return returnValue;
}

intrusive_ptr<Expression> Parser::parseCartesianProductExpression() {
	intrusive_ptr<Expression> exp = this->parseListExpression();
	const Token t = this->getToken();
	if (t.getType()!=TT_CART_PROD) {
		this->ungetToken(t);
		return exp;
	}
	vector<intrusive_ptr<Expression> > operands;
	operands.push_back(exp);
	for(;;) {
		exp = this->parseListExpression();
		operands.push_back(exp);
		const Token t = this->getToken();
		if (t.getType()!=TT_CART_PROD) {
			this->ungetToken(t);
			break;
		}
	}
	return new CartesianProductExpression(operands);
}

intrusive_ptr<Expression> Parser::parseListExpression() {
	const intrusive_ptr<Expression> exp = this->parseAdditiveExpression();
	const Token t = this->getToken();
	if (t.getType()!=TT_DOTDOT) {
		this->ungetToken(t);
		return exp;
	}
	return new DotDotExpression(exp, t.getBegin(), t.getEnd(), this->parseAdditiveExpression());
}

intrusive_ptr<Expression> Parser::parseListPrimary(const Token &tokOpenSquare) {
	vector<intrusive_ptr<Expression> > exps;
	Token t = this->getToken();
	if (t.getType()==TT_CLOSED_SQUARE)
		return new ListExpression(tokOpenSquare.getBegin(), t.getEnd(), exps);
	if (t.getType()==TT_IDENTIFIER) {
		const Token tokId = t;
		t = this->getToken();
		if (t.getType()==TT_IN) {
			const intrusive_ptr<Expression> exp1 = this->parseExpression();
			expecting(TT_PIPE);
			const intrusive_ptr<Expression> exp2 = this->parseExpression();
			const Token tokClosedSquare = this->expecting(TT_CLOSED_SQUARE);
			return new IdInExpSuchThatExp(tokOpenSquare.getBegin(), tokClosedSquare.getEnd(), tokId, exp1, exp2);
		}
		this->ungetToken(t);
		this->ungetToken(tokId);
	} else
		this->ungetToken(t);
	const intrusive_ptr<Expression> exp1 = this->parseExpression();
	t = this->getToken();
	if (t.getType()==TT_COMMA || t.getType()==TT_CLOSED_SQUARE) {
		exps.push_back(exp1);
		while (t.getType()==TT_COMMA) {
			exps.push_back(this->parseExpression());
			t = this->getToken();
		}
		this->ungetToken(t);
		const Token tokClosedSquare = this->expecting(TT_CLOSED_SQUARE);
		return new ListExpression(tokOpenSquare.getBegin(), tokClosedSquare.getEnd(), exps);
	} else
		this->ungetToken(t);
	this->expecting(TT_PIPE);
	const Token tokId = this->expecting(TT_IDENTIFIER);
	this->expecting(TT_IN);
	const intrusive_ptr<Expression> exp2 = this->parseEqualityExpression();
	intrusive_ptr<Expression> exp3;
	t = this->getToken();
	if (t.getType()==TT_AND) {
		exp3 = this->parseExpression();
	} else
		this->ungetToken(t);
	const Token tokClosedSquare = this->expecting(TT_CLOSED_SQUARE);
	return new ExpSuchThatIdInExpAndExp(tokOpenSquare.getBegin(), tokClosedSquare.getEnd(), exp1, tokId, exp2, exp3);
}

intrusive_ptr<Expression> Parser::parseAdditiveExpression() {
	const intrusive_ptr<Expression> firstExp = this->parseMultiplicativeExpression();
	vector<intrusive_ptr<Expression> > exps;
	vector<Token> ops;
	exps.push_back(firstExp);
	for(;;) {
		const Token t = this->getToken();
		const TokenType type = t.getType();
		if (type!=TT_PLUS && type!=TT_MINUS) {
			this->ungetToken(t);
			break;
		}
		ops.push_back(t);
		exps.push_back(this->parseMultiplicativeExpression());
	}
	const size_t size = exps.size();
	if (size==1)
		return firstExp;
	if (size==2) {
		const Token t = ops.front();
		switch (t.getType()) {
		case TT_PLUS:
			return new SumExpression(exps[0], t.getBegin(), t.getEnd(), exps[1]);
		case TT_MINUS:
			return new SubtractionExpression(exps[0], t.getBegin(), t.getEnd(), exps[1]);
		default:
			assert(false);
		}
	}
	return new SummationExpression(exps, ops);
}

intrusive_ptr<Expression> Parser::parseMultiplicativeExpression() {
	static intrusive_ptr<Line> fakeLineForStar(new Line(0, "*", intrusive_ptr<LineProvider>()));
	static CharPointer cpStar(fakeLineForStar, 0);
	const intrusive_ptr<Expression> firstExp = this->parseUnaryExpression();
	bool usedStarOrModAfterSlashOrColon=false;
	int nSlashesOrColons=0;
	vector<intrusive_ptr<Expression> > exps;
	vector<Token> ops;
	exps.push_back(firstExp);
	for(;;) {
		const Token t = this->getToken();
		const TokenType type = t.getType();
		if (type==TT_IDENTIFIER) {
			if (!this->lexer->getStatus().isInCoCoA4CompatibilityMode())
				throw UnexpectedTokenException("Unexpected identifier; are you forgetting a \"*\" or a \";\" ?", t);
			this->ungetToken(t);
			ops.push_back(Token(cpStar, cpStar, TT_STAR));
			exps.push_back(this->parseUnaryExpression());
			usedStarOrModAfterSlashOrColon = nSlashesOrColons>0;
			continue;
		}
		switch(type) {
		case TT_STAR:
			usedStarOrModAfterSlashOrColon = nSlashesOrColons>0;
			break;
		case TT_SLASH:
			++nSlashesOrColons;
			break;
		case TT_MOD:
			usedStarOrModAfterSlashOrColon = nSlashesOrColons>0;
			break;
		case TT_COLON:
			++nSlashesOrColons;
			break;
		default:
			this->ungetToken(t);
			goto longBreak; // I'm neither ashamed nor afraid of using an "evil" goto in situations like this one
		}
		ops.push_back(t);
        const Token AfterOp = this->getToken();
        const bool SignAfterOp = (AfterOp.getType() == TT_PLUS || AfterOp.getType() == TT_MINUS);
        this->ungetToken(AfterOp);
		exps.push_back(this->parseUnaryExpression());
        if (SignAfterOp)
          this->reportWarning("Using + or - right after * or / or : is ugly; please use brackets", exps.back()->getBegin(), exps.back()->getEnd(), WS_LOW);

	}
longBreak:
	const size_t size = exps.size();
	if (size==1)
		return firstExp;
	if (nSlashesOrColons>1)
		this->reportWarning("Using two or more operators \"/\" or \":\" at the same level is potentially ambiguous; please parenthesize", firstExp->getBegin(), exps.back()->getEnd(), WS_LOW);
	if ( usedStarOrModAfterSlashOrColon )
		this->reportWarning("Using operators \"*\" or \"%\" after \"/\" or \":\" is potentially ambiguous; please parenthesize", firstExp->getBegin(), exps.back()->getEnd(), WS_LOW);
	if (size==2) {
		const Token t = ops.front();
		switch (t.getType()) {
		case TT_STAR:
			return new ProductExpression(exps[0], t.getBegin(), t.getEnd(), exps[1]);
		case TT_SLASH:
			return new DivisionExpression(exps[0], t.getBegin(), t.getEnd(), exps[1]);
		case TT_MOD:
			return new ModuloExpression(exps[0], t.getBegin(), t.getEnd(), exps[1]);
		case TT_COLON:
			return new ColonExpression(exps[0], t.getBegin(), t.getEnd(), exps[1]);
		default:
			assert(false);
		}
	}
	return new MultiplicationExpression(exps, ops);
}

  // 2021-02-17 somewhat hacked for redmine 1540
  intrusive_ptr<Expression> Parser::parsePowerExpression(const Token &primaryTok, bool IsExponent)
  {
	const intrusive_ptr<Expression> base = this->parseSelectors(this->parsePrimary(primaryTok));
	const Token t = this->getToken();
	if (t.getType()!=TT_POWER)
        {
		this->ungetToken(t);
		return base;
	}
        const Token AfterHat = this->getToken();
        const bool SignAfterHat = (AfterHat.getType() == TT_PLUS || AfterHat.getType() == TT_MINUS);
        this->ungetToken(AfterHat);
        if (SignAfterHat)
          throw UnexpectedTokenException("Exponent may not start with plus or minus; use brackets", AfterHat);

	const intrusive_ptr<Expression> exponent = this->parsePowerExpression(this->getToken(), true); // tell parser this is an exponent
        if (IsExponent)
          this->reportWarning("\"^\" is right associative; use brackets to avoid this warning", base->getBegin(), exponent->getEnd(), WS_LOW);

	return new PowerExpression(base, t.getBegin(), t.getEnd(), exponent);
}

intrusive_ptr<Expression> Parser::parseSelectors(intrusive_ptr<Expression> targetExp) {
	intrusive_ptr<Expression> returnValue = targetExp;
	for(;;) {
		const Token t = this->getToken();
		const TokenType type = t.getType();
		if (type!=TT_DOT && type!=TT_OPEN_SQUARE && type!=TT_OPEN_ROUND) {
			this->ungetToken(t);
			return returnValue;
		}
		switch(type) {
		case TT_DOT:
			returnValue = new FieldAccessExpression(returnValue, this->expecting(TT_IDENTIFIER));
			break;
		case TT_OPEN_SQUARE:
			{
				const vector<intrusive_ptr<Expression> > eList = this->parseExpressionList();
				const Token tokClosedBracket = this->expecting(TT_CLOSED_SQUARE);
				returnValue = new IndexedAccessExpression(returnValue, eList, tokClosedBracket.getEnd());
			}
			break;
		case TT_OPEN_ROUND:
			{
				ScopeData &sd = this->currentScopeData();
				const bool insideFun = sd.type==ScopeData::ST_DEFINE || sd.type==ScopeData::ST_LAMBDA;
				if (insideFun)
					if (intrusive_ptr<Identifier> id = dynamic_pointer_cast<Identifier>(returnValue))
						sd.maybeImplicitImports.insert(id->identifier);
				vector<Argument> args;
				const Token maybeEllipsisOrClosedRound = this->getToken();
				if (maybeEllipsisOrClosedRound.getType()==TT_ELLIPSIS) {
					const Token tokClosedRound = this->expecting(TT_CLOSED_ROUND);
					if (!this->currentScopeData().thereIsEllipsisParameter)
						this->reportError("Ellipsis can be used only inside fn-procs with a variable number of arguments", maybeEllipsisOrClosedRound.getBegin(), maybeEllipsisOrClosedRound.getEnd());
					returnValue = new InvocationExpression(returnValue, boost::shared_ptr<Token>(new Token(t)), args, tokClosedRound.getEnd(), this->currentPackage);
				} else {
					this->ungetToken(maybeEllipsisOrClosedRound);
					if (maybeEllipsisOrClosedRound.getType()!=TT_CLOSED_ROUND)
						args = this->parseArgumentList();
					const Token tokClosedRound = this->expecting(TT_CLOSED_ROUND);
					if (intrusive_ptr<Identifier> id = dynamic_pointer_cast<Identifier>(returnValue))
						id->arity = args.size();
					returnValue = new InvocationExpression(returnValue, boost::shared_ptr<Token>(), args, tokClosedRound.getEnd(), this->currentPackage);
				}
			}
			break;
		default:
			assert(false);
			throw FailedAssertionException("Parser::parseSelectors unknown token type");
		}
	}
}

intrusive_ptr<RecordExpression> Parser::parseRecord(const Token &tokRecord) {
	vector<RecordField> fields;
	Token t = this->expecting(TT_OPEN_SQUARE);
	t = this->getToken();
	if (t.getType()!=TT_CLOSED_SQUARE) {
		this->ungetToken(t);
		for(;;) {
			RecordField f;
			f.name = new Identifier(this->expecting(TT_IDENTIFIER));
			t = this->getToken();
			TokenType tt = t.getType();
			if (tt==TT_EQUAL)
				this->reportWarning("You should use \":=\" instead of \"=\" inside records", t.getBegin(), t.getEnd(), WS_LOW);
			else if (tt!=TT_ASSIGN)
				throw UnexpectedTokenException("Expected \":=\" before the field initialization expression", t);
			f.initExp = this->parseExpression();
			fields.push_back(f);
			t = this->getToken();
			tt = t.getType();
			if (tt==TT_CLOSED_SQUARE)
				break;
			if (tt!=TT_COMMA)
                        {
                          this->ungetToken(t);
                          throw UnexpectedTokenException("Expected a comma or a close-square-bracket to end the record definition", t);
                        }
		}
	}
	return new RecordExpression(tokRecord.getBegin(), fields, t.getEnd());
}

namespace {
	void decodeStringLiteral(const Token &token, string &unescapedString) {
		const string escapedString = token.getBegin().stringTo(token.getEnd());
		const string::size_type escapedLen = escapedString.length();
		bool multiline = escapedLen>3 && escapedString[0]=='\"' && escapedString[1]=='\"' && escapedString[2]=='\"';
		const int lenMinusClosingQuote = escapedLen-(multiline ? 3 : 1);
		assert(lenMinusClosingQuote>=(multiline ? 3 : 1));
		assert(escapedString[0]=='"' || escapedString[0]=='\'');
		for(int a=multiline ? 3 : 1; a<lenMinusClosingQuote; ++a) {
			char c = escapedString[a];
			if (c!='\\')
				unescapedString.append(1, c);
			else {
				c = escapedString[++a];
				assert(a<lenMinusClosingQuote);
				switch (c) {
				case 'n':
					unescapedString.append(1, '\n');
					break;
				case 'r':
					unescapedString.append(1, '\r');
					break;
				case 't':
					unescapedString.append(1, '\t');
					break;
				case 'a':
					unescapedString.append(1, '\a');
					break;
				case '\'':
					unescapedString.append(1, '\'');
					break;
				case '\"':
					unescapedString.append(1, '"');
					break;
				case '\\':
					unescapedString.append(1, '\\');
					break;
				case 'x': {
					char c = escapedString[++a];
					int d1 = isdigit(c) ? c-'0' : (tolower(c)-'a'+10);
					c = escapedString[++a];
					int d2 = isdigit(c) ? c-'0' : (tolower(c)-'a'+10);
					unescapedString.append(1, static_cast<char>((d1<<4)|d2));
				}
				default:
					;// nop (the lexer should have already reported the problem)
				}
			}
		}
	}
}

string Parser::normalizePackageName(const Token &t) {
	string packageName(t.lexeme());
	if (packageName=="$") {
		if (!this->currentPackage.length())
			this->reportError("Cannot use the special package-identifier \"$\" outside a package", t.getBegin(), t.getEnd());
		else
			packageName = this->currentPackage;
	}
	return packageName;
}

intrusive_ptr<Expression> Parser::parsePrimary(Token t) {
	Token tokPackagename = Token::EndOfFile;
	bool thereIsPackagename = false;
	string packageName;
	switch (const TokenType tt = t.getType()) {
	case TT_LAMBDA:
		return this->parseLambdaExpression(t);
	//case TT_OPEN_COMPAT:
	case TT_COMPATIBILITY_MARKER:
		{
			const intrusive_ptr<Expression> exp = this->parseExpression();
			//const Token t2 = this->expecting(TT_CLOSE_COMPAT);
			const Token t2 = this->expecting(TT_COMPATIBILITY_MARKER);
			return new BackwardCompatibleExp(exp, t.getBegin(), t2.getEnd());
		}
	case TT_RECORD:
		return this->parseRecord(t);
	case TT_OPEN_SQUARE:
		return this->parseListPrimary(t);
	case TT_OPEN_ROUND:
		{
			const intrusive_ptr<Expression> exp = this->parseExpression();
			const Token tokClosedRound = this->expecting(TT_CLOSED_ROUND);
			return exp;
		}
	case TT_ISDEFINED:
		{
			const Token tokOpenRound = this->expecting(TT_OPEN_ROUND);
			const Token tokId = this->expecting(TT_IDENTIFIER);
			const Token tokClosedRound = this->expecting(TT_CLOSED_ROUND);
			return new IsDefinedExp(new Identifier(tokId), tokOpenRound.getBegin(), tokClosedRound.getEnd());
		}
	case TT_PACKAGENAME:
		tokPackagename = t;
		packageName = normalizePackageName(tokPackagename);
		thereIsPackagename = true;
		{
			const Token shouldBeDot = this->getToken();
			if (shouldBeDot.getType()!=TT_DOT) {
				this->ungetToken(shouldBeDot);
				throw UnexpectedTokenException("A package name must be followed by a member access (that is, a dot followed by the member name)", shouldBeDot);
			}
		}
		t = this->expecting(TT_IDENTIFIER);
		if (packageName==TopLevelPackageName)
			this->reportError(TopLevelPackageName+" is not a real package, use \"TopLevel\" declarations to import top-level names", tokPackagename.getBegin(), tokPackagename.getEnd());
		// fall through
	case TT_IDENTIFIER:
		{
			intrusive_ptr<Expression> identifier(thereIsPackagename ?
					static_cast<Expression *>(new FullyQualifiedIdentifier(tokPackagename, packageName, t)) :
					static_cast<Expression *>(new Identifier(t)));
			const Token maybeColonColon = this->getToken();
			if (maybeColonColon.getType()!=TT_COLONCOLON) {
				this->ungetToken(maybeColonColon);
				return identifier;
			}
			Token nextToken = this->getToken();
			return new ScopedExpression(identifier, maybeColonColon.getBegin(), maybeColonColon.getEnd(), parsePrimary(nextToken));
		}
	case TT_INT_LITERAL:
		return new IntLiteral(t);
	case TT_FLOAT_LITERAL:
		return new FloatLiteral(t);
	case TT_TRUE:
	case TT_FALSE:
		return new BoolLiteral(t.getBegin(), t.getEnd(), tt);
	case TT_STRING_LITERAL: {
		string unescapedString;
		decodeStringLiteral(t, unescapedString);
		CharPointer begin = t.getBegin();
		CharPointer end = t.getEnd();
                // JAA 20140626 commented out code to allow juxtaposed string literals
                // (seems quite unimportant for CoCoA-5)
		// for(;;) {
		// 	t = this->getToken();
		// 	if (t.getType()!=TT_STRING_LITERAL) {
		// 		this->ungetToken(t);
		// 		break;
		// 	}
		// 	decodeStringLiteral(t, unescapedString);
		// 	end = t.getEnd();
		// }
		return new StringLiteral(begin, end, unescapedString);
	}
	case TT_EOF:
		throw UnexpectedTokenException("Found End-of-Input while looking for the start of an expression", t);
	default:
		this->ungetToken(t);
		throw UnexpectedTokenException("Invalid start of expression", t);
	}
}

void Parser::checkThereIsntAnotherUnaryPlusOrMinus(const Token &tokFirstOp) {
	const Token t = this->getToken();
	this->ungetToken(t);
	const TokenType tt = t.getType();
	if (tt==TT_MINUS && tokFirstOp.getType()==TT_MINUS)
          throw UnexpectedTokenException("Using two unary minus operators in a row is not allowed because it looks suspiciously like as a mistyped single-line comment; please parenthesize the inner expression or remove the blank(s) to make your intentions clear", t);
	else if (tt==TT_PLUS || tt==TT_MINUS)
          throw UnexpectedTokenException("Using two unary operators in a row is not allowed; please check your expression and parenthesize its inner sub-expression if this is really what you want", t);
}

intrusive_ptr<Expression> Parser::parseUnaryExpression() {
	const Token t = this->getToken();
	switch (t.getType()) {
	case TT_PLUS:
		this->checkThereIsntAnotherUnaryPlusOrMinus(t);
		return new UnaryPlusExpression(t.getBegin(), parseUnaryExpression());
	case TT_MINUS:
		this->checkThereIsntAnotherUnaryPlusOrMinus(t);
		return new UnaryMinusExpression(t.getBegin(), parseUnaryExpression());
	default:
		return this->parsePowerExpression(t);
	}
}

} // namespace ParserNS
} // namespace CoCoA
