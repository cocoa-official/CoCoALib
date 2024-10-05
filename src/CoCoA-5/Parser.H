#ifndef CoCoA5_Parser_H
#define CoCoA5_Parser_H

//   Copyright (c) 2009-2019  John Abbott,  Anna M. Bigatti
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

#include <vector>
#include <set>
#include <algorithm>
#include "AST.H"

namespace CoCoA {
namespace ParserNS {

class Parser : public LexerNS::ReferenceCountedObject {
	ParserStatus status;
	std::string currentPackage;
	struct ScopeData {
		enum ScopeType { ST_TOPLEVEL, ST_PACKAGE, ST_DEFINE, ST_LAMBDA };
		ScopeType type;
		int nNestingLoops;
		bool thereIsEllipsisParameter;
		bool thereAreReturnsWithExpr;
		bool thereAreReturnsWithoutExpr;
		std::set<std::string> maybeLocals; // names used as assignment-targets or in ref parameter-passing
		std::set<std::string> maybeImplicitImports; // names used for calling fn-procs
		std::vector<std::string> iterationVars;
		std::vector<std::string> loopLabels;
		ScopeData(ScopeType type) :
			type(type),
			nNestingLoops(0),
			thereIsEllipsisParameter(false),
			thereAreReturnsWithExpr(false),
			thereAreReturnsWithoutExpr(false)
		{}
		bool isIterationVar(const std::string &identifier) {
			return std::find(this->iterationVars.begin(), this->iterationVars.end(), identifier)!=this->iterationVars.end();
		}
	};
	std::stack<ScopeData> scopeDataStack;
	ScopeData &currentScopeData() {
		assert(this->scopeDataStack.size()>=1);
		return this->scopeDataStack.top();
	}
	ScopeData::ScopeType currentScopeType() {
		return this->currentScopeData().type;
	}
	LexerNS::CharPointer lastRecoveredPosition;
	boost::intrusive_ptr<LexerNS::Lexer> lexer;
	LexerNS::Token expecting(LexerNS::TokenType);
	std::string normalizePackageName(const LexerNS::Token &t);
	boost::intrusive_ptr<AST::RingDefinition> parseRingDefinition();
	boost::intrusive_ptr<AST::IndeterminateDeclaration> parseIndeterminateDecl(bool isInsideElim);
	std::vector<boost::intrusive_ptr<AST::Expression> > parseExpressionList();
	std::vector<AST::Argument> parseArgumentList();
	boost::intrusive_ptr<AST::Expression> parseExpression();
	boost::intrusive_ptr<AST::Expression> parseConditionalOrExpression();
	boost::intrusive_ptr<AST::Expression> parseConditionalAndExpression();
	boost::intrusive_ptr<AST::Expression> parseEqualityExpression();
	boost::intrusive_ptr<AST::Expression> parseRelationalExpression();
	boost::intrusive_ptr<AST::Expression> parseCartesianProductExpression();
	boost::intrusive_ptr<AST::Expression> parseListExpression();
	boost::intrusive_ptr<AST::Expression> parseAdditiveExpression();
	boost::intrusive_ptr<AST::Expression> parseMultiplicativeExpression();
  boost::intrusive_ptr<AST::Expression> parsePowerExpression(const LexerNS::Token &primaryTok, bool IsExponent = false);
	boost::intrusive_ptr<AST::Expression> parseUnaryExpression();
	boost::intrusive_ptr<AST::Expression> parsePrimary(LexerNS::Token t);
	boost::intrusive_ptr<AST::RecordExpression> parseRecord(const LexerNS::Token &tokRecord);
	boost::intrusive_ptr<AST::Expression> parseSelectors(boost::intrusive_ptr<AST::Expression> targetExp);
	boost::intrusive_ptr<AST::Expression> parseListPrimary(const LexerNS::Token &tokOpenSquare);
	boost::intrusive_ptr<AST::UseStatement> parseUseStatement(const LexerNS::Token &tokUse);
	void parseSetStatement(const LexerNS::Token &tokSet);
	void parseUnsetStatement(const LexerNS::Token &tokUnset);
	boost::intrusive_ptr<AST::ForStatement> parseForStatement(const LexerNS::Token &tokFor, const std::string &label);
	boost::intrusive_ptr<AST::ForeachStatement> parseForeachStatement(const LexerNS::Token &tokForeach, const std::string &label);
	boost::intrusive_ptr<AST::IfStatement> parseIfStatement(const LexerNS::Token &tokIf);
	boost::intrusive_ptr<AST::ReturnStatement> parseReturnStatement(const LexerNS::Token &tokReturn);
	boost::intrusive_ptr<AST::UsingStatement> parseUsingStatement(const LexerNS::Token &tokUsing);
	boost::intrusive_ptr<AST::PrintStatement> parsePrintStatement(const LexerNS::Token &tokPrintOrPrintln);
	boost::intrusive_ptr<AST::Statement> parseBreakOrContinueStatement(const LexerNS::Token &tokBreakOrContinue);
	boost::intrusive_ptr<AST::SkipStatement> parseSkipStatement(const LexerNS::Token &tokSkip);
	boost::intrusive_ptr<AST::ProtectStatement> parseProtectStatement(const LexerNS::Token &tokProtect);
	boost::intrusive_ptr<AST::UnprotectStatement> parseUnprotectStatement(const LexerNS::Token &tokUnProtect);
	boost::intrusive_ptr<AST::CiaoOrQuitStatement> parseCiaoOrQuitStatement(const LexerNS::Token &tokCiaoOrQuit);
	boost::intrusive_ptr<AST::DescribeStatement> parseDescribeStatement(const LexerNS::Token &tokDescribe);
	boost::intrusive_ptr<AST::RepeatUntilStatement> parseRepeatUntilStatement(const LexerNS::Token &tokRepeat, const std::string &label);
	boost::intrusive_ptr<AST::WhileStatement> parseWhileStatement(const LexerNS::Token &tokWhile, const std::string &label);
	void parseClearStatement(const LexerNS::Token &tokClear);
	void parseDeleteStatement(const LexerNS::Token &tokDelete);
	void parseDestroyStatement(const LexerNS::Token &tokDestroy);
	boost::intrusive_ptr<AST::AliasStatement> parseAliasStatement(const LexerNS::Token &tokAlias);
	boost::intrusive_ptr<AST::BlockStatement> parseBlockStatement(const LexerNS::Token &tokBlock);
	void parseCatchStatement(const LexerNS::Token &tokCatch);
	boost::intrusive_ptr<AST::PackageStatement> parsePackageDefinition(const LexerNS::Token &tokPackage);
	boost::intrusive_ptr<AST::SourceStatement> parseSourceStatement(const LexerNS::Token &tokSource);
	boost::intrusive_ptr<AST::SourceRegionStatement> parseSourceRegionStatement(const LexerNS::Token &tokSourceRegion);
	boost::intrusive_ptr<AST::DefineStatement> parseDefineStatement(const LexerNS::Token &tokDefine);
	boost::intrusive_ptr<AST::TryStatement> parseTryStatement(const LexerNS::Token &tokTry);
	boost::intrusive_ptr<AST::LambdaExpression> parseLambdaExpression(const LexerNS::Token &tokLambda);
	boost::intrusive_ptr<AST::Statements> parseStatements(const LexerNS::Token &tokBeforeStatements);
	boost::intrusive_ptr<AST::Statement> parseFunBodyStatement();
	void parseImports(const std::string &fnName, const std::set<std::string> &paramNames, std::vector<AST::Import> &imports);
	boost::intrusive_ptr<AST::FunctionDeclaration> parseFunctionDeclaration(const std::string &fnName);
	void checkThereIsntAnotherUnaryPlusOrMinus(const LexerNS::Token &tokFirstOp);
	std::vector<AST::Param> parseParamList(int &nMandatoryParameters, std::set<std::string> &paramNames);
	boost::intrusive_ptr<AST::Expression> parseExpOrNewLine();
	LexerNS::Token getTokenWithVirtualSemicolon(const LexerNS::CharPointer &position);
	std::vector<AST::Binding> parseBindingList(LexerNS::CharPointer &position);
	std::vector<boost::intrusive_ptr<AST::Identifier> > parseIdList(LexerNS::CharPointer &position);
	LexerNS::Token expectingSemicolon(const LexerNS::CharPointer &position, const char * const errorMsg);
	LexerNS::Token expectingSemicolon(const LexerNS::CharPointer &position) { return expectingSemicolon(position, 0); }
	void enteringBlock(const std::string &blockName, bool isLoop, const std::string label = "");
	void exitingBlock(bool isLoop);
	void maybeLocal(const std::string &identifier);
	void reportSkippedInput(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	boost::intrusive_ptr<AST::Statement> tryParseFunBodyStatement();
	inline LexerNS::Token getToken() { return this->lexer->getToken(this->status); }
	inline void ungetToken(const LexerNS::Token &t) { this->lexer->ungetToken(t); }
	void reportError(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void reportWarning(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to, LexerNS::WarningSeverity severity);
	void reportWarningMayBeRemoved(const std::string &prefixMsg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void findLocalsImportsAndLambdas(boost::intrusive_ptr<AST::FunctionDeclaration> funDecl, const std::string &name);
	void checkLoopLabel(const LexerNS::Token &t);
public:
	explicit Parser(boost::intrusive_ptr<LexerNS::Lexer> lexer) :
		lastRecoveredPosition(LexerNS::CharPointer::Null),
		lexer(lexer)
	{
		assert(lexer);
		this->scopeDataStack.push(ScopeData(ScopeData::ST_TOPLEVEL));
	}
	inline const ParserStatus & getStatus() const { return status; }
	void tryToRecover(const LexerNS::CharPointer &beginPosition);
	boost::intrusive_ptr<AST::Statement> parseTopLevelStatement();
	static const std::string ARGV;
	static const std::string TopLevelPackageName;
};

class UnexpectedTokenException : public LexerNS::ExceptionWithSourcePosition {
	static std::string humanReadable(LexerNS::TokenType tt);
public:
	const LexerNS::Token found;
	UnexpectedTokenException(const LexerNS::TokenType expected, const LexerNS::Token & found) :
            // See redmine 673 for help understanding the next complicated line
            ExceptionWithSourcePosition("I was expecting "+humanReadable(expected)+" but I've found the "+found.LexemeClass()+" \""+found.lexeme()+"\"", found, true),
		found(found)
	{}
	UnexpectedTokenException(const std::string msg, const LexerNS::Token & found) :
		ExceptionWithSourcePosition(msg, found, true),
		found(found)
	{}
	~UnexpectedTokenException() throw () {}
};

} // namespace ParserNS
} // namespace CoCoA


// Close include guard:
#endif
