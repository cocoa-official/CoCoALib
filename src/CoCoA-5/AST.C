//   Copyright (c) 2009  John Abbott,  Anna M Bigatti
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

#include <cctype>
#include <algorithm>
#include <iterator>
//#include <boost/foreach.hpp>
//#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "Interpreter.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"

#include "C5Utils.H"

namespace CoCoA {
namespace AST {

using namespace std;
using namespace boost;
using namespace CoCoA::LexerNS;

#define IMPLEMENT_ACCEPT_METHOD(CLASS) void CLASS::accept(ParsedObjectVisitor *v) { v->visit(*this); }

IMPLEMENT_ACCEPT_METHOD(BinaryExpression)
IMPLEMENT_ACCEPT_METHOD(UnaryExpression)
IMPLEMENT_ACCEPT_METHOD(IntLiteral)
IMPLEMENT_ACCEPT_METHOD(FloatLiteral)
IMPLEMENT_ACCEPT_METHOD(StringLiteral)
IMPLEMENT_ACCEPT_METHOD(BoolLiteral)
IMPLEMENT_ACCEPT_METHOD(IdInExpSuchThatExp)
IMPLEMENT_ACCEPT_METHOD(ExpSuchThatIdInExpAndExp)
IMPLEMENT_ACCEPT_METHOD(ListExpression)
IMPLEMENT_ACCEPT_METHOD(InvocationExpression)
IMPLEMENT_ACCEPT_METHOD(UnaryMinusExpression)
IMPLEMENT_ACCEPT_METHOD(UnaryPlusExpression)
IMPLEMENT_ACCEPT_METHOD(IndexedAccessExpression)
IMPLEMENT_ACCEPT_METHOD(CartesianProductExpression)
IMPLEMENT_ACCEPT_METHOD(SummationExpression)
IMPLEMENT_ACCEPT_METHOD(MultiplicationExpression)
IMPLEMENT_ACCEPT_METHOD(OrExpression)
IMPLEMENT_ACCEPT_METHOD(PowerExpression)
IMPLEMENT_ACCEPT_METHOD(DotDotExpression)
IMPLEMENT_ACCEPT_METHOD(SumExpression)
IMPLEMENT_ACCEPT_METHOD(SubtractionExpression)
IMPLEMENT_ACCEPT_METHOD(ProductExpression)
IMPLEMENT_ACCEPT_METHOD(DivisionExpression)
IMPLEMENT_ACCEPT_METHOD(ColonExpression)
IMPLEMENT_ACCEPT_METHOD(ModuloExpression)
IMPLEMENT_ACCEPT_METHOD(AndExpression)
IMPLEMENT_ACCEPT_METHOD(EqualExpression)
IMPLEMENT_ACCEPT_METHOD(NotEqualExpression)
IMPLEMENT_ACCEPT_METHOD(IsInExpression)
IMPLEMENT_ACCEPT_METHOD(LessThanExpression)
IMPLEMENT_ACCEPT_METHOD(LessOrEqualExpression)
IMPLEMENT_ACCEPT_METHOD(GreaterThanExpression)
IMPLEMENT_ACCEPT_METHOD(GreaterOrEqualExpression)
IMPLEMENT_ACCEPT_METHOD(ScopedExpression)
IMPLEMENT_ACCEPT_METHOD(Statements)
IMPLEMENT_ACCEPT_METHOD(EvalStatement)
IMPLEMENT_ACCEPT_METHOD(ReturnStatement)
IMPLEMENT_ACCEPT_METHOD(SourceStatement)
IMPLEMENT_ACCEPT_METHOD(SourceRegionStatement)
IMPLEMENT_ACCEPT_METHOD(DescribeStatement)
IMPLEMENT_ACCEPT_METHOD(EmptyStatement)
IMPLEMENT_ACCEPT_METHOD(BreakStatement)
IMPLEMENT_ACCEPT_METHOD(ContinueStatement)
IMPLEMENT_ACCEPT_METHOD(SkipStatement)
IMPLEMENT_ACCEPT_METHOD(ProtectStatement)
IMPLEMENT_ACCEPT_METHOD(UnprotectStatement)
IMPLEMENT_ACCEPT_METHOD(CiaoOrQuitStatement)
IMPLEMENT_ACCEPT_METHOD(AssignmentStatement)
IMPLEMENT_ACCEPT_METHOD(Identifier)
IMPLEMENT_ACCEPT_METHOD(FullyQualifiedIdentifier)
IMPLEMENT_ACCEPT_METHOD(IndeterminateDeclaration)
IMPLEMENT_ACCEPT_METHOD(RingDefinition)
IMPLEMENT_ACCEPT_METHOD(RingAssignStatement)
IMPLEMENT_ACCEPT_METHOD(FieldAccessExpression)
IMPLEMENT_ACCEPT_METHOD(RecordExpression)
IMPLEMENT_ACCEPT_METHOD(UseStatement)
IMPLEMENT_ACCEPT_METHOD(TimeStatement)
//IMPLEMENT_ACCEPT_METHOD(SetStatement)
//IMPLEMENT_ACCEPT_METHOD(UnsetStatement)
IMPLEMENT_ACCEPT_METHOD(ForStatement)
IMPLEMENT_ACCEPT_METHOD(ForeachStatement)
IMPLEMENT_ACCEPT_METHOD(WhileStatement)
IMPLEMENT_ACCEPT_METHOD(RepeatUntilStatement)
IMPLEMENT_ACCEPT_METHOD(IfStatement)
IMPLEMENT_ACCEPT_METHOD(UsingStatement)
IMPLEMENT_ACCEPT_METHOD(AliasStatement)
IMPLEMENT_ACCEPT_METHOD(BlockStatement)
IMPLEMENT_ACCEPT_METHOD(DefineStatement)
IMPLEMENT_ACCEPT_METHOD(LambdaExpression)
IMPLEMENT_ACCEPT_METHOD(PackageStatement)
IMPLEMENT_ACCEPT_METHOD(PrintStatement)
IMPLEMENT_ACCEPT_METHOD(HelpStatement)
IMPLEMENT_ACCEPT_METHOD(BackwardCompatibleExp)
IMPLEMENT_ACCEPT_METHOD(IsDefinedExp)
IMPLEMENT_ACCEPT_METHOD(FunctionDeclaration)
IMPLEMENT_ACCEPT_METHOD(TryStatement)

#undef IMPLEMENT_ACCEPT_METHOD

namespace
{
	bool IsBlankOrSemicolon(char c) { return IsBlank(c) || c==';'; }
}

HelpStatement::HelpStatement(const LexerNS::Token &tokHelp) :
	Statement(tokHelp.getBegin(), tokHelp.getEnd()),
	lexeme(tokHelp.lexeme()),
	topic(algorithm::trim_left_copy(algorithm::trim_right_copy_if(this->lexeme.substr(1, string::npos), IsBlankOrSemicolon)))
{
	assert(lexeme.length() && lexeme[0]=='?');
}

InvocationExpression::InvocationExpression(const intrusive_ptr<Expression> targetExp, const boost::shared_ptr<const Token> tokEllipsis, const vector<Argument> &args, const CharPointer &closedParenthesisPosition, const string packageName) :
	Expression(targetExp->getBegin(), closedParenthesisPosition),
	targetExp(targetExp),
	tokEllipsis(tokEllipsis),
	ellipsis(tokEllipsis),
	args(args),
	packageName(packageName.length() ? packageName : ParserNS::Parser::TopLevelPackageName)
{
	assert(targetExp);
}

bool Identifier::isLeftValue() const {
	return true;
}

bool FullyQualifiedIdentifier::isLeftValue() const {
	return true;
}

BigRat FloatLiteral::buildTheRational(const Token &t) {
	std::string l = t.lexeme();
	const size_t dotPos = l.find('.');
	assert(dotPos!=l.npos);
	assert(l.length()-dotPos>0);
	l.erase(dotPos, 1);
	return BigRat(BigIntFromString(l), power(BigInt(10), BigInt((l.length()-dotPos))));
}

bool Expression::isLeftValue() const {
	return false;
}

bool FieldAccessExpression::isLeftValue() const {
	return targetExp->isLeftValue();
}

bool IndexedAccessExpression::isLeftValue() const {
	return targetExp->isLeftValue();
}

void DumpAsTreeVisitor::indent() const {
	for(int a=0; a<this->indentationLevel; ++a)
		this->out << "  ";
}

void DumpAsTreeVisitor::dumpIndented(const intrusive_ptr<ParsedObject> o) {
	++this->indentationLevel;
	o->accept(this);
	--this->indentationLevel;
}

ostream &operator<<(ostream &out, const ParsedObject * const parsedObject) {
	return parsedObject->dumpAsString(out);
}

namespace {
	void outVarData(ostream &out, const StaticEnv::VarData &vd) {
		out << "{" << vd.depth << ", " << vd.index << ", isCap=" << vd.isCapturedValue << ", isIt=" << vd.isIterationVar << ", isII=" << vd.isImplictlyImported << "}";
	}
}

ostream &Identifier::dumpAsString(ostream &out) const {
	out << this->identifier;
	outVarData(out, this->varData);
	return out;
}

void DumpAsTreeVisitor::visit(Identifier &id) {
	this->indent();
	this->out << id.identifier;
	outVarData(this->out, id.varData);
	this->out << "\n";
}

ostream &FullyQualifiedIdentifier::dumpAsString(ostream &out) const {
	return out << this->pkgName << '.' << this->id;
}

void DumpAsTreeVisitor::visit(FullyQualifiedIdentifier &fullyQualifiedIdentifier) {
	this->indent();
	this->out << fullyQualifiedIdentifier.pkgName << '.' << fullyQualifiedIdentifier.id << '\n';
}

ostream &EvalStatement::dumpAsString(ostream &out) const {
	return out << exp << "; ";
}

void DumpAsTreeVisitor::visit(EvalStatement &evalStatement) {
	this->indent();
	out << "; (eval-statement)\n";
	this->dumpIndented(evalStatement.exp);
}

ostream &ReturnStatement::dumpAsString(ostream &out) const {
	out << "Return ";
	if (exp)
		out << exp;
	return out << "; ";
}

void DumpAsTreeVisitor::visit(ReturnStatement &returnStatement) {
	this->indent();
	out << "Return\n";
	if (returnStatement.exp)
		this->dumpIndented(returnStatement.exp);
}

ostream &SourceStatement::dumpAsString(ostream &out) const {
	return out << (this->ttype==TT_SOURCE ? "Source ": (this->ttype==TT_LOAD ? "Load" : "<< ")) << this->exp << "; ";
}

void DumpAsTreeVisitor::visit(SourceStatement &sourceStatement) {
	this->indent();
	this->out << (sourceStatement.ttype==TT_SOURCE ? "Source\n": (sourceStatement.ttype==TT_LOAD ? "Load\n": "<<\n"));
	this->dumpIndented(sourceStatement.exp);
}

ostream &SourceRegionStatement::dumpAsString(ostream &out) const {
  return out << "SourceRegion " << this->expFromLine << ", " << this->expFromChar << " to " << this->expToLine << ", " << this->expToChar << " in " << this->expFileName << "; ";
}

void DumpAsTreeVisitor::visit(SourceRegionStatement &sourceregionStatement) {
	this->indent();
	this->out << "SourceRegion\n";
	this->dumpIndented(sourceregionStatement.expFromLine);
	this->dumpIndented(sourceregionStatement.expFromChar);
	this->dumpIndented(sourceregionStatement.expToLine);
	this->dumpIndented(sourceregionStatement.expToChar);
	this->dumpIndented(sourceregionStatement.expFileName);
}

ostream &DescribeStatement::dumpAsString(ostream &out) const {
	return out << "Describe " << exp << "; ";
}

void DumpAsTreeVisitor::visit(DescribeStatement &describeStatement) {
	this->indent();
	out << "Describe\n";
	this->dumpIndented(describeStatement.exp);
}

ostream &UseStatement::dumpAsString(ostream &out) const {
	out << "Use ";
	if (identifier)
		out << identifier;
	if (ringDefinition) {
		if (identifier)
			out << " ::= ";
		out << ringDefinition;
	}
	return out << "; ";
}

void DumpAsTreeVisitor::visit(UseStatement &useStatement) {
	this->indent();
	out << "Use\n";
	if (useStatement.identifier)
		this->dumpIndented(useStatement.identifier);
	if (useStatement.ringDefinition) {
		this->indent();
		out << "::=\n";
		this->dumpIndented(useStatement.ringDefinition);
	}
}

/*ostream &SetStatement::dumpAsString(ostream &out) const {
	out << "Set " << identifier;
	if (exp)
		out << " := " << exp;
	return out << "; ";
}

void DumpAsTreeVisitor::visit(SetStatement &setStatement) {
	this->indent();
	out << "Set\n";
	this->dumpIndented(setStatement.identifier);
	if (setStatement.exp) {
		this->indent();
		out << ":=\n";
		this->dumpIndented(setStatement.exp);
	}
}

ostream &UnsetStatement::dumpAsString(ostream &out) const {
	return out << "Unset " << identifier << "; ";
}

void DumpAsTreeVisitor::visit(UnsetStatement &unsetStatement) {
	this->indent();
	this->out << "Unset\n";
	this->dumpIndented(unsetStatement.identifier);
}*/

ostream &BinaryExpression::dumpAsString(ostream &out) const {
	const string op = this->beginOperatorSourcePosition.stringTo(this->endOperatorSourcePosition);
	if (isalpha(op[0]))
		return out << leftExp << ' ' << op << ' ' << rightExp;
	return out << leftExp << op << rightExp;
}

void DumpAsTreeVisitor::visit(BinaryExpression &binaryExpression) {
	this->indent();
	this->out << binaryExpression.beginOperatorSourcePosition.stringTo(binaryExpression.endOperatorSourcePosition) << "\n";
	this->dumpIndented(binaryExpression.leftExp);
	this->dumpIndented(binaryExpression.rightExp);
}

ostream &UnaryExpression::dumpAsString(ostream &out) const {
	if (isalpha(op[0]))
		return out << op << ' ' << exp;
	return out << op << exp;
}

void DumpAsTreeVisitor::visit(UnaryExpression &unaryExpression) {
	this->indent();
	this->out << unaryExpression.op << "\n";
	this->dumpIndented(unaryExpression.exp);
}

ostream &LiteralExpression::dumpAsString(ostream &out) const {
	return out << beginSourcePosition.stringTo(endSourcePosition);
}

void DumpAsTreeVisitor::visit(LiteralExpression &literalExpression) {
	this->indent();
	this->out << literalExpression.getBegin().stringTo(literalExpression.getEnd()) << '\n';
}

ostream &AssignmentStatement::dumpAsString(ostream &out) const {
	return out << leftExp << " := " << rightExp << "; ";
}

void DumpAsTreeVisitor::visit(AssignmentStatement &assignmentStatement) {
	this->indent();
	this->out << ":=\n";
	this->dumpIndented(assignmentStatement.leftExp);
	this->dumpIndented(assignmentStatement.rightExp);
}

ostream &RingAssignStatement::dumpAsString(ostream &out) const {
	return out << this->leftExp << " ::= " << this->ringDef << "; ";
}

void DumpAsTreeVisitor::visit(RingAssignStatement &ringAssignStatement) {
	this->indent();
	out << "::=\n";
	this->dumpIndented(ringAssignStatement.leftExp);
	this->dumpIndented(ringAssignStatement.ringDef);
}

ostream &IsDefinedExp::dumpAsString(ostream &out) const {
	return out << "IsDefined(" << this->expId << ")";
}

void DumpAsTreeVisitor::visit(IsDefinedExp &isDefinedExp) {
	this->indent();
	this->out << "IsDefined(" << isDefinedExp.expId << ")\n";
}

ostream &TimeStatement::dumpAsString(ostream &out) const {
	return out << "Time " << stmt;
}

void DumpAsTreeVisitor::visit(TimeStatement &timeStatement) {
	this->indent();
	out << "Time\n";
	this->dumpIndented(timeStatement.stmt);
}

ostream &BackwardCompatibleExp::dumpAsString(ostream &out) const {
	return out << "${" << innerExp << "}$";
}

void DumpAsTreeVisitor::visit(BackwardCompatibleExp &backwardCompatibleExp) {
	this->indent();
	out << "${\n";
	this->dumpIndented(backwardCompatibleExp.innerExp);
	this->indent();
	out << "}$\n";
}

namespace {
	template <typename T>
	void outputExpList(ostream &out, const vector<T> &list, const string &separator) {
		bool first=true;
		for (const T& e: list)
                {
                  if (first)
                    first=false;
                  else
                    out << separator;
                  out << e;
		}
	}

	void outputArgumentList(ostream &out, const vector<Argument> &list) {
		bool first=true;
		for (const Argument &arg: list)
                {
                  if (first)
                    first=false;
                  else
                    out << ", ";
                  if (arg.byRef)
                    out << "Ref ";
                  out << arg.exp;
		}
	}

	template <typename T>
	void outputCommaSeparatedExpList(ostream &out, const vector<T> &list) {
		outputExpList(out, list, ", ");
	}
}

ostream &CartesianProductExpression::dumpAsString(ostream &out) const {
	outputExpList(out, operands, "><");
	return out;
}

void DumpAsTreeVisitor::visit(CartesianProductExpression &cartesianProductExpression) {
	this->indent();
	out << "><\n";
	for (intrusive_ptr<Expression> e: cartesianProductExpression.operands)
		this->dumpIndented(e);
}

ostream &IdInExpSuchThatExp::dumpAsString(ostream &out) const {
	return out << "[" << identifier << " In " << exp1 << " | " << exp2 << "]";
}

void DumpAsTreeVisitor::visit(IdInExpSuchThatExp &idInExpSuchThatExp) {
	this->indent();
	this->out << "[ " << idInExpSuchThatExp.identifier << " In ... | ...]\n";
	this->dumpIndented(idInExpSuchThatExp.exp1);
	this->dumpIndented(idInExpSuchThatExp.exp2);
}

ostream &ExpSuchThatIdInExpAndExp::dumpAsString(ostream &out) const {
	out << "[" << exp1 << " | " << identifier << " In " << exp2;
	if (optionalExp3)
		out << " And " << optionalExp3;
	return out << "]";
}

void DumpAsTreeVisitor::visit(ExpSuchThatIdInExpAndExp &expSuchThatIdInExpAndExp) {
	this->indent();
	out << "[ ... | " << expSuchThatIdInExpAndExp.identifier << " In ... ";
	if (expSuchThatIdInExpAndExp.optionalExp3)
		out << "And ... ";
	out << "]\n";
	this->dumpIndented(expSuchThatIdInExpAndExp.exp1);
	this->dumpIndented(expSuchThatIdInExpAndExp.exp2);
	if (expSuchThatIdInExpAndExp.optionalExp3)
		this->dumpIndented(expSuchThatIdInExpAndExp.optionalExp3);
}

ostream &ListExpression::dumpAsString(ostream &out) const {
	out << "[";
	outputCommaSeparatedExpList(out, exps);
	return out << "]";
}

void DumpAsTreeVisitor::visit(ListExpression &listExpression) {
	this->indent();
	out << "[]\n";
	for (intrusive_ptr<Expression> e: listExpression.exps)
		this->dumpIndented(e);
}

ostream &IndeterminateDeclaration::dumpAsString(ostream &out) const {
	out << identifier;
	if (!ranges.empty()) {
		out << "[";
		bool first = true;
		for(vector<pair <intrusive_ptr<Expression>, intrusive_ptr<Expression> > >::const_iterator it=ranges.begin(); it!=ranges.end(); ++it) {
			if (first)
				first = false;
			else
				out << ", ";
			out << it->first << ".." << it->second;
		}
		return out << "]";
	}
	return out;
}

void DumpAsTreeVisitor::visit(IndeterminateDeclaration &indeterminateDeclaration) {
	indeterminateDeclaration.identifier->accept(this);
	if (indeterminateDeclaration.ranges.empty())
		return;
	++this->indentationLevel;
	this->indent();
	this->out << "[]\n";
	for(vector<pair <intrusive_ptr<Expression>, intrusive_ptr<Expression> > >::const_iterator it=indeterminateDeclaration.ranges.begin(); it!=indeterminateDeclaration.ranges.end(); ++it) {
		this->indent();
		out << ".." << "\n";
		this->dumpIndented(it->first);
		this->dumpIndented(it->second);
	}
	--this->indentationLevel;
}

ostream &RingDefinition::dumpAsString(ostream &out) const {
	out << identifier;
	if (optionalExp)
		out << "/(" << optionalExp << ")";
	if (indeterminates.size()) {
		out << "[";
		outputCommaSeparatedExpList(out, indeterminates);
		out << "]";
	}
	return out;
}

void DumpAsTreeVisitor::visit(RingDefinition &ringDefinition) {
	this->indent();
	this->out << ringDefinition.identifier << '\n';
	if (ringDefinition.optionalExp) {
		this->indent();
		this->out << "/()\n";
		this->dumpIndented(ringDefinition.optionalExp);
	}
	if (ringDefinition.indeterminates.size()) {
		this->indent();
		this->out << "[]\n";
		for (intrusive_ptr<IndeterminateDeclaration> id: ringDefinition.indeterminates)
			this->dumpIndented(id);
	}
}

ostream &FieldAccessExpression::dumpAsString(ostream &out) const {
	return out << targetExp << '.' << name;
}

void DumpAsTreeVisitor::visit(FieldAccessExpression &fieldAccessExpression) {
	this->indent();
	this->out << "." << fieldAccessExpression.name << "\n";
	this->dumpIndented(fieldAccessExpression.targetExp);
}

ostream &IndexedAccessExpression::dumpAsString(ostream &out) const {
	out << targetExp << "[";
	outputCommaSeparatedExpList(out, indexes);
	return out << "]";
}

void DumpAsTreeVisitor::visit(IndexedAccessExpression &indexedAccessExpression) {
	this->indent();
	out << "[]\n";
	this->dumpIndented(indexedAccessExpression.targetExp);
	for (intrusive_ptr<Expression> e: indexedAccessExpression.indexes)
		this->dumpIndented(e);
}

ostream &InvocationExpression::dumpAsString(ostream &out) const {
	out << targetExp << "(";
	if (ellipsis)
		out << "...";
	else
		outputArgumentList(out, args);
	return out << ")";
}

void DumpAsTreeVisitor::visit(InvocationExpression &invocationExpression) {
	this->indent();
	this->out << "()\n";
	this->dumpIndented(invocationExpression.targetExp);
	++this->indentationLevel;
	if (invocationExpression.ellipsis) {
		++this->indentationLevel;
		this->indent();
		--this->indentationLevel;
		this->out << "...\n";
	} else
          for (const Argument &arg: invocationExpression.args)
          {
            if (arg.byRef) {
              this->indent();
              this->out << "Ref\n";
            }
            this->dumpIndented(arg.exp);
          }
	--this->indentationLevel;
}

ostream &PrintStatement::dumpAsString(ostream &out) const {
	out << (ttype==TT_PRINT ? "Print " : "PrintLn ");
	outputCommaSeparatedExpList(out, exps);
	if (onExp)
		out << " On " << onExp;
	return out << "; ";
}

void DumpAsTreeVisitor::visit(PrintStatement &printStatement) {
	this->indent();
	this->out << (printStatement.ttype==TT_PRINT ? "Print\n" : "PrintLn\n");
	for (intrusive_ptr<Expression> e: printStatement.exps)
		this->dumpIndented(e);
	if (printStatement.onExp) {
		this->indent();
		this->out << "On\n";
		this->dumpIndented(printStatement.onExp);
	}
}

ostream &AliasStatement::dumpAsString(ostream &out) const {
	out << "alias ";
	bool first=true;
	for (const Binding &b: bindings)
        {
          if (first)
            first=false;
          else
            out << ", ";
          out << b.identifier << " := " << b.packageName;
	}
	/* if (this->thereIsIn)
		return out << "In "<< this->statements << "EndAlias "; */
	return out << "; ";
}

void DumpAsTreeVisitor::visit(AliasStatement &aliasStatement) {
	this->indent();
	out << "alias\n";
	for (const Binding &b: aliasStatement.bindings)
        {
          ++this->indentationLevel;
          this->indent();
          out << b.identifier << " := " << b.packageName << "\n";
          --this->indentationLevel;
	}
	/* if (aliasStatement.thereIsIn) {
		this->indent();
		this->out << "In\n";
		this->dumpIndented(aliasStatement.statements);
		this->indent();
		this->out << "EndAlias\n";
	} */
}

ostream &BlockStatement::dumpAsString(ostream &out) const {
	return out << "Block " << this->statements << "EndBlock ";
}

void DumpAsTreeVisitor::visit(BlockStatement &blockStatement) {
	this->indent();
	out << "Block\n";
	this->dumpIndented(blockStatement.statements);
	this->indent();
	out << "EndBlock\n";
}

namespace {
	void outParams(bool thereIsEllipsis, const vector<Param> &params, ostream &out) {
		if (thereIsEllipsis) {
			out << "...";
			return;
		}
		bool first=true;
		for (const Param &p: params)
                {
                  if (first)
                    first=false;
                  else
                    out << ", ";
                  out << (p.byRef ? "Ref ":"") << (p.opt ? "Opt ":"") << p.name;
		}
	}
}

ostream &DefineStatement::dumpAsString(ostream &out) const {
	return out << "Define " << this->name << this->funDecl << "EndDefine";
}

namespace {

	string importTypeToString(const Import::ImportType tt) {
		switch (tt) {
		case Import::IT_BYREF: return "ImportByRef ";
		case Import::IT_BYVALUE: return "ImportByValue ";
		case Import::IT_TOPLEVEL: return "TopLevel ";
		}
		assert(false);
		return "";
	}

}

ostream &Import::dumpAsString(ostream &out) const {
	out << importTypeToString(this->type) << (this->implicit ? "/* (inferred) */ ":"") << this->name << "; ";
	if (this->type==IT_BYVALUE) {
		assert(this->expId);
		out << "/* byValueIndex=" << this->byValueIndex << " source=" << this->expId << " */";
	}
	return out;
}

void DumpAsTreeVisitor::visit(Import &import) {
	this->indent();
	import.dumpAsString(out);
	out << "\n";
}

void DumpAsTreeVisitor::visit(DefineStatement &defineStatement) {
	this->indent();
	out << "Define " << defineStatement.name;
	++this->indentationLevel;
	defineStatement.funDecl->accept(this);
	--this->indentationLevel;
	out << "EndDefine\n";
}

ostream &LambdaExpression::dumpAsString(ostream &out) const {
	return out << "Func " << this->funDecl << "EndFunc";
}

void DumpAsTreeVisitor::visit(LambdaExpression &lambdaExpression) {
	this->indent();
	out << "Func ";
	this->dumpIndented(lambdaExpression.funDecl);
	this->indent();
	out << "EndFunc\n";
}

ostream &FunctionDeclaration::dumpAsString(ostream &out) const {
	out << "(";
	outParams(this->thereIsEllipsis, this->params, out);
	out << ") ";
	for (const Import &i: this->imports)
		i.dumpAsString(out);
	return out << this->statements;
}

void DumpAsTreeVisitor::visit(FunctionDeclaration &fnDecl) {
	out << '(';
	outParams(fnDecl.thereIsEllipsis, fnDecl.params, out);
	out << ")\n";
	++this->indentationLevel;
	this->indent();
	out << "// Locals: ";
	for(map<string, int>::const_iterator it=fnDecl.localNames.begin(); it!=fnDecl.localNames.end(); ++it)
		out << it->first << "(" << it->second << ") ";
	out << "\n";
	for (const Import &i: fnDecl.imports)
        {
          this->indent();
          i.dumpAsString(this->out);
          this->out << '\n';
	}
	--this->indentationLevel;
	this->dumpIndented(fnDecl.statements);
}

ostream &RecordExpression::dumpAsString(ostream &out) const {
	out << "record[";
	bool first=true;
	for (const RecordField &f: fields)
        {
          if (first)
            first=false;
          else
            out << ", ";
          out << f.name << " := " << f.initExp;
	}
	return out << "]";
}

void DumpAsTreeVisitor::visit(RecordExpression &recordExpression) {
	this->indent();
	out << "record[\n";
	for (const RecordField &f: recordExpression.fields)
        {
          this->indent();
          out << ":=\n";
          this->dumpIndented(f.name);
          this->dumpIndented(f.initExp);
	}
	this->indent();
	out << "]\n";
}

ostream &PackageStatement::dumpAsString(ostream &out) const {
	return out << "Package " << name << ' ' << this->statements << "EndPackage ";
}

void DumpAsTreeVisitor::visit(PackageStatement &packageStatement) {
	this->indent();
	out << "Package " << packageStatement.name << '\n';
	this->dumpIndented(packageStatement.statements);
	this->indent();
	out << "EndPackage\n";
}

ostream &HelpStatement::dumpAsString(ostream &out) const {
	return out << lexeme << '\n';
}

void DumpAsTreeVisitor::visit(HelpStatement &helpStatement) {
	this->indent();
	out << helpStatement.lexeme << '\n';
}

char operandToChar(TokenType tt) {
	switch (tt) {
	case TT_PLUS:
		return '+';
	case TT_MINUS:
		return '-';
	case TT_STAR:
		return '*';
	case TT_SLASH:
		return '/';
	case TT_MOD:
		return '%';
	case TT_COLON:
		return ':';
	default:
		assert(false);
		return '?';
	}
}

namespace {
	ostream &dumpSequenceOfOpsAsString(ostream &out, vector<Token>::const_iterator opIt, vector<intrusive_ptr<Expression> >::const_iterator expIt, vector<intrusive_ptr<Expression> >::const_iterator expEnd) {
		for(bool isFirst=true; expIt!=expEnd; ++expIt) {
			if (isFirst)
				isFirst = false;
			else
				out << (*opIt++).lexeme();
			out << *expIt;
		}
		return out;
	}

}

void DumpAsTreeVisitor::dumpSequenceOfOpsAsTree(vector<Token>::const_iterator opIt, vector<intrusive_ptr<Expression> >::const_iterator expIt, vector<intrusive_ptr<Expression> >::const_iterator expEnd) {
	++this->indentationLevel;
	for(bool isFirst=true; expIt!=expEnd; ++expIt) {
		if (isFirst) {
			isFirst = false;
		} else {
			this->indent();
			out << (*opIt++).lexeme() << '\n';
		}
		this->dumpIndented(*expIt);
	}
	--this->indentationLevel;
}

ostream &SummationExpression::dumpAsString(ostream &out) const {
	return dumpSequenceOfOpsAsString(out, operators.begin(), operands.begin(), operands.end());
}

void DumpAsTreeVisitor::visit(SummationExpression &summationExpression) {
	this->indent();
	out << "summation\n";
	this->dumpSequenceOfOpsAsTree(summationExpression.operators.begin(), summationExpression.operands.begin(), summationExpression.operands.end());
}

ostream &MultiplicationExpression::dumpAsString(ostream &out) const {
	return dumpSequenceOfOpsAsString(out, operators.begin(), operands.begin(), operands.end());
}

void DumpAsTreeVisitor::visit(MultiplicationExpression &multiplicationExpression) {
	this->indent();
	out << "multiplication\n";
	this->dumpSequenceOfOpsAsTree(multiplicationExpression.operators.begin(), multiplicationExpression.operands.begin(), multiplicationExpression.operands.end());
}

ostream &ForStatement::dumpAsString(ostream &out) const {
	out << "For " << this->identifier << ":=" << this->beginExp << " To " << this->endExp;
	if (this->stepExp)
		out << " Step " << stepExp;
	return out << " Do " << this->statements << "EndFor ";
}

void DumpAsTreeVisitor::visit(ForStatement &forStatement) {
	this->indent();
	out << "For " << forStatement.identifier << ":=\n";
	this->dumpIndented(forStatement.beginExp);
	this->indent();
	out << "To\n";
	this->dumpIndented(forStatement.endExp);
	if (forStatement.stepExp) {
		this->indent();
		out << "Step\n";
		this->dumpIndented(forStatement.stepExp);
	}
	this->indent();
	out << "Do\n";
	this->dumpIndented(forStatement.statements);
	this->indent();
	out << "EndFor\n";
}

ostream &UsingStatement::dumpAsString(ostream &out) const {
	return out << "Using " << identifier << " Do " << this->statements << "EndUsing ";
}

void DumpAsTreeVisitor::visit(UsingStatement &usingStatement) {
	this->indent();
	out << "Using " << usingStatement.identifier << " Do\n";
	this->dumpIndented(usingStatement.statements);
	this->indent();
	out << "EndUsing\n";
}

ostream &ForeachStatement::dumpAsString(ostream &out) const {
	return out << "Foreach " << this->identifier << " In " << inExp << " Do " << this->statements << "EndForeach ";
}

void DumpAsTreeVisitor::visit(ForeachStatement &foreachStatement) {
	this->indent();
	out << "Foreach " << foreachStatement.identifier << " In\n";
	this->dumpIndented(foreachStatement.inExp);
	this->indent();
	out << "Do\n";
	this->dumpIndented(foreachStatement.statements);
	this->indent();
	out << "EndForeach\n";
}

ostream &TryStatement::dumpAsString(ostream &out) const {
	return out << "Try " << this->tryStatements << "UponError " << this->identifier << " Do " << this->uponErrorStatements << "EndTry ";
}

void DumpAsTreeVisitor::visit(TryStatement &tryStatement) {
	this->indent();
	out << "Try\n";
	this->dumpIndented(tryStatement.tryStatements);
	this->indent();
	out << "UponError " << tryStatement.identifier << " Do\n";
	this->dumpIndented(tryStatement.uponErrorStatements);
	this->indent();
	out << "EndTry\n";
}

ostream &Statements::dumpAsString(ostream &out) const {
  for (intrusive_ptr<Statement> s: this->statements)
		out << s << " ";
	return out;
}

void DumpAsTreeVisitor::visit(Statements &statements) {
  for (intrusive_ptr<Statement> s: statements.statements)
		s->accept(this);
}

ostream &WhileStatement::dumpAsString(ostream &out) const {
	return out << "While " << exp << " Do "<< this->statements << "EndWhile ";
}

void DumpAsTreeVisitor::visit(WhileStatement &whileStatement) {
	this->indent();
	out << "While\n";
	this->dumpIndented(whileStatement.exp);
	this->indent();
	out << "Do\n";
	this->dumpIndented(whileStatement.statements);
	this->indent();
	out << "EndWhile\n";
}

ostream &RepeatUntilStatement::dumpAsString(ostream &out) const {
	out << "Repeat " << this->statements;
	if (!this->optExp)
		return out << "EndRepeat";
	return out << "Until " << this->optExp << "; ";
}

void DumpAsTreeVisitor::visit(RepeatUntilStatement &repeatUntilStatement) {
	this->indent();
	out << "Repeat\n";
	this->dumpIndented(repeatUntilStatement.statements);
	this->indent();
	if (repeatUntilStatement.optExp) {
		out << "Until\n";
		this->dumpIndented(repeatUntilStatement.optExp);
	} else
		out << "EndRepeat\n";
}

ostream &EmptyStatement::dumpAsString(ostream &out) const {
	return out << "; ";
}

void DumpAsTreeVisitor::visit(EmptyStatement &) {
	this->indent();
	this->out << "; (empty-statement)\n";
}

ostream &IfStatement::dumpAsString(ostream &out) const {
	bool first = true;
	for (const IfBranch &ib: branches)
        {
          if (ib.optExp) {
            if (first) {
              first = false;
              out << "If ";
            } else
              out << "Elif ";
            out << ib.optExp << " Then ";
          } else
            out << "Else ";
          out << ib.statements;
	}
	return out << "EndIf ";
}

void DumpAsTreeVisitor::visit(IfStatement &ifStatement) {
	bool first = true;
	for (const IfBranch &ib: ifStatement.branches)
        {
          this->indent();
          if (ib.optExp) {
            if (first) {
              first = false;
              out << "If\n";
            } else
              out << "Elif\n";
            this->dumpIndented(ib.optExp);
            this->indent();
            out << "Then\n";
          } else
            out << "Else\n";
          this->dumpIndented(ib.statements);
	}
	this->indent();
	out << "EndIf\n";
}

ostream &BreakStatement::dumpAsString(ostream &out) const {
	return out << "Break; ";
}

void DumpAsTreeVisitor::visit(BreakStatement &) {
	this->indent();
	this->out << "Break\n";
}

ostream &ContinueStatement::dumpAsString(ostream &out) const {
	return out << "Continue; ";
}

void DumpAsTreeVisitor::visit(ContinueStatement &) {
	this->indent();
	out << "Continue\n";
}

ostream &SkipStatement::dumpAsString(ostream &out) const {
	return out << "Skip; ";
}

void DumpAsTreeVisitor::visit(SkipStatement &) {
	this->indent();
	out << "Skip\n";
}

ostream &ProtectStatement::dumpAsString(ostream &out) const {
	out << "Protect " << this->expId;
	if (this->optExp)
		out << " : " << this->optExp;
	return out << "; ";
}

void DumpAsTreeVisitor::visit(ProtectStatement &protectStatement) {
	this->indent();
	out << "Protect " << protectStatement.expId;
	if (protectStatement.optExp) {
		this->out << " :\n";
		this->dumpIndented(protectStatement.optExp);
	} else
		out << "\n";
}

ostream &UnprotectStatement::dumpAsString(ostream &out) const {
	return out << "Unprotect " << this->expId << "; ";
}

void DumpAsTreeVisitor::visit(UnprotectStatement &unprotectStatement) {
	this->indent();
	out << "Unprotect " << unprotectStatement.expId << "\n";
}

ostream &CiaoOrQuitStatement::dumpAsString(ostream &out) const {
  return out << (ttype==TT_CIAO?"Ciao; ":(ttype==TT_EXIT?"Exit; ":"Quit; "));
}

void DumpAsTreeVisitor::visit(CiaoOrQuitStatement &ciaoOrQuitStatement) {
	this->indent();
	this->out << (ciaoOrQuitStatement.ttype==TT_CIAO?"Ciao\n":(ciaoOrQuitStatement.ttype==TT_EXIT?"Exit\n":"Quit\n"));
}

void ParsedObjectVisitor::visit(BinaryExpression &binaryExp) {
	binaryExp.leftExp->accept(this);
	binaryExp.rightExp->accept(this);
}

void ParsedObjectVisitor::visit(UnaryExpression &unaryExp) {
	unaryExp.exp->accept(this);
}

void ParsedObjectVisitor::visit(LiteralExpression &) {}

void ParsedObjectVisitor::visit(IntLiteral &intLiteral) {
	this->visit(static_cast<LiteralExpression &>(intLiteral));
}

void ParsedObjectVisitor::visit(FloatLiteral &floatLiteral) {
	this->visit(static_cast<LiteralExpression &>(floatLiteral));
}

void ParsedObjectVisitor::visit(StringLiteral &stringLiteral) {
	this->visit(static_cast<LiteralExpression &>(stringLiteral));
}

void ParsedObjectVisitor::visit(BoolLiteral &boolLiteral) {
	this->visit(static_cast<LiteralExpression &>(boolLiteral));
}

void ParsedObjectVisitor::visit(IdInExpSuchThatExp &idInExpSuchThatExp) {
	idInExpSuchThatExp.exp1->accept(this);
	idInExpSuchThatExp.exp2->accept(this);
}

void ParsedObjectVisitor::visit(ExpSuchThatIdInExpAndExp &expSuchThatIdInExpAndExp) {
	expSuchThatIdInExpAndExp.exp1->accept(this);
	expSuchThatIdInExpAndExp.exp2->accept(this);
	if (expSuchThatIdInExpAndExp.optionalExp3)
		expSuchThatIdInExpAndExp.optionalExp3->accept(this);
}

void ParsedObjectVisitor::visit(ListExpression &listExpression) {
  for (intrusive_ptr<Expression> e: listExpression.exps)
			e->accept(this);
}

void ParsedObjectVisitor::visit(InvocationExpression &invocationExp) {
	invocationExp.targetExp->accept(this);
	for (const Argument &arg: invocationExp.args)
		arg.exp->accept(this);
}

#define VISIT_UNARY(CLASS) void ParsedObjectVisitor::visit(CLASS &e) { this->visit(static_cast<UnaryExpression &>(e)); }
	VISIT_UNARY(UnaryMinusExpression)
	VISIT_UNARY(UnaryPlusExpression)
#undef VISIT_UNARY

void ParsedObjectVisitor::visit(IndexedAccessExpression &indexedAccessExpression) {
	indexedAccessExpression.targetExp->accept(this);
	for (intrusive_ptr<Expression> e: indexedAccessExpression.indexes)
		e->accept(this);
}

void ParsedObjectVisitor::visit(CartesianProductExpression &cartesianProductExpression) {
  for (intrusive_ptr<Expression> e: cartesianProductExpression.operands)
		e->accept(this);
}

void ParsedObjectVisitor::visit(SummationExpression &summationExpression) {
  for (intrusive_ptr<Expression> e: summationExpression.operands)
		e->accept(this);
}

void ParsedObjectVisitor::visit(MultiplicationExpression &multiplicationExpression) {
  for (intrusive_ptr<Expression> e: multiplicationExpression.operands)
		e->accept(this);
}

#define VISIT_BINARY(CLASS) void ParsedObjectVisitor::visit(CLASS &e) { this->visit(static_cast<BinaryExpression &>(e)); }
	VISIT_BINARY(OrExpression)
	VISIT_BINARY(PowerExpression)
	VISIT_BINARY(DotDotExpression)
	VISIT_BINARY(SumExpression)
	VISIT_BINARY(SubtractionExpression)
	VISIT_BINARY(ProductExpression)
	VISIT_BINARY(DivisionExpression)
	VISIT_BINARY(ColonExpression)
	VISIT_BINARY(ModuloExpression)
	VISIT_BINARY(AndExpression)
	VISIT_BINARY(EqualExpression)
	VISIT_BINARY(NotEqualExpression)
	VISIT_BINARY(IsInExpression)
	VISIT_BINARY(LessThanExpression)
	VISIT_BINARY(LessOrEqualExpression)
	VISIT_BINARY(GreaterThanExpression)
	VISIT_BINARY(GreaterOrEqualExpression)
	VISIT_BINARY(ScopedExpression)
#undef VISIT_BINARY

void ParsedObjectVisitor::visit(Statements &stmts) {
          for (intrusive_ptr<Statement> s: stmts.statements)
		s->accept(this);
}

void ParsedObjectVisitor::visit(EvalStatement &evalStmt) {
	evalStmt.exp->accept(this);
}

void ParsedObjectVisitor::visit(ReturnStatement &returnStmt) {
	if (returnStmt.exp)
		returnStmt.exp->accept(this);
}

void ParsedObjectVisitor::visit(SourceStatement &sourceStmt) {
	sourceStmt.exp->accept(this);
}

void ParsedObjectVisitor::visit(SourceRegionStatement &sourceregionStmt) {
	sourceregionStmt.expFromLine->accept(this);
	sourceregionStmt.expFromChar->accept(this);
	sourceregionStmt.expToLine->accept(this);
	sourceregionStmt.expToChar->accept(this);
	sourceregionStmt.expFileName->accept(this);
}

void ParsedObjectVisitor::visit(DescribeStatement &describeStmt) {
	describeStmt.exp->accept(this);
}

void ParsedObjectVisitor::visit(EmptyStatement &) {}

void ParsedObjectVisitor::visit(BreakStatement &) {}

void ParsedObjectVisitor::visit(ContinueStatement &) {}

void ParsedObjectVisitor::visit(SkipStatement &) {}

void ParsedObjectVisitor::visit(ProtectStatement &protectStmt) {
	protectStmt.expId->accept(this);
	if (protectStmt.optExp)
		protectStmt.optExp->accept(this);
}

void ParsedObjectVisitor::visit(UnprotectStatement &unprotectStmt) {
	unprotectStmt.expId->accept(this);
}

void ParsedObjectVisitor::visit(CiaoOrQuitStatement &) {}

void ParsedObjectVisitor::visit(AssignmentStatement &assignmentStmt) {
	assignmentStmt.leftExp->accept(this);
	assignmentStmt.rightExp->accept(this);
}

void ParsedObjectVisitor::visit(Identifier &) {}

void ParsedObjectVisitor::visit(FullyQualifiedIdentifier &) {}

void ParsedObjectVisitor::visit(IndeterminateDeclaration &) {
}

void ParsedObjectVisitor::visit(RingDefinition &rd) {
	rd.identifier->accept(this);
	if (rd.optionalExp)
		rd.optionalExp->accept(this);
	for (boost::intrusive_ptr<IndeterminateDeclaration> ind: rd.indeterminates)
        {
          typedef pair<intrusive_ptr<Expression>, intrusive_ptr<Expression> > ExpPair;
          for (const ExpPair &p: ind->ranges)
          {
            p.first->accept(this);
            if (p.second)
              p.second->accept(this);
          }
	}
}

void ParsedObjectVisitor::visit(RingAssignStatement &ra) {
	ra.leftExp->accept(this);
	ra.ringDef->accept(this);
}

void ParsedObjectVisitor::visit(FieldAccessExpression &fieldAccessExpression) {
	fieldAccessExpression.targetExp->accept(this);
}

void ParsedObjectVisitor::visit(RecordExpression &recordExpression) {
  for (const RecordField &f: recordExpression.fields)
		f.initExp->accept(this);
}

void ParsedObjectVisitor::visit(UseStatement &useStatement) {
	if (useStatement.identifier)
		useStatement.identifier->accept(this);
	useStatement.ringDefinition->accept(this);
}

void ParsedObjectVisitor::visit(TimeStatement &timeStatement) {
	timeStatement.stmt->accept(this);
}

/*void ParsedObjectVisitor::visit(SetStatement &setStatement) {
	setStatement.identifier->accept(this);
	setStatement.exp->accept(this);
}

void ParsedObjectVisitor::visit(UnsetStatement &unsetStatement) {
	unsetStatement.identifier->accept(this);
}*/

void ParsedObjectVisitor::visit(ForStatement &forStmt) {
	forStmt.beginExp->accept(this);
	forStmt.endExp->accept(this);
	if (forStmt.stepExp)
		forStmt.stepExp->accept(this);
	forStmt.statements->accept(this);
}

void ParsedObjectVisitor::visit(ForeachStatement &foreachStmt) {
	foreachStmt.inExp->accept(this);
	foreachStmt.statements->accept(this);
}

void ParsedObjectVisitor::visit(TryStatement &tryStatement) {
	tryStatement.tryStatements->accept(this);
	tryStatement.uponErrorStatements->accept(this);
}

void ParsedObjectVisitor::visit(WhileStatement &whileStmt) {
	whileStmt.exp->accept(this);
	whileStmt.statements->accept(this);
}

void ParsedObjectVisitor::visit(RepeatUntilStatement &repeatStmt) {
	repeatStmt.statements->accept(this);
	if (repeatStmt.optExp)
		repeatStmt.optExp->accept(this);
}

void ParsedObjectVisitor::visit(IfStatement &ifStatement) {
  for (const IfBranch &ib: ifStatement.branches)
  {
    if (ib.optExp)
      ib.optExp->accept(this);
    ib.statements->accept(this);
  }
}

void ParsedObjectVisitor::visit(UsingStatement &usingStmt) {
	usingStmt.identifier->accept(this);
	usingStmt.statements->accept(this);
}

void ParsedObjectVisitor::visit(AliasStatement &) {
	/* if (aliasStmt.statements)
		aliasStmt.statements->accept(this); */
}

void ParsedObjectVisitor::visit(BlockStatement &blockStmt) {
	blockStmt.statements->accept(this);
}

void ParsedObjectVisitor::visit(FunctionDeclaration &funDecl) {
	funDecl.statements->accept(this);
}

void ParsedObjectVisitor::visit(DefineStatement &defineStmt) {
	defineStmt.expId->accept(this);
	defineStmt.funDecl->accept(this);
}

void ParsedObjectVisitor::visit(LambdaExpression &lambda) {
	lambda.funDecl->accept(this);
}

void ParsedObjectVisitor::visit(PackageStatement &packageStmt) {
	packageStmt.statements->accept(this);
}

void ParsedObjectVisitor::visit(PrintStatement &printStmt) {
  for (intrusive_ptr<Expression> e: printStmt.exps)
		e->accept(this);
	if (printStmt.onExp)
		printStmt.onExp->accept(this);
}

void ParsedObjectVisitor::visit(HelpStatement &) {}

void ParsedObjectVisitor::visit(BackwardCompatibleExp &backwardCompatibleExp) {
	backwardCompatibleExp.innerExp->accept(this);
}

void ParsedObjectVisitor::visit(IsDefinedExp &isDefinedExp) {
	isDefinedExp.expId->accept(this);
}

StaticEnv::VarData StaticEnv::lookup(const string &identifier) const {
	if (!parent)
		return VarData(IN_THE_WILD, -1, false, false, false); // must be considered IN_THE_WILD, even if it might exist at top-level, because it has not been imported
	IdMap::const_iterator it = this->identifierMap.find(identifier);
	if (it==this->identifierMap.end()) {
		VarData p = parent->lookup(identifier);
		if (p.depth>=0)
			++p.depth;
		return p;
	}
	return it->second;
}

void StaticEnv::add(const string &identifier, int depth, int slot, bool isIterationVar, bool isCapturedValue, bool isImplictlyImported) {
	this->identifierMap.insert(make_pair(identifier, StaticEnv::VarData(depth, slot, isIterationVar, isCapturedValue, isImplictlyImported)));
}

void AliasEnv::dump(intrusive_ptr<InterpreterNS::OSTREAM> out) const {
	if (this->aliases.size()) {
		out->print("alias ");
		bool first=true;
		for(AliasMap::const_iterator it = this->aliases.begin(); it!=this->aliases.end(); ++it) {
			if (first)
				first = false;
			else
				out->print(",\n      ");
			out->print(it->first)->print(" := ")->print(it->second);
		}
		out->print(";\n");
	} else
		out->print("/* There are no aliases */\n");
}

void AliasEnv::add(const string &from, const string &to) {
	this->aliases.erase(from);
	this->aliases.insert(make_pair(from, to));
}

bool AliasEnv::lookup(const string &identifier, string &alias) const {
	AliasMap::const_iterator it = this->aliases.find(identifier);
	if (it==this->aliases.end()) {
		if (parent)
			return parent->lookup(identifier, alias);
		return false;
	}
	alias = it->second;
	return true;
}

} // namespace AST
} // namespace CoCoA
