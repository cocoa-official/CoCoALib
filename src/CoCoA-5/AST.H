#ifndef CoCoA5_AST_H
#define CoCoA5_AST_H

//   Copyright (c) 2009  John Abbott,  Anna M. Bigatti
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

#include <ostream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <boost/intrusive_ptr.hpp>

#ifdef C5IDE
class QTreeWidget;
class QTreeWidgetItem;
#endif // #ifdef C5IDE

#include "CoCoA/BigInt.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/QuasiPoly.H"
#include "CoCoA/SparsePolyOps-hilbert.H"
#include "CoCoA/degree.H"
#include "CoCoA/factorization.H"
#include "Lexer.H"


namespace CoCoA {

  class ConstMatrixView; // fwd decl -- defined in matrix.H
  class ModuleElem;      // fwd decl -- defined in module.H
  class RingElem;        // fwd decl -- defined in ring.H
  class RingHom;         // fwd decl -- defined in RingHom.H
  class ideal;           // fwd decl -- defined in ideal.H
  class matrix;          // fwd decl -- defined in matrix.H
  class module;          // fwd decl -- defined in module.H
  class ring;            // fwd decl -- defined in ring.H
  class JBMill;

#ifdef C5IDE
namespace IDE {
	class Debugger;
} // namespace IDE
#endif // #ifdef C5IDE

template <typename T, typename U> inline boost::intrusive_ptr<T> intrusive_ptr_cast(boost::intrusive_ptr<U> p) {
	assert(boost::dynamic_pointer_cast<T>(p));
	return boost::static_pointer_cast<T>(p);
}

namespace AST {
	class ParsedObject;
	class InvocationExpression;
}

namespace InterpreterNS {
	class LeftValue;
	class RightValue;
	class BOOL;
	class INT;
	class RAT;
	class RECORD;
	class STRING;
	class LIST;
	class RuntimeEnvironment;
	class Frame;
	class VariableSlot;
	class ERROR_MESG;
	class OSTREAM;
	class ISTREAM;
	class Dispatcher;
	class TYPE;
	class RING;
	class RINGELEM;
	class MAT;
	class RINGHOM;
	class IDEAL;
	class MODULE;
	class MODULEELEM;
	class JBMillValue;

	class Value : public LexerNS::ReferenceCountedObject {
		friend class Dispatcher;
	public:
		virtual std::ostream &dumpAsString(std::ostream &out) const = 0;
		virtual ~Value() {}
		virtual boost::intrusive_ptr<RightValue> asRightValue() = 0; // important: there must be no side-effects in this conversion (the conversion may fail, though)
		virtual boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd) = 0;
		virtual boost::intrusive_ptr<Value> indexedByString(boost::intrusive_ptr<STRING> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd) = 0;
		virtual bool canBeIndexedByBigInt() { return false; }
		virtual bool canBeIndexedByString() { return false; }
    // inline -- defined in Interpreter.H
		static boost::intrusive_ptr<RINGELEM> from(const RingElem &theRingElem);
		static boost::intrusive_ptr<RING> from(const ring &theRing);
		static boost::intrusive_ptr<RAT> from(const BigRat &qq);
		static boost::intrusive_ptr<INT> from(const BigInt &N);
		static boost::intrusive_ptr<INT> from(const long n);
//		static boost::intrusive_ptr<INT> from(const size_t n);
//		static boost::intrusive_ptr<INT> from(const CoCoA::MachineInt &n);
		template<typename T> static boost::intrusive_ptr<LIST> from(const std::vector<T> &v);
		static boost::intrusive_ptr<MAT> from(const matrix &theMatrix);
		static boost::intrusive_ptr<MAT> from(ConstMatrixView theMatrix);
		static boost::intrusive_ptr<BOOL> from(const bool b);
		static boost::intrusive_ptr<RINGHOM> from(const RingHom &theRingHom);
		static boost::intrusive_ptr<STRING> from(const std::string &str);
		template <typename T>
			static inline boost::intrusive_ptr<RightValue> simplifiedValueFrom(const T &v) {
				return Value::from(v);
			}
		static boost::intrusive_ptr<IDEAL> from(const ideal &theIdeal);
		static boost::intrusive_ptr<MODULE> from(const module &theModule);
		static boost::intrusive_ptr<MODULEELEM> from(const ModuleElem &theModuleElem);
    // not inline -- defined in Interpreter.C
		static boost::intrusive_ptr<JBMillValue> from(const JBMill &theJBMill);
    //    template <typename T> ???
		static boost::intrusive_ptr<RECORD> from(const factorization<BigInt> &f);
		static boost::intrusive_ptr<RECORD> from(const factorization<RingElem> &f);
		static boost::intrusive_ptr<RECORD> from(const HPSeries &s);
		static boost::intrusive_ptr<LIST> from(const degree &d);
		static boost::intrusive_ptr<LIST> from(const QuasiPoly &f);
	};

	class VariableName;

	class LeftValue : public Value {
	public:
		virtual boost::intrusive_ptr<VariableName> getBase() = 0;
		virtual bool assignmentNeedsOwnership() const = 0;
		virtual void assign(boost::intrusive_ptr<RightValue> value, const LexerNS::CharPointer &valueExpBegin, const LexerNS::CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv) = 0;
                virtual void checkAssignable(RuntimeEnvironment *runtimeEnv) const = 0; // see redmine 435
		void obtainOwnership();
		boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
		boost::intrusive_ptr<Value> indexedByString(boost::intrusive_ptr<STRING> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
		bool canBeIndexedByBigInt() { return true; }
		bool canBeIndexedByString() { return true; }
	};

	struct SnapshotFrame {
		boost::intrusive_ptr<const AST::InvocationExpression> invocationExp;
		boost::intrusive_ptr<const AST::ParsedObject> block;
		SnapshotFrame(const Frame &frame);
	};

	class RuntimeException : public LexerNS::ExceptionWithSourcePosition {
	public:
		RuntimeException(const std::string &reason, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to, const Frame * const = nullptr);
		RuntimeException(const std::string &reason, boost::intrusive_ptr<const AST::ParsedObject> po, const Frame * const = nullptr);
		RuntimeException(const std::string &reason, const LexerNS::Token &token, const Frame * const = nullptr);
		std::vector<SnapshotFrame> snapshot;
          ~RuntimeException() throw() {}
	};

	class CannotBeIndexedException : public RuntimeException
        {
          static std::string errorMessage(const std::string &objectType)
            {
              return "The type "+objectType+" cannot be indexed, but there is an index";
            }
	public:
          CannotBeIndexedException(const std::string &objectType, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) :
              RuntimeException(errorMessage(objectType), from, to)
            {}
        };

	class WrongTypeException : public RuntimeException {
		static std::string errorMessage(const std::string &expectedType, const std::string &foundType) {
			return "Expecting type "+expectedType+", but found type "+foundType;
		}
	public:
		const std::string expectedType, foundType;
		WrongTypeException(const std::string &expectedType, const std::string &foundType, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) :
			RuntimeException(errorMessage(expectedType, foundType), from, to),
			expectedType(expectedType),
			foundType(foundType)
		{}
		WrongTypeException(const std::string &expectedType, const std::string &foundType,  boost::intrusive_ptr<const AST::ParsedObject> po) :
			RuntimeException(errorMessage(expectedType, foundType), po),
			expectedType(expectedType),
			foundType(foundType)
		{}
		WrongTypeException(const std::string &expectedType, const std::string &foundType, const LexerNS::Token &token) :
			RuntimeException(errorMessage(expectedType, foundType), token),
			expectedType(expectedType),
			foundType(foundType)
		{}
		~WrongTypeException() throw () {}
	};

	class VoidIsNotValueException : public RuntimeException {
	public:
		VoidIsNotValueException(boost::intrusive_ptr<const AST::ParsedObject> po) :
			RuntimeException("This is not a value", po)
		{}
		~VoidIsNotValueException() throw () {}
	};

	class RightValue : public Value {
	public:
		virtual boost::intrusive_ptr<RightValue> untagged() { return this; }
		virtual boost::intrusive_ptr<const RightValue> untagged() const { return this; }
		boost::intrusive_ptr<RightValue> asRightValue() { return this; }
		virtual boost::intrusive_ptr<RightValue> unaryMinus(const LexerNS::CharPointer &opPosition, RuntimeEnvironment * const runtimeEnv);
		virtual boost::intrusive_ptr<RightValue> clone() = 0;
		virtual bool needsToBeCopiedBeforeChanges() const = 0;
		virtual boost::intrusive_ptr<TYPE> getType() const = 0;
		virtual void describe(boost::intrusive_ptr<OSTREAM> out) const;
		boost::intrusive_ptr<Value> indexedByBigInt(boost::intrusive_ptr<INT> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
		boost::intrusive_ptr<Value> indexedByString(boost::intrusive_ptr<STRING> index, const LexerNS::CharPointer &targetExpBegin, const LexerNS::CharPointer &targetExpEnd, const LexerNS::CharPointer &indexExpBegin, const LexerNS::CharPointer &indexExpEnd);
	#ifdef C5IDE
		virtual void initTreeWidget(QTreeWidgetItem *) const;
	#endif // #ifdef C5IDE
	};

	class ImmutableRightValue : public RightValue {
	public:
		boost::intrusive_ptr<RightValue> clone() { return this; }
		bool needsToBeCopiedBeforeChanges() const { return false; }
	};

	class VoidValue : public ImmutableRightValue {
		VoidValue() {}
	public:
		static boost::intrusive_ptr<TYPE> type;
		std::ostream &dumpAsString(std::ostream &out) const;
		static boost::intrusive_ptr<VoidValue> theInstance;
		boost::intrusive_ptr<TYPE> getType() const;
		virtual ~VoidValue() /*throw ()*/ {}
	#ifdef C5IDE
		void initTreeWidget(QTreeWidgetItem *) const;
	#endif // #ifdef C5IDE
	};

	class INT : public ImmutableRightValue {
	public:
		static boost::intrusive_ptr<TYPE> type;
		explicit INT(const CoCoA::BigInt &theBigInt) : theBigInt(theBigInt) {}
		explicit INT(const CoCoA::MachineInt &n) : theBigInt(BigInt(n)) {}
		static boost::intrusive_ptr<INT> zero, one, minusOne;
		std::ostream &dumpAsString(std::ostream &out) const;
		const CoCoA::BigInt theBigInt;
		boost::intrusive_ptr<RightValue> unaryMinus(const LexerNS::CharPointer &opPosition, RuntimeEnvironment * const runtimeEnv);
		boost::intrusive_ptr<TYPE> getType() const;
		void describe(boost::intrusive_ptr<OSTREAM> out) const;
		static inline boost::intrusive_ptr<INT> fromInt(const int i) {
			switch(i) {
			case 0: return zero;
			case 1: return one;
			case -1: return minusOne;
			default:
				return new INT(BigInt(i));
			}
		}
	};

	class RAT : public ImmutableRightValue {
	public:
		static boost::intrusive_ptr<TYPE> type;
		explicit RAT(const CoCoA::BigRat &theBigRat) : theBigRat(theBigRat) {}
		std::ostream &dumpAsString(std::ostream &out) const;
		const CoCoA::BigRat theBigRat;
		boost::intrusive_ptr<RightValue> unaryMinus(const LexerNS::CharPointer &opPosition, RuntimeEnvironment * const runtimeEnv);
		boost::intrusive_ptr<TYPE> getType() const;
	};

	class TYPE : public ImmutableRightValue {
		static std::map<std::string, boost::intrusive_ptr<TYPE> > taggedTypes;
		friend class InterpreterNS::RuntimeEnvironment;
	public:
		static const int NUM_OF_DISPATCHABLE_TYPES = 24;
		static const int DISPATCH_INDEX_OF_BIGINT = 0;
		static const int DISPATCH_INDEX_OF_BIGRAT = 1;
		static const int DISPATCH_INDEX_OF_RECORD = 2;
		static const int DISPATCH_INDEX_OF_BOOL = 3;
		static const int DISPATCH_INDEX_OF_VOID = 4;
		static const int DISPATCH_INDEX_OF_STRING = 5;
		static const int DISPATCH_INDEX_OF_LIST = 6;
		static const int DISPATCH_INDEX_OF_FUNCTION = 7;
		static const int DISPATCH_INDEX_OF_TYPE = 8;
		static const int DISPATCH_INDEX_OF_ERROR_MESG = 9;
		static const int DISPATCH_INDEX_OF_OSTREAM = 10;
		static const int DISPATCH_INDEX_OF_ISTREAM = 11;
//JAA 20140901		static const int DISPATCH_INDEX_OF_ZMOD = 11;
		static const int DISPATCH_INDEX_OF_RINGELEM = 12;
		static const int DISPATCH_INDEX_OF_RATFUN = 13;
		static const int DISPATCH_INDEX_OF_MODULEELEM = 14;
		static const int DISPATCH_INDEX_OF_IDEAL = 15;
		static const int DISPATCH_INDEX_OF_MODULE = 16;
		static const int DISPATCH_INDEX_OF_MAT = 17;
		static const int DISPATCH_INDEX_OF_RING = 18;
		static const int DISPATCH_INDEX_OF_PACKAGE = 19;
		static const int DISPATCH_INDEX_OF_INTMAP = 20;
		static const int DISPATCH_INDEX_OF_RINGHOM = 21;
		static const int DISPATCH_INDEX_OF_MATRIXROW = 22;
		static const int DISPATCH_INDEX_OF_JBMILL = 24;
		TYPE(const std::string name, int dispatchIndex, const boost::intrusive_ptr<TYPE> parent = 0) :
			name(name),
			parent(parent),
			dispatchIndex(dispatchIndex)
		{}
		static boost::intrusive_ptr<TYPE> type;
		const std::string name;
		const boost::intrusive_ptr<TYPE> parent;
		const int dispatchIndex;
		std::ostream &dumpAsString(std::ostream &out) const;
		boost::intrusive_ptr<TYPE> getType() const;
		static boost::intrusive_ptr<TYPE> tagType(const std::string &tag);
		virtual ~TYPE() /*throw ()*/ {}
		bool isProperSubtypeOf(int otherDI) const;
		static boost::intrusive_ptr<const RightValue> convert(boost::intrusive_ptr<const RightValue> from, boost::intrusive_ptr<const RightValue> to);
	#ifdef C5IDE
		void initTreeWidget(QTreeWidgetItem *) const;
	#endif // #ifdef C5IDE
	};

} // namespace InterpreterNS

namespace ParserNS {
	class Parser;
} // namespace ParserNS

namespace AST {

class StaticEnv : public LexerNS::ReferenceCountedObject {
	InterpreterNS::RuntimeEnvironment *runtimeEnvironment;
#ifdef C5IDE
	friend class IDE::Debugger;
#endif // C5IDE
public:
	static const int IN_THE_WILD = -2; // non-imported, (presumably) top-level
	static const int TOP_LEVEL = -1;
	class VarData {
        public: // all members are public
		int depth;
		int index; // can be in the var-slots (when !isCapturedValue) or in the capturedValue-slots (when isCapturedValue)
		bool isIterationVar;
		bool isCapturedValue;
		bool isImplictlyImported;
		VarData(int depth, int index, bool isIterationVar, bool isCapturedValue, bool isImplictlyImported) :
			depth(depth),
			index(index),
			isIterationVar(isIterationVar),
			isCapturedValue(isCapturedValue),
			isImplictlyImported(isImplictlyImported)
		{
			assert(!(this->isCapturedValue && this->isIterationVar));
			assert(!this->isIterationVar || this->index==0);
		}
		VarData() :
			depth(IN_THE_WILD),
			index(-1),
			isIterationVar(false),
			isCapturedValue(false),
			isImplictlyImported(false)
		{}
		InterpreterNS::Frame *tryToFindFrame(InterpreterNS::RuntimeEnvironment *runtimeEnv) const;
#ifdef C5IDE
		const InterpreterNS::Frame *debuggerTryToFindFrame(const InterpreterNS::Frame *topLevelFrame, const InterpreterNS::Frame *currentFrame) const;
#endif // #ifdef C5IDE
	};
	const boost::intrusive_ptr<StaticEnv> parent;
	StaticEnv(InterpreterNS::RuntimeEnvironment *runtimeEnvironment) :
		runtimeEnvironment(runtimeEnvironment)
	{
		assert(this->runtimeEnvironment);
	}
	explicit StaticEnv(boost::intrusive_ptr<StaticEnv> parent) :
		runtimeEnvironment(0),
		parent(parent)
	{}
	void add(const std::string &identifier, int depth, int slot, bool isIterationVar, bool isCapturedValue, bool isImplictlyImported);
	void collectSimilarlyNamedIdentifiers(const std::string &id, int arity, std::vector<std::string> &NearMatches, bool &thereIsAnExactMatch) const;
	VarData lookup(const std::string &identifier) const;
private:
	typedef std::map<std::string, VarData> IdMap;
	IdMap identifierMap;
};

class AliasEnv : public LexerNS::ReferenceCountedObject {
public:
	explicit AliasEnv(boost::intrusive_ptr<AliasEnv> parent) :
		parent(parent)
	{}
	void add(const std::string &from, const std::string &to);
	bool lookup(const std::string &identifier, std::string &alias) const;
	const boost::intrusive_ptr<AliasEnv> parent;
	void dump(boost::intrusive_ptr<InterpreterNS::OSTREAM> out) const;
private:
	typedef std::map<std::string, std::string> AliasMap;
	AliasMap aliases;
};

class ParsedObjectVisitor;

class ParsedObject : public LexerNS::ReferenceCountedObject {
protected:
	LexerNS::CharPointer beginSourcePosition, endSourcePosition;
public:
	ParsedObject(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition) :
		beginSourcePosition(beginSourcePosition),
		endSourcePosition(endSourcePosition)
	{
		assert(*beginSourcePosition);
		assert(*endSourcePosition);
	}
	boost::intrusive_ptr<StaticEnv> staticEnv;
	inline const LexerNS::CharPointer &getBegin() const { return this->beginSourcePosition; }
	inline const LexerNS::CharPointer &getEnd() const { return this->endSourcePosition; }
	virtual std::ostream &dumpAsString(std::ostream &out) const = 0;
	std::string toSourceString() const { return beginSourcePosition.stringTo(endSourcePosition); }
	virtual void accept(ParsedObjectVisitor *) = 0;
	virtual ~ParsedObject() {}
#ifdef C5IDE
	virtual bool skipDebugging() const;
#endif
};

#define DECLARE_ACCEPT_METHOD void accept(ParsedObjectVisitor *v);

std::ostream &operator<<(std::ostream &out, const ParsedObject * const parsedObject);

class Expression : public ParsedObject {
	virtual boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const = 0;
public:
	Expression(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition) :
		ParsedObject(beginSourcePosition, endSourcePosition)
	{}
	virtual bool isLeftValue() const;
#ifdef C5IDE
	boost::intrusive_ptr<InterpreterNS::Value> eval(InterpreterNS::RuntimeEnvironment *re) const;
#else // #ifdef C5IDE
	inline boost::intrusive_ptr<InterpreterNS::Value> eval(InterpreterNS::RuntimeEnvironment *re) const {
		return this->implEval(re);
	}
#endif // #ifdef C5IDE
	template <typename T> boost::intrusive_ptr<T> evalAs(InterpreterNS::RuntimeEnvironment *) const;
	inline boost::intrusive_ptr<InterpreterNS::RightValue> evalAsRightValueVoidIsOk(InterpreterNS::RuntimeEnvironment *runtimeEnv) const  {
		return this->eval(runtimeEnv)->asRightValue();
	}
};

template <> inline boost::intrusive_ptr<InterpreterNS::RightValue> Expression::evalAs<InterpreterNS::RightValue>(InterpreterNS::RuntimeEnvironment * const runtimeEnv) const {
	const boost::intrusive_ptr<InterpreterNS::RightValue> rv = this->eval(runtimeEnv)->asRightValue();
	if (boost::dynamic_pointer_cast<InterpreterNS::VoidValue>(rv))
		throw InterpreterNS::VoidIsNotValueException(this);
	return rv;
}

template <typename T> boost::intrusive_ptr<T> Expression::evalAs(InterpreterNS::RuntimeEnvironment *runtimeEnv) const {
	const boost::intrusive_ptr<InterpreterNS::RightValue> rv = this->evalAs<InterpreterNS::RightValue>(runtimeEnv);
	const boost::intrusive_ptr<T> v = boost::dynamic_pointer_cast<T>(rv);
	if (!v)
		throw InterpreterNS::WrongTypeException(T::type->name, rv->getType()->name,	this);
	return v;
}

template <> inline boost::intrusive_ptr<InterpreterNS::RAT> Expression::evalAs<InterpreterNS::RAT>(InterpreterNS::RuntimeEnvironment *runtimeEnv) const {
	const boost::intrusive_ptr<InterpreterNS::RightValue> v = this->evalAs<InterpreterNS::RightValue>(runtimeEnv);
	const boost::intrusive_ptr<InterpreterNS::RAT> qq = boost::dynamic_pointer_cast<InterpreterNS::RAT>(v);
	if (!qq) {
		if (const boost::intrusive_ptr<InterpreterNS::INT> N = boost::dynamic_pointer_cast<InterpreterNS::INT>(v))
                  return new InterpreterNS::RAT(BigRat(N->theBigInt,1));
		throw InterpreterNS::WrongTypeException(InterpreterNS::RAT::type->name, v->getType()->name, this);
	}
	return qq;
}

template <> inline boost::intrusive_ptr<InterpreterNS::LeftValue> Expression::evalAs<InterpreterNS::LeftValue>(InterpreterNS::RuntimeEnvironment * const runtimeEnv) const {
	const boost::intrusive_ptr<InterpreterNS::Value> v = this->eval(runtimeEnv);
	if (const boost::intrusive_ptr<InterpreterNS::LeftValue> lv = boost::dynamic_pointer_cast<InterpreterNS::LeftValue>(v))
		return lv;
	throw InterpreterNS::RuntimeException("Left-value expected", this);
}

class BinaryExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> leftExp, rightExp;
	const LexerNS::CharPointer beginOperatorSourcePosition, endOperatorSourcePosition;
	LexerNS::TokenType operatorTT;
	std::ostream &dumpAsString(std::ostream &out) const;
	BinaryExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const LexerNS::TokenType operatorTT, const boost::intrusive_ptr<Expression> rightExp) :
		Expression(leftExp->getBegin(), rightExp->getEnd()),
		leftExp(leftExp),
		rightExp(rightExp),
		beginOperatorSourcePosition(beginOperatorSourcePosition),
		endOperatorSourcePosition(endOperatorSourcePosition),
		operatorTT(operatorTT)
	{
		assert(leftExp);
		assert(rightExp);
		assert(beginOperatorSourcePosition.stringTo(endOperatorSourcePosition).length()>0);
	}
};

class UnaryExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> exp;
	const std::string op;
	std::ostream &dumpAsString(std::ostream &out) const;
	UnaryExpression(const LexerNS::CharPointer &beginSourcePosition, const std::string &op, const boost::intrusive_ptr<Expression> exp) :
		Expression(beginSourcePosition, exp->getEnd()),
		exp(exp),
		op(op)
	{
		assert(exp);
		assert(op.length()>0);
	}
};

class LiteralExpression : public Expression {
public:
	std::ostream &dumpAsString(std::ostream &out) const;
	LiteralExpression(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition) :
		Expression(beginSourcePosition, endSourcePosition)
	{}
};

class IntLiteral : public LiteralExpression {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<InterpreterNS::INT> theBigInt;
	IntLiteral(const LexerNS::Token &t);
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class FloatLiteral : public LiteralExpression {
	static CoCoA::BigRat buildTheRational(const LexerNS::Token &t);
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<InterpreterNS::RightValue> theBigRat; // just RightValue (not RAT) because it might be simplified to a INT
	FloatLiteral(const LexerNS::Token &t);
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class StringLiteral : public LiteralExpression {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<InterpreterNS::STRING> theString;
	StringLiteral(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition, const std::string &unescapedString);
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class BoolLiteral : public LiteralExpression {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<InterpreterNS::BOOL> theBool;
	BoolLiteral(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition, LexerNS::TokenType tokenType);
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class IdInExpSuchThatExp : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	const LexerNS::Token tokIdentifier;
	const std::string identifier;
	const boost::intrusive_ptr<Expression> exp1, exp2;
	IdInExpSuchThatExp(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition, const LexerNS::Token &tokIdentifier, const boost::intrusive_ptr<Expression> exp1, const boost::intrusive_ptr<Expression> exp2) :
		Expression(beginSourcePosition, endSourcePosition),
		tokIdentifier(tokIdentifier),
		identifier(tokIdentifier.lexeme()),
		exp1(exp1),
		exp2(exp2)
	{
		assert(identifier.length()>0);
		assert(exp1);
		assert(exp2);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class ExpSuchThatIdInExpAndExp : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	const LexerNS::Token tokIdentifier;
	const std::string identifier;
	const boost::intrusive_ptr<Expression> exp1, exp2, optionalExp3;
	ExpSuchThatIdInExpAndExp(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition, const boost::intrusive_ptr<Expression> exp1, const LexerNS::Token &tokIdentifier, const boost::intrusive_ptr<Expression> exp2, const boost::intrusive_ptr<Expression> optionalExp3) :
		Expression(beginSourcePosition, endSourcePosition),
		tokIdentifier(tokIdentifier),
		identifier(tokIdentifier.lexeme()),
		exp1(exp1),
		exp2(exp2),
		optionalExp3(optionalExp3)
	{
		assert(identifier.length()>0);
		assert(exp1);
		assert(exp2);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class ListExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::vector<boost::intrusive_ptr<Expression> > exps;
	ListExpression(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition, const std::vector<boost::intrusive_ptr<Expression> > &exps) :
		Expression(beginSourcePosition, endSourcePosition),
		exps(exps)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class UnaryMinusExpression : public UnaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	UnaryMinusExpression(const LexerNS::CharPointer &beginSourcePosition, const boost::intrusive_ptr<Expression> exp) :
		UnaryExpression(beginSourcePosition, "-", exp)
	{
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class UnaryPlusExpression : public UnaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	UnaryPlusExpression(const LexerNS::CharPointer &beginSourcePosition, const boost::intrusive_ptr<Expression> exp) :
		UnaryExpression(beginSourcePosition, "+", exp)
	{
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class IndexedAccessExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> targetExp;
	const std::vector<boost::intrusive_ptr<Expression> > indices;
	std::ostream &dumpAsString(std::ostream &out) const;
	bool isLeftValue() const;
	IndexedAccessExpression(const boost::intrusive_ptr<Expression> targetExp, const std::vector<boost::intrusive_ptr<Expression> > &indices, const LexerNS::CharPointer &closedBracketPosition) :
		Expression(targetExp->getBegin(), closedBracketPosition),
		targetExp(targetExp),
		indices(indices)
	{
		assert(targetExp);
		assert(indices.size()>0);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class Argument {
public:
	bool byRef;
	boost::intrusive_ptr<Expression> exp;
	bool synthetized;
	Argument(const bool byRef, const boost::intrusive_ptr<Expression> exp, const bool synthetized=false) :
		byRef(byRef),
		exp(exp),
		synthetized(synthetized)
	{}
};

class InvocationExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	boost::intrusive_ptr<Expression> targetExp;
	const boost::shared_ptr<const LexerNS::Token> tokEllipsis;
	const bool ellipsis;
	std::vector<Argument> args;
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string packageName;
	InvocationExpression(const boost::intrusive_ptr<Expression> targetExp, const boost::shared_ptr<const LexerNS::Token> tokEllipsis, const std::vector<Argument> &args, const LexerNS::CharPointer &closedParenthesisPosition, const std::string packageName);
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
	int checkNumberOfArgs(int n) const;
	int checkNumberOfArgs(int nMin, int nMax) const;
#ifdef C5IDE
	bool skipDebugging() const;
#endif
};

class CartesianProductExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const std::vector<boost::intrusive_ptr<Expression> > operands;
	std::ostream &dumpAsString(std::ostream &out) const;
	CartesianProductExpression(const std::vector<boost::intrusive_ptr<Expression> > &operands) :
		Expression(operands.front()->getBegin(), operands.back()->getEnd()),
		operands(operands)
	{
		assert(operands.size()>0);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class SummationExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const std::vector<boost::intrusive_ptr<Expression> > operands;
	const std::vector<LexerNS::Token> operators;
	std::ostream &dumpAsString(std::ostream &out) const;
	SummationExpression(const std::vector<boost::intrusive_ptr<Expression> > &operands, const std::vector<LexerNS::Token> &operators) :
		Expression(operands.front()->getBegin(), operands.back()->getEnd()),
		operands(operands),
		operators(operators)
	{
		assert(operands.size()>2);
		assert(operators.size()==(operands.size()-1));
#ifndef NDEBUG
		for (const LexerNS::Token &t: operators)
                {
                  LexerNS::TokenType tt = t.getType();
                  assert(tt==LexerNS::TT_PLUS || tt==LexerNS::TT_MINUS);
		}
#endif // #ifndef NDEBUG
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class MultiplicationExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const std::vector<boost::intrusive_ptr<Expression> > operands;
	const std::vector<LexerNS::Token> operators;
	std::ostream &dumpAsString(std::ostream &out) const;
	MultiplicationExpression(const std::vector<boost::intrusive_ptr<Expression> > &operands, const std::vector<LexerNS::Token> &operators) :
		Expression(operands.front()->getBegin(), operands.back()->getEnd()),
		operands(operands),
		operators(operators)
	{
		assert(operands.size()>2);
		assert(operators.size()==(operands.size()-1));
#ifndef NDEBUG
		for (const LexerNS::Token &t: operators)
                {
                  LexerNS::TokenType tt = t.getType();
                  assert(tt==LexerNS::TT_STAR || tt==LexerNS::TT_SLASH ||
                         tt==LexerNS::TT_MOD || tt==LexerNS::TT_COLON);
		}
#endif // #ifndef NDEBUG
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class OrExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	OrExpression(const boost::intrusive_ptr<Expression> leftExp,const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_OR, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class PowerExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	PowerExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_POWER, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class DotDotExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	DotDotExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_DOTDOT, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class SumExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	SumExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_PLUS, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class SubtractionExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	SubtractionExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_MINUS, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class ProductExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	ProductExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_STAR, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class DivisionExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	DivisionExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_SLASH, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class ColonExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	ColonExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_COLON, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class ModuloExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	ModuloExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_MOD, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class AndExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	AndExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_AND, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class EqualExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	EqualExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_EQUAL, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class NotEqualExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	NotEqualExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_NOTEQUAL, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class IsInExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	IsInExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_ISIN, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class LessThanExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	LessThanExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_LT, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class LessOrEqualExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	LessOrEqualExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_LE, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class GreaterThanExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	GreaterThanExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_GT, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class GreaterOrEqualExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	GreaterOrEqualExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_GE, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class ScopedExpression : public BinaryExpression {
public:
	DECLARE_ACCEPT_METHOD
	ScopedExpression(const boost::intrusive_ptr<Expression> leftExp, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition, const boost::intrusive_ptr<Expression> rightExp) :
		BinaryExpression(leftExp, beginOperatorSourcePosition, endOperatorSourcePosition, LexerNS::TT_COLONCOLON, rightExp)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class Statement : public ParsedObject {
	virtual void implExecute(InterpreterNS::RuntimeEnvironment *) = 0;
public:
	Statement(const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition) :
		ParsedObject(beginSourcePosition, endSourcePosition)
	{}
  virtual bool IamEmptyStmt() const { return false; }
#ifdef C5IDE
	void execute(InterpreterNS::RuntimeEnvironment *re);
#else // #ifdef C5IDE
	inline void execute(InterpreterNS::RuntimeEnvironment *re) {
		this->implExecute(re);
	}
#endif // #ifdef C5IDE
};

class Statements : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	Statements(const LexerNS::Token &tokBeforeStatements, const std::vector<boost::intrusive_ptr<Statement> > &statements) :
			Statement(statements.size()==0 ? tokBeforeStatements.getEnd() : statements.front()->getBegin(),
					  statements.size()==0 ? tokBeforeStatements.getEnd() : statements.back()->getEnd()),
			statements(statements)
		{
		}
	std::ostream &dumpAsString(std::ostream &out) const;
	void implExecute(InterpreterNS::RuntimeEnvironment *);
	const std::vector<boost::intrusive_ptr<Statement> > statements;
#ifdef C5IDE
	bool skipDebugging() const;
#endif
};

class EvalStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> exp;
	std::ostream &dumpAsString(std::ostream &out) const;
	EvalStatement(const boost::intrusive_ptr<Expression> exp, const LexerNS::CharPointer &semicolonPosition) :
		Statement(exp->getBegin(), semicolonPosition),
		exp(exp)
	{
		assert(exp);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
#ifdef C5IDE
	bool skipDebugging() const;
#endif
};

class ReturnStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> exp;
	std::ostream &dumpAsString(std::ostream &out) const;
	ReturnStatement(const LexerNS::CharPointer &returnPosition, const boost::intrusive_ptr<Expression> exp, const LexerNS::CharPointer &semicolonPosition) :
		Statement(returnPosition, semicolonPosition),
		exp(exp)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class SourceStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> exp;
	const LexerNS::TokenType ttype;
	std::ostream &dumpAsString(std::ostream &out) const;
	SourceStatement(const LexerNS::Token &tokSource, const boost::intrusive_ptr<Expression> exp, const LexerNS::CharPointer &semicolonPosition) :
		Statement(tokSource.getBegin(), semicolonPosition),
		exp(exp),
		ttype(tokSource.getType())
	{
		assert(exp);
		assert(this->ttype==LexerNS::TT_SOURCE || this->ttype==LexerNS::TT_SOURCE_AS_LSHIFT || this->ttype==LexerNS::TT_LOAD);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class SourceRegionStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> expFromLine;
	const boost::intrusive_ptr<Expression> expFromChar;
	const boost::intrusive_ptr<Expression> expToLine;
	const boost::intrusive_ptr<Expression> expToChar;
	const boost::intrusive_ptr<Expression> expFileName;
	const LexerNS::TokenType ttype;
	std::ostream &dumpAsString(std::ostream &out) const;
	SourceRegionStatement(const LexerNS::Token &tokSource, const boost::intrusive_ptr<Expression> exp1, const boost::intrusive_ptr<Expression> exp2, const boost::intrusive_ptr<Expression> exp3, const boost::intrusive_ptr<Expression> exp4, const boost::intrusive_ptr<Expression> exp5, const LexerNS::CharPointer &semicolonPosition) :
		Statement(tokSource.getBegin(), semicolonPosition),
		expFromLine(exp1),
		expFromChar(exp2),
		expToLine(exp3),
		expToChar(exp4),
		expFileName(exp5),
		ttype(tokSource.getType())
	{
		assert(expFromLine);
		assert(expFromChar);
		assert(expToLine);
		assert(expToChar);
		assert(expFileName);
		assert(this->ttype==LexerNS::TT_SOURCEREGION);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class DescribeStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> exp;
	std::ostream &dumpAsString(std::ostream &out) const;
	DescribeStatement(const LexerNS::CharPointer &describePosition, const boost::intrusive_ptr<Expression> exp, const LexerNS::CharPointer &semicolonPosition) :
		Statement(describePosition, semicolonPosition),
		exp(exp)
	{
		assert(exp);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class EmptyStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	EmptyStatement(const LexerNS::CharPointer &semicolonPosition) :
		Statement(semicolonPosition, semicolonPosition)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
  virtual bool IamEmptyStmt() const { return true; }
#ifdef C5IDE
	virtual bool skipDebugging() const;
#endif
};

class BreakStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string label;
	BreakStatement(const LexerNS::CharPointer &breakPosition, const LexerNS::CharPointer &semicolonPosition, const std::string &label) :
		Statement(breakPosition, semicolonPosition),
		label(label)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class ContinueStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string label;
	ContinueStatement(const LexerNS::CharPointer &continuePosition, const LexerNS::CharPointer &semicolonPosition, const std::string &label) :
		Statement(continuePosition, semicolonPosition),
		label(label)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class SkipStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	SkipStatement(const LexerNS::CharPointer &skipPosition, const LexerNS::CharPointer &semicolonPosition) :
		Statement(skipPosition, semicolonPosition)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
//???  virtual bool IamEmptyStmt() const { return true; }
};

class Identifier;

class ProtectStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<Identifier> expId;
	boost::intrusive_ptr<Expression> optExp;
	ProtectStatement(const LexerNS::CharPointer &protectPosition, boost::intrusive_ptr<Identifier> expId, boost::intrusive_ptr<Expression> optExp, const LexerNS::CharPointer &semicolonPosition) :
		Statement(protectPosition, semicolonPosition),
		expId(expId),
		optExp(optExp)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class UnprotectStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	boost::intrusive_ptr<Identifier> expId;
	UnprotectStatement(const LexerNS::CharPointer &unprotectPosition, boost::intrusive_ptr<Identifier> expId, const LexerNS::CharPointer &semicolonPosition) :
		Statement(unprotectPosition, semicolonPosition),
		expId(expId)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class CiaoOrQuitStatement : public Statement {
public:
	const LexerNS::TokenType ttype;
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	CiaoOrQuitStatement(const LexerNS::Token &tokCiaoOrQuit, const LexerNS::CharPointer &semicolonPosition) :
		Statement(tokCiaoOrQuit.getBegin(), semicolonPosition),
		ttype(tokCiaoOrQuit.getType())
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class AssignmentStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> leftExp, rightExp;
	std::ostream &dumpAsString(std::ostream &out) const;
	AssignmentStatement(const boost::intrusive_ptr<Expression> leftExp, const boost::intrusive_ptr<Expression> rightExp, const LexerNS::CharPointer &semicolonPosition) :
		Statement(leftExp->getBegin(), semicolonPosition),
		leftExp(leftExp),
		rightExp(rightExp)
	{
		assert(leftExp);
		assert(rightExp);
		// note: we cannot assert leftExp->isLeftValue() here; when leftExp is not a LV, the parser reports an error but builds the tree anyway
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class Identifier : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const LexerNS::Token tokenId;
	std::string identifier;
	StaticEnv::VarData varData;
	int arity;
	std::ostream &dumpAsString(std::ostream &out) const;
	explicit Identifier(const LexerNS::Token &tokenId) :
		Expression(tokenId.getBegin(), tokenId.getEnd()),
		tokenId(tokenId),
		identifier(tokenId.lexeme()),
		arity(-1)
	{
		// note: tokenId may be of type TT_ORD or TT_ELIM (this allows to call those specials functions)
		assert(this->identifier.length()>0);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
	bool isLeftValue() const;
	InterpreterNS::VariableSlot *getVariableSlot(InterpreterNS::RuntimeEnvironment *);
};

class IndeterminateDeclaration : public ParsedObject {
public:
	boost::intrusive_ptr<Identifier> identifier;
	const std::vector<std::pair <boost::intrusive_ptr<Expression>, boost::intrusive_ptr<Expression> > > ranges;
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	IndeterminateDeclaration(const boost::intrusive_ptr<Identifier> identifier, const std::vector<std::pair <boost::intrusive_ptr<Expression>, boost::intrusive_ptr<Expression> > > & ranges, const LexerNS::CharPointer & endSourcePosition) :
		ParsedObject(identifier->getBegin(), endSourcePosition),
		identifier(identifier),
		ranges(ranges)
	{}
};

class RingDefinition : public ParsedObject {
public:
	const boost::intrusive_ptr<Identifier> identifier;
	const boost::intrusive_ptr<Expression> optionalExp;
	const std::vector<boost::intrusive_ptr<IndeterminateDeclaration> > indeterminates;
	const LexerNS::Token maybeOrderToken; // this part is incomplete...
	const LexerNS::TokenType maybeOrderTT;
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	RingDefinition(const boost::intrusive_ptr<Identifier> identifier, const boost::intrusive_ptr<Expression> optionalExp, const std::vector<boost::intrusive_ptr<IndeterminateDeclaration> > & indeterminates, const LexerNS::Token &maybeOrderToken, const LexerNS::CharPointer &endSourcePosition) :
		ParsedObject(identifier->getBegin(), endSourcePosition),
		identifier(identifier),
		optionalExp(optionalExp),
		indeterminates(indeterminates),
		maybeOrderToken(maybeOrderToken),
		maybeOrderTT(maybeOrderToken.getType())
	{
	}
	boost::intrusive_ptr<InterpreterNS::RING> eval(InterpreterNS::RuntimeEnvironment *) const;
};

class RingAssignStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> leftExp;
	const boost::intrusive_ptr<RingDefinition> ringDef;
	std::ostream &dumpAsString(std::ostream &out) const;
	RingAssignStatement(const boost::intrusive_ptr<Expression> leftExp, const boost::intrusive_ptr<RingDefinition> ringDef, const LexerNS::CharPointer &semicolonPosition) :
		Statement(leftExp->getBegin(), semicolonPosition),
		leftExp(leftExp),
		ringDef(ringDef)
	{
		assert(leftExp);
		assert(ringDef);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class FullyQualifiedIdentifier : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const std::string pkgName, id;
	std::ostream &dumpAsString(std::ostream &out) const;
	const LexerNS::Token tokPkgname;
	const LexerNS::Token tokId;
	explicit FullyQualifiedIdentifier(const LexerNS::Token &tokPkgname, const std::string &pkgName, const LexerNS::Token &tokId) :
		Expression(tokPkgname.getBegin(), tokId.getEnd()),
		pkgName(pkgName),
		id(tokId.lexeme()),
		tokPkgname(tokPkgname),
		tokId(tokId)
	{
		assert(this->pkgName.length()>0);
		assert(this->id.length()>0);
	}
	bool isLeftValue() const;
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class FieldAccessExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> targetExp;
	const LexerNS::Token tokName;
	const std::string name;
	std::ostream &dumpAsString(std::ostream &out) const;
	bool isLeftValue() const;
	FieldAccessExpression(const boost::intrusive_ptr<Expression> targetExp, const LexerNS::Token &tokName) :
		Expression(targetExp->getBegin(), tokName.getEnd()),
		targetExp(targetExp),
		tokName(tokName),
		name(tokName.lexeme())
	{
		assert(targetExp);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

struct RecordField {
	boost::intrusive_ptr<Identifier> name;
	boost::intrusive_ptr<Expression> initExp;
};

class RecordExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	const std::vector<RecordField> fields;
	std::ostream &dumpAsString(std::ostream &out) const;
	RecordExpression(const LexerNS::CharPointer &recordPosition, const std::vector<RecordField> &fields, const LexerNS::CharPointer &closedSquarePosition) :
		Expression(recordPosition, closedSquarePosition),
		fields(fields)
	{
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class UseStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Identifier> identifier;
	const boost::intrusive_ptr<RingDefinition> ringDefinition;
	std::ostream &dumpAsString(std::ostream &out) const;
	UseStatement(const LexerNS::CharPointer &usePosition, const boost::intrusive_ptr<Identifier> identifier, const boost::intrusive_ptr<RingDefinition> ringDefinition, const LexerNS::CharPointer &semicolonPosition) :
		Statement(usePosition, semicolonPosition),
		identifier(identifier),
		ringDefinition(ringDefinition)
	{
		assert(identifier || ringDefinition);
	}
	virtual void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class TimeStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Statement> stmt;
	std::ostream &dumpAsString(std::ostream &out) const;
	TimeStatement(boost::intrusive_ptr<Statement> stmt, const LexerNS::CharPointer &beginSourcePosition) :
		Statement(beginSourcePosition, stmt->getEnd()),
		stmt(stmt)
	{
		assert(stmt);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

/*class SetStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Identifier> identifier;
	const boost::intrusive_ptr<Expression> exp;
	std::ostream &dumpAsString(std::ostream &out) const;
	SetStatement(const LexerNS::CharPointer &setPosition, const boost::intrusive_ptr<Identifier> identifier, const boost::intrusive_ptr<Expression> exp, const LexerNS::CharPointer &semicolonPosition) :
		Statement(setPosition, semicolonPosition),
		identifier(identifier),
		exp(exp)
	{
		assert(identifier);
	}
};

class UnsetStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Identifier> identifier;
	std::ostream &dumpAsString(std::ostream &out) const;
	UnsetStatement(const LexerNS::CharPointer &setPosition, const boost::intrusive_ptr<Identifier> identifier, const LexerNS::CharPointer &semicolonPosition) :
		Statement(setPosition, semicolonPosition),
		identifier(identifier)
	{
		assert(identifier);
	}
};*/

class ForStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const LexerNS::Token tokIdentifier;
	const std::string identifier;
	const boost::intrusive_ptr<Expression> beginExp;
	const boost::intrusive_ptr<Expression> endExp;
	const boost::intrusive_ptr<Expression> stepExp;
	const boost::intrusive_ptr<Statements> statements;
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string label;
	ForStatement(const LexerNS::CharPointer &forPosition, const LexerNS::Token &tokIdentifier, boost::intrusive_ptr<Expression> beginExp, boost::intrusive_ptr<Expression> endExp, boost::intrusive_ptr<Expression> stepExp, const boost::intrusive_ptr<Statements> statements, const LexerNS::CharPointer &endforPosition, const std::string &label) :
		Statement(forPosition, endforPosition),
		tokIdentifier(tokIdentifier),
		identifier(tokIdentifier.lexeme()),
		beginExp(beginExp),
		endExp(endExp),
		stepExp(stepExp),
		statements(statements),
		label(label)
	{
		assert(beginExp);
		assert(endExp);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class UsingStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Identifier> identifier;
	const boost::intrusive_ptr<Statements> statements;
	std::ostream &dumpAsString(std::ostream &out) const;
	UsingStatement(const LexerNS::CharPointer &usingPosition, const boost::intrusive_ptr<Identifier> identifier, const boost::intrusive_ptr<Statements> statements, const LexerNS::CharPointer &endusingPosition) :
		Statement(usingPosition, endusingPosition),
		identifier(identifier),
		statements(statements)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class ForeachStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const LexerNS::Token tokIdentifier;
	const std::string identifier;
	const boost::intrusive_ptr<Expression> inExp;
	const boost::intrusive_ptr<Statements> statements;
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string label;
	ForeachStatement(const LexerNS::CharPointer &foreachPosition, const LexerNS::Token &tokIdentifier, boost::intrusive_ptr<Expression> inExp, const boost::intrusive_ptr<Statements> statements, const LexerNS::CharPointer &endforPosition, const std::string &label) :
		Statement(foreachPosition, endforPosition),
		tokIdentifier(tokIdentifier),
		identifier(tokIdentifier.lexeme()),
		inExp(inExp),
		statements(statements),
		label(label)
	{
		assert(inExp);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class TryStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const LexerNS::Token tokIdentifier;
	const std::string identifier;
	const boost::intrusive_ptr<Statements> tryStatements;
	const boost::intrusive_ptr<Statements> uponErrorStatements;
	std::ostream &dumpAsString(std::ostream &out) const;
	TryStatement(const LexerNS::CharPointer &tryPosition, const boost::intrusive_ptr<Statements> tryStatements, const LexerNS::Token &tokIdentifier, const boost::intrusive_ptr<Statements> uponErrorStatements, const LexerNS::CharPointer &endTryPosition) :
		Statement(tryPosition, endTryPosition),
		tokIdentifier(tokIdentifier),
		identifier(tokIdentifier.lexeme()),
		tryStatements(tryStatements),
		uponErrorStatements(uponErrorStatements)
	{
		assert(this->tryStatements);
		assert(this->uponErrorStatements);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class RepeatUntilStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Statements> statements;
	const boost::intrusive_ptr<Expression> optExp;
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string label;
	RepeatUntilStatement(const LexerNS::CharPointer &repeatPosition, const boost::intrusive_ptr<Statements> statements, boost::intrusive_ptr<Expression> optExp, const LexerNS::CharPointer &semicolonPosition, const std::string &label) :
		Statement(repeatPosition, semicolonPosition),
		statements(statements),
		optExp(optExp),
		label(label)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class WhileStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Expression> exp;
	const boost::intrusive_ptr<Statements> statements;
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::string label;
	WhileStatement(const LexerNS::CharPointer &whilePosition, boost::intrusive_ptr<Expression> exp, const boost::intrusive_ptr<Statements> statements, const LexerNS::CharPointer &endwhilePosition, const std::string &label) :
        Statement(whilePosition, endwhilePosition),
		exp(exp),
		statements(statements),
		label(label)
	{
		assert(exp);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

struct IfBranch {
	boost::intrusive_ptr<Expression> optExp; // no exp means it's the "else" branch
	boost::intrusive_ptr<Statements> statements;
};

class IfStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const std::vector<IfBranch> branches;
	std::ostream &dumpAsString(std::ostream &out) const;
	IfStatement(const LexerNS::CharPointer &ifPosition, const std::vector<IfBranch> branches, const LexerNS::CharPointer &endifPosition) :
		Statement(ifPosition, endifPosition),
		branches(branches)
	{
		assert(branches.size()>=1);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

struct Binding {
	boost::intrusive_ptr<Identifier> identifier;
	std::string packageName;
};

class AliasStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const std::vector<Binding> bindings;
	//const bool thereIsIn;
	//const boost::intrusive_ptr<Statements> statements;
	std::ostream &dumpAsString(std::ostream &out) const;
	AliasStatement(const LexerNS::CharPointer &aliasPosition, const std::vector<Binding> &bindings,
			//const bool thereIsIn, const boost::intrusive_ptr<Statements> statements,
			const LexerNS::CharPointer &endaliasOrSemicolonPosition) :
		Statement(aliasPosition, endaliasOrSemicolonPosition),
		bindings(bindings)
		//thereIsIn(thereIsIn),
		//statements(statements)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class BlockStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Statements> statements;
	std::ostream &dumpAsString(std::ostream &out) const;
	BlockStatement(const LexerNS::CharPointer &blockPosition, const boost::intrusive_ptr<Statements> statements, const LexerNS::CharPointer &endblockPosition) :
		Statement(blockPosition, endblockPosition),
		statements(statements)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class Param {
public: // all members are public
	std::string name;
	bool byRef;
	bool opt;
	Param() : byRef(false), opt(false) {}
};

class Import {
public: // all members are public
	enum ImportType { IT_TOPLEVEL, IT_BYREF, IT_BYVALUE };
	boost::intrusive_ptr<AST::Identifier> expId;
	std::string name;
	ImportType type;
	bool implicit;
	int byValueIndex;
	Import(const boost::intrusive_ptr<AST::Identifier> expId, const std::string &name, ImportType type, bool implicit) :
		expId(expId),
		name(name),
		type(type),
		implicit(implicit),
		byValueIndex(-1)
	{
		assert(!this->implicit || (!this->expId && this->type==IT_TOPLEVEL));
		assert(this->implicit || this->expId);
		assert(!this->expId || this->expId->identifier==this->name);
	}
	std::ostream &dumpAsString(std::ostream &out) const;
};

class LambdaExpression;

class FunctionDeclaration : public ParsedObject {
public:
	DECLARE_ACCEPT_METHOD
	const bool thereAreReturnsWithExpr;
	const bool thereIsEllipsis;
	const int nMandatoryParameters;
	const std::vector<Param> params;
	std::vector<Import> imports;
	std::map<std::string, int> localNames;
	const boost::intrusive_ptr<Statements> statements;
	void setSourcePosition(ParsedObject *po) {
		this->beginSourcePosition = po->getBegin();
		this->endSourcePosition = po->getEnd();
	}
	std::string fnName;
	FunctionDeclaration(const LexerNS::CharPointer &openRoundPosition, const bool thereAreReturnsWithExpr, const bool thereIsEllipsis, const int nMandatoryParameters, const std::vector<Param> &params, const std::vector<Import> &imports, const boost::intrusive_ptr<Statements> statements) :
		ParsedObject(openRoundPosition, statements->getEnd()),
		thereAreReturnsWithExpr(thereAreReturnsWithExpr),
		thereIsEllipsis(thereIsEllipsis),
		nMandatoryParameters(nMandatoryParameters),
		params(params),
		imports(imports),
		statements(statements)
	{
#ifndef NDEBUG
		if (thereIsEllipsis) {
			assert(params.size()==0);
			assert(nMandatoryParameters==0);
		} else {
			assert(nMandatoryParameters<=static_cast<int>(params.size()));
			for(std::vector<Param>::size_type a=0; static_cast<int>(a)<nMandatoryParameters; ++a)
				assert(!params[a].opt);
			for(std::vector<Param>::size_type a=nMandatoryParameters; a<params.size(); ++a) {
				assert(params[a].opt);
				assert(!params[a].byRef);
			}
		}
#endif // #ifndef NDEBUG
	}
	std::ostream &dumpAsString(std::ostream &out) const;
};

class DefineStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const boost::intrusive_ptr<Identifier> expId;
	const std::string name;
	boost::intrusive_ptr<FunctionDeclaration> funDecl;
	std::ostream &dumpAsString(std::ostream &out) const;
	DefineStatement(const LexerNS::CharPointer &definePosition, const boost::intrusive_ptr<Identifier> expId, boost::intrusive_ptr<FunctionDeclaration> funDecl, const LexerNS::CharPointer &endDefinePosition) :
		Statement(definePosition, endDefinePosition),
		expId(expId),
		name(expId->identifier),
		funDecl(funDecl)
	{
		this->funDecl->setSourcePosition(this);
		this->funDecl->fnName = this->name;
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class LambdaExpression : public Expression {
public:
	DECLARE_ACCEPT_METHOD
	boost::intrusive_ptr<FunctionDeclaration> funDecl;
	std::ostream &dumpAsString(std::ostream &out) const;
	LambdaExpression(const LexerNS::CharPointer &lambdaPosition, boost::intrusive_ptr<FunctionDeclaration> funDecl, const LexerNS::CharPointer &endLambdaPosition) :
		Expression(lambdaPosition, endLambdaPosition),
		funDecl(funDecl)
	{
		this->funDecl->setSourcePosition(this);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class PackageStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const LexerNS::Token tokName;
	const std::string name;
	const boost::intrusive_ptr<Statements> statements;
	std::ostream &dumpAsString(std::ostream &out) const;
	const std::vector<LexerNS::Token> exportedNames;
	const std::set<std::string> memberNames;
	PackageStatement(const LexerNS::CharPointer &packagePosition, const LexerNS::Token &tokName, const boost::intrusive_ptr<Statements> statements, const LexerNS::CharPointer &endPackagePosition, const std::vector<LexerNS::Token> &exportedNames, const std::set<std::string> &memberNames) :
		Statement(packagePosition, endPackagePosition),
		tokName(tokName),
		name(tokName.lexeme()),
		statements(statements),
		exportedNames(exportedNames),
		memberNames(memberNames)
	{
		assert(statements);
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class HelpStatement : public Statement {
public:
	DECLARE_ACCEPT_METHOD
	const std::string lexeme;
	const std::string topic;
	std::ostream &dumpAsString(std::ostream &out) const;
	HelpStatement(const LexerNS::Token &tokHelp);
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class PrintStatement : public Statement {
public:
	const LexerNS::TokenType ttype;
	DECLARE_ACCEPT_METHOD
	const std::vector<boost::intrusive_ptr<Expression> > exps;
	const boost::intrusive_ptr<Expression> onExp;
	std::ostream &dumpAsString(std::ostream &out) const;
	PrintStatement(const LexerNS::Token tokPrintOrPrintLn, const std::vector<boost::intrusive_ptr<Expression> > &exps, const boost::intrusive_ptr<Expression> onExp, const LexerNS::CharPointer &endifPosition) :
		Statement(tokPrintOrPrintLn.getBegin(), endifPosition),
		ttype(tokPrintOrPrintLn.getType()),
		exps(exps),
		onExp(onExp)
	{
	}
	void implExecute(InterpreterNS::RuntimeEnvironment *);
};

class IsDefinedExp : public Expression {
public:
	const boost::intrusive_ptr<Identifier> expId;
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	IsDefinedExp(const boost::intrusive_ptr<Identifier> expId, const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition) :
		Expression(beginSourcePosition, endSourcePosition),
		expId(expId)
	{}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class BackwardCompatibleExp : public Expression {
public:
	const boost::intrusive_ptr<Expression> innerExp;
	DECLARE_ACCEPT_METHOD
	std::ostream &dumpAsString(std::ostream &out) const;
	BackwardCompatibleExp(boost::intrusive_ptr<Expression> innerExp, const LexerNS::CharPointer &beginSourcePosition, const LexerNS::CharPointer & endSourcePosition) :
		Expression(beginSourcePosition, endSourcePosition),
		innerExp(innerExp)
	{
		assert(innerExp);
	}
	boost::intrusive_ptr<InterpreterNS::Value> implEval(InterpreterNS::RuntimeEnvironment *) const;
};

class ParsedObjectVisitor {
public:
	virtual void visit(BinaryExpression &binaryExp);
	virtual void visit(UnaryExpression &unaryExp);
	virtual void visit(LiteralExpression &);
	virtual void visit(IntLiteral &);
	virtual void visit(FloatLiteral &);
	virtual void visit(StringLiteral &);
	virtual void visit(BoolLiteral &);
	virtual void visit(IdInExpSuchThatExp &idInExpSuchThatExp);
	virtual void visit(ExpSuchThatIdInExpAndExp &expSuchThatIdInExpAndExp);
	virtual void visit(ListExpression &listExpression);
	virtual void visit(InvocationExpression &invocationExp);
	#define VISIT_UNARY(CLASS) virtual void visit(CLASS &e);
		VISIT_UNARY(UnaryMinusExpression)
		VISIT_UNARY(UnaryPlusExpression)
	#undef VISIT_UNARY
	virtual void visit(IndexedAccessExpression &indexedAccessExpression);
	virtual void visit(CartesianProductExpression &cartesianProductExpression);
	virtual void visit(SummationExpression &summationExpression);
	virtual void visit(MultiplicationExpression &multiplicationExpression);
	#define VISIT_BINARY(CLASS) virtual void visit(CLASS &e);
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
	virtual void visit(Statements &stmts);
	virtual void visit(EvalStatement &evalStmt);
	virtual void visit(ReturnStatement &returnStmt);
	virtual void visit(SourceStatement &sourceStmt);
	virtual void visit(SourceRegionStatement &sourceregionStmt);
	virtual void visit(DescribeStatement &describeStmt);
	virtual void visit(EmptyStatement &);
	virtual void visit(BreakStatement &);
	virtual void visit(ContinueStatement &);
	virtual void visit(SkipStatement &);
	virtual void visit(ProtectStatement &protectStmt);
	virtual void visit(UnprotectStatement &);
	virtual void visit(CiaoOrQuitStatement &);
	virtual void visit(AssignmentStatement &assignmentStmt);
	virtual void visit(TryStatement &);
	virtual void visit(Identifier &);
	virtual void visit(FullyQualifiedIdentifier &);
	virtual void visit(IndeterminateDeclaration &);
	virtual void visit(RingDefinition &);
	virtual void visit(RingAssignStatement &);
	virtual void visit(FieldAccessExpression &fieldAccessExpression);
	virtual void visit(RecordExpression &recordExpression);
	virtual void visit(UseStatement &useStatement);
	virtual void visit(TimeStatement &timeStatement);
	//virtual void visit(SetStatement &setStatement);
	//virtual void visit(UnsetStatement &unsetStatement);
	virtual void visit(ForStatement &forStmt);
	virtual void visit(ForeachStatement &foreachStmt);
	virtual void visit(WhileStatement &whileStmt);
	virtual void visit(RepeatUntilStatement &repeatUntilStatement);
	virtual void visit(IfStatement &ifStatement);
	virtual void visit(UsingStatement &usingStmt);
	virtual void visit(AliasStatement &aliasStmt);
	virtual void visit(BlockStatement &blockStmt);
	virtual void visit(FunctionDeclaration &funDecl);
	virtual void visit(DefineStatement &defineStmt);
	virtual void visit(LambdaExpression &lambda);
	virtual void visit(PackageStatement &packageStmt);
	virtual void visit(PrintStatement &printStmt);
	virtual void visit(HelpStatement &);
	virtual void visit(BackwardCompatibleExp &backwardCompatibleExp);
	virtual void visit(IsDefinedExp &);
	virtual ~ParsedObjectVisitor() /*throw ()*/ {}
};

class DumpAsTreeVisitor : public ParsedObjectVisitor {
	std::ostream &out;
	int indentationLevel;
	void indent() const;
	void dumpIndented(const boost::intrusive_ptr<ParsedObject> o);
	void dumpSequenceOfOpsAsTree(std::vector<LexerNS::Token>::const_iterator opIt, std::vector<boost::intrusive_ptr<Expression> >::const_iterator expIt, std::vector<boost::intrusive_ptr<Expression> >::const_iterator expEnd);
public:
	explicit DumpAsTreeVisitor(std::ostream &out) :
		out(out),
		indentationLevel(0)
	{}
  using ParsedObjectVisitor::visit; // disables warnings of overloading
	void visit(Import &import);
	void visit(BinaryExpression &binaryExp);
	void visit(UnaryExpression &unaryExp);
	void visit(LiteralExpression &);
	void visit(IdInExpSuchThatExp &idInExpSuchThatExp);
	void visit(ExpSuchThatIdInExpAndExp &expSuchThatIdInExpAndExp);
	void visit(ListExpression &listExpression);
	void visit(InvocationExpression &invocationExp);
	void visit(IndexedAccessExpression &indexedAccessExpression);
	void visit(CartesianProductExpression &cartesianProductExpression);
	void visit(SummationExpression &summationExpression);
	void visit(MultiplicationExpression &multiplicationExpression);
	void visit(Statements &stmts);
	void visit(EvalStatement &evalStmt);
	void visit(ReturnStatement &returnStmt);
	void visit(SourceStatement &sourceStmt);
	void visit(SourceRegionStatement &sourceregionStmt);
	void visit(DescribeStatement &describeStmt);
	void visit(EmptyStatement &);
	void visit(BreakStatement &);
	void visit(ContinueStatement &);
	void visit(SkipStatement &);
	void visit(ProtectStatement &protectStmt);
	void visit(UnprotectStatement &);
	void visit(CiaoOrQuitStatement &);
	void visit(AssignmentStatement &assignmentStmt);
	void visit(Identifier &);
	void visit(FullyQualifiedIdentifier &);
	void visit(IndeterminateDeclaration &);
	void visit(RingDefinition &);
	void visit(RingAssignStatement &);
	void visit(FieldAccessExpression &fieldAccessExpression);
	void visit(RecordExpression &recordExpression);
	void visit(UseStatement &useStatement);
	void visit(TimeStatement &timeStatement);
	//void visit(SetStatement &setStatement);
	//void visit(UnsetStatement &unsetStatement);
	void visit(ForStatement &forStmt);
	void visit(ForeachStatement &foreachStmt);
	void visit(WhileStatement &whileStmt);
	void visit(RepeatUntilStatement &repeatUntilStatement);
	void visit(TryStatement &);
	void visit(IfStatement &ifStatement);
	void visit(UsingStatement &usingStmt);
	void visit(AliasStatement &aliasStmt);
	void visit(BlockStatement &blockStmt);
	void visit(FunctionDeclaration &funDecl);
	void visit(DefineStatement &defineStmt);
	void visit(LambdaExpression &lambda);
	void visit(PackageStatement &packageStmt);
	void visit(PrintStatement &printStmt);
	void visit(HelpStatement &);
	void visit(BackwardCompatibleExp &backwardCompatibleExp);
	void visit(IsDefinedExp &);
	~DumpAsTreeVisitor() /*throw ()*/ {}
};

} // namespace AST
} // namespace CoCoA

// Close include guard:
#endif
