//   Copyright (c)  2010  John Abbott,  Anna Bigatti
//   Orig author: 2010 Giovanni Lagorio (lagorio@disi.unige.it)
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

#include "CoCoA/BigIntOps.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"

////#include <limits>
#include <iostream>
#include <exception>
#include <cassert>
#include <cerrno>
// used just for ENOENT (after call to stat -- replace when changing to c++17)
#include <set>
#include <boost/scope_exit.hpp>
#include <boost/algorithm/string/case_conv.hpp>

//#include "CoCoA/library.H" // already included by AST->Parser->Interpreter
#include "Interpreter.H"
#include "globals.H"
#include "BuiltInFunctions.H"
#include "CoCoALibSupplement.H"
// CoCoALibSupplement needed by impls of value::from  (see around line 3220)
#include "OnlineHelp.H"

namespace CoCoA {

  // Get a helpful string from a CoCoA::exception
  std::string message_from_cocoalib(const exception& exc)
  {
    return message(exc) + "\n--> [CoCoALib fn] " + context(exc);
  }


namespace InterpreterNS {

using namespace std;
using namespace boost;
using namespace CoCoA::AST;
using namespace CoCoA::LexerNS;
using namespace CoCoA::ParserNS;

intrusive_ptr<BOOL> BOOL::trueValue(new BOOL(true)), BOOL::falseValue(new BOOL(false));
intrusive_ptr<INT> INT::zero(new INT(BigInt(0))), INT::one(new INT(BigInt(1))), INT::minusOne(new INT(BigInt(-1)));
intrusive_ptr<STRING> STRING::empty(new STRING(""));

namespace // anonymous
{
ring ringForList(const intrusive_ptr<const LIST> list) {
  bool foundBigRat = false;
  const LIST::ContainerType::size_type size = list->size();
  for (LIST::ContainerType::size_type a=0; a<size; ++a) {
    const intrusive_ptr<RightValue> e = list->getValue(a);
    if (const intrusive_ptr<RINGELEM> ringElem = dynamic_pointer_cast<RINGELEM>(e))
      return owner(ringElem->theRingElem);
    else if (dynamic_pointer_cast<RAT>(e))
      foundBigRat = true;
  }
  if (foundBigRat)
    return RingQQ();
  return RingZZ();
}


vector<RingElem> toRingElemList(const intrusive_ptr<const LIST> list, const intrusive_ptr<const Expression> exp, const ring &ring) {
  const LIST::ContainerType::size_type size = list->size();
  vector<RingElem> result;
  for (LIST::ContainerType::size_type a=0; a<size; ++a)
  {
    const intrusive_ptr<RightValue> e = list->getValue(a);
    if (const intrusive_ptr<RINGELEM> ringElem = dynamic_pointer_cast<RINGELEM>(e))
      result.push_back(ringElem->theRingElem);
    else if (const intrusive_ptr<INT> N = dynamic_pointer_cast<INT>(e))
      result.push_back(RingElem(ring, N->theBigInt));
    else if (const intrusive_ptr<RAT> qq = dynamic_pointer_cast<RAT>(e))
      result.push_back(RingElem(ring, qq->theBigRat));
    else
      throw RuntimeException("A list containing RINGELEM, INT or RAT values is required", exp);
  }
  return result;
 }

}  // anonymous namespace


vector<RingElem> RuntimeEnvironment::evalArgAsListOfRingElem(const ring &R, const Argument &arg)
{ return toRingElemList(this->evalArgAs<LIST>(arg), arg.exp, R); }


vector<RingElem> RuntimeEnvironment::evalArgAsListOfRingElem(const Argument &arg)
{
	const intrusive_ptr<const LIST> list = this->evalArgAs<LIST>(arg);
	return toRingElemList(list, arg.exp, ringForList(list));
}


vector<RingElem> RuntimeEnvironment::SpecializeToListOfRingElem(const boost::intrusive_ptr<const LIST> list, const Argument &arg)
{ return toRingElemList(list, arg.exp, ringForList(list)); }


vector<RingElem> RuntimeEnvironment::evalRVAsListOfRingElem(const ring &R, const boost::intrusive_ptr<const RightValue> RV, const Argument &arg)
{
	const intrusive_ptr<const LIST> list = boost::dynamic_pointer_cast<const LIST>(RV);
	return toRingElemList(list, arg.exp, R);
}


vector<RingElem> RuntimeEnvironment::evalRVAsListOfRingElem(const boost::intrusive_ptr<const RightValue> RV, const Argument &arg)
{
	const intrusive_ptr<const LIST> list = boost::dynamic_pointer_cast<const LIST>(RV);
	return toRingElemList(list, arg.exp, ringForList(list));
}


vector<symbol> RuntimeEnvironment::evalArgAsListOfSymbols(const Argument &arg)
{
	int which;
  auto v = this->evalArgAsT1orT2<STRING, LIST>(arg, which);
  switch (which)
  {
  case 1: return symbols(RefTo<string>(v));
  case 2:
  {
    vector<symbol> syms;
    intrusive_ptr<LIST> l = this->evalArgAs<LIST>(arg);
    LIST::ContainerType::size_type size = l->size();
    for (unsigned long i=0; i<size; ++i)
    {
      if (const boost::intrusive_ptr<STRING> s = boost::dynamic_pointer_cast<STRING>(l->getValue(i)))
      { syms.push_back(symbol(s->theString)); continue; }
      if (const boost::intrusive_ptr<RECORD> r = boost::dynamic_pointer_cast<RECORD>(l->getValue(i)))
      {
        if (!r->getFieldNoCheck("head"))
          throw RuntimeException("record must have \"head\" field (\"indices\" field is optional)", arg.exp);

        const boost::intrusive_ptr<STRING> s = boost::dynamic_pointer_cast<STRING>(r->getFieldNoCheck("head"));
        if (s == nullptr) throw RuntimeException("field \"head\" must be a string", arg.exp);

        if (!r->getFieldNoCheck("indices"))
        { syms.push_back(symbol(s->theString)); continue;}

        vector<long> indices;
        const boost::intrusive_ptr<LIST> inds = boost::dynamic_pointer_cast<LIST>(r->getFieldNoCheck("indices"));
        if (inds == nullptr) throw RuntimeException("field \"indices\" must be a list (of integers)", arg.exp);
        long tmp;
        LIST::ContainerType::size_type NumIndices = inds->size();
        for (unsigned long j=0; j<NumIndices; ++j)
        {
          if (boost::dynamic_pointer_cast<INT>(inds->getValue(j)) == nullptr)
            throw RuntimeException("Indices must be integers", arg.exp);

          if (IsConvertible(tmp, (boost::dynamic_pointer_cast<INT>(inds->getValue(j)))->theBigInt))
            indices.push_back(tmp);
          else
            throw RuntimeException("All indices must fit into a machine-integer", arg.exp);
        }
        syms.push_back(symbol(s->theString, indices));
      }
      else throw RuntimeException("List of symbols must have strings or records", arg.exp);
    };
    return syms;
  };  // LIST
  default: throw RuntimeException(ERRORMissingCode(v), arg.exp);
  }
}


long RuntimeEnvironment::evalArgAsLong(const Argument &arg) {
	const intrusive_ptr<INT> N = this->evalArgAs<INT>(arg);
  long n;
  if (!IsConvertible(n, N->theBigInt))
    throw RuntimeException("number is too big", arg.exp);
  return n;
}

intrusive_ptr<Value> LeftValue::indexedByBigInt(intrusive_ptr<INT> index, const CharPointer &targetExpBegin, const CharPointer &targetExpEnd, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	return new IntegerIndexedAccess(this, index->theBigInt, targetExpBegin, targetExpEnd, indexExpBegin, indexExpEnd);
}

intrusive_ptr<Value> LeftValue::indexedByString(intrusive_ptr<STRING> index, const CharPointer &targetExpBegin, const CharPointer &targetExpEnd, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	return new StringIndexedAccess(this, index->theString, targetExpBegin, targetExpEnd, indexExpBegin, indexExpEnd);
}

intrusive_ptr<Value> RightValue::indexedByBigInt(intrusive_ptr<INT> /* index */, const CharPointer &targetExpBegin, const CharPointer &targetExpEnd, const CharPointer & /* indexExpBegin */, const CharPointer & /* indexExpEnd */) {
	throw NonIntegerIndexableException(this, targetExpBegin, targetExpEnd);
}

intrusive_ptr<Value> RightValue::indexedByString(intrusive_ptr<STRING> /* index */, const CharPointer &targetExpBegin, const CharPointer &targetExpEnd, const CharPointer & /* indexExpBegin */, const CharPointer & /* indexExpEnd */) {
	throw NonStringIndexableException(this, targetExpBegin, targetExpEnd);
}

intrusive_ptr<Value> RECORD::indexedByString(intrusive_ptr<STRING> str, const CharPointer & /* targetExpBegin */, const CharPointer & /* targetExpEnd */, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
  const intrusive_ptr<Value> v = this->getField(str->theString, indexExpBegin, indexExpBegin);
	if (!v)
		throw FieldNotFoundException(str->theString, indexExpBegin, indexExpEnd, this);
	return v;
}

intrusive_ptr<Value> LIST::indexedByBigInt(intrusive_ptr<INT> N, const CharPointer & /* targetExpBegin */, const CharPointer & /* targetExpEnd */, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	long l;
	if (!IsConvertible(l, N->theBigInt) || l<=0 || static_cast<ContainerType::size_type>(l)>this->size())
		throw IndexOutOfRangeException(N->theBigInt, this->size(), indexExpBegin, indexExpEnd);
	return this->getValue(l-1);
}

intrusive_ptr<Value> STRING::indexedByBigInt(intrusive_ptr<INT> N, const CharPointer & /* targetExpBegin */, const CharPointer & /* targetExpEnd */, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	long l;
	if (!IsConvertible(l, N->theBigInt) || l<=0 || static_cast<string::size_type>(l)>this->theString.length())
		throw IndexOutOfRangeException(N->theBigInt, this->theString.length(), indexExpBegin, indexExpEnd);
	return new STRING(string(1, this->theString[l-1]));
}

intrusive_ptr<Value> IntMapValue::indexedByBigInt(intrusive_ptr<INT> N, const CharPointer & /* targetExpBegin */, const CharPointer & /* targetExpEnd */, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	long l;
	if (!IsConvertible(l, N->theBigInt))
		throw RuntimeException("Index is not a machine-integer", indexExpBegin, indexExpEnd);
	if (const intrusive_ptr<RightValue> result = this->getValue(l))
		return result;
	throw RuntimeException("Invalid index", indexExpBegin, indexExpEnd);
}

intrusive_ptr<Value> MAT::indexedByBigInt(intrusive_ptr<INT> N, const CharPointer & /* targetExpBegin */, const CharPointer & /* targetExpEnd */, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	long l;
	if (!IsConvertible(l, N->theBigInt))
		throw RuntimeException("Index is not a machine-integer", indexExpBegin, indexExpEnd);
	const size_t row = static_cast<size_t>(l);
	if (l>=1 && row<=this->numRows()) // note: cannot use row>=1 because size_t is unsigned
		return new MatrixRowValue(this, row-1);
	throw RuntimeException("Invalid row index", indexExpBegin, indexExpEnd);
}

intrusive_ptr<Value> MatrixRowValue::indexedByBigInt(intrusive_ptr<INT> N, const CharPointer & /* targetExpBegin */, const CharPointer & /* targetExpEnd */, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	long l;
	if (!IsConvertible(l, N->theBigInt))
		throw RuntimeException("Index is not a machine-integer", indexExpBegin, indexExpEnd);
	const size_t col = static_cast<size_t>(l);
	if (l>=1 && col<=this->matrix->numColumns()) // note: cannot use col>=1 because size_t is unsigned
		return new RINGELEM(this->matrix->theMatrix(this->nRow, col-1));
	throw RuntimeException("Invalid column index", indexExpBegin, indexExpEnd);
}

intrusive_ptr<Value> MODULEELEM::indexedByBigInt(intrusive_ptr<INT> N, const CharPointer & /* targetExpBegin */, const CharPointer & /* targetExpEnd */, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) {
	long l;
  long nc = NumCompts(this->theModuleElem);
	if (!IsConvertible(l, N->theBigInt) || l<=0 || l>nc)
		throw IndexOutOfRangeException(N->theBigInt, nc, indexExpBegin, indexExpEnd);
	return Value::from((this->theModuleElem)[l-1]);
}

IntegerIndexedAccess::IntegerIndexedAccess(const intrusive_ptr<LeftValue> targetLV, const BigInt &index,  const CharPointer &targetExpBegin, const CharPointer &targetExpEnd, const CharPointer &indexExpBegin, const CharPointer &indexExpEnd) :
	IndexedAccess(targetLV, targetExpBegin, targetExpEnd, indexExpBegin, indexExpEnd),
	index(index)
{
	intrusive_ptr<RightValue> target = this->targetLV->asRightValue();
	if (!target->canBeIndexedByBigInt())
		throw NonIntegerIndexableException(target, targetExpBegin, targetExpEnd);
}

SnapshotFrame::SnapshotFrame(const Frame &frame) :
	invocationExp(frame.invocationExp),
	block(frame.block)
{
	assert(this->invocationExp);
	assert(this->block);
}

RuntimeException::RuntimeException(const string &reason, const CharPointer &from, const CharPointer &to, const Frame * const) :
	ExceptionWithSourcePosition(reason, from, to, false)
{}

RuntimeException::RuntimeException(const string &reason, intrusive_ptr<const ParsedObject> po, const Frame * const) :
	ExceptionWithSourcePosition(reason, po->getBegin(), po->getEnd(), false)
{}

RuntimeException::RuntimeException(const string &reason, const Token &token, const Frame * const) :
	ExceptionWithSourcePosition(reason, token, false)
{}

PackageValue::PackageValue(const intrusive_ptr<const PackageStatement> pkgDecl) :
		pkgDecl(pkgDecl),
		pkgName(pkgDecl->name),
		prefix(pkgDecl->name+".")
{}

void RightValue::describe(intrusive_ptr<OSTREAM> out) const {
	out->print("A value of type ")->println(this->getType());
}

void INT::describe(intrusive_ptr<OSTREAM> out) const {
	if (this->theBigInt==BigInt(42))
		out->println("The Answer to the Ultimate Question of Life, the Universe and Everything");
	else
		this->RightValue::describe(out);
}

void RING::describe(intrusive_ptr<OSTREAM> out) const {
  ostringstream os;
  this->theRing->myOutputSelfLong(os);
  out->print(os.str())->newline()->flush();
}

void RINGHOM::describe(intrusive_ptr<OSTREAM> out) const {
	out->println(this);
}

void BuiltInFunction::describe(intrusive_ptr<OSTREAM> out) const {
	out->println("A built-in function");
}

void UserDefinedFunction::describe(intrusive_ptr<OSTREAM> out) const {
  //	out->println(this->fnDecl->getBegin().stringTo(this->fnDecl->getEnd()));
	out->println("A user defined function.\n\
For more information use print, starting(\"..\"), or PackageOf(\"..\")");
}

void PackageValue::describe(intrusive_ptr<OSTREAM> out) const
{
  const intrusive_ptr<const PackageStatement> pkgDecl = this->pkgDecl;
  const string &pkgName = pkgDecl->name;
  if (this->pkgDecl->memberNames.empty())
  {
    out->print("The package ")->print(pkgName)->println(" is empty.\n");
    return;
  }
  set<string> exports;
  if (pkgDecl->exportedNames.empty())
  { out->print("The package ")->print(pkgName)->print(" has no exported names.\n"); }
  else
  {
    out->print("The package ")->print(pkgName)->print(" exports the following names:\n");
    // Put names into a set so that they are printed in alphabetical order.
    for (const Token &t: pkgDecl->exportedNames)
      exports.insert(t.lexeme());
    for (const string &s: exports)
      out->print("* ")->println(s);
  }
  if (pkgDecl->memberNames.size() > exports.size())
  {
    out->print("\nThe package ")->print(pkgName)->print(" also has the following non-exported members:\n");
    for (const string &s: pkgDecl->memberNames)
      if (exports.find(s) == exports.end()) // print s only if it not an export
        out->print("* ")->print(pkgName)->print(".")->println(s);
  }
}

OutputFileStreamValue::OutputFileStreamValue(intrusive_ptr<STRING> filename, intrusive_ptr<const Expression> sourceExp) :
	CppOSTREAM(stream, true, false),
	stream(filename->theString.c_str())
{
	if (!stream.is_open())
		throw RuntimeException("Cannot open file \""+filename->theString+"\" for writing", sourceExp);
}

intrusive_ptr<OSTREAM> RuntimeEnvironment::getOutputStream() {
	return this->standardOutput;
}

TaggedValue::TaggedValue(boost::intrusive_ptr<RightValue> value, const std::string &tag) :
	type(TYPE::tagType(tag)),
	value(value),
	tag(tag)
{}

map<string, intrusive_ptr<TYPE> > TYPE::taggedTypes;

intrusive_ptr<TYPE> TYPE::tagType(const string &tag) {
	map<string, intrusive_ptr<TYPE> >::const_iterator it = taggedTypes.find(tag);
	if (it!=taggedTypes.end())
		return it->second;
	intrusive_ptr<TYPE> newTaggedType(new TYPE("TAGGED(\""+tag+"\")", -1));
	taggedTypes.insert(make_pair(tag, newTaggedType));
	return newTaggedType;
}

ostream &TYPE::dumpAsString(ostream &out) const {
	return out << this->name;
}

ostream &RINGELEM::dumpAsString(ostream &out) const {
	return out << this->theRingElem;
}

ostream &IDEAL::dumpAsString(ostream &out) const {
	return out << this->theIdeal;
}

ostream &MODULE::dumpAsString(ostream &out) const {
	return out << this->theModule;
}

ostream &MODULEELEM::dumpAsString(ostream &out) const {
	return out << this->theModuleElem;
}

ostream &RING::dumpAsString(ostream &out) const {
	return out << this->theRing;
}

ostream &IntMapValue::dumpAsString(ostream &out) const {
	out << "{";
	bool first=true;
	for (MapType::const_iterator pos = this->map.begin(); pos!=this->map.end(); ++pos) {
		if (first)
			first = false;
		else
			out << ", ";
		out << pos->first << " -> " << pos->second;
		pos->second->dumpRefCountAsString(out);
	}
	return out << "}";
}

intrusive_ptr<TYPE> BOOL::type(new TYPE("BOOL", TYPE::DISPATCH_INDEX_OF_BOOL));
intrusive_ptr<TYPE> FUNCTION::type(new TYPE("FUNCTION", TYPE::DISPATCH_INDEX_OF_FUNCTION));
intrusive_ptr<TYPE> LIST::type(new TYPE("LIST", TYPE::DISPATCH_INDEX_OF_LIST));
intrusive_ptr<TYPE> INT::type(new TYPE("INT", TYPE::DISPATCH_INDEX_OF_BIGINT));
intrusive_ptr<TYPE> RAT::type(new TYPE("RAT", TYPE::DISPATCH_INDEX_OF_BIGRAT));
intrusive_ptr<TYPE> RECORD::type(new TYPE("RECORD", TYPE::DISPATCH_INDEX_OF_RECORD));
intrusive_ptr<TYPE> TYPE::type(new TYPE("TYPE", TYPE::DISPATCH_INDEX_OF_TYPE));
intrusive_ptr<TYPE> STRING::type(new TYPE("STRING", TYPE::DISPATCH_INDEX_OF_STRING));
intrusive_ptr<TYPE> VoidValue::type(new TYPE("VOID", TYPE::DISPATCH_INDEX_OF_VOID));
intrusive_ptr<TYPE> ERROR_MESG::type(new TYPE("ERROR_MESG", TYPE::DISPATCH_INDEX_OF_ERROR_MESG));
intrusive_ptr<TYPE> OSTREAM::type(new TYPE("OSTREAM", TYPE::DISPATCH_INDEX_OF_OSTREAM));
intrusive_ptr<TYPE> ISTREAM::type(new TYPE("ISTREAM", TYPE::DISPATCH_INDEX_OF_ISTREAM));
//JAA 20140901 intrusive_ptr<TYPE> ZZModValue::type(new TYPE("ZMOD", TYPE::DISPATCH_INDEX_OF_ZMOD));
intrusive_ptr<TYPE> RINGELEM::type(new TYPE("RINGELEM", TYPE::DISPATCH_INDEX_OF_RINGELEM));
intrusive_ptr<TYPE> RatFunValue::type(new TYPE("RATFUN", TYPE::DISPATCH_INDEX_OF_RATFUN));
intrusive_ptr<TYPE> IDEAL::type(new TYPE("IDEAL", TYPE::DISPATCH_INDEX_OF_IDEAL));
intrusive_ptr<TYPE> MODULE::type(new TYPE("MODULE", TYPE::DISPATCH_INDEX_OF_MODULE));
intrusive_ptr<TYPE> MODULEELEM::type(new TYPE("MODULEELEM", TYPE::DISPATCH_INDEX_OF_MODULEELEM));
intrusive_ptr<TYPE> MAT::type(new TYPE("MAT", TYPE::DISPATCH_INDEX_OF_MAT));
intrusive_ptr<TYPE> RING::type(new TYPE("RING", TYPE::DISPATCH_INDEX_OF_RING));
intrusive_ptr<TYPE> PackageValue::type(new TYPE("PACKAGE", TYPE::DISPATCH_INDEX_OF_PACKAGE));
intrusive_ptr<TYPE> IntMapValue::type(new TYPE("INTMAP", TYPE::DISPATCH_INDEX_OF_INTMAP));
intrusive_ptr<TYPE> RINGHOM::type(new TYPE("RINGHOM", TYPE::DISPATCH_INDEX_OF_RINGHOM));
intrusive_ptr<TYPE> MatrixRowValue::type(new TYPE("MATRIXROW", TYPE::DISPATCH_INDEX_OF_MATRIXROW));


intrusive_ptr<TYPE> TaggedValue::getType()  const {return this->type;}
intrusive_ptr<TYPE> BOOL::getType()         const {return type;}
intrusive_ptr<TYPE> FUNCTION::getType() const {return type;}
intrusive_ptr<TYPE> LIST::getType()         const {return type;}
intrusive_ptr<TYPE> INT::getType()          const {return type;}
intrusive_ptr<TYPE> RAT::getType()          const {return type;}
intrusive_ptr<TYPE> RECORD::getType()       const {return type;}
intrusive_ptr<TYPE> STRING::getType()       const {return type;}
intrusive_ptr<TYPE> VoidValue::getType()    const {return type;}
intrusive_ptr<TYPE> ERROR_MESG::getType()   const {return type;}
intrusive_ptr<TYPE> OSTREAM::getType()      const {return type;}
intrusive_ptr<TYPE> ISTREAM::getType()      const {return type;}
intrusive_ptr<TYPE> TYPE::getType()         const {return type;}
intrusive_ptr<TYPE> PackageValue::getType() const {return type;}
intrusive_ptr<TYPE> RINGELEM::getType()     const {return type;}
intrusive_ptr<TYPE> IDEAL::getType()        const {return type;}
intrusive_ptr<TYPE> MODULE::getType()       const {return type;}
intrusive_ptr<TYPE> MODULEELEM::getType()   const {return type;}
intrusive_ptr<TYPE> RING::getType()         const {return type;}
intrusive_ptr<TYPE> IntMapValue::getType()  const {return type;}
intrusive_ptr<TYPE> RINGHOM::getType()      const {return type;}
intrusive_ptr<TYPE> MAT::getType()          const {return type;}
intrusive_ptr<TYPE> MatrixRowValue::getType() const {return type;}

vector<RING::SymbolPair> RING::allSymbolValues() {
	vector<SymbolPair> result;
	createAllSymbolValues(this->theRing, result);
	return result;
}

void RING::createAllSymbolValues(const ring &R, std::vector<SymbolPair> &allIndets)
{
  vector<symbol> syms;
  if (IsPolyRing(R))
    syms = symbols(CoeffRing(R));
  else
    syms = symbols(R);
  for (const symbol &s: syms)
    allIndets.push_back(make_pair(s, new RINGELEM(RingElem(R, s))));
  if (!IsPolyRing(R)) return;
  const long nvars = NumIndets(R);
  allIndets.reserve(allIndets.size() + nvars);
  for (long i=0; i < nvars; ++i)
    allIndets.push_back(make_pair(IndetSymbol(R,i), new RINGELEM(indet(R,i))));
}

void RING::removeInjectedIndeterminates(RuntimeEnvironment *runtimeEnv) {
	Frame * const tlFrame = runtimeEnv->getTopLevelFrame();
	for (int slot: this->injectedSlots) {
		VariableSlot &vs = tlFrame->varSlots[slot];
		vs.value = 0;
		vs.unprotect();
	}
	this->injectedSlots.clear();
}

void RING::injectIndeterminates(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RingDefinition> ringDefinition, const string &ringName) {
	assert(this->injectedSlots.empty());
	string protectedNames, msgPrefix, msgSuffix;
	Frame * const tlFrame = runtimeEnv->getTopLevelFrame();
	vector<SymbolPair> allIndets(this->allSymbolValues());
	for (const SymbolPair &pair: allIndets)
        {
          const symbol &sym = pair.first;
          const string hd = head(sym);
          if (hd==ringName)
            runtimeEnv->interpreter->errorReporter->reportWarning("The name of an indeterminate corresponds to the name of the ring", ringDefinition->getBegin(), ringDefinition->getEnd(), WS_NORMAL);
          const int slot = runtimeEnv->slotFor(hd);
          if (slot<0)
            continue;
          if (tlFrame->varSlots[slot].isProtected()) {
            if (protectedNames.length()==0) {
              protectedNames = hd;
              msgPrefix = " ";
              msgSuffix = " is ";
            } else {
              protectedNames += ", ";
              protectedNames += hd;
              msgPrefix = "s ";
              msgSuffix = " are ";
            }
          }
	}
	if (protectedNames.length())
		throw RuntimeException("Cannot use the ring because the name"+msgPrefix+protectedNames+msgSuffix+"protected", ringDefinition);
// UNUSED?	int indetIndex=0;
	map<string, intrusive_ptr<IntMapValue> > head2map;
	for (const SymbolPair &pair: allIndets)
        {
          const symbol &sym = pair.first;
          const int nSubscripts = NumSubscripts(sym);
          const intrusive_ptr<RINGELEM> indeterminate(pair.second);
          const string hd(head(sym));
          //cout << "sym=" << sym << ", hd=" << hd << ", nSubscripts = " << nSubscripts << ", indeterminate=" << indeterminate << endl;
          if (nSubscripts==0)
            this->injectedSlots.push_back(runtimeEnv->setTopLevelVar(hd, indeterminate, VariableSlot::VSF_SystemProtected));
          else {
            map<string, intrusive_ptr<IntMapValue> >::const_iterator it = head2map.find(hd);
            intrusive_ptr<IntMapValue> map;
            if (it==head2map.end()) {
              map = new IntMapValue();
              head2map.insert(make_pair(hd, map));
            } else
              map = it->second;
            this->injectedSlots.push_back(runtimeEnv->setTopLevelVar(hd, map, VariableSlot::VSF_SystemProtected));
            for (int a=0; a<(nSubscripts-1); ++a) {
              IntMapValue::KeyType i = subscript(sym, a);
              intrusive_ptr<IntMapValue> v = dynamic_pointer_cast<IntMapValue>(map->getValue(i));
              if (!v) {
                v = new IntMapValue();
                map->addValue(i, v);
              }
              map = v;
            }
            map->addValue(subscript(sym, nSubscripts-1), indeterminate);
          }
// UNUSED?          ++indetIndex;
	}
}

namespace {
	class FakeExpression : public Expression {
	private:
		const intrusive_ptr<RightValue> value;
		intrusive_ptr<Value> implEval(RuntimeEnvironment *) const {
			return this->value;
		}
		bool skipDebugging() const {
			return true;
		}
	public:
		explicit FakeExpression(const intrusive_ptr<RightValue> value, const CharPointer &beginSourcePosition, const CharPointer & endSourcePosition) :
			Expression(beginSourcePosition, endSourcePosition),
			value(value)
		{}
		bool isLeftValue() const { return false; }
		ostream &dumpAsString(ostream &out) const { assert(false); return out; }
		void accept(ParsedObjectVisitor *) { assert(false); }
	};

}

intrusive_ptr<RightValue> TaggedValue::unaryMinus(const CharPointer &opPosition, RuntimeEnvironment * const runtimeEnv) {
	runtimeEnv->interpreter->reportWarning("Performing this operation requires an implicit untagging; use untagged() to avoid this warning", opPosition, opPosition);
	return this->value->unaryMinus(opPosition, runtimeEnv);
}

intrusive_ptr<RightValue> TaggedValue::clone() {
	return new TaggedValue(this->value->clone(), this->tag);
}

ostream &TaggedValue::dumpAsString(ostream &out) const {
	out << "tagged(";
	this->value->dumpAsString(out);
	intrusive_ptr<STRING> s = new STRING(this->tag);
	return out << ", " << s << ")";
}

ostream &ERROR_MESG::dumpAsString(ostream &out) const {
	out << "error(";
	intrusive_ptr<STRING> s = new STRING(this->message);
	return out << s << ")";
}

ostream &OSTREAM::dumpAsString(ostream &out) const {
	return out << "<out-stream>";
}

ostream &ISTREAM::dumpAsString(ostream &out) const {
	return out << "<in-stream>";
}

intrusive_ptr<RightValue> OSTREAM::close(RuntimeEnvironment *, intrusive_ptr<const Expression> sourceExp) {
	if (!this->canBeClosed)
		throw RuntimeException("This out-stream cannot be closed", sourceExp);
	if (this->isClosed)
		throw RuntimeException("This out-stream is already closed", sourceExp);
	this->isClosed = true;
	// we don't actually close anything here, because that's a subclass responsibility
	return VoidValue::theInstance;
}

intrusive_ptr<RightValue> OutputStringStreamValue::close(RuntimeEnvironment *runtimeEnv, intrusive_ptr<const Expression> sourceExp) {
	OSTREAM::close(runtimeEnv, sourceExp);
	return new STRING(this->ss.str());
}

intrusive_ptr<RightValue> OutputFileStreamValue::close(RuntimeEnvironment *runtimeEnv, intrusive_ptr<const Expression> sourceExp) {
	OSTREAM::close(runtimeEnv, sourceExp);
	this->stream.close();
	return VoidValue::theInstance;
}

InputStringStreamValue::InputStringStreamValue(boost::intrusive_ptr<STRING> str) :
            CppISTREAM(ss, true, false),
            str(str),
            ss(str->theString)
	{}


InputFileStreamValue::InputFileStreamValue(intrusive_ptr<STRING> filename, intrusive_ptr<const Expression> sourceExp) :
	CppISTREAM(stream, true, false),
	stream(filename->theString.c_str())
{
	if (!stream.is_open())
		throw RuntimeException("Cannot open file \""+filename->theString+"\" for reading", sourceExp);
}

boost::intrusive_ptr<RightValue> InputFileStreamValue::close(RuntimeEnvironment *runtimeEnv, intrusive_ptr<const Expression> sourceExp)
{
	ISTREAM::close(runtimeEnv, sourceExp);
	this->stream.close();
	return VoidValue::theInstance;
}


boost::intrusive_ptr<RightValue> ISTREAM::close(RuntimeEnvironment *, intrusive_ptr<const Expression> sourceExp) {
	if (!this->canBeClosed)
		throw RuntimeException("This in-stream cannot be closed", sourceExp);
	if (this->isClosed)
		throw RuntimeException("This in-stream is already closed", sourceExp);
	this->isClosed = true;
	// we don't actually close anything here, because that's a subclass responsibility
	return VoidValue::theInstance;
}


boost::intrusive_ptr<OSTREAM> OutputBoostSocketStreamValue::print(const std::string &str) { (*out) << str; return this; }
boost::intrusive_ptr<OSTREAM> OutputBoostSocketStreamValue::print(boost::intrusive_ptr<const RightValue> v) { (*out) << v; return this; }
boost::intrusive_ptr<OSTREAM> OutputBoostSocketStreamValue::newline() { (*out) << '\n' << std::flush; return this; }
void OutputBoostSocketStreamValue::flush() { out->flush(); }



namespace {
	string stripPkgName(const string &s) {
		string::size_type dotPos = s.find('.', 0);
		assert(dotPos!=string::npos);
		return s.substr(dotPos+1);
	}

	string getPkgNameWithDot(const string &s) {
		static const string::size_type l = ParserNS::Parser::TopLevelPackageName.length();
		if (s.substr(0, l)==ParserNS::Parser::TopLevelPackageName)
			return "";
		string::size_type dotPos = s.find('.', 0);
		assert(dotPos!=string::npos);
		return s.substr(0, dotPos+1);
	}
}

void printTaggedValue(RuntimeEnvironment *runtimeEnv, intrusive_ptr<TaggedValue> v, intrusive_ptr<OSTREAM> out, intrusive_ptr<const Expression> sourceExp) {
	const string recordName(getPkgNameWithDot(v->tag)+"PrintTagged");
	int slot = runtimeEnv->slotFor(recordName);
	if (slot<0) {
		runtimeEnv->interpreter->reportWarning("Cannot find \""+recordName+"\", so I'm implicitly untagging the value", sourceExp);
implicit_untagging:
		out->print(runtimeEnv, v->untagged(), sourceExp);
		return;
	}
	intrusive_ptr<RECORD> record = dynamic_pointer_cast<RECORD>(runtimeEnv->getTopLevelFrame()->varSlots[slot].value);
	if (!record) {
		runtimeEnv->interpreter->reportWarning("The variable \""+recordName+"\" is not a record, so I'm implicitly untagging the value", sourceExp);
		goto implicit_untagging;
	}
	const string fieldName(stripPkgName(v->tag));
	intrusive_ptr<RightValue> fieldValue = record->getFieldNoCheck(fieldName);
	if (!fieldValue) {
		runtimeEnv->interpreter->reportWarning("The variable \""+recordName+"\" does not contain a field named \""+fieldName+"\", so I'm implicitly untagging the value", sourceExp);
		goto implicit_untagging;
	}
	intrusive_ptr<FUNCTION> printingFun = dynamic_pointer_cast<FUNCTION>(fieldValue);
	const string funName(recordName+"."+fieldName);
	if (!printingFun) {
		runtimeEnv->interpreter->reportWarning(funName+" does not contain a function, so I'm implicitly untagging the value", sourceExp);
		goto implicit_untagging;
	}
	if (!printingFun->canBeCalledWith(2)) {
		runtimeEnv->interpreter->reportWarning("The function "+funName+" cannot receive two arguments, so I'm implicitly untagging the value", sourceExp);
		goto implicit_untagging;
	}
	if (const intrusive_ptr<UserDefinedFunction> udf = dynamic_pointer_cast<UserDefinedFunction>(printingFun)) {
		for (int a=0; a<=1; ++a)
			if (udf->fnDecl->params[a].byRef) {
				runtimeEnv->interpreter->reportWarning("The function "+funName+" expects a by-ref argument, so I'm implicitly untagging the value", sourceExp);
				goto implicit_untagging;
			}
	}
	vector<Argument> args;
	args.push_back(Argument(false, new FakeExpression(out, sourceExp->getBegin(), sourceExp->getEnd()), true));
	args.push_back(Argument(false, new FakeExpression(v->untagged(), sourceExp->getBegin(), sourceExp->getEnd()), true));
	intrusive_ptr<InvocationExpression> fakeInvoke =
			new InvocationExpression(
					new FakeExpression(printingFun, sourceExp->getBegin(), sourceExp->getEnd()),
					boost::shared_ptr<Token>(),
					args,
					sourceExp->getEnd(),
					""
			);
	intrusive_ptr<Value> retValue = fakeInvoke->eval(runtimeEnv);
	if (!dynamic_pointer_cast<VoidValue>(retValue))
		throw RuntimeException("Printing-procedures must not return any value", sourceExp);
}

void OSTREAM::print(RuntimeEnvironment *runtimeEnv, intrusive_ptr<RightValue> v,  intrusive_ptr<const Expression> sourceExp) {
	if (this->isClosed)
		throw RuntimeException("Cannot print to a closed out-stream", sourceExp);
	if (const intrusive_ptr<STRING> s = dynamic_pointer_cast<STRING>(v))
		this->print(s->theString);
	else if (const intrusive_ptr<TaggedValue> tv = dynamic_pointer_cast<TaggedValue>(v))
		printTaggedValue(runtimeEnv, tv, this, sourceExp);
	else
		this->print(v);
}

void CppOSTREAM::flush() {
	out.flush();
}

intrusive_ptr<OSTREAM> CppOSTREAM::print(const string &s) {
	this->out << s;
	return this;
}

intrusive_ptr<OSTREAM> CppOSTREAM::print(intrusive_ptr<const RightValue> v) {
	this->out << v;
	return this;
}

intrusive_ptr<OSTREAM> CppOSTREAM::newline() {
	this->out << endl;
	return this;
}


// Returns "\n" to indicate EOF (or any other read failure)
std::string CppISTREAM::getline()
{
  if (!in) return "\n"; // THROWING CAUSES A CRASH: throw IOException("GetLine from closed <in-stream>");
  string str;
  std::getline(in, str);
///???  if (in.eof() && str.empty()) return "\n";
  return str;
}


bool TaggedValue::needsToBeCopiedBeforeChanges() const {
	return this->value->needsToBeCopiedBeforeChanges();
}

FramePointer::FramePointer(Frame *const frame) :
	frame(frame),
	id(frame ? frame->id : 0)
{
}

Frame *FramePointer::toCheckedPointer() const {
#ifdef VERBOSE_RUNTIME_DEBUG
	//cout << "toCheckedPointer this->frame=" << this->frame << endl;
	//if (this->frame)
	//	cout << "expected id=" << this->id << ", found " << this->frame->id << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	return this->frame ? (this->frame->id==this->id ? this->frame : nullptr) : nullptr;
}

Frame *FramePointer::toNonNullCheckedPointer(intrusive_ptr<const Identifier> expId) const {
	Frame *f = this->toCheckedPointer();
	if (!f)
		throw DeadEnviromentException(expId);
	return f;
}

intrusive_ptr<StaticEnv> CheckNamesVisitor::addEnvironmentForIterationVar(const string &identifier, const CharPointer &from, const CharPointer &to) {
	StaticEnv::VarData vd = this->env->lookup(identifier);
	if (vd.depth!=StaticEnv::IN_THE_WILD)
		this->errorReporter->reportWarning("This name hides an outer one", from, to);
	this->env = new StaticEnv(this->env);
	this->env->add(identifier, 0, 0, true, false, false);
	return this->env;
}

void CheckNamesVisitor::visit(IdInExpSuchThatExp &idInExpSuchThatExp) {
	idInExpSuchThatExp.exp1->accept(this);
	idInExpSuchThatExp.staticEnv = this->addEnvironmentForIterationVar(idInExpSuchThatExp.identifier, idInExpSuchThatExp.tokIdentifier.getBegin(), idInExpSuchThatExp.tokIdentifier.getEnd());
	idInExpSuchThatExp.exp2->accept(this);
	this->popEnv();
}

void CheckNamesVisitor::visit(ExpSuchThatIdInExpAndExp &expSuchThatIdInExpAndExp) {
	expSuchThatIdInExpAndExp.exp2->accept(this);
	expSuchThatIdInExpAndExp.staticEnv = this->addEnvironmentForIterationVar(expSuchThatIdInExpAndExp.identifier, expSuchThatIdInExpAndExp.tokIdentifier.getBegin(), expSuchThatIdInExpAndExp.tokIdentifier.getEnd());
	expSuchThatIdInExpAndExp.exp1->accept(this);
	if (expSuchThatIdInExpAndExp.optionalExp3)
		expSuchThatIdInExpAndExp.optionalExp3->accept(this);
	this->popEnv();
}

void CheckNamesVisitor::visit(FieldAccessExpression &fieldAccessExp) {
	if (const intrusive_ptr<Identifier> identifier = dynamic_pointer_cast<Identifier>(fieldAccessExp.targetExp)) {
		StaticEnv::VarData vd = this->env->lookup(identifier->identifier);
		if (vd.depth==StaticEnv::IN_THE_WILD) {
			string alias;
			if (this->aliases->lookup(identifier->identifier, alias)) {
				identifier->identifier = alias;
				assert(vd.index<0);
				assert(!vd.isCapturedValue);
				assert(!vd.isIterationVar);
				assert(!vd.isImplictlyImported);
				vd.depth = StaticEnv::TOP_LEVEL;
			} else {
				if (this->insideFunctions)
					this->reportError(VariableNotFoundException::errorMessage(identifier->identifier, identifier->arity, this->env), identifier->getBegin(), identifier->getEnd());
				else
					vd.depth = StaticEnv::TOP_LEVEL; // outside fn-procs, anything goes...
			}
		} else {
			string alias;
			if (this->aliases->lookup(identifier->identifier, alias))
				this->errorReporter->reportWarning("I'm using the local/imported variable here, but note that there is an alias, with the same name, for the package "+alias, identifier->getBegin(), identifier->getEnd(), WS_PEDANTIC);
		}
		identifier->varData = vd;
	} else
		fieldAccessExp.targetExp->accept(this);
}

void CheckNamesVisitor::visit(Identifier &identifier) {
	StaticEnv::VarData vd = this->env->lookup(identifier.identifier);
	if (vd.depth==StaticEnv::IN_THE_WILD) {
		if (this->runtimeEnvironment->isRegisteredType(identifier.identifier))
			vd.depth = StaticEnv::TOP_LEVEL;
		else {
			if (this->insideFunctions)
				this->reportError(VariableNotFoundException::errorMessage(identifier.identifier, identifier.arity, this->env), identifier.getBegin(), identifier.getEnd());
			else {
				vd.depth = StaticEnv::TOP_LEVEL; // outside fn-procs, anything goes...
				vd.isImplictlyImported = true; // this is important for non-fnproc package member declarations
			}
		}
	}
	identifier.varData = vd;
}

void CheckNamesVisitor::reportError(const string &msg, const CharPointer &from, const CharPointer &to) {
	this->foundErrors = true;
	this->errorReporter->reportError(msg, from, to);
}

void CheckNamesVisitor::visit(ForStatement &forStmt) {
	forStmt.beginExp->accept(this);
	if (forStmt.stepExp)
		forStmt.stepExp->accept(this);
	forStmt.endExp->accept(this);
	forStmt.staticEnv = this->addEnvironmentForIterationVar(forStmt.identifier, forStmt.tokIdentifier.getBegin(), forStmt.tokIdentifier.getEnd());
	forStmt.statements->accept(this);
	this->popEnv();
}

void CheckNamesVisitor::visit(ForeachStatement &foreachStmt) {
	foreachStmt.inExp->accept(this);
	foreachStmt.staticEnv = this->addEnvironmentForIterationVar(foreachStmt.identifier, foreachStmt.tokIdentifier.getBegin(), foreachStmt.tokIdentifier.getEnd());
	foreachStmt.statements->accept(this);
	this->popEnv();
}

void CheckNamesVisitor::visit(TryStatement &tryStmt) {
	tryStmt.tryStatements->accept(this);
	tryStmt.staticEnv = this->addEnvironmentForIterationVar(tryStmt.identifier, tryStmt.tokIdentifier.getBegin(), tryStmt.tokIdentifier.getEnd());
	tryStmt.uponErrorStatements->accept(this);
	this->popEnv();
}

void CheckNamesVisitor::visit(FunctionDeclaration &funDecl) {
	bool oldInsideFunctions = this->insideFunctions;
	this->insideFunctions = true;
	funDecl.staticEnv = this->env = new StaticEnv(this->env);
	for (Import &i: funDecl.imports)
        {
          StaticEnv::VarData vd = (i.type==Import::IT_BYVALUE ? this->env->parent : this->env)->lookup(i.name);
          const int depth = vd.depth;
          if (i.type==Import::IT_TOPLEVEL) {
            if (i.implicit && depth>=0)
              this->errorReporter->reportWarning("The name \""+i.name+"\" has been implicitly imported from TopLevel; note that importing the name by-ref could yield a different result (add an explicit import to avoid this message)", funDecl.getBegin(), funDecl.getEnd());
            this->env->add(i.name, StaticEnv::TOP_LEVEL, -1, false, false, i.implicit);
            continue;
          }
          if (depth==StaticEnv::IN_THE_WILD) {
            assert(i.expId);
            this->reportError("Cannot find "+i.name, i.expId->getBegin(), i.expId->getEnd());
            continue;
          }
          if (i.type==Import::IT_BYVALUE) {
            assert(i.expId);
            i.expId->varData = vd;
            assert(!i.implicit);
            this->env->add(i.name, 0, i.byValueIndex, false, true, false);
          } else {
            assert(i.type==Import::IT_BYREF);
            if (vd.isIterationVar) {
              assert(i.expId);
              this->reportError("Iteration variables cannot be imported by-ref", i.expId->getBegin(), i.expId->getEnd());
            } if (vd.isCapturedValue) {
              assert(i.expId);
              this->reportError("Names imported by-value cannot be re-imported by-ref", i.expId->getBegin(), i.expId->getEnd());
            } else {
              assert(!i.implicit);
              this->env->add(i.name, depth, vd.index, false, false, false);
            }
          }
	}
	for (map<string, int>::const_iterator it=funDecl.localNames.begin(); it!=funDecl.localNames.end(); ++it)
		env->add(it->first, 0, it->second, false, false, false);
	funDecl.statements->accept(this);
	this->popEnv();
	this->insideFunctions = oldInsideFunctions;
}

void CheckNamesVisitor::visit(InvocationExpression &invocationExp) {
	ParsedObjectVisitor::visit(invocationExp);
	for (const Argument &arg: invocationExp.args)
		if (arg.byRef)
			if (const intrusive_ptr<Identifier> id = dynamic_pointer_cast<Identifier>(arg.exp)) {
				if (id->varData.isCapturedValue)
					this->reportError("Cannot pass by-ref an imported by value", arg.exp->getBegin(), arg.exp->getEnd());
			}
}

void CheckNamesVisitor::visit(AssignmentStatement &assignmentStmt) {
	ParsedObjectVisitor::visit(assignmentStmt);
	if (const intrusive_ptr<Identifier> id = dynamic_pointer_cast<Identifier>(assignmentStmt.leftExp)) {
		if (id->varData.isCapturedValue)
			this->reportError("Cannot assign to a variable imported by value", assignmentStmt.leftExp->getBegin(), assignmentStmt.leftExp->getEnd());
	}
}

void CheckNamesVisitor::visit(DefineStatement &defineStmt) {
	defineStmt.funDecl->accept(this);
}

void CheckNamesVisitor::visit(AliasStatement &aliasStmt) {
	assert(this->aliases);
	//if (aliasStmt.statements)
	//	this->aliases = new AliasEnv(this->aliases);
	for (const Binding &b: aliasStmt.bindings)
        {
          const string newAlias(b.identifier->identifier);
          if ( /* aliasStmt.statements || */ this->insidePackage) {
            string useless;
            if (this->aliases->parent->lookup(newAlias, useless))
              this->errorReporter->reportWarning("This alias hides an outer one", b.identifier->getBegin(), b.identifier->getEnd());
          }
          this->aliases->add(newAlias, b.packageName);
	}
	/* if (aliasStmt.statements) {
		aliasStmt.statements->accept(this);
		this->aliases = this->aliases->parent;
		assert(this->aliases);
	} */
}

void CheckNamesVisitor::visit(PackageStatement &p) {
	assert(this->aliases);
	this->insidePackage = true;
	this->aliases = new AliasEnv(this->aliases);
	p.statements->accept(this);
	this->aliases = this->aliases->parent;
	assert(this->aliases);
	this->insidePackage = false;
}

void ResolvePackageNamesVisitor::visit(InvocationExpression &invocationExp) {
	const intrusive_ptr<Identifier> id = dynamic_pointer_cast<Identifier>(invocationExp.targetExp);
	if (!id || !id->varData.isImplictlyImported || this->packageStatement.memberNames.find(id->identifier)==this->packageStatement.memberNames.end())
		invocationExp.targetExp->accept(this);
	else {
		StaticEnv::VarData &vd = id->varData;
		assert(vd.depth == StaticEnv::TOP_LEVEL);
		assert(!vd.isCapturedValue);
		assert(!vd.isIterationVar);
		id->identifier = this->packageStatement.name+"."+id->identifier;
		vd.isImplictlyImported = false;
		//cout << "\n\n### " << id << " --> " << invocationExp.targetExp << " ###\n\n";
	}
	for (const Argument &arg: invocationExp.args)
		arg.exp->accept(this);
}

namespace {
	// based on: http://www.merriampark.com/ld.htm and http://www.merriampark.com/ldcpp.htm
	int levenshteinDistance(string source, string target) {
	  const int n = source.length();
	  const int m = target.length();
	  if (n == 0)
		return m;
	  if (m == 0)
		return n;
	  to_lower(source);
	  to_lower(target);
	  vector<vector<int> > matrix(n+1);
	  for (int i = 0; i <= n; ++i)
		matrix[i].resize(m+1);
	  for (int i = 0; i <= n; ++i)
		matrix[i][0]=i;
	  for (int j = 0; j <= m; ++j)
		matrix[0][j]=j;
	  for (int i = 1; i <= n; ++i) {
		const char s_i = source[i-1];
		for (int j = 1; j <= m; ++j) {
		  const char t_j = target[j-1];
		  const int cost = s_i == t_j ? 0 : 1;
		  const int above = matrix[i-1][j];
		  const int left = matrix[i][j-1];
		  const int diag = matrix[i-1][j-1];
		  int cell = min(above + 1, min(left + 1, diag + cost));
		  // Step 6A: Cover transposition, in addition to deletion,
		  // insertion and substitution. This step is taken from:
		  // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's
		  // Enhanced Dynamic Programming ASM Algorithm"
		  // (http://www.acm.org/~hlb/publications/asm/asm.html)
		  if (i>2 && j>2) {
			int trans=matrix[i-2][j-2]+1;
			if (source[i-2]!=t_j) ++trans;
			if (s_i!=target[j-2]) ++trans;
			if (cell>trans) cell=trans;
		  }
		  matrix[i][j]=cell;
		}
	  }
	  return matrix[n][m];
	}

	inline int maximumDistanceForSimilarIdentifiers(const string &id) {
		const int len = id.length();
		assert(len>0);
		if (len<=1)
			return 0; // use only case-insensitivity when looking for one-char identifiers (otherwise all of them would be considered similar)
    return (len+3)/4; // allow 1-char diff for lengths 2,3,4; 2-char diff for lengths 5-8, etc.
	}
}

void Interpreter::reportInterrupt(const CharPointer &from, const CharPointer &to) {
	this->errorReporter->reportInterrupt(from, to);
}

void Interpreter::reportError(const string &msg) {
	this->errorReporter->reportError(msg);
}

void Interpreter::reportError(const string &msg, const CharPointer &from, const CharPointer &to) {
	this->errorReporter->reportError(msg, from, to);
}

intrusive_ptr<const StaticEnv> RuntimeEnvironment::findStaticEnv() {
	Frame *f = this->currentFrame;
//JAA	for(; f>this->frames; --f)
          for (; f>getTopLevelFrame(); --f)
		if (f->userdefinedFun) {
			assert(f->userdefinedFun->fnDecl);
			assert(f->userdefinedFun->fnDecl->staticEnv);
			return f->userdefinedFun->fnDecl->staticEnv;
		}
	return topLevelStaticEnv;
}

Frame *RuntimeEnvironment::pushFrame(const FramePointer &accessLink, intrusive_ptr<const InvocationExpression> invocationExp, intrusive_ptr<const ParsedObject> block, intrusive_ptr<const UserDefinedFunction> userdefinedFun) {
	assert(block);
#ifdef VERBOSE_RUNTIME_DEBUG
	cout << "pushFrame for " << block << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	if (this->currentFrame==&frames.back())
          throw RuntimeException("Too many nested scopes", invocationExp ? invocationExp : block);
	Frame * const f = ++this->currentFrame;
	f->id = ++this->nextFrameId;
	assert(f->id); // f->id==0 would mean we've just overflowed!
	f->accessLink = accessLink;
	f->invocationExp = invocationExp;
	f->block = block;
	f->userdefinedFun = userdefinedFun;
#ifdef C5IDE
	f->singleStepOnPop = false;
#endif // #ifdef C5IDE
	return f;
}

Frame *RuntimeEnvironment::pushIterationFrame(intrusive_ptr<const ParsedObject> block) {
	assert(block);
	assert(block->staticEnv);
	Frame *f = this->pushFrame(this->currentFrame, 0, block, 0);
	this->currentFrame->varSlots.push_back(VariableSlot(VoidValue::theInstance, static_cast<VariableSlot::Flags>(VariableSlot::VSF_SystemProtected|VariableSlot::VSF_IterationVariable)));
	return f;
}

void RuntimeEnvironment::popFrame() {
	Frame * const f = this->currentFrame--;
	assert(f > getTopLevelFrame());
#ifdef C5IDE
	if (f->singleStepOnPop)
		this->interpreter->singleStepExecution = true;
#endif // #ifdef C5IDE
	f->id = 0;
	f->accessLink.reset();
	f->varSlots.clear();
	f->invocationExp = nullptr;
	f->block = nullptr;
	f->userdefinedFun = nullptr;
#ifdef VERBOSE_RUNTIME_DEBUG
	if (this->currentFrame==&frames.back())
		cout << "popFrame, back to global\n";
	else
		cout << "popFrame, back to: " << this->currentFrame->block << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
}

vector<SnapshotFrame> RuntimeEnvironment::takeSnapshot() {
	vector<SnapshotFrame> result;
	const Frame *f = this->currentFrame;
	do {
		if (f->invocationExp) {
			result.push_back(SnapshotFrame(*f));
		}
	} while (--f != getTopLevelFrame());
	return result;
}

void RuntimeEnvironment::registerType(intrusive_ptr<TYPE> type) {
#ifndef NDEBUG
	bool inserted =
#endif
	this->registeredTypeNames.insert(type->name).second;
	assert(inserted);
	this->setTopLevelVar(type->name, type, VariableSlot::VSF_SystemProtected);
}

intrusive_ptr<LIST> RuntimeEnvironment::currentTypes() {
	intrusive_ptr<LIST> result(new LIST);
	for (const string &s: this->registeredTypeNames)
        {
          const int slot = this->slotFor(s);
          assert(slot>=0 && static_cast<vector<VariableSlot>::size_type>(slot)<getTopLevelFrame()->varSlots.size());
          result->addValue(intrusive_ptr_cast<TYPE>(getTopLevelFrame()->varSlots[slot].value));
	}
	// Below p is more-or-less of type pair<string, intrusive_ptr<TYPE> >;
	for (const auto& p: TYPE::taggedTypes)
		result->addValue(p.second);
	return result;
}

  RuntimeEnvironment::RuntimeEnvironment(Interpreter * const interpreter, intrusive_ptr<OSTREAM> standardOutput, long MaxStackSize) :
    frames(MaxStackSize),
    currentFrame(&(this->frames[0])),
	nextFrameId(1),
	standardOutput(standardOutput),
	topLevelStaticEnv(new StaticEnv(this)),
	topLevelAliases(new AliasEnv(0)),
	itSlot(this->setTopLevelVar(Interpreter::IT, VoidValue::theInstance, VariableSlot::VSF_SystemProtected)),
	currentRingSlot(this->setTopLevelVar(Interpreter::CURRENT_RING, VoidValue::theInstance, VariableSlot::VSF_SystemProtected)),
	interpreter(interpreter)
{
	this->initBuiltInFunctions();
	this->registerType(BOOL::type);
	this->registerType(FUNCTION::type);
	this->registerType(LIST::type);
	this->registerType(INT::type);
	this->registerType(RAT::type);
	this->registerType(RECORD::type);
	this->registerType(TYPE::type);
	this->registerType(STRING::type);
	this->registerType(VoidValue::type);
	this->registerType(ERROR_MESG::type);
	this->registerType(OSTREAM::type);
	this->registerType(ISTREAM::type);
//JAA 20140901	this->registerType(ZZModValue::type);
	this->registerType(RINGELEM::type);
	this->registerType(RatFunValue::type);
	this->registerType(IDEAL::type);
	this->registerType(MODULE::type);
	this->registerType(MODULEELEM::type);
	this->registerType(MAT::type);
	this->registerType(RING::type);
	this->registerType(IntMapValue::type);
	this->registerType(RINGHOM::type);
	this->registerType(MatrixRowValue::type);
	this->setTopLevelVar("ZZ", new RING(RingZZ()), VariableSlot::VSF_SystemProtected);
	this->setTopLevelVar("QQ", new RING(RingQQ()), VariableSlot::VSF_SystemProtected);
	this->setTopLevelVar("P", new RING(NewPolyRing(RingQQ(),symbols("x,y,z"))), VariableSlot::VSF_None); // *** SEE ALSO CoCoA-5/packages/init.cocoa5  and  Banner.C ***
	this->initMaps();
}

intrusive_ptr<LIST> RuntimeEnvironment::TopLevelFunctions() {
	intrusive_ptr<LIST> result(new LIST);
	for (map<string, int>::const_iterator it = this->topLevelIdentifiers.begin(); it!=this->topLevelIdentifiers.end(); ++it) {
		assert(it->first.length());
		assert(it->second>=0 && it->second<static_cast<int>(getTopLevelFrame()->varSlots.size()));
		const VariableSlot &vslot = getTopLevelFrame()->varSlots[it->second];
		if (intrusive_ptr<FUNCTION> f = dynamic_pointer_cast<FUNCTION>(vslot.value)) {
			intrusive_ptr<RECORD> record(new RECORD);
			result->addValue(record);
			record->setFieldNoCheck("name", new STRING(it->first));
			record->setFieldNoCheck("IsExported", Value::from(vslot.hasBeenExported()));
		}
	}
	return result;
}


namespace // anonymous
{

	intrusive_ptr<PackageValue> FindOwnerPkg(RuntimeEnvironment *runtimeEnv, const string& name)
  {
    for (const VariableSlot &vs: runtimeEnv->getTopLevelFrame()->varSlots)
      if (intrusive_ptr<PackageValue> p =dynamic_pointer_cast<PackageValue>(vs.value))
        for (const Token &t: p->pkgDecl->exportedNames)
          if (t.lexeme()==name)  return p;
		return nullptr;
	}

} // namespace // anonymous


intrusive_ptr<STRING> RuntimeEnvironment::PackageOf(const std::string &name)
{
  // should check name is valid as an identifier
  intrusive_ptr<PackageValue> pkgptr = FindOwnerPkg(this, name);
  if (pkgptr != nullptr)
    return intrusive_ptr<STRING>(new STRING(pkgptr->pkgName));
  // should check it exists or not, etc...
  return intrusive_ptr<STRING>(new STRING("<not a package-exported name>"));
}


  long RuntimeEnvironment::ResizeStack(long NewSize, const boost::intrusive_ptr<const AST::Expression> OrigExp)
  {
    if (NewSize < 2)
      throw RuntimeException("Ridiculous stack size", OrigExp);
    const long CurrStackHeight = 1+(currentFrame-&frames[0]);
    if (CurrStackHeight >= NewSize)
      throw RuntimeException("Stack size too small", OrigExp);
    frames.resize(NewSize);
    return CurrStackHeight;
  }

int RuntimeEnvironment::slotFor(const string &id) {
	map<string, int>::const_iterator i=this->topLevelIdentifiers.find(id);
	return i!=this->topLevelIdentifiers.end() ? i->second : -1;
}

int RuntimeEnvironment::setTopLevelVar(const string &id, intrusive_ptr<Value> v, VariableSlot::Flags flags) {
	const int slot = this->slotFor(id);
	if (slot>=0) {
		VariableSlot &vs = this->frames[0].varSlots[slot];
		vs.value = v;
		vs.flags = flags;
		return slot;
	}
	const int newSlot = this->frames[0].varSlots.size();
	this->frames[0].varSlots.push_back(VariableSlot(v, flags));
	this->topLevelIdentifiers.insert(make_pair(id, newSlot));
	return newSlot;
}

Interpreter::Interpreter(bool warnAboutCocoa5, intrusive_ptr<LineProvider> lineProvider, intrusive_ptr<ErrorReporter> errorReporter,  intrusive_ptr<OSTREAM> standardOutput, /* intrusive_ptr<ISTREAM> standardInput,*/ bool fullCoCoALibError, long MaxStackSize) :
	warnAboutCocoa5(warnAboutCocoa5),
	fullCoCoALibError(fullCoCoALibError),
	lineProvider(lineProvider),
	runtimeEnvironment(this, standardOutput, MaxStackSize),
////	controlC(false),
#ifdef C5IDE
	singleStepExecution(false),
	doStepOver(false),
#endif // #ifdef C5IDE
	errorReporter(errorReporter)
#ifdef C5IDE
	, status(IS_WAITING_FOR_COMMAND) // while not 100% true (the call to this->run() will make this happen), it's a safe initialization value
#endif // #ifdef C5IDE
{
	assert(this->lineProvider);
	assert(this->errorReporter);
}

string VariableNotFoundException::errorMessage(const string &varName, int arity, intrusive_ptr<const StaticEnv> env) {
	string message("Cannot find a variable named \"");
	message += varName;
	message += "\" in scope";
	vector<string> NearMatches;
	bool thereIsAnExactMatch=false;
	env->collectSimilarlyNamedIdentifiers(varName, arity, NearMatches, thereIsAnExactMatch);
	if (thereIsAnExactMatch)
		return message+", but there is one in an outside scope.  You're probably missing an import statement";
	const int size = NearMatches.size();
	if (size==0)
		return message;
	if (size==1) {
		message += ".\nA similarly named variable (that you might need to import) is: \"";
		message += NearMatches.front();
		message += '\"';
		return message;
	}
	message += ".\nSimilarly named variables (that you might need to import) are: ";
	bool first = true;
	for (const string &id: NearMatches)
        {
          if (first)
            first = false;
          else
            message += ", ";
          message += '\"';
          message += id;
          message += '\"';
	}
	return message;
}

void RECORD::collectSimilarlyNamedFields(const string &fieldName, set<string> &set) const {
	const int maxDistance = maximumDistanceForSimilarIdentifiers(fieldName);
	for (MapType::const_iterator pos = this->fields.begin(); pos!=this->fields.end(); ++pos)
		if (levenshteinDistance(fieldName, pos->first)<=maxDistance)
			set.insert(pos->first);
}

vector<string> RECORD::myFieldNamesStrings() const {
	vector<string> result;
	for (MapType::const_iterator pos = this->fields.begin(); pos!=this->fields.end(); ++pos)
		result.push_back(pos->first);
	return result;
}

intrusive_ptr<LIST> RECORD::fieldNames() const {
	intrusive_ptr<LIST> result(new LIST);
	for (MapType::const_iterator pos = this->fields.begin(); pos!=this->fields.end(); ++pos)
		result->addValue(new STRING(pos->first));
	return result;
}

string FieldNotFoundException::errorMessage(const string &fieldName, intrusive_ptr<const RECORD> rv) {
	string message("Cannot find a field named \"");
	message += fieldName;
	message += '\"';
	if (rv->numberOfFields()==0)
		message += ". Note: the record is actually empty; that is, it has no fields at all";
	else {
		set<string> s;
		rv->collectSimilarlyNamedFields(fieldName, s);
		const int size = s.size();
		if (size==0)
			return message;
		if (size==1) {
			message += ".\nA similarly named (existing) field is: \"";
			message += *s.begin();
			message += '\"';
			return message;
		}
		message += ".\nSimilarly named (existing) fields are: ";
		bool first = true;
		for (const string &id: s)
                {
                  if (first)
                    first = false;
                  else
                    message += ", ";
                  message += '\"';
                  message += id;
                  message += '\"';
		}
	}
	return message;
}

void Interpreter::reportWarning(const string &msg, const intrusive_ptr<const ParsedObject> po) {
	this->errorReporter->reportWarning(msg, po->getBegin(), po->getEnd());
}

void Interpreter::reportWarning(const string &msg, const CharPointer &from, const CharPointer &to) {
	this->errorReporter->reportWarning(msg, from, to);
}

intrusive_ptr<Value> RuntimeEnvironment::evalArg(const Argument &arg, EvalKind evalKind) {
	const bool byRef = evalKind==EVAL_BY_REF;
	if (arg.byRef && !byRef)
		throw CorrespondingParameterNotReferenceException(arg.exp);
	intrusive_ptr<Value> v;
	if (byRef) {
		if (!arg.byRef && !arg.synthetized)
			this->interpreter->reportWarning("You should use \"ref\" when passing arguments by reference", arg.exp);
		v = arg.exp->eval(this);
		if (!dynamic_pointer_cast<LeftValue>(v))
			throw RuntimeException("Only left-values can be passed by reference", arg.exp);
	} else
		v = arg.exp->evalAs<RightValue>(this);
	return v;
}

const string Interpreter::IT("It");
const string Interpreter::CURRENT_RING("CurrentRing");

namespace // anonymous
{

  void CheckProtection(VariableSlot *vs, const intrusive_ptr<const Identifier> id, RuntimeEnvironment *runtimeEnv, const string& action)
  {
    const string& name = id->identifier;
    if (vs->isProtected())
    {
      const string& reason(vs->getProtectionReason());
      if (!reason.empty())
        throw ProtectedVariableException(name, "Cannot "+action+" \""+name+"\" (protected: "+reason+")", id);
      if (vs->isSystemProtected())
        throw ProtectedVariableException(name, "Cannot "+action+" \""+name+"\" (system-protected variable)", id);
      throw ProtectedVariableException(name, "Cannot "+action+" \""+name+"\" (user-protected variable)", id);
    }
    if (vs->isInPackage())
    {
      /* const string& reason(vs->getProtectionReason()); */
      /* if (!reason.empty()) */
      /*   throw ProtectedVariableException(name, "Cannot set \""+name+"\" (protected: "+reason+")", id); */
      throw ProtectedVariableException(name, "Cannot "+action+" \""+name+"\" (package-exported variable: "+FindOwnerPkg(runtimeEnv, name)->pkgName+")", id);
    }
  }

} // end of anonymous namespace

bool TYPE::isProperSubtypeOf(int otherDI) const {
	switch (this->dispatchIndex) {
	case DISPATCH_INDEX_OF_BIGRAT:
		return otherDI==DISPATCH_INDEX_OF_RINGELEM;
	case DISPATCH_INDEX_OF_BIGINT:
		return otherDI==DISPATCH_INDEX_OF_BIGRAT || otherDI==DISPATCH_INDEX_OF_RINGELEM;
	case DISPATCH_INDEX_OF_RINGHOM:
		return otherDI==DISPATCH_INDEX_OF_FUNCTION;
	case DISPATCH_INDEX_OF_RECORD:
	case DISPATCH_INDEX_OF_BOOL:
	case DISPATCH_INDEX_OF_VOID:
	case DISPATCH_INDEX_OF_STRING:
	case DISPATCH_INDEX_OF_LIST:
	case DISPATCH_INDEX_OF_FUNCTION:
	case DISPATCH_INDEX_OF_TYPE:
	case DISPATCH_INDEX_OF_ERROR_MESG:
	case DISPATCH_INDEX_OF_OSTREAM:
	case DISPATCH_INDEX_OF_ISTREAM:
//JAA 20140901	case DISPATCH_INDEX_OF_ZMOD:
	case DISPATCH_INDEX_OF_RATFUN:
	case DISPATCH_INDEX_OF_MODULEELEM:
	case DISPATCH_INDEX_OF_IDEAL:
	case DISPATCH_INDEX_OF_MODULE:
	case DISPATCH_INDEX_OF_MAT:
	case DISPATCH_INDEX_OF_RING:
	case DISPATCH_INDEX_OF_PACKAGE:
	case DISPATCH_INDEX_OF_RINGELEM:
	case DISPATCH_INDEX_OF_INTMAP:
	case DISPATCH_INDEX_OF_MATRIXROW:
		return false;
	default:
		assert(false);
	}
	return false;
}

intrusive_ptr<const RightValue> TYPE::convert(intrusive_ptr<const RightValue> from, intrusive_ptr<const RightValue> to) {
	assert(from);
	assert(to);
	assert(from->getType()->isProperSubtypeOf(to->getType()->dispatchIndex));
	switch (to->getType()->dispatchIndex) {
		case DISPATCH_INDEX_OF_BIGRAT:
                  return new RAT(BigRat(intrusive_ptr_cast<const INT>(from)->theBigInt));
		case DISPATCH_INDEX_OF_RINGELEM:
			if (from->getType()->dispatchIndex==DISPATCH_INDEX_OF_BIGINT)
				return new RINGELEM(RingElem(owner(intrusive_ptr_cast<const RINGELEM>(to)->theRingElem), intrusive_ptr_cast<const INT>(from)->theBigInt));
			return new RINGELEM(RingElem(owner(intrusive_ptr_cast<const RINGELEM>(to)->theRingElem), intrusive_ptr_cast<const RAT>(from)->theBigRat));
		case DISPATCH_INDEX_OF_FUNCTION:
			return from;
	}
	assert(false);
	return nullptr;
}


intrusive_ptr<RightValue> RuntimeEnvironment::binaryOperatorDispatch(intrusive_ptr<const RightValue> leftOp, intrusive_ptr<const RightValue> rightOp, DispatchMapType map, const LexerNS::CharPointer &beginOperatorSourcePosition, const LexerNS::CharPointer &endOperatorSourcePosition)
{
  this->interpreter->checkForInterrupts(beginOperatorSourcePosition, endOperatorSourcePosition);
  const intrusive_ptr<const TYPE> leftType(leftOp->getType());
  assert(leftType);
  int leftIndex = leftType->dispatchIndex;
  if (leftIndex<0)
  {
    const intrusive_ptr<const TaggedValue> tv = dynamic_pointer_cast<const TaggedValue>(leftOp);
    assert(tv);
    this->interpreter->reportWarning("I'm implicitly untagging the left operand; use untagged() to avoid this warning", beginOperatorSourcePosition, endOperatorSourcePosition);
    leftOp = tv->untagged();
    leftIndex = leftOp->getType()->dispatchIndex;
  }
  assert(leftIndex>=0 && leftIndex<TYPE::NUM_OF_DISPATCHABLE_TYPES);
  const intrusive_ptr<const TYPE> rightType(rightOp->getType());
  assert(rightType);
  int rightIndex = rightType->dispatchIndex;
  if (rightIndex<0)
  {
    const intrusive_ptr<const TaggedValue> tv = dynamic_pointer_cast<const TaggedValue>(rightOp);
    assert(tv);
    this->interpreter->reportWarning("I'm implicitly untagging the right operand; use untagged() to avoid this warning", beginOperatorSourcePosition, endOperatorSourcePosition);
    rightOp = tv->untagged();
    rightIndex = rightOp->getType()->dispatchIndex;
  }
  assert(rightIndex>=0 && rightIndex<TYPE::NUM_OF_DISPATCHABLE_TYPES);
  for (;;)
  {
    BinaryOpFunc f = map[leftIndex][rightIndex];
    if (f)
      try
      {
        return f(this, leftOp, rightOp, beginOperatorSourcePosition, endOperatorSourcePosition);
      }
      catch (const ErrorInfo& err)
      {
        this->announceCLE(err);
        throw RuntimeException(message_from_cocoalib(err), beginOperatorSourcePosition, endOperatorSourcePosition);
      }
      catch (const InterruptReceived& intr)
      {
        throw InterruptException(beginOperatorSourcePosition, endOperatorSourcePosition);
      }
    if (leftOp->getType()->isProperSubtypeOf(rightIndex) && map[rightIndex][rightIndex])
    {
      try
      {
        leftOp = TYPE::convert(leftOp, rightOp);
      }
      catch (const ErrorInfo& err)
      {
        this->announceCLE(err);
        throw RuntimeException(message_from_cocoalib(err), beginOperatorSourcePosition, endOperatorSourcePosition);
      }
      catch (const InterruptReceived& /*intr*/)
      {
        throw InterruptException(beginOperatorSourcePosition, endOperatorSourcePosition);
      }

      assert(leftOp->getType()->dispatchIndex == rightIndex);
      leftIndex = rightIndex;
      continue;
    }
    if (rightOp->getType()->isProperSubtypeOf(leftIndex) && map[leftIndex][leftIndex])
    {
      try
      {
        rightOp = TYPE::convert(rightOp, leftOp);
      }
      catch (const ErrorInfo& err)
      {
        this->announceCLE(err);
        throw RuntimeException(message_from_cocoalib(err), beginOperatorSourcePosition, endOperatorSourcePosition);
      }
      catch (const InterruptReceived& /*intr*/)
      {
        throw InterruptException(beginOperatorSourcePosition, endOperatorSourcePosition);
      }

      assert(rightOp->getType()->dispatchIndex == leftIndex);
      rightIndex = leftIndex;
      continue;
    }
    const string leftType = leftOp->getType()->name;
    const string rightType = rightOp->getType()->name;
    throw RuntimeException("I don't know how to evaluate operator "+beginOperatorSourcePosition.stringTo(endOperatorSourcePosition)+" between "+leftType+" and "+rightType, beginOperatorSourcePosition, endOperatorSourcePosition);
  }
}

void RuntimeEnvironment::announceCLE(const ErrorInfo& err) {
	if (this->interpreter->fullCoCoALibError) {
		ostringstream os;
		ANNOUNCE(os, err);
		this->standardOutput->print(os.str());
	}
}

//----- operators ==, !=, <, >, <=, >=, +, -, *, /  -----------------------
// template functions for implementation in CoCoALib

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opEqual(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::from(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) == theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opNotEqual(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::from(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) != theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opLessThan(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::from(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) < theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opLessOrEqual(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::from(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) <= theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opGreaterOrEqual(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::from(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) >= theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opGreaterThan(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::from(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) > theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opPlus(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::simplifiedValueFrom(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) + theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opMinus(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::simplifiedValueFrom(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) - theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opStar(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::simplifiedValueFrom(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) * theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

template<typename LeftType, typename RightType>
intrusive_ptr<RightValue> opSlash(RuntimeEnvironment *, const intrusive_ptr<const RightValue> LeftTypeOp, const intrusive_ptr<const RightValue> RightTypeOp, const CharPointer &, const CharPointer &) {
	return Value::simplifiedValueFrom(theValue(intrusive_ptr_cast<const LeftType>(LeftTypeOp)) / theValue(intrusive_ptr_cast<const RightType>(RightTypeOp)));
}

//----- operators <, >, <=, >= --------------------------------------------
// not implemented as such in CoCoALib

intrusive_ptr<RightValue> opLessThan_RINGELEM_RINGELEM(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const RINGELEM> a = intrusive_ptr_cast<const RINGELEM>(leftOp);
	intrusive_ptr<const RINGELEM> b = intrusive_ptr_cast<const RINGELEM>(rightOp);
  if (IsOrderedDomain(owner(a->theRingElem)))
    return Value::from(a->theRingElem < b->theRingElem);
  if (IsMonomial(a->theRingElem) && IsMonomial(b->theRingElem)
      && IsOne(LC(a->theRingElem)) && IsOne(LC(b->theRingElem)))
    return Value::from(LPP(a->theRingElem) < LPP(b->theRingElem));
  throw RuntimeException("RINGELEM can be compared only if in ordered domain or if power-products", beginOperatorSourcePosition, endOperatorSourcePosition);
}


intrusive_ptr<RightValue> opLessOrEqual_RINGELEM_RINGELEM(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const RINGELEM> a = intrusive_ptr_cast<const RINGELEM>(leftOp);
	intrusive_ptr<const RINGELEM> b = intrusive_ptr_cast<const RINGELEM>(rightOp);
  if (IsOrderedDomain(owner(a->theRingElem)))
    return Value::from(a->theRingElem <= b->theRingElem);
  if (IsMonomial(a->theRingElem) && IsMonomial(b->theRingElem)
      && IsOne(LC(a->theRingElem)) && IsOne(LC(b->theRingElem)))
    return Value::from(LPP(a->theRingElem) <= LPP(b->theRingElem));
  throw RuntimeException("RINGELEM can be compared only if in ordered domain or if power-products", beginOperatorSourcePosition, endOperatorSourcePosition);
}


intrusive_ptr<RightValue> opGreaterThan_RINGELEM_RINGELEM(RuntimeEnvironment *rv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  return opLessThan_RINGELEM_RINGELEM(rv, rightOp, leftOp, beginOperatorSourcePosition, endOperatorSourcePosition);
}


intrusive_ptr<RightValue> opGreaterOrEqual_RINGELEM_RINGELEM(RuntimeEnvironment *rv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  return opLessOrEqual_RINGELEM_RINGELEM(rv, rightOp, leftOp, beginOperatorSourcePosition, endOperatorSourcePosition);
}


  intrusive_ptr<RightValue> opLessThan_IDEAL_IDEAL(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("IDEAL<IDEAL is obsolete: use IsContained(IDEAL,IDEAL) instead", beginOperatorSourcePosition, endOperatorSourcePosition);
}

  intrusive_ptr<RightValue> opLessOrEqual_IDEAL_IDEAL(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("IDEAL<=IDEAL is obsolete: use IsContained(IDEAL,IDEAL) instead", beginOperatorSourcePosition, endOperatorSourcePosition);
}

  intrusive_ptr<RightValue> opGreaterOrEqual_IDEAL_IDEAL(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("IDEAL>=IDEAL is obsolete: use IsContained(IDEAL,IDEAL) instead", beginOperatorSourcePosition, endOperatorSourcePosition);
}

  intrusive_ptr<RightValue> opGreaterThan_IDEAL_IDEAL(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("IDEAL>IDEAL is obsolete: use IsContained instead(IDEAL,IDEAL)", beginOperatorSourcePosition, endOperatorSourcePosition);
}

  intrusive_ptr<RightValue> opLessThan_MODULE_MODULE(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("MODULE<MODULE is obsolete: use IsContained(MODULE,MODULE) instead", beginOperatorSourcePosition, endOperatorSourcePosition);
}

  intrusive_ptr<RightValue> opLessOrEqual_MODULE_MODULE(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("MODULE<=MODULE is obsolete: use IsContained(MODULE,MODULE) instead", beginOperatorSourcePosition, endOperatorSourcePosition);
}

  intrusive_ptr<RightValue> opGreaterOrEqual_MODULE_MODULE(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("MODULE>=MODULE is obsolete: use IsContained(MODULE,MODULE) instead", beginOperatorSourcePosition, endOperatorSourcePosition);
}

  intrusive_ptr<RightValue> opGreaterThan_MODULE_MODULE(RuntimeEnvironment *, const intrusive_ptr<const RightValue> /*leftOp*/, const intrusive_ptr<const RightValue> /*rightOp*/, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
  throw RuntimeException("MODULE>MODULE is obsolete: use IsContained(MODULE,MODULE) instead", beginOperatorSourcePosition, endOperatorSourcePosition);
}



//----- opPower_ ---------------------------------------------------------
intrusive_ptr<RightValue> opPower_RINGELEM_INT(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &, const CharPointer &) {
	intrusive_ptr<const INT> exponent = intrusive_ptr_cast<const INT>(rightOp);
	return new RINGELEM(power(intrusive_ptr_cast<const RINGELEM>(leftOp)->theRingElem, exponent->theBigInt));
}

intrusive_ptr<RightValue> opPower_INT_INT(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &, const CharPointer &) {
	intrusive_ptr<const INT> exponent = intrusive_ptr_cast<const INT>(rightOp);
	if (exponent->theBigInt<0)
          return Value::simplifiedValueFrom(power(BigRat(intrusive_ptr_cast<const INT>(leftOp)->theBigInt), exponent->theBigInt));
	return new INT(power(intrusive_ptr_cast<const INT>(leftOp)->theBigInt, exponent->theBigInt));
}

intrusive_ptr<RightValue> opPower_RAT_INT(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &, const CharPointer &) {
	return new RAT(power(intrusive_ptr_cast<const RAT>(leftOp)->theBigRat, intrusive_ptr_cast<const INT>(rightOp)->theBigInt));
}

intrusive_ptr<RightValue> opPower_MAT_INT(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &, const CharPointer &) {
	intrusive_ptr<const INT> exponent = intrusive_ptr_cast<const INT>(rightOp);
	return Value::from(power(intrusive_ptr_cast<const MAT>(leftOp)->theMatrix, exponent->theBigInt));
}

intrusive_ptr<RightValue> opPower_IDEAL_INT(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &, const CharPointer &) {
	intrusive_ptr<const INT> exponent = intrusive_ptr_cast<const INT>(rightOp);
	return Value::from(power(intrusive_ptr_cast<const IDEAL>(leftOp)->theIdeal, exponent->theBigInt));
}

// intrusive_ptr<RightValue> opPower_RING_INT(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
// 	intrusive_ptr<const INT> exponent = intrusive_ptr_cast<const INT>(rightOp);
//   long n;
//   if (!IsConvertible(n, exponent->theBigInt) || n<0)
//     throw RuntimeException("Function values cannot be compared", beginOperatorSourcePosition, endOperatorSourcePosition);
// 	return Value::from(NewFreeModule(intrusive_ptr_cast<const RING>(leftOp)->theRing, n));
// }

//----- opSlash_ ---------------------------------------------------------

intrusive_ptr<RightValue> opSlash_INT_INT(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &, const CharPointer &) {
	return Value::simplifiedValueFrom(BigRat(intrusive_ptr_cast<const INT>(leftOp)->theBigInt, intrusive_ptr_cast<const INT>(rightOp)->theBigInt));
}

intrusive_ptr<RightValue> opSlash_RING_IDEAL(RuntimeEnvironment *, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &, const CharPointer &) {
	return Value::from(NewQuotientRing(intrusive_ptr_cast<const RING>(leftOp)->theRing, intrusive_ptr_cast<const IDEAL>(rightOp)->theIdeal));
}

// intrusive_ptr<RightValue> opEqual_Function_Function(RuntimeEnvironment *, const intrusive_ptr<const RightValue>, const intrusive_ptr<const RightValue>, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
// 	throw RuntimeException("Function values cannot be compared", beginOperatorSourcePosition, endOperatorSourcePosition);
// }

// intrusive_ptr<RightValue> opNotEqual_Function_Function(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
// 	return Value::from(!intrusive_ptr_cast<BOOL>(opEqual_Function_Function(runtimeEnv, leftOp, rightOp, beginOperatorSourcePosition, endOperatorSourcePosition))->theBool);
// }

intrusive_ptr<RightValue> opEqual_LIST_LIST(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const LIST> l1 = intrusive_ptr_cast<const LIST>(leftOp);
	intrusive_ptr<const LIST> l2 = intrusive_ptr_cast<const LIST>(rightOp);
	const LIST::ContainerType::size_type size = l1->size();
	if (size!=l2->size())
		return BOOL::falseValue;
	for (LIST::ContainerType::size_type a=0; a<size; ++a) {
		intrusive_ptr<BOOL> b = intrusive_ptr_cast<BOOL>(runtimeEnv->binaryOperatorDispatch(l1->getValue(a), l2->getValue(a), RuntimeEnvironment::opEqualMap, beginOperatorSourcePosition, endOperatorSourcePosition));
		if (!b->theBool)
			return b;
	}
	return BOOL::trueValue;
}

intrusive_ptr<RightValue> opPlus_LIST_LIST(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const LIST> l1 = intrusive_ptr_cast<const LIST>(leftOp);
	intrusive_ptr<const LIST> l2 = intrusive_ptr_cast<const LIST>(rightOp);
	intrusive_ptr<LIST> result = new LIST();
	const LIST::ContainerType::size_type size = l1->size();
	if (size!=l2->size())
		throw RuntimeException("Cannot sum lists of different sizes", beginOperatorSourcePosition, endOperatorSourcePosition);
	for (LIST::ContainerType::size_type a=0; a<size; ++a)
		result->addValue(runtimeEnv->binaryOperatorDispatch(
								l1->getValue(a),
								l2->getValue(a),
								RuntimeEnvironment::opPlusMap,
								beginOperatorSourcePosition,
								endOperatorSourcePosition));
	return result;
}

intrusive_ptr<RightValue> opMinus_LIST_LIST(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const LIST> l1 = intrusive_ptr_cast<const LIST>(leftOp);
	intrusive_ptr<const LIST> l2 = intrusive_ptr_cast<const LIST>(rightOp);
	intrusive_ptr<LIST> result = new LIST();
	const LIST::ContainerType::size_type size = l1->size();
	if (size!=l2->size())
		throw RuntimeException("Cannot subtract lists of different sizes", beginOperatorSourcePosition, endOperatorSourcePosition);
	for (LIST::ContainerType::size_type a=0; a<size; ++a)
		result->addValue(runtimeEnv->binaryOperatorDispatch(
								l1->getValue(a),
								l2->getValue(a),
								RuntimeEnvironment::opMinusMap,
								beginOperatorSourcePosition,
								endOperatorSourcePosition));
	return result;
}

intrusive_ptr<RightValue> opNotEqual_LIST_LIST(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	return Value::from(!intrusive_ptr_cast<BOOL>(opEqual_LIST_LIST(runtimeEnv, leftOp, rightOp, beginOperatorSourcePosition, endOperatorSourcePosition))->theBool);
}

intrusive_ptr<RightValue> opEqual_RECORD_RECORD(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const RECORD> r1 = intrusive_ptr_cast<const RECORD>(leftOp);
	intrusive_ptr<const RECORD> r2 = intrusive_ptr_cast<const RECORD>(rightOp);
	const RECORD::MapType::size_type numberOfFields = r1->numberOfFields();
	if (numberOfFields!=r2->numberOfFields())
		return BOOL::falseValue;
	RECORD::MapType::const_iterator end = r1->end();
	int a=0;
	for (RECORD::MapType::const_iterator begin = r1->begin(); begin!=end; ++begin, ++a) {
		intrusive_ptr<const RightValue> v2 = r2->getFieldNoCheck((*begin).first);
		if (!v2)
			return BOOL::falseValue;
		intrusive_ptr<BOOL> b = intrusive_ptr_cast<BOOL>(runtimeEnv->binaryOperatorDispatch((*begin).second, v2, RuntimeEnvironment::opEqualMap, beginOperatorSourcePosition, endOperatorSourcePosition));
		if (!b->theBool)
			return b;
	}
	return BOOL::trueValue;
}

intrusive_ptr<RightValue> opNotEqual_RECORD_RECORD(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	return Value::from(!intrusive_ptr_cast<BOOL>(opEqual_RECORD_RECORD(runtimeEnv, leftOp, rightOp, beginOperatorSourcePosition, endOperatorSourcePosition))->theBool);
}

intrusive_ptr<RightValue> opStar_IntRatRE_LIST(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const LIST> list = intrusive_ptr_cast<const LIST>(rightOp);
	intrusive_ptr<LIST> result = new LIST();
	const LIST::ContainerType::size_type size = list->size();
	for (LIST::ContainerType::size_type a=0; a<size; ++a)
		result->addValue(runtimeEnv->binaryOperatorDispatch(
								leftOp,
								list->getValue(a),
								RuntimeEnvironment::opStarMap,
								beginOperatorSourcePosition,
								endOperatorSourcePosition));
	return result;
}

intrusive_ptr<RightValue> opStar_LIST_IntRatRE(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const LIST> list = intrusive_ptr_cast<const LIST>(leftOp);
	intrusive_ptr<LIST> result = new LIST();
	const LIST::ContainerType::size_type size = list->size();
	for (LIST::ContainerType::size_type a=0; a<size; ++a)
		result->addValue(runtimeEnv->binaryOperatorDispatch(
								list->getValue(a),
								rightOp,
								RuntimeEnvironment::opStarMap,
								beginOperatorSourcePosition,
								endOperatorSourcePosition));
	return result;
}

intrusive_ptr<RightValue> opSlash_LIST_IntRatRE(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer &beginOperatorSourcePosition, const CharPointer &endOperatorSourcePosition) {
	intrusive_ptr<const LIST> list = intrusive_ptr_cast<const LIST>(leftOp);
	intrusive_ptr<LIST> result = new LIST();
	const LIST::ContainerType::size_type size = list->size();
	for (LIST::ContainerType::size_type a=0; a<size; ++a)
		result->addValue(runtimeEnv->binaryOperatorDispatch(
								list->getValue(a),
								rightOp,
								RuntimeEnvironment::opSlashMap,
								beginOperatorSourcePosition,
								endOperatorSourcePosition));
	return result;
}


intrusive_ptr<RightValue> opColon_IDEAL_IDEAL(RuntimeEnvironment * /*runtimeEnv*/, const intrusive_ptr<const RightValue> leftOp, const intrusive_ptr<const RightValue> rightOp, const CharPointer & /*beginOperatorSourcePosition*/, const CharPointer & /*endOperatorSourcePosition*/) {
	const intrusive_ptr<const IDEAL> l1 = intrusive_ptr_cast<const IDEAL>(leftOp);
	const intrusive_ptr<const IDEAL> l2 = intrusive_ptr_cast<const IDEAL>(rightOp);
	return Value::from(colon(l1->theIdeal, l2->theIdeal));
}

RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opPlusMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opMinusMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opStarMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opSlashMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opLessThanMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opLessOrEqualMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opEqualMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opNotEqualMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opGreaterOrEqualMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opGreaterThanMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opPowerMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opModMap;
RuntimeEnvironment::DispatchMapType RuntimeEnvironment::opColonMap;

void RuntimeEnvironment::initMaps() {
	static bool initialized=false;
	if (initialized)
		return;
	initialized = true;

  //------ STRING ----------------------------------------
	opPlusMap[TYPE::DISPATCH_INDEX_OF_STRING][TYPE::DISPATCH_INDEX_OF_STRING] = opPlus<STRING, STRING>;
	opLessThanMap[TYPE::DISPATCH_INDEX_OF_STRING][TYPE::DISPATCH_INDEX_OF_STRING] = opLessThan<STRING, STRING>;
	opLessOrEqualMap[TYPE::DISPATCH_INDEX_OF_STRING][TYPE::DISPATCH_INDEX_OF_STRING] = opLessOrEqual<STRING, STRING>;
	opEqualMap[TYPE::DISPATCH_INDEX_OF_STRING][TYPE::DISPATCH_INDEX_OF_STRING] = opEqual<STRING, STRING>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_STRING][TYPE::DISPATCH_INDEX_OF_STRING] = opNotEqual<STRING, STRING>;
	opGreaterOrEqualMap[TYPE::DISPATCH_INDEX_OF_STRING][TYPE::DISPATCH_INDEX_OF_STRING] = opGreaterOrEqual<STRING, STRING>;
	opGreaterThanMap[TYPE::DISPATCH_INDEX_OF_STRING][TYPE::DISPATCH_INDEX_OF_STRING] = opGreaterThan<STRING, STRING>;

  //------ BIGINT ----------------------------------------
	opPlusMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opPlus<INT, INT>;
	opMinusMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opMinus<INT, INT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opStar<INT, INT>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opSlash_INT_INT;
	opLessThanMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opLessThan<INT, INT>;
	opLessOrEqualMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opLessOrEqual<INT, INT>;
	opEqualMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opEqual<INT, INT>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opNotEqual<INT, INT>;
	opGreaterOrEqualMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opGreaterOrEqual<INT, INT>;
	opGreaterThanMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opGreaterThan<INT, INT>;
	opPowerMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opPower_INT_INT;

  //------ BIGRAT ----------------------------------------
	opPlusMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opPlus<RAT, RAT>;
	opMinusMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opMinus<RAT, RAT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opStar<RAT, RAT>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opSlash<RAT, RAT>;
	opLessThanMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opLessThan<RAT, RAT>;
	opLessOrEqualMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opLessOrEqual<RAT, RAT>;
	opEqualMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opEqual<RAT, RAT>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opNotEqual<RAT, RAT>;
	opGreaterOrEqualMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opGreaterOrEqual<RAT, RAT>;
	opGreaterThanMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opGreaterThan<RAT, RAT>;
	opPowerMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opPower_RAT_INT;

///???	opEqualMap[TYPE::DISPATCH_INDEX_OF_FUNCTION][TYPE::DISPATCH_INDEX_OF_FUNCTION] = opEqual_Function_Function;
///???	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_FUNCTION][TYPE::DISPATCH_INDEX_OF_FUNCTION] = opNotEqual_Function_Function;

  //------ LIST ----------------------------------------
	opEqualMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_LIST] = opEqual_LIST_LIST;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_LIST] = opNotEqual_LIST_LIST;
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_LIST] = opStar_IntRatRE_LIST;
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_LIST] = opStar_IntRatRE_LIST;
	opStarMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_LIST] = opStar_IntRatRE_LIST;
	opStarMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_BIGINT] = opStar_LIST_IntRatRE;
	opStarMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opStar_LIST_IntRatRE;
	opStarMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opStar_LIST_IntRatRE;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_BIGINT] = opSlash_LIST_IntRatRE;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opSlash_LIST_IntRatRE;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opSlash_LIST_IntRatRE;
	opPlusMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_LIST] = opPlus_LIST_LIST;
	opMinusMap[TYPE::DISPATCH_INDEX_OF_LIST][TYPE::DISPATCH_INDEX_OF_LIST] = opMinus_LIST_LIST;

  //------ MAT ----------------------------------------
	opStarMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_MAT] = opStar<MAT, MAT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_MAT] = opStar<INT, MAT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_MAT] = opStar<RAT, MAT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_MAT] = opStar<RINGELEM, MAT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opStar<MAT, INT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opStar<MAT, RAT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opStar<MAT, RINGELEM>;
	opPlusMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_MAT] = opPlus<MAT, MAT>;
	opMinusMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_MAT] = opMinus<MAT, MAT>;
	opEqualMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_MAT] = opEqual<MAT, MAT>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_MAT] = opNotEqual<MAT, MAT>;
	opPowerMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opPower_MAT_INT;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_BIGINT] = opSlash<MAT, INT>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opSlash<MAT, RAT>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_MAT][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opSlash<MAT, RINGELEM>;

  //------ RECORD ----------------------------------------
	opEqualMap[TYPE::DISPATCH_INDEX_OF_RECORD][TYPE::DISPATCH_INDEX_OF_RECORD] = opEqual_RECORD_RECORD;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_RECORD][TYPE::DISPATCH_INDEX_OF_RECORD] = opNotEqual_RECORD_RECORD;

  //------ BOOL ----------------------------------------
	opEqualMap[TYPE::DISPATCH_INDEX_OF_BOOL][TYPE::DISPATCH_INDEX_OF_BOOL] = opEqual<BOOL, BOOL>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_BOOL][TYPE::DISPATCH_INDEX_OF_BOOL] = opNotEqual<BOOL, BOOL>;

  //------ TYPE ----------------------------------------
	opEqualMap[TYPE::DISPATCH_INDEX_OF_TYPE][TYPE::DISPATCH_INDEX_OF_TYPE] = opEqual<TYPE, TYPE>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_TYPE][TYPE::DISPATCH_INDEX_OF_TYPE] = opNotEqual<TYPE, TYPE>;

  //------ RINGELEM ----------------------------------------
	opPlusMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opPlus<RINGELEM, RINGELEM>;
	opMinusMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opMinus<RINGELEM, RINGELEM>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opStar<RINGELEM, RINGELEM>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opSlash<RINGELEM, RINGELEM>;
  //	opLessThanMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opLessThan<RINGELEM, RINGELEM>;
  //	opLessOrEqualMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opLessOrEqual<RINGELEM, RINGELEM>;
	opLessThanMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opLessThan_RINGELEM_RINGELEM;
	opLessOrEqualMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opLessOrEqual_RINGELEM_RINGELEM;
	opEqualMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opEqual<RINGELEM, RINGELEM>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opNotEqual<RINGELEM, RINGELEM>;
  //	opGreaterOrEqualMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opGreaterOrEqual<RINGELEM, RINGELEM>;
  //	opGreaterThanMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opGreaterThan<RINGELEM, RINGELEM>;
	opGreaterOrEqualMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opGreaterOrEqual_RINGELEM_RINGELEM;
	opGreaterThanMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opGreaterThan_RINGELEM_RINGELEM;
	opPowerMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_BIGINT] = opPower_RINGELEM_INT;

  //------ MODULEELEM ----------------------------------------
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGINT][TYPE::DISPATCH_INDEX_OF_MODULEELEM] = opStar<INT, MODULEELEM>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_BIGRAT][TYPE::DISPATCH_INDEX_OF_MODULEELEM] = opStar<RAT, MODULEELEM>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_MODULEELEM] = opStar<RINGELEM, MODULEELEM>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_BIGINT] = opStar<MODULEELEM, INT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opStar<MODULEELEM, RAT>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opStar<MODULEELEM, RINGELEM>;
	opPlusMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_MODULEELEM] = opPlus<MODULEELEM, MODULEELEM>;
	opMinusMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_MODULEELEM] = opMinus<MODULEELEM, MODULEELEM>;
	opEqualMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_MODULEELEM] = opEqual<MODULEELEM, MODULEELEM>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_MODULEELEM] = opNotEqual<MODULEELEM, MODULEELEM>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_BIGINT] = opSlash<MODULEELEM, INT>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_BIGRAT] = opSlash<MODULEELEM, RAT>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_MODULEELEM][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opSlash<MODULEELEM, RINGELEM>;

  //------ RING ----------------------------------------
	opEqualMap[TYPE::DISPATCH_INDEX_OF_RING][TYPE::DISPATCH_INDEX_OF_RING] = opEqual<RING, RING>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_RING][TYPE::DISPATCH_INDEX_OF_RING] = opNotEqual<RING, RING>;
	opSlashMap[TYPE::DISPATCH_INDEX_OF_RING][TYPE::DISPATCH_INDEX_OF_IDEAL] = opSlash_RING_IDEAL;
  //	opPowerMap[TYPE::DISPATCH_INDEX_OF_RING][TYPE::DISPATCH_INDEX_OF_BIGINT] = opPower_RING_INT;

  //------ MODULE ----------------------------------------
	opEqualMap[TYPE::DISPATCH_INDEX_OF_MODULE][TYPE::DISPATCH_INDEX_OF_MODULE] = opEqual<MODULE, MODULE>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_MODULE][TYPE::DISPATCH_INDEX_OF_MODULE] = opNotEqual<MODULE, MODULE>;
  //	opSlashMap[TYPE::DISPATCH_INDEX_OF_RING][TYPE::DISPATCH_INDEX_OF_IDEAL] = opSlash_RING_IDEAL;
	opGreaterThanMap[TYPE::DISPATCH_INDEX_OF_MODULE][TYPE::DISPATCH_INDEX_OF_MODULE] = opGreaterThan_MODULE_MODULE;
	opGreaterOrEqualMap[TYPE::DISPATCH_INDEX_OF_MODULE][TYPE::DISPATCH_INDEX_OF_MODULE] = opGreaterOrEqual_MODULE_MODULE;
	opLessThanMap[TYPE::DISPATCH_INDEX_OF_MODULE][TYPE::DISPATCH_INDEX_OF_MODULE] = opLessThan_MODULE_MODULE;
	opLessOrEqualMap[TYPE::DISPATCH_INDEX_OF_MODULE][TYPE::DISPATCH_INDEX_OF_MODULE] = opLessOrEqual_MODULE_MODULE;

  //------ IDEAL ----------------------------------------
	opPlusMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opPlus<IDEAL, IDEAL>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_RINGELEM] = opStar<IDEAL, RINGELEM>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_RINGELEM][TYPE::DISPATCH_INDEX_OF_IDEAL] = opStar<RINGELEM, IDEAL>;
	opStarMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opStar<IDEAL, IDEAL>;
	opPowerMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_BIGINT] = opPower_IDEAL_INT;
	opEqualMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opEqual<IDEAL, IDEAL>;
	opNotEqualMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opNotEqual<IDEAL, IDEAL>;
	opColonMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opColon_IDEAL_IDEAL;
	opGreaterThanMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opGreaterThan_IDEAL_IDEAL;
	opGreaterOrEqualMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opGreaterOrEqual_IDEAL_IDEAL;
	opLessThanMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opLessThan_IDEAL_IDEAL;
	opLessOrEqualMap[TYPE::DISPATCH_INDEX_OF_IDEAL][TYPE::DISPATCH_INDEX_OF_IDEAL] = opLessOrEqual_IDEAL_IDEAL;
}

ostream &VariableName::dumpAsString(ostream &out) const {
	return out << "<" << this->expId->identifier << " @ " <<
			static_cast<void *>(this->framePointer.toCheckedPointer()) << "[" << this->slotIndex << "]>";
}

ostream &FieldAccess::dumpAsString(std::ostream &out) const {
	return out << "<access to field " << this->fieldName << " of " << this->targetLV << ">";
}

ostream &IntegerIndexedAccess::dumpAsString(std::ostream &out) const {
	return out << "<access to index " << this->index << " of " << this->targetLV << ">";
}

ostream &StringIndexedAccess::dumpAsString(std::ostream &out) const {
	return out << "<indexed-access to field " << this->fieldName << " of " << this->targetLV << ">";
}

intrusive_ptr<RECORD> FieldAccess::targetAsRECORD() {
	intrusive_ptr<RightValue> target = this->targetLV->asRightValue();
	intrusive_ptr<RECORD> rv = dynamic_pointer_cast<RECORD>(target);
	if (!rv)
		throw WrongTypeException(RECORD::type->name, target->getType()->name, this->fieldAccessExp->targetExp);
	return rv;
}

intrusive_ptr<RECORD> StringIndexedAccess::targetAsRECORD() {
	intrusive_ptr<RightValue> target = this->targetLV->asRightValue();
	intrusive_ptr<RECORD> rv = dynamic_pointer_cast<RECORD>(target);
	if (!rv)
		throw WrongTypeException(
				RECORD::type->name,
				target->getType()->name,
				this->targetExpBegin,
				this->targetExpEnd);
	return rv;
}

intrusive_ptr<RightValue> IntegerIndexedAccess::targetAsListOrMatrixOrMatrixRow(int &which) {
	const intrusive_ptr<RightValue> target = this->targetLV->asRightValue();
	if (dynamic_pointer_cast<LIST>(target)) {
		which = 1;
		return target;
	}
	if (dynamic_pointer_cast<MAT>(target)) {
		which = 2;
		return target;
	}
	if (dynamic_pointer_cast<MatrixRowValue>(target)) {
		which = 3;
		return target;
	}
	throw WrongTypeException(
			LIST::type->name + ", " + MAT::type->name + " or " + MatrixRowValue::type->name,
			target->getType()->name,
			this->targetExpBegin,
			this->targetExpEnd);
}

LIST::ContainerType::size_type IntegerIndexedAccess::getIndexFromBigInt(const BigInt &N, intrusive_ptr<const LIST> list) const {
	long l;
	if (!IsConvertible(l, N) || l<=0 || static_cast<LIST::ContainerType::size_type>(l)>list->size())
		throw IndexOutOfRangeException(N, list->size(), this->indexExpBegin, this->indexExpEnd);
	return l-1;
}

string::size_type IntegerIndexedAccess::getIndexFromBigInt(const BigInt &N, intrusive_ptr<const STRING> str) const {
	long l;
	if (!IsConvertible(l, N) || l<=0 || static_cast<string::size_type>(l)>str->theString.size())
		throw IndexOutOfRangeException(N, str->theString.size(), this->indexExpBegin, this->indexExpEnd);
	return l-1;
}

size_t IntegerIndexedAccess::getIndexFromBigInt(const BigInt &N, intrusive_ptr<const MatrixRowValue> mrv) const {
	long l;
	if (!IsConvertible(l, N) || l<=0 || static_cast<size_t>(l)>mrv->matrix->numColumns())
		throw IndexOutOfRangeException(N, mrv->matrix->numColumns(), this->indexExpBegin, this->indexExpEnd);
	return l-1;
}

intrusive_ptr<RightValue> IntegerIndexedAccess::asRightValue() {
	return this->targetLV->asRightValue()->indexedByBigInt(new INT(this->index), this->targetExpBegin, this->targetExpEnd, this->indexExpBegin, this->indexExpEnd)->asRightValue();
}

void IntegerIndexedAccess::assign(boost::intrusive_ptr<RightValue> value, const CharPointer &valueExpBegin, const CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv) {
	int which;
	intrusive_ptr<RightValue> rv = this->targetAsListOrMatrixOrMatrixRow(which);
	switch (which) {
	case 1 /* LIST */ : {
			const intrusive_ptr<LIST> list = intrusive_ptr_cast<LIST>(rv);
			list->setValue(this->getIndexFromBigInt(index, list), value);
		}
		break;
	case 2 /* MAT */ : {
			throw RuntimeException("Matrix rows cannot be assigned", this->targetExpBegin, this->indexExpEnd);
		}
		break;
	case 3 /* MatrixRowValue */ : {
			try
                        {
				const intrusive_ptr<MatrixRowValue> mrv = intrusive_ptr_cast<MatrixRowValue>(rv);
				if (const intrusive_ptr<RINGELEM> ringElem = dynamic_pointer_cast<RINGELEM>(value))
					SetEntry(mrv->matrix->theMatrix, mrv->nRow, this->getIndexFromBigInt(index, mrv), ringElem->theRingElem);
				else if (const intrusive_ptr<INT> N = dynamic_pointer_cast<INT>(value))
					SetEntry(mrv->matrix->theMatrix, mrv->nRow, this->getIndexFromBigInt(index, mrv), N->theBigInt);
				else if (const intrusive_ptr<RAT> qq = dynamic_pointer_cast<RAT>(value))
					SetEntry(mrv->matrix->theMatrix, mrv->nRow, this->getIndexFromBigInt(index, mrv), qq->theBigRat);
				else
					throw RuntimeException(value->getType()->name+" cannot be assigned to matrix cells", valueExpBegin, valueExpEnd);
			}
                        catch (const ErrorInfo& err)
                        {
                          runtimeEnv->announceCLE(err);
                          throw RuntimeException(message_from_cocoalib(err), valueExpBegin, valueExpEnd);
			}
		}
		break;
	default:
		assert(false);
	}
}

void IntegerIndexedAccess::checkAssignable(RuntimeEnvironment* /*runtimeEnv*/) const
{ /* always possible */ }


intrusive_ptr<VariableName> IndexedAccess::getBase() {
	return this->targetLV->getBase();
}

bool IndexedAccess::assignmentNeedsOwnership() const {
  return true;  // 2019-04-29  cannot simply change to false (see redmine 1228, comment 18)
}

intrusive_ptr<RightValue> FieldAccess::asRightValue() {
	intrusive_ptr<RECORD> rv = this->targetAsRECORD();
	intrusive_ptr<RightValue> fieldValue = rv->getField(this->fieldName, this->fieldAccessExp->tokName.getBegin(), this->fieldAccessExp->tokName.getEnd());
	if (!fieldValue)
		throw FieldNotFoundException(this->fieldName, this->fieldAccessExp->tokName.getBegin(), this->fieldAccessExp->tokName.getEnd(), rv);
	return fieldValue;
}

intrusive_ptr<RightValue> StringIndexedAccess::asRightValue() {
	intrusive_ptr<RECORD> rv = this->targetAsRECORD();
	intrusive_ptr<RightValue> fieldValue = rv->getField(this->fieldName, this->indexExpBegin, this->indexExpEnd);
	if (!fieldValue)
		throw FieldNotFoundException(this->fieldName,  this->indexExpBegin, this->indexExpEnd, rv);
	return fieldValue;
}

void StringIndexedAccess::assign(boost::intrusive_ptr<RightValue> value, const CharPointer & valueExpBegin, const CharPointer & valueExpEnd, RuntimeEnvironment * /*runtimeEnv*/) {
	intrusive_ptr<RECORD> rv = this->targetAsRECORD();
	rv->setField(this->fieldName, value, valueExpBegin, valueExpEnd);
}

void StringIndexedAccess::checkAssignable(RuntimeEnvironment* /*runtimeEnv*/) const
{ /* always possible */ }


void FieldAccess::assign(boost::intrusive_ptr<RightValue> value, const CharPointer & valueExpBegin, const CharPointer & valueExpEnd, RuntimeEnvironment * /*runtimeEnv*/) {
	intrusive_ptr<RECORD> rv = this->targetAsRECORD();
	rv->setField(this->fieldName, value, valueExpBegin, valueExpEnd);
}

void FieldAccess::checkAssignable(RuntimeEnvironment* /*runtimeEnv*/) const
{ /* always possible */ }

intrusive_ptr<VariableName> FieldAccess::getBase() {
	return this->targetLV->getBase();
}

bool FieldAccess::assignmentNeedsOwnership() const {
	return true;
}

ostream &ReferenceVariable::dumpAsString(ostream &out) const {
	return out << "<reference to " << this->framePointer.toCheckedPointer() << "[" << this->slotIndex << "]>";
}

VariableSlot *ReferenceVariable::referencedSlot() const {
	return &(this->framePointer.toNonNullCheckedPointer(this->expId)->varSlots[this->slotIndex]);
}

intrusive_ptr<RightValue> ReferenceVariable::asRightValue() {
#ifdef VERBOSE_RUNTIME_DEBUG
	cout << "ReferenceVariable::asRightValue expId=" << this->expId << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	try {
		return this->referencedSlot()->value->asRightValue();
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		this->throwIncompatibleChanges();
		return VoidValue::theInstance; // useless, but it keeps the compiler happy ;-)
	}
}

void ReferenceVariable::throwIncompatibleChanges() const {
	throw RuntimeException(
			"Reference variable \""+this->expId->identifier+"\" references an undefined variable or its underlying structure has changed in an incompatible way",
			this->expId);
}

void ReferenceVariable::assign(intrusive_ptr<RightValue> value, const CharPointer &valueExpBegin, const CharPointer &valueExpEnd, RuntimeEnvironment *runtimeEnv) {
	try {
		VariableSlot *vs = this->referencedSlot();
		CheckProtection(vs, this->expId, runtimeEnv, "set");
		return intrusive_ptr_cast<LeftValue>(vs->value)->assign(value, valueExpBegin, valueExpEnd, runtimeEnv);
	} catch (const ProtectedVariableException &pe) {
		throw ProtectedVariableException(pe.protectedVarName,
				"Reference variable "+this->expId->identifier+" cannot be assigned because "+pe.protectedVarName +
								", its (directly or indirectly) referenced variable, is protected",
				this->expId);
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		this->throwIncompatibleChanges();
	}
}

void ReferenceVariable::checkAssignable(RuntimeEnvironment* runtimeEnv) const
{
	try {
		VariableSlot *vs = this->referencedSlot();
		CheckProtection(vs, this->expId, runtimeEnv, "set");
	} catch (const ProtectedVariableException &pe) {
		throw ProtectedVariableException(pe.protectedVarName,
				"Reference variable "+this->expId->identifier+" cannot be assigned because "+pe.protectedVarName +
								", its (directly or indirectly) referenced variable, is protected",
				this->expId);
	} catch (const InterruptException &) {
		throw;
	} catch (const RuntimeException &) {
		this->throwIncompatibleChanges();
	}
}


intrusive_ptr<VariableName> ReferenceVariable::getBase() {
	return this->referencedLV()->getBase();
}

intrusive_ptr<VariableName> VariableName::getBase() {
	return this;
}

intrusive_ptr<RightValue> VariableName::asRightValue() {
	Frame * const f = this->framePointer.toNonNullCheckedPointer(this->expId);
	if (this->slotIndex<0) {
		assert(f==this->runtimeEnvironment->getTopLevelFrame());
		this->slotIndex = this->runtimeEnvironment->slotFor(this->expId->identifier);
		if (this->slotIndex<0)
			throw VariableNotFoundException(this->expId, this->runtimeEnvironment);
	}
#ifdef VERBOSE_RUNTIME_DEBUG
	cout << "VariableName::asRightValue for " << this->expId->identifier << ", this->slotIndex=" << this->slotIndex
			<< ", f->varSlots.size()=" << f->varSlots.size() << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	assert(static_cast<vector<VariableSlot>::size_type>(this->slotIndex)<f->varSlots.size());
	VariableSlot &vs = f->varSlots[this->slotIndex];
	if (!vs.value) { // this can only happen when a package has been reloaded (and some of the previous version members are not defined anymore) or an indeterminate has been removed (by Use-ing another ring)
		assert(f == this->runtimeEnvironment->getTopLevelFrame());
		throw VariableNotFoundException(this->expId, this->runtimeEnvironment);
	}
#ifdef VERBOSE_RUNTIME_DEBUG
	cout << "f->varSlots[this->slotIndex].value=" << f->varSlots[this->slotIndex].value
			<< ", f->varSlots[this->slotIndex].value->asRightValue()=" << f->varSlots[this->slotIndex].value->asRightValue() << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	return f->varSlots[this->slotIndex].value->asRightValue();
}

bool VariableName::assignmentNeedsOwnership() const {
	return false;
}

bool ReferenceVariable::assignmentNeedsOwnership() const {
	return this->referencedLV()->assignmentNeedsOwnership();
}

void VariableName::assign(boost::intrusive_ptr<RightValue> value, const CharPointer & /*valueExpBegin*/, const CharPointer & /*valueExpEnd*/, RuntimeEnvironment * /*runtimeEnv*/) {
	Frame *f = this->framePointer.toNonNullCheckedPointer(this->expId);
	if (this->slotIndex<0) {
		assert(f==this->runtimeEnvironment->getTopLevelFrame());
		this->slotIndex = this->runtimeEnvironment->slotFor(this->expId->identifier);
		if (this->slotIndex<0) {
			this->slotIndex = this->runtimeEnvironment->setTopLevelVar(this->expId->identifier, value);
			return;
		}
	}
	VariableSlot *vs = &(f->varSlots[this->slotIndex]);
	CheckProtection(vs, this->expId, this->runtimeEnvironment, "set");
	vs->value = value;
}

void VariableName::checkAssignable(RuntimeEnvironment* /*runtimeEnv*/) const
{
	Frame *f = this->framePointer.toNonNullCheckedPointer(this->expId);
        auto index = this->slotIndex;
	if (index < 0)
        {
		assert(f==this->runtimeEnvironment->getTopLevelFrame());
		index = this->runtimeEnvironment->slotFor(this->expId->identifier);
		if (index < 0) return;
	}
	VariableSlot *vs = &(f->varSlots[index]);
	CheckProtection(vs, this->expId, this->runtimeEnvironment, "set");
}


ostream &operator<<(ostream &out, const Value * const value) {
	return value->dumpAsString(out);
}

ostream &INT::dumpAsString(std::ostream &out) const {
	return out << this->theBigInt;
}

intrusive_ptr<RightValue> INT::unaryMinus(const CharPointer &, RuntimeEnvironment * const) {
	return intrusive_ptr<INT>(new INT(-this->theBigInt));
}

intrusive_ptr<RightValue> RAT::unaryMinus(const CharPointer &, RuntimeEnvironment * const) {
	return intrusive_ptr<RAT>(new RAT(-this->theBigRat));
}

intrusive_ptr<RightValue> RINGELEM::unaryMinus(const CharPointer &, RuntimeEnvironment * const) {
	return intrusive_ptr<RINGELEM>(new RINGELEM(-this->theRingElem));
}

intrusive_ptr<RightValue> MAT::unaryMinus(const CharPointer &, RuntimeEnvironment * const) {
	return intrusive_ptr<MAT>(new MAT(-this->theMatrix));
}

// ??? JAA cannot figure this one out currently :-(   ???
// intrusive_ptr<RightValue> LIST::unaryMinus(const CharPointer &, RuntimeEnvironment * const) {
//   intrusive_ptr<const LIST> l1 = this->theList;
// 	intrusive_ptr<LIST> result = new LIST();
//         typedef LIST::ContainerType::size_type INDEX;
// 	const INDEX size = l1->size();
// 	for(INDEX i=0; i<size; ++i)
// 		result->addValue(runtimeEnv->binaryOperatorDispatch(
// 								l1->getValue(i),
// 								RuntimeEnvironment::opMinusMap,
// 								beginOperatorSourcePosition,
// 								endOperatorSourcePosition));
// 	return result;
// }

intrusive_ptr<RightValue> RightValue::unaryMinus(const CharPointer &op, RuntimeEnvironment * const) {
	throw RuntimeException("I don't know how to \"negate\" a "+this->getType()->name, op, op);
}

ostream &RAT::dumpAsString(std::ostream &out) const {
	return out << this->theBigRat;
}

ostream &RECORD::dumpAsString(std::ostream &out) const {
	out << "record[";
	bool first=true;
	for (MapType::const_iterator pos = this->fields.begin(); pos!=this->fields.end(); ++pos) {
		if (first)
			first = false;
		else
			out << ", ";
		out << pos->first << " := " << pos->second;
		pos->second->dumpRefCountAsString(out);
	}
	return out << "]";
}

RECORD::MapType::size_type RECORD::numberOfFields() const {
	return this->fields.size();
}

RECORD::MapType::const_iterator RECORD::begin() const {
	return this->fields.begin();
}

RECORD::MapType::const_iterator RECORD::end() const {
	return this->fields.end();
}

ostream &BOOL::dumpAsString(std::ostream &out) const {
	return out << (this->theBool ? "true":"false");
}

ostream &STRING::dumpAsString(std::ostream &out) const {
	out << '\"';
	for (char c: this->theString) {
		switch (c) {
		case '\n':
			out << "\\n";
			break;
		case '\r':
			out << "\\r";
			break;
		case '\t':
			out << "\\t";
			break;
		case '\a':
			out << "\\a";
			break;
		case '\'':
			out << "\\\'";
			break;
		case '\"':
			out << "\\\"";
			break;
		case '\\':
			out << "\\\\";
			break;
		default:
			out << c;
		}
	}
	return out << '\"';
}

ostream &VoidValue::dumpAsString(std::ostream &out) const {
	return out << "<void>";
}

intrusive_ptr<VoidValue> VoidValue::theInstance(new VoidValue());

intrusive_ptr<RightValue> RECORD::clone() {
	intrusive_ptr<RECORD> copy(new RECORD);
	for (MapType::const_iterator pos = this->fields.begin(); pos!=this->fields.end(); ++pos)
		copy->fields.insert(make_pair(pos->first, pos->second->clone()));
	return copy;
}

intrusive_ptr<RightValue> LIST::clone() {
	intrusive_ptr<LIST> copy(new LIST(this->size()));
	for (intrusive_ptr<RightValue> v: this->container)
        {
          copy->addValue(v->clone());
	}
	return copy;
}

intrusive_ptr<RightValue> MAT::clone() {
	const matrix &m = this->theMatrix;
	const size_t nRows = NumRows(m);
	const size_t nCols = NumCols(m);
	matrix clone = matrix(m->myZeroClone(RingOf(m), nRows, nCols));
	for (size_t a=0; a<nRows; ++a)
		for (size_t b=0; b<nCols; ++b)
			SetEntry(clone, a, b, m(a, b));
	return new MAT(clone);
}

bool LIST::needsToBeCopiedBeforeChanges() const {
	if (this->getRefCount()>1)
		return true;
	for (ContainerType::const_iterator i = this->container.begin(); i!=this->container.end(); ++i)
		if ((*i)->needsToBeCopiedBeforeChanges())
			return true;
	return false;
}

bool MAT::needsToBeCopiedBeforeChanges() const {
	return this->getRefCount()>1;
}

intrusive_ptr<RightValue> LIST::getValue(ContainerType::size_type index) const {
	assert(index<=this->container.size());
	return this->container[index];
}

void LIST::setValue(ContainerType::size_type index, intrusive_ptr<RightValue> newValue) {
	assert(index<=this->container.size());
	this->container[index] = newValue;
}

void LIST::addValue(intrusive_ptr<RightValue> newValue) {
	this->container.push_back(newValue);
}

void LIST::insertValue(ContainerType::size_type position, intrusive_ptr<RightValue> newValue) {
	this->container.insert(this->container.begin()+position, newValue);
}

void LIST::removeValue(ContainerType::size_type position) {
	this->container.erase(this->container.begin()+position);
}

LIST::ContainerType::size_type LIST::size() const {
	return this->container.size();
}

ostream &MAT::dumpAsString(ostream &out) const {
	return out << this->theMatrix;
}

ostream &MatrixRowValue::dumpAsString(ostream &out) const {
	return out << "<matrix-row>";
}

ostream &LIST::dumpAsString(ostream &out) const {
	out << "[";
	bool first = true;
	for (ContainerType::const_iterator i = this->container.begin(); i!=this->container.end(); ++i) {
/// does not work properly????          CoCoA::CheckForInterrupt("Printing LIST");
		if (first)
			first = false;
		else
			out << ",  ";
		out << *i;
		(*i)->dumpRefCountAsString(out);
	}
	return out << "]";
}

void RECORD::CheckFieldName(const std::string &fieldName, const boost::intrusive_ptr<const AST::Expression> expr)
{
  if (!LexerNS::Lexer::isValidIdentifier(fieldName))
  {
    const bool IsKeyword = LexerNS::Lexer::isKeyword(fieldName);
    if (!IsKeyword)
      throw RuntimeException("Invalid field name \""+fieldName+"\"", expr);
    else
      throw RuntimeException("Invalid field name \""+fieldName+"\" - must not be a keyword", expr);
  }
}

  void RECORD::CheckFieldName(const std::string &fieldName, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to)
{
  if (!LexerNS::Lexer::isValidIdentifier(fieldName))
  {
    const bool IsKeyword = LexerNS::Lexer::isKeyword(fieldName);
    if (!IsKeyword)
      throw RuntimeException("Invalid field name \""+fieldName+"\"", from,to);
    else
      throw RuntimeException("Invalid field name \""+fieldName+"\" - must not be a keyword", from,to);
  }
}


intrusive_ptr<RightValue> RECORD::getFieldNoCheck(const string &fieldName) const {
  CheckFieldName(fieldName, CharPointer::Null, CharPointer::Null);
	MapType::const_iterator pos = this->fields.find(fieldName);
	return (pos!=this->fields.end()) ? pos->second : nullptr;
}

  intrusive_ptr<RightValue> RECORD::getField(const string &fieldName, const boost::intrusive_ptr<const AST::Expression> expr) const
  {
  CheckFieldName(fieldName, expr);
  return getFieldNoCheck(fieldName);
}

  intrusive_ptr<RightValue> RECORD::getField(const string &fieldName, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) const
  {
    CheckFieldName(fieldName, from,to);
    return getFieldNoCheck(fieldName);
}

void RECORD::setFieldNoCheck(const string &fieldName, intrusive_ptr<RightValue> newValue) {
	MapType::iterator pos = this->fields.find(fieldName);
	if (pos!=this->fields.end())
		pos->second = newValue;
	else
		this->fields.insert(make_pair(fieldName, newValue));
}

  void RECORD::setField(const string &fieldName, intrusive_ptr<RightValue> newValue, const boost::intrusive_ptr<const AST::Expression> expr)
{
  CheckFieldName(fieldName, expr);
  setFieldNoCheck(fieldName, newValue);
}


  void RECORD::setField(const string &fieldName, intrusive_ptr<RightValue> newValue, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to)
{
  CheckFieldName(fieldName, from, to);
  setFieldNoCheck(fieldName, newValue);
}

bool RECORD::needsToBeCopiedBeforeChanges() const {
	if (this->getRefCount()>1)
		return true;
	for (MapType::const_iterator pos = this->fields.begin(); pos!=this->fields.end(); ++pos)
		if (pos->second->needsToBeCopiedBeforeChanges())
			return true;
	return false;
}

/* ostream &FUNCTION::dumpAsString(ostream &out) const { */
/* 	return out << "<function>"; */
/* } */

// 2020-12 Design 1546
ostream &BuiltInFunction::dumpAsString(std::ostream &out) const {
	return out << "<built-in function>";
}

// 2020-12 Design 1546
ostream &UserDefinedFunction::dumpAsString(std::ostream &out) const {
	return out << this->fnDecl->getBegin().stringTo(this->fnDecl->getEnd());
}

ostream &RINGHOM::dumpAsString(ostream &out) const {
	return out << this->theRingHom;
}

string WrongNumberOfArgumentsException::formatExpectedArgs(int found, int min, int max) {
	assert(min<=max);
	string msg("Wrong number of arguments; found: "+boost::lexical_cast<std::string>(found)+", expecting: ");
	if (min==max)
		return msg+lexical_cast<string>(min);
	return msg+lexical_cast<string>(min) + ".." + lexical_cast<string>(max);
}

UserDefinedFunction::UserDefinedFunction(RuntimeEnvironment *runtimeEnv, const intrusive_ptr<const FunctionDeclaration> fnDecl) :
	fnDecl(fnDecl),
	fp(runtimeEnv->getCurrentFrame())
{
	assert(this->fnDecl);
	for (const Import &import: fnDecl->imports)
		if (import.type==Import::IT_BYVALUE) {
			assert(this->capturedValues.size()==static_cast<std::vector<boost::intrusive_ptr<RightValue> >::size_type>(import.byValueIndex));
			this->capturedValues.push_back(import.expId->evalAs<RightValue>(runtimeEnv));
		}
}

Frame *RuntimeEnvironment::getCurrentFunctionFrame() {
	Frame *f = this->currentFrame;
	assert(f>getTopLevelFrame());
	for (;;) {
		if (f->userdefinedFun)
			return f;
		--f;
		assert(f>getTopLevelFrame());
	}
}

Frame *RuntimeEnvironment::getCurrentFunctionFrameOrNull() {
	Frame *f = this->currentFrame;
	for (;;) {
          if (f==getTopLevelFrame())
			return nullptr;
		if (f->userdefinedFun)
			return f;
		--f;
	}
}

intrusive_ptr<RightValue> UserDefinedFunction::eval(const intrusive_ptr<const InvocationExpression> invocationExpression, RuntimeEnvironment *runtimeEnv) const {
	const vector<Param>::size_type nParameters = this->fnDecl->params.size();
	const vector<Param>::size_type nMandatoryParameters = static_cast<vector<Param>::size_type>(this->fnDecl->nMandatoryParameters);
	intrusive_ptr<LIST> ellipsisValues(new LIST());
	const bool calledWithEllipsis = invocationExpression->ellipsis;
	const bool receivesEllipsis = this->fnDecl->thereIsEllipsis;
	if (calledWithEllipsis) {
		const Frame * const f = runtimeEnv->getCurrentFunctionFrame();
		assert(f->varSlots.size());
		const VariableSlot &vs = f->varSlots.front();
		assert(vs.isArgs());
		ellipsisValues = intrusive_ptr_cast<LIST>(vs.value);
		if (!receivesEllipsis) {
			const LIST::ContainerType::size_type nArgs = ellipsisValues->size();
			if (nArgs<nMandatoryParameters || nArgs>nParameters)
				throw WrongNumberOfArgumentsException(invocationExpression, nArgs, nMandatoryParameters, nParameters);
			for (const Param &param: this->fnDecl->params)
				if (param.byRef)
					throw RuntimeException("Cannot use ellipsis on functions/procedures expecting by-ref parameters", invocationExpression);
		}
	} else { // called with actual arguments (that is, not "...")
		const vector<Argument>::size_type nArgs = invocationExpression->args.size();
		if (receivesEllipsis) {
                  for (const Argument &arg: invocationExpression->args)
				if (arg.byRef)
					throw RuntimeException("Cannot pass by-ref arguments to functions/procedures expecting a variable number of arguments", invocationExpression);
		} else
			if (nArgs<nMandatoryParameters || nArgs>nParameters)
				throw WrongNumberOfArgumentsException(invocationExpression, nArgs, nMandatoryParameters, nParameters);
	}
	vector<VariableSlot> slots;
	if (!receivesEllipsis) {
		if (calledWithEllipsis) {
			const LIST::ContainerType::size_type nArgs = ellipsisValues->size();
			for (vector<Param>::size_type a=0; a<nParameters; ++a) {
				if (a<nArgs) {
					const intrusive_ptr<RightValue> v = ellipsisValues->getValue(a);
					assert(v);
					slots.push_back(VariableSlot(v, VariableSlot::VSF_None));
				} else {
					assert(this->fnDecl->params[a].opt);
					slots.push_back(VariableSlot(VoidValue::theInstance, VariableSlot::VSF_None));
				}
			}
		} else { // don't receive ellipsis and is not called with ellipsis
			const vector<Argument>::size_type nArgs = invocationExpression->args.size();
			for (vector<Param>::size_type a=0; a<nParameters; ++a) {
				const Param &param = this->fnDecl->params[a];
				if (a<nArgs)
					slots.push_back(VariableSlot(runtimeEnv->evalArg(invocationExpression->args[a], param.byRef ? RuntimeEnvironment::EVAL_BY_REF:RuntimeEnvironment::EVAL_BY_VALUE), VariableSlot::VSF_None));
				else {
					assert(param.opt);
					slots.push_back(VariableSlot(VoidValue::theInstance, VariableSlot::VSF_None));
				}
			}
		}
	} else { // receives ellipsis
		slots.push_back(VariableSlot(ellipsisValues, VariableSlot::VSF_Args));
		if (!calledWithEllipsis) {
			const vector<Argument>::size_type nArgs = invocationExpression->args.size();
			for (vector<Argument>::size_type a=0; a<nArgs; ++a)
				ellipsisValues->addValue(runtimeEnv->evalArg(invocationExpression->args[a], RuntimeEnvironment::EVAL_BY_VALUE)->asRightValue());
		}
	}
	Frame *frame = runtimeEnv->pushFrame(this->fp, invocationExpression, this->fnDecl, this);
	BOOST_SCOPE_EXIT( (runtimeEnv) ) {
		runtimeEnv->popFrame();
	} BOOST_SCOPE_EXIT_END
	const std::map<std::string, int>::size_type localSize = this->fnDecl->localNames.size();
	assert(slots.size()<=localSize);
	while (slots.size()<localSize)
		slots.push_back(VariableSlot(VoidValue::theInstance, VariableSlot::VSF_None));
	frame->varSlots.swap(slots);
#ifdef VERBOSE_RUNTIME_DEBUG
	{
		cout << "\n======\nUserDefinedFunction::eval, frames->varSlots:";
		int i = 0;
		for (const VariableSlot &vs: frame->varSlots)
			cout << "[" << i++ << "].value = " << vs.value << endl;
		cout << "======\n";
	}
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	try {
		this->fnDecl->statements->execute(runtimeEnv);
		if (this->fnDecl->thereAreReturnsWithExpr) {
			const CharPointer p = this->fnDecl->statements->getEnd();
			throw RuntimeException("A function must explicitly Return a value", p, p);
		}
	} catch (const Return &ret) {
		return ret.value;
	} catch (const Break &) {
		assert(false);
	} catch (const Continue &) {
		assert(false);
	} catch (RuntimeException &re) {
		if (re.snapshot.size()==0)
			re.snapshot = runtimeEnv->takeSnapshot();
		throw;
	}
	return VoidValue::theInstance;
}

intrusive_ptr<RightValue> RINGHOM::eval(const intrusive_ptr<const InvocationExpression> invocationExpression, RuntimeEnvironment *runtimeEnv) const
{
	intrusive_ptr<RightValue> argValue;
	CharPointer argBegin(invocationExpression->getBegin());
	CharPointer argEnd(invocationExpression->getEnd());
	if (invocationExpression->ellipsis)
  {
		argBegin = invocationExpression->tokEllipsis->getBegin();
		argEnd = invocationExpression->tokEllipsis->getEnd();
		const Frame * const f = runtimeEnv->getCurrentFunctionFrame();
		assert(f->varSlots.size());
		const VariableSlot &vs = f->varSlots.front();
		assert(vs.isArgs());
		intrusive_ptr<LIST> ellipsisValues(intrusive_ptr_cast<LIST>(vs.value));
		const LIST::ContainerType::size_type nArgs = ellipsisValues->size();
		if (nArgs!=1)
			throw WrongNumberOfArgumentsException(invocationExpression, nArgs, 1);
		argValue = ellipsisValues->getValue(0)->asRightValue();
	}
  else
  {
		const Argument &arg = invocationExpression->args.front();
		argBegin = arg.exp->getBegin();
		argEnd = arg.exp->getEnd();
		const vector<Argument>::size_type nArgs = invocationExpression->args.size();
		if (nArgs!=1)
			throw WrongNumberOfArgumentsException(invocationExpression, nArgs, 1);
	/* 	argValue = intrusive_ptr_cast<RightValue>(runtimeEnv->evalArg(arg, RuntimeEnvironment::EVAL_BY_VALUE)); */
	/* } */
	/* try */
	/* { */
	/* 	if (intrusive_ptr<RINGELEM> poly = dynamic_pointer_cast<RINGELEM>(argValue)) */
	/* 		return new RINGELEM(this->theRingHom(poly->theRingElem)); */
	/* 	if (intrusive_ptr<INT> N = dynamic_pointer_cast<INT>(argValue)) */
	/* 		return new RINGELEM(this->theRingHom(N->theBigInt)); */
	/* 	if (intrusive_ptr<RAT> q = dynamic_pointer_cast<RAT>(argValue)) */
	/* 		return new RINGELEM(this->theRingHom(q->theBigRat)); */
	/* 	if (intrusive_ptr<MAT> q = dynamic_pointer_cast<MAT>(argValue)) */
	/* 		return new MAT(this->theRingHom(q->theMatrix)); */
	/* 	if (intrusive_ptr<LIST> q = dynamic_pointer_cast<LIST>(argValue)) */
	/* 		return new LIST(this->theRingHom(runtimeEnv->evalRVAsListOfRingElem(argValue, arg(0)))); ///// anna */
	/* 	else if (intrusive_ptr<RINGHOM> hom = dynamic_pointer_cast<RINGHOM>(argValue)) */
	/* 		return new RINGHOM(this->theRingHom(hom->theRingHom)); */
	/* 	throw WrongTypeException( */
	/* 			RINGELEM::type->name + " or " + */
	/* 			INT::type->name + " or " + */
	/* 			RAT::type->name + " or " + */
	/* 			MAT::type->name + " or " + */
	/* 			LIST::type->name + " or " + */
  /*       RINGHOM::type->name, */
	/* 			argValue->getType()->name, */
	/* 			argBegin, */
	/* 			argEnd); */
	/* } */
  /* catch (const ErrorInfo& err) */
  /* { */
  /*   runtimeEnv->announceCLE(err); */
  /*   throw RuntimeException(message_from_cocoalib(err), invocationExpression->getBegin(), invocationExpression->getEnd()); */
	/* } */
  /* catch (const InterruptReceived& /\*intr*\/) */
  /* { */
  /*   throw InterruptException(invocationExpression->getBegin(), invocationExpression->getEnd()); */
    int which;
    intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5orT6<RINGELEM, INT, RAT, MAT, LIST, RINGHOM>(ARG(0), which);
    switch (which) {
    case 1: return Value::from(this->theRingHom(RefTo<RingElem>(v)));
    case 2: return Value::from(this->theRingHom(RefTo<BigInt>(v)));
    case 3: return Value::from(this->theRingHom(RefTo<BigRat>(v)));
    case 4: return Value::from(this->theRingHom(RefTo<matrix>(v)));
    case 5: return Value::from(this->theRingHom(runtimeEnv->evalRVAsListOfRingElem(v, ARG(0))));
    case 6: return Value::from(this->theRingHom(RefTo<RingHom>(v)));
    default: return Value::from(false); // just to keep the compiler quiet
    }
  }
  return Value::from(false); // just to keep the compiler quiet
}


intrusive_ptr<RightValue> BuiltInFunction::eval(const intrusive_ptr<const InvocationExpression> invocationExpression, RuntimeEnvironment *runtimeEnv) const {
#ifdef C5IDE
	Interpreter::InterpreterStatus prevStatus = runtimeEnv->interpreter->UpdateStatusStartingBuiltIn();
	BOOST_SCOPE_EXIT( (runtimeEnv)(prevStatus) ) {
		runtimeEnv->interpreter->UpdateStatusEndingBuiltIn(prevStatus);
	} BOOST_SCOPE_EXIT_END
#endif // #ifdef C5IDE
	try {
		if (invocationExpression->ellipsis) {
			const Frame * const f = runtimeEnv->getCurrentFunctionFrame();
			assert(f->varSlots.size());
			const VariableSlot &vs = f->varSlots.front();
			assert(vs.isArgs());
			intrusive_ptr<LIST> ellipsisValues(intrusive_ptr_cast<LIST>(vs.value));
			LIST::ContainerType::size_type nArgs = ellipsisValues->size();
			vector<Argument> args;
			CharPointer begin(invocationExpression->tokEllipsis->getBegin());
			CharPointer end(invocationExpression->tokEllipsis->getEnd());
			for (LIST::ContainerType::size_type a=0; a<nArgs; ++a)
				args.push_back(Argument(false, new FakeExpression(ellipsisValues->getValue(a), begin, end), true));
      boost::shared_ptr<Token> noEllipsis;
			const intrusive_ptr<const InvocationExpression> fakeInvocationExpression(new InvocationExpression(invocationExpression->targetExp, noEllipsis, args, invocationExpression->getEnd(), invocationExpression->packageName));
			return this->function(fakeInvocationExpression, runtimeEnv);
		}
		return this->function(invocationExpression, runtimeEnv);
	}
        catch (const ErrorInfo& err)
        {
          runtimeEnv->announceCLE(err);
          throw RuntimeException(message_from_cocoalib(err), invocationExpression);
	}
        catch (const InterruptReceived& /*intr*/)
        {
          throw InterruptException(invocationExpression);
        }
}

void LeftValue::obtainOwnership() {
	const intrusive_ptr<VariableName> base = this->getBase();
	Frame * const f = base->framePointer.toNonNullCheckedPointer(base->expId);
	if (base->slotIndex<0) {
		assert(f==base->runtimeEnvironment->getTopLevelFrame());
		base->slotIndex = base->runtimeEnvironment->slotFor(base->expId->identifier);
		if (base->slotIndex<0)
			return;
	}
	assert(static_cast<vector<VariableSlot>::size_type>(base->slotIndex)<f->varSlots.size());
	assert(f->varSlots[base->slotIndex].value);
	VariableSlot &varSlot = f->varSlots[base->slotIndex];
	RightValue * const target = intrusive_ptr_cast<RightValue>(varSlot.value).get();
	if (target->needsToBeCopiedBeforeChanges())
		varSlot.value = target->clone();
}


LIST::ContainerType::size_type LIST::checkIndex(const intrusive_ptr<const RightValue> v, const intrusive_ptr<const Expression> originalExp) const
{
	if (const intrusive_ptr<const INT> NV = dynamic_pointer_cast<const INT>(v))
  {
		const BigInt N(NV->theBigInt);
		long l;
		if (!IsConvertible(l, N) || l<=0 || static_cast<ContainerType::size_type>(l)>this->size())
			throw IndexOutOfRangeException(N, this->size(), originalExp);
		return l-1;
	}
	throw WrongTypeException(INT::type->name, v->getType()->name, originalExp);
}


LIST::ContainerType::size_type LIST::checkIndexOrEnd(const intrusive_ptr<const RightValue> v, const intrusive_ptr<const Expression> originalExp) const
{
	if (const intrusive_ptr<const INT> NV = dynamic_pointer_cast<const INT>(v))
  {
		const BigInt N(NV->theBigInt);
		long l;
		if (!IsConvertible(l, N) || l<=0 || static_cast<ContainerType::size_type>(l)>this->size()+1)
			throw IndexOutOfRangeException(N, this->size(), originalExp);
		return l-1;
	}
	throw WrongTypeException(INT::type->name, v->getType()->name, originalExp);
}


namespace {
	void printErrors(intrusive_ptr<ErrorReporter> errorReporter, intrusive_ptr<OSTREAM> output) {
		assert(errorReporter);
		const int nErrors = errorReporter->getErrorCount();
		const int nWarnings = errorReporter->getWarningCount();
		if (nErrors || nWarnings)
			output->print("Got ")->print(lexical_cast<string>(nErrors))->print(" error(s) and ")->print(lexical_cast<string>(nWarnings))->print(" warning(s).\n");
	}

	void execute(intrusive_ptr<Statement> s, RuntimeEnvironment *re, intrusive_ptr<ErrorReporter> errorReporter) {
		assert(re);
#ifdef VERBOSE_RUNTIME_DEBUG
		cout << "Parsed statement = " << s << "\n";
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
		CheckNamesVisitor cnv(re, errorReporter);
		s->accept(&cnv);
#ifdef VERBOSE_RUNTIME_DEBUG
		DumpAsTreeVisitor v(cout);
		s->accept(&v);
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
		if (!errorReporter->getErrorCount())
			s->execute(re);
	}
}

#ifdef C5IDE
int Interpreter::run(IDE::Console *console) {
#else // #ifdef C5IDE
int Interpreter::run() {
#endif // #ifdef C5IDE
#ifdef C5IDE
	assert(this->status==IS_WAITING_FOR_COMMAND);
	loadPackages(console);
#endif // #ifdef C5IDE
	try
        {
          intrusive_ptr<Lexer> lexer = intrusive_ptr<Lexer>(new Lexer(this->errorReporter, this->lineProvider, this->warnAboutCocoa5, true));
          intrusive_ptr<Parser> parser(new Parser(lexer));
          double CmdStartTime = 0.0;
          bool StmtWasNOP = false; // set to "false" for correct first-time behaviour
          for (;;)
          {
#ifdef C5IDE
            this->status = IS_WAITING_FOR_COMMAND;
#endif // #ifdef C5IDE
            if (!StmtWasNOP && GlobalStatusLogStream)
            {
              const double CpuNow=CpuTime();
              GlobalStatusLogStream << 'W' << std::round(10*(CpuNow-CmdStartTime))/10 << std::endl;
            }

            intrusive_ptr<Statement> s;
            try	{
              this->errorReporter->resetErrorCounts();
              lexer->startingTopLevelStatement();
              Token t = lexer->getToken(parser->getStatus());
              if (t.getType()==TT_EOF)
              { cout<<"\n[EOF]\nCoCoA-5 quitting... (cleaning up)"<<endl; return EXIT_SUCCESS; }
              lexer->ungetToken(t);
              s = parser->parseTopLevelStatement();
              if (!s || this->errorReporter->getErrorCount())
                continue;
            } catch (const UnexpectedTokenException &ete) {
              lexer->reportError(ete.reason, ete.found.getBegin(), ete.found.getEnd(), false);
              if (ete.found.getType()==TT_EOF)
              { cout<<"\n[ERROR]\nCoCoA-5 quitting... (cleaning up)"<<endl; return EXIT_FAILURE; }
              if (ete.needsRecovery)
                parser->tryToRecover(ete.from);
              continue;
            } catch (const ExceptionWithSourcePosition &ewsp) {
              lexer->reportError(ewsp.reason, ewsp.from, ewsp.to, false);
              if (ewsp.needsRecovery)
                parser->tryToRecover(ewsp.from);
              continue;
            } catch (const BaseException &e) {
              lexer->reportError(e.reason);
              continue;
            }
            assert(s && this->errorReporter->getErrorCount()==0);

            StmtWasNOP = s->IamEmptyStmt();
            if (!StmtWasNOP)  // if NOP, do not wastefully update status to "running" (for a microsecond?)
            {
            try {
#ifdef C5IDE
              this->status = IS_RUNNING;
              this->runtimeEnvironment.getOutputStream()->flush(); // this allows the GUI thread to update the status label
#endif // #ifdef C5IDE
              if (GlobalStatusLogStream)
              {
                GlobalStatusLogStream << 'R' << std::endl;
                CmdStartTime = CpuTime();
              }
              CoCoA::GetAndResetSignalReceived(); // discard any pending signals
              try {
              execute(s, &(this->runtimeEnvironment), this->errorReporter);
              } catch (const ErrorInfo& err)
              {
                this->runtimeEnvironment.announceCLE(err);
                throw RuntimeException(message_from_cocoalib(err), s);
              }
            } catch (const Ciao &) {
#ifdef C5IDE
              this->status = IS_ENDED;
#endif // #ifdef C5IDE
              cout<<"\nCoCoA-5 quitting... (cleaning up)"<<endl;
              return EXIT_SUCCESS;
            } catch (const InterruptException& intr) {
              this->errorReporter->reportInterrupt(intr);
            } catch (const RuntimeException &re) {
              this->errorReporter->reportError(re);
            }
// catch (const CoCoA::InterruptReceived&) { /* don't know how to handle this ==> ignore it */ }
//???                        {  this->errorReporter->reportError(InterruptException(0)); } // BUG BUG BUG ??? what should the arg be???
            }
          }
	} catch (const BaseException &e) {
		cout << "***ERROR*** UNCAUGHT Interpreter BaseException exc-reason=" << e.reason << endl;
	} catch (const ErrorInfo& err) {
		cout << "***ERROR***  UNCAUGHT CoCoA error" << endl;
		ANNOUNCE(cout, err);
	} catch (const std::exception& exc) {
		cout << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
	} catch(...) {
		cout << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
	}
	return EXIT_FAILURE;
}

void Interpreter::readAndExecute(const std::string &filename, bool calledFromMain, bool immediateExecution, long FromLine, long FromChar, long ToLine, long ToChar) {
	struct stat statbuf;
	if (stat(filename.c_str(), &statbuf)) {
		if (errno==ENOENT)
			throw BaseException("Cannot find a file named \""+filename+"\"");
		throw BaseException("Cannot read \""+filename+"\" information; system-error: "+strerror(errno));
	}
	FilePair pair(statbuf.st_dev, statbuf.st_ino);
        // MinGW BUG??? st_ino is always 0 in MinGW. Most sources thus suggest ignoring it.
        // But this means that we cannot anymore detect recursive sourcing then in MinGW
	if (pair.second != 0 && !this->sourcedFiles.insert(pair).second)
		throw BaseException("Files cannot be recursively Source-d");
	Interpreter &This = *this;
	BOOST_SCOPE_EXIT( (&This)(&pair) ) {
		This.sourcedFiles.erase(pair);
	} BOOST_SCOPE_EXIT_END
	try {
          boost::intrusive_ptr<LineProvider> SourceText;
          // Check if FromLine,FromChar and ToLine,ToChar have default values: if so, just read whole file.
          if (FromLine==1 && FromChar==1 && ToLine==LONG_MAX && ToChar==LONG_MAX) SourceText = new FileLineProvider(filename);
          else SourceText = new FileRegionLineProvider(filename, FromLine,FromChar, ToLine,ToChar);
		intrusive_ptr<Lexer> lexer = new Lexer(this->errorReporter, SourceText, this->warnAboutCocoa5, false);
		intrusive_ptr<Parser> parser(new Parser(lexer));
		vector<intrusive_ptr<Statement> > statements;
		for (;;) {
			intrusive_ptr<Statement> s;
			try	{
				lexer->startingTopLevelStatement();
				Token t = lexer->getToken(parser->getStatus());
				if (t.getType()==TT_EOF)
					break;
				lexer->ungetToken(t);
				s = parser->parseTopLevelStatement();
			} catch (const UnexpectedTokenException &ete) {
				lexer->reportError(ete.reason, ete.found.getBegin(), ete.found.getEnd(), false);
				this->errorReporter->outputStream->newline();
				if (ete.found.getType()==TT_EOF) {
					if (calledFromMain) {
						printErrors(this->errorReporter, this->runtimeEnvironment.getOutputStream());
						return;
					}
					throw BaseException("Unexpected End-Of-Input while reading \""+filename+"\"");
				}
				if (ete.needsRecovery)
					parser->tryToRecover(ete.from);
			} catch (const ExceptionWithSourcePosition &ewsp) {
				lexer->reportError(ewsp.reason, ewsp.from, ewsp.to, false);
				this->errorReporter->outputStream->newline();
				if (ewsp.needsRecovery)
					parser->tryToRecover(ewsp.from);
			}
			if (s) {
				if (immediateExecution && this->errorReporter->getErrorCount()==0) {
					try {
						execute(s, &(this->runtimeEnvironment), this->errorReporter);
					} catch (const InterruptException &) {
						throw;
					} catch (const RuntimeException &re) {
						this->errorReporter->reportError(re);
						this->errorReporter->outputStream->newline();
					}
				} else
					statements.push_back(s);
			}
		}
		try {
                  for (intrusive_ptr<Statement> s: statements) {
				if (this->errorReporter->getErrorCount())
					break;
				execute(s, &(this->runtimeEnvironment), this->errorReporter);
			}
		} catch (const InterruptException &) {
			throw;
		} catch (const RuntimeException &re) {
			this->errorReporter->reportError(re);
			this->errorReporter->outputStream->newline();
		}
	} catch (BaseException &) {
		throw;
	} catch (const std::exception &e) {
		throw BaseException(string(e.what())+"\nWhile trying to read and execute \""+filename+"\"");
	}
	if (calledFromMain)
		printErrors(this->errorReporter, this->runtimeEnvironment.getOutputStream());
	else if (this->errorReporter->getErrorCount())
		throw BaseException("Read and execution of the source file \""+filename+"\" failed");
}

bool UserDefinedFunction::canBeCalledWith(int nArg) {
	const intrusive_ptr<const FunctionDeclaration> fnDecl = this->fnDecl;
	return fnDecl->thereIsEllipsis || (nArg>=fnDecl->nMandatoryParameters && static_cast<vector<Param>::size_type>(nArg)<=fnDecl->params.size());
}

bool BuiltInFunction::canBeCalledWith(int nArg) {
	return this->arityCheck(nArg);
}

bool RINGHOM::canBeCalledWith(int nArg) {
	return nArg==1;
}

intrusive_ptr<VariableName> PackageValue::toVariableName(RuntimeEnvironment *runtimeEnv, const string &memberName, const Token &tokMemberName) {
	const int index = runtimeEnv->slotFor(this->prefix+memberName);
	if (index<0)
		throw RuntimeException("Package "+this->pkgName+" doesn't have a member named "+memberName, tokMemberName);
	return new VariableName(runtimeEnv, runtimeEnv->getTopLevelFrame(), index, new Identifier(tokMemberName));
}


bool RuntimeEnvironment::CheckArity(int arity, map<string, int>::const_iterator pos) const
{
  if (arity < 0) return true;

  const VariableSlot &vs = getTopLevelFrame()->varSlots[pos->second];
  intrusive_ptr<Value> v = vs.value;
  if (!v) return false;
  if (vs.isInPackage())
  {
    const intrusive_ptr<const VariableName> vn = intrusive_ptr_cast<VariableName>(v);
    const int index = vn->expId->varData.index;
    assert(vn->expId->varData.depth==StaticEnv::TOP_LEVEL);
    assert(index>=0 && static_cast<vector<VariableSlot>::size_type>(index)<getTopLevelFrame()->varSlots.size());
    v = getTopLevelFrame()->varSlots[index].value;
    if (!v) return false;
  }
  if (const intrusive_ptr<FUNCTION> fv = dynamic_pointer_cast<FUNCTION>(v))
  {
    if (!fv->canBeCalledWith(arity))
      return false;
  } else
    return false;
  return true;
}

void RuntimeEnvironment::collectSimilarlyNamedIdentifiers(const string &id, int arity, vector<string> &NearMatches, bool &thereIsAnExactMatch) const
{
  const int maxDistance = maximumDistanceForSimilarIdentifiers(id);
  vector<string> dist0, dist1, dist2, dist3plus, DifferentArity;
  for (map<string, int>::const_iterator pos = this->topLevelIdentifiers.begin(); pos!=this->topLevelIdentifiers.end(); ++pos)
  {
    if (pos->first==id)
    {
      if (!getTopLevelFrame()->varSlots[pos->second].value)
        continue;
      thereIsAnExactMatch = true;
      return;
    }
    const int ldist = levenshteinDistance(id, pos->first);
    if (ldist>maxDistance) continue;

    if (!CheckArity(arity, pos)) { DifferentArity.push_back(pos->first); continue; }
    if (ldist == 0) { dist0.push_back(pos->first); continue; }
    if (ldist == 1) { dist1.push_back(pos->first); continue; }
    if (ldist == 2) { dist2.push_back(pos->first); continue; }
    dist3plus.push_back(pos->first);
  }
  // Sort each category of names into alphabetical order
  sort(dist0.begin(), dist0.end());
  sort(dist1.begin(), dist1.end());
  sort(dist2.begin(), dist2.end());
  sort(dist3plus.begin(), dist3plus.end());
  sort(DifferentArity.begin(), DifferentArity.end());
  // Now output the names found in increasing order of remoteness
  NearMatches.insert(NearMatches.end(), dist0.begin(), dist0.end());
  NearMatches.insert(NearMatches.end(), dist1.begin(), dist1.end());
  NearMatches.insert(NearMatches.end(), dist2.begin(), dist2.end());
  NearMatches.insert(NearMatches.end(), dist3plus.begin(), dist3plus.end());
  NearMatches.insert(NearMatches.end(), DifferentArity.begin(), DifferentArity.end());
  if (id=="Z") NearMatches.push_back("ZZ");
  if (id=="Q") NearMatches.push_back("QQ");
}

//////////////  Value::from implementations ///////////////////////////////

//??? template<typename T>  ???
 boost::intrusive_ptr<RECORD> Value::from(const factorization<RingElem>& f)
 {
   intrusive_ptr<RECORD> facsC5(new RECORD);
   facsC5->setFieldNoCheck("factors", Value::from(f.myFactors()));
   facsC5->setFieldNoCheck("multiplicities", Value::from(f.myMultiplicities()));
   facsC5->setFieldNoCheck("RemainingFactor", Value::from(f.myRemainingFactor()));
   return facsC5;
 }

 boost::intrusive_ptr<RECORD> Value::from(const factorization<BigInt>& f)
 {
   intrusive_ptr<RECORD> facsC5(new RECORD);
   facsC5->setFieldNoCheck("factors", Value::from(f.myFactors()));
   facsC5->setFieldNoCheck("multiplicities", Value::from(f.myMultiplicities()));
   facsC5->setFieldNoCheck("RemainingFactor", Value::from(f.myRemainingFactor()));
   return facsC5;
 }

 boost::intrusive_ptr<RECORD> Value::from(const HPSeries &s)
 {
   intrusive_ptr<RECORD> r(new RECORD);
   r->setFieldNoCheck("num", Value::from(num(s)));
   r->setFieldNoCheck("DenFactors", Value::from(DenFactors(s)));
   return r;
 }

 boost::intrusive_ptr<LIST> Value::from(const degree &x)
 { return Value::from(DegreeToVec(x)); }

 boost::intrusive_ptr<LIST> Value::from(const QuasiPoly &x)
 { return Value::from(constituents(x)); }

//////////////  Value::from implementations (end) ////////////////////////

}
 // namespace InterpreterNS

namespace AST {

using namespace std;
using namespace boost;
using namespace CoCoA::InterpreterNS;

IntLiteral::IntLiteral(const Token &t) :
	LiteralExpression(t.getBegin(), t.getEnd()),
	theBigInt(INT::from(BigIntFromString(t.lexeme())))
{
}

FloatLiteral::FloatLiteral(const Token &t) :
		LiteralExpression(t.getBegin(), t.getEnd()),
		theBigRat(RAT::simplifiedValueFrom(buildTheRational(t)))
{
}

StringLiteral::StringLiteral(const CharPointer &beginSourcePosition, const CharPointer & endSourcePosition, const string &unescapedString) :
	LiteralExpression(beginSourcePosition, endSourcePosition),
	theString(STRING::from(unescapedString))
{
}

BoolLiteral::BoolLiteral(const CharPointer &beginSourcePosition, const CharPointer & endSourcePosition, TokenType tokenType) :
	LiteralExpression(beginSourcePosition, endSourcePosition),
	theBool(BOOL::from(tokenType==LexerNS::TT_TRUE))
{
	assert(tokenType==LexerNS::TT_TRUE || tokenType==LexerNS::TT_FALSE);
}

int InvocationExpression::checkNumberOfArgs(int n) const {
	const int nArgs = this->args.size();
	if (nArgs!=n)
		throw WrongNumberOfArgumentsException(this, nArgs, n);
	return nArgs;
}

int InvocationExpression::checkNumberOfArgs(int nMin, int nMax) const {
	const int nArgs = this->args.size();
	if (nArgs<nMin || nArgs>nMax)
		throw WrongNumberOfArgumentsException(this, nArgs, nMin, nMax);
	return nArgs;
}

void StaticEnv::collectSimilarlyNamedIdentifiers(const string &id, int arity, vector<string> &NearMatches, bool &thereIsAnExactMatch) const {
	if (this->runtimeEnvironment) {
		assert(!this->parent);
		this->runtimeEnvironment->collectSimilarlyNamedIdentifiers(id, arity, NearMatches, thereIsAnExactMatch);
		return;
	}
	assert(this->parent);
	const int maxDistance = maximumDistanceForSimilarIdentifiers(id);
	for (IdMap::const_iterator pos = this->identifierMap.begin(); pos!=this->identifierMap.end(); ++pos) {
		if (pos->second.depth<0)
			continue;
		if (pos->first==id) {
			thereIsAnExactMatch = true;
			return;
		}
		if (levenshteinDistance(id, pos->first)<=maxDistance)
			NearMatches.push_back(pos->first);
	}
	if (!thereIsAnExactMatch)
		this->parent->collectSimilarlyNamedIdentifiers(id, arity, NearMatches, thereIsAnExactMatch);
}

namespace {
	intrusive_ptr<Value> dontKnow(const intrusive_ptr<const Expression> e) {
		throw RuntimeException("I don't know (yet) how to evaluate this expression", e->getBegin(), e->getEnd());
	}
}

intrusive_ptr<Value> ColonExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
				this->leftExp->evalAs<RightValue>(runtimeEnv),
				this->rightExp->evalAs<RightValue>(runtimeEnv),
				RuntimeEnvironment::opColonMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> ModuloExpression::implEval(RuntimeEnvironment * /*runtimeEnv*/) const {
	return dontKnow(this);
}

intrusive_ptr<Value> BackwardCompatibleExp::implEval(RuntimeEnvironment *runtimeEnv) const {
	return this->innerExp->eval(runtimeEnv);
}

intrusive_ptr<Value> ScopedExpression::implEval(RuntimeEnvironment *) const {
	return dontKnow(this);
}

intrusive_ptr<Value> UnaryPlusExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return this->exp->evalAs<RightValue>(runtimeEnv);
}

intrusive_ptr<Value> UnaryMinusExpression::implEval(RuntimeEnvironment * const runtimeEnv) const {
	return this->exp->evalAs<RightValue>(runtimeEnv)->unaryMinus(this->getBegin(), runtimeEnv);
}

intrusive_ptr<Value> IntLiteral::implEval(RuntimeEnvironment *) const {
	return this->theBigInt;
}

intrusive_ptr<Value> FloatLiteral::implEval(RuntimeEnvironment *) const {
	return this->theBigRat;
}

intrusive_ptr<Value> BoolLiteral::implEval(RuntimeEnvironment *) const {
	return this->theBool;
}

intrusive_ptr<Value> StringLiteral::implEval(RuntimeEnvironment *) const {
	return this->theString;
}

intrusive_ptr<Value> SumExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opPlusMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> PowerExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opPowerMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> IsInExpression::implEval(RuntimeEnvironment *runtimeEnv) const
{
  intrusive_ptr<RightValue> v = this->leftExp->evalAs<RightValue>(runtimeEnv);
  intrusive_ptr<RightValue> container = this->rightExp->evalAs<RightValue>(runtimeEnv);
  // First case: sthg ISIN list
  if (const intrusive_ptr<LIST> list = dynamic_pointer_cast<LIST>(container))
  {
    LIST::ContainerType::size_type size = list->size();
    for (LIST::ContainerType::size_type a=0; a<size; ++a)
    {
      try
      {
        intrusive_ptr<BOOL> bv = intrusive_ptr_cast<BOOL>(
          runtimeEnv->binaryOperatorDispatch(
            v,
            list->getValue(a),
            RuntimeEnvironment::opEqualMap,
            this->getBegin(),
            this->getEnd()));
        if (bv->theBool)
          return bv;
      }
      catch (const InterruptException &) { throw; }
      catch (const RuntimeException &) { /* ignore errors from opEqual */ }
    }
    return BOOL::falseValue;
  }
  // Second case: substr ISIN string
  if (const intrusive_ptr<STRING> str = dynamic_pointer_cast<STRING>(container))
  {
    intrusive_ptr<STRING> leftStr = dynamic_pointer_cast<STRING>(v);
    if (!leftStr)
      throw WrongTypeException(STRING::type->name, v->getType()->name, this->leftExp);
    return Value::from(str->theString.find(leftStr->theString)!=string::npos);
  }
  // Third case: sthg ISIN ideal
  if (const intrusive_ptr<IDEAL> I = dynamic_pointer_cast<IDEAL>(container))
  {
    try
    {
      // Handle separately different types for sthg: RINGELEM, INT, RAT
      if (intrusive_ptr<RINGELEM> f = dynamic_pointer_cast<RINGELEM>(v))
        return Value::from(IsElem(f->theRingElem, I->theIdeal));
      if (intrusive_ptr<INT> f = dynamic_pointer_cast<INT>(v))
        return Value::from(IsElem(f->theBigInt, I->theIdeal));
      if (intrusive_ptr<RAT> f = dynamic_pointer_cast<RAT>(v))
        return Value::from(IsElem(f->theBigRat, I->theIdeal));
    }
    catch (const ErrorInfo& err) { throw RuntimeException(message_from_cocoalib(err), beginOperatorSourcePosition, endOperatorSourcePosition); }
    catch (const InterruptReceived& /*intr*/) { throw InterruptException(beginOperatorSourcePosition, endOperatorSourcePosition); }

    throw WrongTypeException(RINGELEM::type->name  + " or " +
                             INT::type->name  + " or " +
                             RAT::type->name, v->getType()->name, this->leftExp);
  }
  // Fourth case: moduleelem ISIN module
  if (const intrusive_ptr<MODULE> M = dynamic_pointer_cast<MODULE>(container))
  {
    if (intrusive_ptr<MODULEELEM> mv = dynamic_pointer_cast<MODULEELEM>(v))
      try { return Value::from(IsElem(mv->theModuleElem, M->theModule)); }
      catch (const ErrorInfo& err) { throw RuntimeException(message_from_cocoalib(err), beginOperatorSourcePosition, endOperatorSourcePosition); }
      catch (const InterruptReceived& /*intr*/) { throw InterruptException(beginOperatorSourcePosition, endOperatorSourcePosition); }
    throw WrongTypeException(MODULEELEM::type->name,
                             v->getType()->name, this->leftExp);
  }
  // Not a recognised container type, so give error
  throw WrongTypeException(
    LIST::type->name + " or " + STRING::type->name + " or " +
    IDEAL::type->name + " or " + MODULE::type->name,
    container->getType()->name,
    this->rightExp);
}

intrusive_ptr<Value> AndExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<BOOL> left = this->leftExp->evalAs<BOOL>(runtimeEnv);
	if (!left->theBool)
		return left;
	return this->rightExp->evalAs<BOOL>(runtimeEnv);
}

intrusive_ptr<Value> OrExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<BOOL> left = this->leftExp->evalAs<BOOL>(runtimeEnv);
	if (left->theBool)
		return left;
	return this->rightExp->evalAs<BOOL>(runtimeEnv);
}

intrusive_ptr<Value> LessThanExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opLessThanMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> LessOrEqualExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opLessOrEqualMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> EqualExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opEqualMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> NotEqualExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opNotEqualMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> GreaterOrEqualExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opGreaterOrEqualMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> GreaterThanExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opGreaterThanMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> SubtractionExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opMinusMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> ProductExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opStarMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> DivisionExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return runtimeEnv->binaryOperatorDispatch(
			this->leftExp->evalAs<RightValue>(runtimeEnv),
			this->rightExp->evalAs<RightValue>(runtimeEnv),
			RuntimeEnvironment::opSlashMap,
			this->beginOperatorSourcePosition,
			this->endOperatorSourcePosition);
}

intrusive_ptr<Value> SummationExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	vector<Token>::const_iterator opIt = this->operators.begin();
	vector<intrusive_ptr<Expression> >::const_iterator expIt = this->operands.begin();
	vector<intrusive_ptr<Expression> >::const_iterator expEnd = this->operands.end();
	assert(expIt!=expEnd);
	intrusive_ptr<RightValue> currValue = (*expIt++)->evalAs<RightValue>(runtimeEnv);
	while (expIt!=expEnd) {
		assert(opIt!=operators.end());
		const Token op = *opIt++;
		const TokenType tt = op.getType();
		assert(tt==TT_PLUS || tt==TT_MINUS);
		currValue = runtimeEnv->binaryOperatorDispatch(
					currValue,
					(*expIt++)->evalAs<RightValue>(runtimeEnv),
					tt == TT_PLUS ? RuntimeEnvironment::opPlusMap : RuntimeEnvironment::opMinusMap,
					op.getBegin(),
					op.getEnd());
	}
	return currValue;
}

intrusive_ptr<Value> MultiplicationExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	vector<Token>::const_iterator opIt = this->operators.begin();
	vector<intrusive_ptr<Expression> >::const_iterator expIt = this->operands.begin();
	vector<intrusive_ptr<Expression> >::const_iterator expEnd = this->operands.end();
	assert(expIt!=expEnd);
	intrusive_ptr<RightValue> currValue = (*expIt++)->evalAs<RightValue>(runtimeEnv);
	while (expIt!=expEnd) {
		assert(opIt!=operators.end());
		const Token op = *opIt++;
		RuntimeEnvironment::DispatchMapType *map;
		switch(op.getType()) {
		case TT_STAR:
			map = &RuntimeEnvironment::opStarMap;
			break;
		case TT_SLASH:
			map = &RuntimeEnvironment::opSlashMap;
			break;
		case TT_MOD:
			map = &RuntimeEnvironment::opModMap;
			break;
		case TT_COLON:
			map = &RuntimeEnvironment::opColonMap;
			break;
		default:
			assert(false);
			throw RuntimeException("Something rather unpleasant has happened (op.getType() returned an unknown value) and I'm way too confused to continue", op.getBegin(), op.getEnd());
		}
		currValue = runtimeEnv->binaryOperatorDispatch(currValue, (*expIt++)->evalAs<RightValue>(runtimeEnv), *map, op.getBegin(), op.getEnd());
	}
	return currValue;
}

Frame *StaticEnv::VarData::tryToFindFrame(RuntimeEnvironment *runtimeEnv) const {
	assert(this->depth>StaticEnv::IN_THE_WILD);
	Frame *f;
	if (this->depth==StaticEnv::TOP_LEVEL) {
#ifdef VERBOSE_RUNTIME_DEBUG
		cout << "StaticEnv::VarData::tryToFindFrame f = toplevel\n";
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
		f = runtimeEnv->getTopLevelFrame();
		assert(f);
	} else {
		f = runtimeEnv->getCurrentFrame();
		assert(f);
#ifdef VERBOSE_RUNTIME_DEBUG
		cout << "StaticEnv::VarData::tryToFindFrame f=" << f << ", depth=" << this->depth << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
		for (int a=0; a<this->depth; ++a) {
#ifdef VERBOSE_RUNTIME_DEBUG
			cout << "StaticEnv::VarData::tryToFindFrame a=" << a << ", f=";
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
			f = f->accessLink.toCheckedPointer();
#ifdef VERBOSE_RUNTIME_DEBUG
			cout << "StaticEnv::VarData::tryToFindFrame final f=" << f << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
			if (!f)
				break;
		}
	}
	return f;
}

#ifdef C5IDE
const Frame *StaticEnv::VarData::debuggerTryToFindFrame(const Frame *topLevelFrame, const Frame *currentFrame) const {
	assert(this->depth>StaticEnv::IN_THE_WILD);
	const Frame *f;
	if (this->depth==StaticEnv::TOP_LEVEL) {
		f = topLevelFrame;
		assert(f);
	} else {
		f = currentFrame;
		assert(f);
		for (int a=0; a<this->depth; ++a) {
			f = f->accessLink.toCheckedPointer();
			if (!f)
				break;
		}
	}
	return f;
}
#endif // #ifdef C5IDE

intrusive_ptr<Value> Identifier::implEval(RuntimeEnvironment *runtimeEnv) const {
#ifdef VERBOSE_RUNTIME_DEBUG
	cout << "Identifier::eval Finding frame for " << this->identifier << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	Frame * const f = this->varData.tryToFindFrame(runtimeEnv);
	if (!f)
		throw DeadEnviromentException(this);
	if (this->varData.isCapturedValue) {
		assert(f->userdefinedFun);
		assert(static_cast<vector<intrusive_ptr<RightValue> >::size_type>(this->varData.index) < f->userdefinedFun->capturedValues.size());
		const intrusive_ptr<Value> v(f->userdefinedFun->capturedValues[this->varData.index]);
#ifdef VERBOSE_RUNTIME_DEBUG
		cout << "Identifier::eval: it's the captured value v=" << v << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
		return v;
	}
	const int index = this->varData.index;
#ifdef VERBOSE_RUNTIME_DEBUG
	cout << "Identifier::eval index=" << index << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	if (index>=0) {
		const intrusive_ptr<Value> v = f->varSlots[index].value;
#ifdef VERBOSE_RUNTIME_DEBUG
		cout << "Identifier::eval v=" << v << endl;
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
		if (dynamic_pointer_cast<LeftValue>(v)) {
#ifdef VERBOSE_RUNTIME_DEBUG
			cout << "Identifier::eval returning a ReferenceVariable\n";
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
			return new ReferenceVariable(this, f, index);
		}
	}
#ifdef VERBOSE_RUNTIME_DEBUG
		cout << "Identifier::eval returning a VariableName\n";
#endif // #ifdef VERBOSE_RUNTIME_DEBUG
	return new VariableName(runtimeEnv, f, index, this);
}

intrusive_ptr<Value> IsDefinedExp::implEval(RuntimeEnvironment *runtimeEnv) const {
	if (this->expId->varData.isCapturedValue) {
		runtimeEnv->interpreter->reportWarning("This expression will be always true, since imported values are always defined", this);
		return BOOL::trueValue;
	}
	VariableSlot *vs;
	try {
		vs = this->expId->getVariableSlot(runtimeEnv);
//JAA                vs->asRightValue(); // throws if not defined; otherwise ignore result
	} catch (VariableNotFoundException &) {
		return BOOL::falseValue;
	}
//JAA        catch (const RunTimeException&)
//JAA        { return BOOL::falseValue; }
	assert(vs->value); // getVariableSlot should throw an exception if vs->value==0 (this might happen when a package is reloaded or an indeterminate removed)
	if (dynamic_pointer_cast<VoidValue>(vs->value)) // void means optional-arg without value
		return BOOL::falseValue;
	if (const intrusive_ptr<LeftValue> lv = dynamic_pointer_cast<LeftValue>(vs->value))
		try {
			lv->asRightValue();
		} catch (const InterruptException &) {
			throw;
		} catch (const RuntimeException &) {
			return BOOL::falseValue;
		}
	return BOOL::trueValue;
}

intrusive_ptr<Value> FieldAccessExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<Value> v = this->targetExp->eval(runtimeEnv);
	if (const intrusive_ptr<LeftValue> lv = dynamic_pointer_cast<LeftValue>(v)) {
		if (const intrusive_ptr<VariableName> vn = dynamic_pointer_cast<VariableName>(lv))
			if (vn->expId->identifier[0]=='$') {
				intrusive_ptr<PackageValue> pv = intrusive_ptr_cast<PackageValue>(vn->asRightValue());
				return pv->toVariableName(runtimeEnv, this->name, this->tokName);
			}
		return new FieldAccess(lv, this->name, intrusive_ptr<const FieldAccessExpression>(this));
	}
	intrusive_ptr<RECORD> rv = dynamic_pointer_cast<RECORD>(v);
	if (!rv)
		throw WrongTypeException(
				RECORD::type->name+" or "+PackageValue::type->name,
				intrusive_ptr_cast<RightValue>(v)->getType()->name,
				this->targetExp->getBegin(),
				this->targetExp->getEnd());
	intrusive_ptr<RightValue> fieldValue = rv->getField(this->name, this->tokName.getBegin(), this->tokName.getEnd());
	if (!fieldValue)
		throw FieldNotFoundException(this->name, this->tokName.getBegin(), this->tokName.getEnd(), rv);
	return fieldValue;
}

intrusive_ptr<Value> IndexedAccessExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<Value> v = this->targetExp->eval(runtimeEnv);
	const CharPointer targetExpBegin = this->targetExp->getBegin();
	CharPointer targetExpEnd = this->targetExp->getEnd();
        // 2020-11-04 Loop below is a little complicated so that it can give helpful err mesgs (redmine 1527)
        for (const intrusive_ptr<Expression>& indexExp: this->indices) // added reference to avoid compiler warning
        {
          const CharPointer indexExpBegin = indexExp->getBegin();
          const CharPointer indexExpEnd = indexExp->getEnd();
          if (!v->asRightValue()->canBeIndexedByBigInt() && !v->asRightValue()->canBeIndexedByString())
            throw CannotBeIndexedException(v->asRightValue()->getType()->name, indexExpBegin, indexExpEnd);
          const intrusive_ptr<RightValue> indexValue = indexExp->evalAs<RightValue>(runtimeEnv);
          if (intrusive_ptr<INT> N = dynamic_pointer_cast<INT>(indexValue))
            v = v->indexedByBigInt(N, targetExpBegin, targetExpEnd, indexExpBegin, indexExpEnd);
          else if (intrusive_ptr<STRING> str = dynamic_pointer_cast<STRING>(indexValue))
            v = v->indexedByString(str, targetExpBegin, targetExpEnd, indexExpBegin, indexExpEnd);
          else
          {
            // Index was neither STRING nor INT, so throw error -- err mesg depends on whether INT or STRING index was expected.
            if (v->asRightValue()->canBeIndexedByBigInt())
              throw WrongTypeException(INT::type->name, indexValue->getType()->name, indexExpBegin, indexExpEnd);
            else
              throw WrongTypeException(STRING::type->name, indexValue->getType()->name, indexExpBegin, indexExpEnd);
          }
          targetExpEnd = indexExpEnd;
	}
	return v;
}

intrusive_ptr<Value> RecordExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<RECORD> r(new RECORD());
	for (const RecordField &rf: this->fields) {
		intrusive_ptr<const Identifier> idExp = intrusive_ptr_cast<const Identifier>(rf.name);
                RECORD::CheckFieldName(idExp->identifier, idExp);
		if (r->getFieldNoCheck(idExp->identifier))
			throw RuntimeException("Duplicate field", idExp);
		r->setFieldNoCheck(idExp->identifier, rf.initExp->evalAs<RightValue>(runtimeEnv));
	}
	return r;
}

intrusive_ptr<Value> ListExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<LIST> l(new LIST());
	for (intrusive_ptr<const Expression> e: this->exps) {
		l->addValue(e->evalAs<RightValue>(runtimeEnv));
	}
	return l;
}

intrusive_ptr<Value> CartesianProductExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<LIST> previous(new LIST());
	previous->addValue(new LIST(this->operands.size()));
	intrusive_ptr<LIST> current;
	for (intrusive_ptr<const Expression> e: this->operands)
        {
          intrusive_ptr<LIST> l = e->evalAs<LIST>(runtimeEnv);
          const LIST::ContainerType::size_type lSize = l->size();
          const LIST::ContainerType::size_type previousSize = previous->size();
          current = new LIST();
          for (LIST::ContainerType::size_type a=0; a<previousSize; ++a) {
            for (LIST::ContainerType::size_type b=0; b<lSize; ++b) {
              intrusive_ptr<LIST> cloned = intrusive_ptr_cast<LIST>(previous->getValue(a)->clone());
              cloned->addValue(l->getValue(b));
              current->addValue(cloned);
            }
          }
          previous = current;
	}
	return current;
}

intrusive_ptr<Value> ExpSuchThatIdInExpAndExp::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<LIST> list = this->exp2->evalAs<LIST>(runtimeEnv);
	intrusive_ptr<LIST> result(new LIST());
	Frame *f = runtimeEnv->pushIterationFrame(this);
	BOOST_SCOPE_EXIT( (runtimeEnv) ) {
		runtimeEnv->popFrame();
	} BOOST_SCOPE_EXIT_END
	VariableSlot &varSlot = f->varSlots.front();
	const int size = list->size();
	for (int a=0; a<size; ++a) {
		intrusive_ptr<RightValue> sourceElem = list->getValue(a);
		assert(sourceElem);
		varSlot.value = sourceElem;
		if (!this->optionalExp3 || this->optionalExp3->evalAs<BOOL>(runtimeEnv)->theBool)
			result->addValue(this->exp1->evalAs<RightValue>(runtimeEnv));
	}
	return result;
}

intrusive_ptr<Value> IdInExpSuchThatExp::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<LIST> list = this->exp1->evalAs<LIST>(runtimeEnv);
	Frame *f = runtimeEnv->pushIterationFrame(this);
	BOOST_SCOPE_EXIT( (runtimeEnv) ) {
		runtimeEnv->popFrame();
	} BOOST_SCOPE_EXIT_END
	VariableSlot &varSlot = f->varSlots.front();
	intrusive_ptr<LIST> result(new LIST());
	const int size = list->size();
	for (int a=0; a<size; ++a) {
		intrusive_ptr<RightValue> elem = list->getValue(a);
		assert(elem);
		varSlot.value = elem;
		if (this->exp2->evalAs<BOOL>(runtimeEnv)->theBool)
			result->addValue(elem);
	}
	return result;
}

intrusive_ptr<Value> DotDotExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	const intrusive_ptr<LIST> l(new LIST());
  const intrusive_ptr<RightValue>
    lo = this->leftExp->evalAs<RightValue>(runtimeEnv),
    hi = this->rightExp->evalAs<RightValue>(runtimeEnv);
  std::string ErrMesg = "range operator \"dot-dot\": only between INTs or between indets of the same ring";
  intrusive_ptr<INT> NL, NR;
  intrusive_ptr<RINGELEM> xL, xR;
	if ((NL=dynamic_pointer_cast<INT>(lo)) && (NR=dynamic_pointer_cast<INT>(hi)))
  {
		BigInt lower(NL->theBigInt), upper(NR->theBigInt);
                if (upper-lower > 10000000)  // BUG??? arbitrary max size
                  throw RuntimeException("dot-dot range limited to 10000000 values",this);
		for (; lower<=upper; ++lower) l->addValue(new INT(lower));
	}
  else
    if ((xL=dynamic_pointer_cast<RINGELEM>(lo)) && (xR=dynamic_pointer_cast<RINGELEM>(hi)))
    {
      const ring theRing(owner(xL->theRingElem));
      if (theRing!=owner(xR->theRingElem)) throw RuntimeException(ErrMesg,this);
      if (!IsPolyRing(theRing)) throw RuntimeException(ErrMesg, this);
      const PolyRing polyRing = theRing;
      long lower, upper;
      if (!IsIndet(lower,xL->theRingElem)) throw RuntimeException(ErrMesg,this->leftExp);
      if (!IsIndet(upper,xR->theRingElem)) throw RuntimeException(ErrMesg,this->rightExp);
      for (; lower<=upper; ++lower)
        l->addValue(new RINGELEM(indet(polyRing, lower)));
    }
    else
      throw RuntimeException(ErrMesg, this);
	return l;
}


intrusive_ptr<Value> InvocationExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	intrusive_ptr<RightValue> target = this->targetExp->evalAs<RightValue>(runtimeEnv);
	intrusive_ptr<FUNCTION> fun = dynamic_pointer_cast<FUNCTION>(target);
	if (!fun)
		throw WrongTypeException(FUNCTION::type->name,
                             target->getType()->name + " (maybe you forgot \"*\"?)",
                             this->targetExp);
	return fun->eval(this, runtimeEnv);
}


namespace {
	intrusive_ptr<Value> dontKnow(const intrusive_ptr<const Statement> s) {
		throw RuntimeException("I don't know (yet) how to execute this statement", s->getBegin(), s->getEnd());
	}
}

void UseStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	const intrusive_ptr<RING> ring = this->ringDefinition->eval(runtimeEnv);
	string ringName;
	if (this->identifier) {
		ringName = this->identifier->identifier;
		const intrusive_ptr<LeftValue> leftValue = this->identifier->evalAs<LeftValue>(runtimeEnv);
		if (leftValue->assignmentNeedsOwnership())
			leftValue->obtainOwnership();
		leftValue->assign(ring, this->ringDefinition->getBegin(), this->ringDefinition->getEnd(), runtimeEnv);
	}
	if (runtimeEnv->currentRing)
		runtimeEnv->currentRing->removeInjectedIndeterminates(runtimeEnv);
	runtimeEnv->currentRing = ring;
	runtimeEnv->getTopLevelFrame()->varSlots[runtimeEnv->currentRingSlot].value = ring;
	ring->injectIndeterminates(runtimeEnv, this->ringDefinition, ringName);
}

void UsingStatement::implExecute(RuntimeEnvironment *) {
	dontKnow(this);
}

void DescribeStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	this->exp->evalAs<RightValue>(runtimeEnv)->describe(runtimeEnv->getOutputStream());
}

void HelpStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	try {
		ostringstream os;
		OnlineHelp::PrintMan(os, this->topic);
		runtimeEnv->getOutputStream()->print(os.str())->flush();
	} catch (const std::exception& err) {
		throw RuntimeException(err.what(), this);
	}
}

intrusive_ptr<RING> RingDefinition::eval(InterpreterNS::RuntimeEnvironment *runtimeEnv) const {
	try {
		intrusive_ptr<RING> ring = this->identifier->evalAs<RING>(runtimeEnv);
		if (this->optionalExp) {
			const intrusive_ptr<RightValue> quotient = this->optionalExp->evalAs<RightValue>(runtimeEnv);
			if (const intrusive_ptr<INT> N = dynamic_pointer_cast<INT>(quotient))
                          ring = new RING(NewQuotientRing(ring->theRing, ideal(RingElem(ring->theRing, N->theBigInt))));
			else if (const intrusive_ptr<RINGELEM> poly = dynamic_pointer_cast<RINGELEM>(quotient))
				ring = new RING(NewQuotientRing(ring->theRing, ideal(poly->theRingElem)));
			else if (const intrusive_ptr<IDEAL> ideal = dynamic_pointer_cast<IDEAL>(quotient))
				ring = new RING(NewQuotientRing(ring->theRing, ideal->theIdeal));
			else
				throw WrongTypeException(INT::type->name + ", " + RINGELEM::type->name + " or " + IDEAL::type->name, quotient->getType()->name, this->optionalExp);
		}
		if (!this->indeterminates.empty()) {
			vector<symbol> symbols;
			for (intrusive_ptr<IndeterminateDeclaration> ind: this->indeterminates)
                        {
                          if (ind->ranges.empty()) {
                            symbols.push_back(symbol(ind->identifier->identifier));
                            continue;
                          }
                          vector<pair<long, long> > ranges;
                          typedef pair<intrusive_ptr<Expression>, intrusive_ptr<Expression> > ExpPair;
                          for (const ExpPair &p: ind->ranges)
                          {
                            intrusive_ptr<INT> lower = p.first->evalAs<INT>(runtimeEnv);
                            long l1;
                            if (!IsConvertible(l1, lower->theBigInt))
                              throw RuntimeException("The index does not fit a machine-integer", p.first);
                            intrusive_ptr<INT> upper = p.second ? p.second->evalAs<INT>(runtimeEnv) : lower;
                            long l2;
                            if (!IsConvertible(l2, upper->theBigInt))
                              throw RuntimeException("The index does not fit a machine-integer", p.second);
                            if (l2<l1)
                              throw RuntimeException("The upper bound is smaller than the lower bound", p.first->getBegin(), p.second->getEnd());
                            ranges.push_back(make_pair(l1, l2));
                          }
                          vector<long> indices;
                          for (vector<pair<long, long> >::size_type a=0; a<ranges.size(); ++a)
                          indices.push_back(ranges[a].first);
                          for (;;) {
                            symbols.push_back(symbol(ind->identifier->identifier, indices));
                            int incIndex = indices.size()-1;
                            while (incIndex>=0)
                              if ( ++indices[incIndex]>ranges[incIndex].second ) {
                                indices[incIndex] = ranges[incIndex].first;
                                --incIndex;
                              } else
                                break;
                            if (incIndex==-1)
                              break;
                          }
			}
			/*ostringstream os;
			for (const symbol &s: symbols)
				os << ' ' << s;
			runtimeEnv->getOutputStream()->print("Symbols:")->println(os.str());*/
			switch (this->maybeOrderTT) {
			case TT_LEX:
				ring = new RING(NewPolyRing(ring->theRing, symbols, lex));
				break;
			case TT_DEGLEX:
				ring = new RING(NewPolyRing(ring->theRing, symbols, StdDegLex));
				break;
			case TT_DEGREVLEX:
				ring = new RING(NewPolyRing(ring->theRing, symbols, StdDegRevLex));
				break;
			case TT_XEL:
				ring = new RING(NewPolyRing(ring->theRing, symbols, xel));
				break;
			default:
				ring = new RING(NewPolyRing(ring->theRing, symbols));
			}
		}
		return ring;
	}
        catch (const ErrorInfo& err)
        {
          runtimeEnv->announceCLE(err);
          throw RuntimeException(message_from_cocoalib(err), this);
	}
        catch (const InterruptReceived& /*intr*/)
        {
          throw InterruptException(this);
        }
}

void RingAssignStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	const intrusive_ptr<LeftValue> leftValue = this->leftExp->evalAs<LeftValue>(runtimeEnv);
	runtimeEnv->interpreter->checkForInterrupts(this->leftExp);
	const intrusive_ptr<RING> ring = this->ringDef->eval(runtimeEnv);
	if (leftValue->assignmentNeedsOwnership())
		leftValue->obtainOwnership();
	leftValue->assign(ring, this->ringDef->getBegin(), this->ringDef->getEnd(), runtimeEnv);
}

void TimeStatement::implExecute(RuntimeEnvironment *) {
	dontKnow(this);
}

void EmptyStatement::implExecute(RuntimeEnvironment *) {
	// nothing to do here, what a lovely statement ;-)
}

void SkipStatement::implExecute(RuntimeEnvironment *) {
	// nothing to do here, what a lovely statement ;-)
}

VariableSlot *Identifier::getVariableSlot(RuntimeEnvironment *runtimeEnv) {
	Frame * const f = this->varData.tryToFindFrame(runtimeEnv);
	if (!f)
		throw DeadEnviromentException(this);
	int index = this->varData.index;
	if (index<0) {
		assert(f==runtimeEnv->getTopLevelFrame());
		index = runtimeEnv->slotFor(this->identifier);
		if (index<0)
			throw VariableNotFoundException(this, runtimeEnv);
		this->varData.index = index;
	}
	VariableSlot &vs = f->varSlots[index];
	if (!vs.value) { // this can only happen when a package has been reloaded (and some of the previous version members are not defined anymore) or an indeterminate has been removed (by Use-ing another ring)
		assert(f==runtimeEnv->getTopLevelFrame());
		throw VariableNotFoundException(this, runtimeEnv);
	}
	return &vs;
}

void ProtectStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	if (this->expId->varData.isCapturedValue)
		throw RuntimeException("Cannot protect a(n imported) value", this->expId);
	VariableSlot *vs = this->expId->getVariableSlot(runtimeEnv);
	if (vs->isSystemProtected())
		throw RuntimeException("This variable is already system-protected", this->expId);
	if (vs->isProtected())
        {
          string reason;
          if (!vs->getProtectionReason().empty())
            reason = " ("+vs->getProtectionReason()+")";
		throw RuntimeException("Variable already protected"+reason, this->expId);
        }
	string reason;
	if (this->optExp)
		reason = this->optExp->evalAs<STRING>(runtimeEnv)->theString;
	vs->protect(reason);
}

void UnprotectStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	if (this->expId->varData.isCapturedValue)
		throw RuntimeException("Cannot unprotect a(n imported) value", this->expId);
	VariableSlot *vs = this->expId->getVariableSlot(runtimeEnv);
	if (vs->isSystemProtected())
    CheckProtection(vs, this->expId, runtimeEnv, "unprotect");
///		throw RuntimeException("Cannot unprotect a system-protected variable", this->expId);
	if (vs->isInPackage())
///		throw RuntimeException("Cannot unprotect a package-exported variable", this->expId);
    CheckProtection(vs, this->expId, runtimeEnv, "unprotect");
	if (!vs->isProtected())
		runtimeEnv->interpreter->reportWarning("This variable was not protected", this->expId);
	vs->unprotect();
}

void IfStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
  for (const IfBranch &branch: this->branches)
		if (!branch.optExp || branch.optExp->evalAs<BOOL>(runtimeEnv)->theBool) {
			branch.statements->execute(runtimeEnv);
			return;
		}
}

void CiaoOrQuitStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	runtimeEnv->getOutputStream()->println("Ciao!");
	throw Ciao();
}

void SourceStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	const string filename = this->exp->evalAs<STRING>(runtimeEnv)->theString;
	try {
		runtimeEnv->interpreter->readAndExecute(filename, false, this->ttype!=TT_LOAD);
	} catch (const RuntimeException &) {
		throw;
	} catch (const BaseException &be) {
		throw RuntimeException(be.reason, this->getBegin(), this->getEnd());
	}
}

void SourceRegionStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
  try
  {
    const long FromLine = ConvertTo<long>(this->expFromLine->evalAs<INT>(runtimeEnv)->theBigInt);
    const long FromChar = ConvertTo<long>(this->expFromChar->evalAs<INT>(runtimeEnv)->theBigInt);
    const long ToLine = ConvertTo<long>(this->expToLine->evalAs<INT>(runtimeEnv)->theBigInt);
    const long ToChar = ConvertTo<long>(this->expToChar->evalAs<INT>(runtimeEnv)->theBigInt);

    const string filename = this->expFileName->evalAs<STRING>(runtimeEnv)->theString;

    runtimeEnv->interpreter->readAndExecute(filename, false, true, FromLine,FromChar, ToLine,ToChar);
  }
  catch (const RuntimeException &)
  {
    throw;
  }
  catch (const BaseException &be)
  {
    throw RuntimeException(be.reason, this->getBegin(), this->getEnd());
  }
  catch (const InterruptReceived& /*intr*/)
  {
    throw InterruptException(this->getBegin(), this->getEnd());
  }
}

void EvalStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	intrusive_ptr<RightValue> v = this->exp->evalAsRightValueVoidIsOk(runtimeEnv);
	assert(v);
	if (!dynamic_pointer_cast<VoidValue>(v)) {
		Frame *currentFrame = runtimeEnv->getCurrentFrame();
		if (currentFrame->invocationExp)
			throw RuntimeException("Ignored value inside fn-proc!  All values must be assigned"
					" or operated upon by a command/procedure (maybe you forgot to assign, print, or return it?)", this->getBegin(), this->getEnd());
		else
			runtimeEnv->getTopLevelFrame()->varSlots[runtimeEnv->itSlot].value = v;
		const intrusive_ptr<OSTREAM> output(runtimeEnv->getOutputStream());
		output->print(runtimeEnv, v, this->exp);
		output->newline();
#ifdef PRINT_REFCOUNTS
		output->print(v->refCountAsString());
#endif // #ifdef PRINT_REFCOUNTS
	}
}

void PrintStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	intrusive_ptr<OSTREAM> out = this->onExp ?
			this->onExp->evalAs<OSTREAM>(runtimeEnv) :
			runtimeEnv->getOutputStream();
	for (intrusive_ptr<const Expression> e: this->exps)
        {
          intrusive_ptr<RightValue> v = e->evalAs<RightValue>(runtimeEnv);
          out->print(runtimeEnv, v, e);
	}
	if (this->ttype==TT_PRINTLN)
		out->newline();
	out->flush();
}

void ForStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	BigInt begin(this->beginExp->evalAs<INT>(runtimeEnv)->theBigInt);
	BigInt end(this->endExp->evalAs<INT>(runtimeEnv)->theBigInt);
	BigInt step(1);
	if (this->stepExp) {
		step = this->stepExp->evalAs<INT>(runtimeEnv)->theBigInt;
		if (IsZero(step))
			throw RuntimeException("Step value cannot be zero", this->stepExp);
	}
	bool goingUp = step>0;
	Frame *f = runtimeEnv->pushIterationFrame(this);
	BOOST_SCOPE_EXIT( (runtimeEnv) ) {
		runtimeEnv->popFrame();
	} BOOST_SCOPE_EXIT_END
	VariableSlot &varSlot = f->varSlots.front();
	while( (goingUp && begin<=end) || (!goingUp && begin>=end)) {
		runtimeEnv->interpreter->checkForInterrupts(this);
		varSlot.value = new INT(begin);
		try {
			this->statements->execute(runtimeEnv);
		} catch (const Break &b) {
			if (b.label.length() && this->label!=b.label)
				throw;
			return;
		} catch (const Continue &c) {
			if (c.label.length() && this->label!=c.label)
				throw;
		}
		begin += step;
	}
}

void Statements::implExecute(RuntimeEnvironment *runtimeEnv) {
  for (intrusive_ptr<Statement> s: this->statements)
		s->execute(runtimeEnv);
}

void ForeachStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	intrusive_ptr<LIST> list = this->inExp->evalAs<LIST>(runtimeEnv);
	Frame *f = runtimeEnv->pushIterationFrame(this);
	BOOST_SCOPE_EXIT( (runtimeEnv) ) {
		runtimeEnv->popFrame();
	} BOOST_SCOPE_EXIT_END
	VariableSlot &varSlot = f->varSlots.front();
	const int size = list->size();
	for (int a=0; a<size; ++a) {
		runtimeEnv->interpreter->checkForInterrupts(this);
		intrusive_ptr<RightValue> sourceElem = list->getValue(a);
		assert(sourceElem);
		varSlot.value = sourceElem;
		try {
			this->statements->execute(runtimeEnv);
		} catch (const Break &b) {
			if (b.label.length() && this->label!=b.label)
				throw;
			return;
		} catch (const Continue &c) {
			if (c.label.length() && this->label!=c.label)
				throw;
			continue;
		}
	}
}

void TryStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	try	{
		this->tryStatements->execute(runtimeEnv);
	} catch (InterruptException &) {
		throw;
	} catch (RuntimeException &exc) {
		Frame *f = runtimeEnv->pushIterationFrame(this);
		BOOST_SCOPE_EXIT( (runtimeEnv) ) {
			runtimeEnv->popFrame();
		} BOOST_SCOPE_EXIT_END
		VariableSlot &varSlot = f->varSlots.front();
		varSlot.value = new ERROR_MESG(exc.reason);
		this->uponErrorStatements->execute(runtimeEnv);
	}
}

void BlockStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	this->statements->execute(runtimeEnv);
}

void AliasStatement::implExecute(RuntimeEnvironment *) {
	//if (this->statements)
	//	this->statements->execute(runtimeEnv);
}

void WhileStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	while(this->exp->evalAs<BOOL>(runtimeEnv)->theBool) {
		runtimeEnv->interpreter->checkForInterrupts(this);
		try {
			this->statements->execute(runtimeEnv);
		} catch (const Break &b) {
			if (b.label.length() && this->label!=b.label)
				throw;
			return;
		} catch (const Continue &c) {
			if (c.label.length() && this->label!=c.label)
				throw;
		}
	}
}

void RepeatUntilStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	const bool forever = !this->optExp;
	do {
		runtimeEnv->interpreter->checkForInterrupts(this);
		try {
			this->statements->execute(runtimeEnv);
		} catch (const Break &b) {
			if (b.label.length() && this->label!=b.label)
				throw;
			return;
		} catch (const Continue &c) {
			if (c.label.length() && this->label!=c.label)
				throw;
		}
	} while (forever || !this->optExp->evalAs<BOOL>(runtimeEnv)->theBool);
}


void PackageStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	Frame * const tlFrame = runtimeEnv->getTopLevelFrame();
	int tlIndex = runtimeEnv->slotFor(this->name);
	if (tlIndex>=0) {
		assert(tlFrame->varSlots[tlIndex].value);
		intrusive_ptr<PackageValue> oldPkg = intrusive_ptr_cast<PackageValue>(tlFrame->varSlots[tlIndex].value);
		runtimeEnv->interpreter->reportWarning("Package "+this->name+" has been redefined", this->beginSourcePosition, this->tokName.getEnd());
		for (int index: oldPkg->ownedTopLevelIndices)
                {
                  assert(index>=0 && static_cast<vector<VariableSlot>::size_type >(index)<tlFrame->varSlots.size());
                  VariableSlot &vs = tlFrame->varSlots[index];
                  assert(vs.value);
                  vs.resetFlags();
                  vs.value = nullptr;
		}
	}
	const string protectionReason = "Owned by package "+this->name;
	intrusive_ptr<PackageValue> newPackage(new PackageValue(this));
	const int pkgIndex = runtimeEnv->setTopLevelVar(this->name, newPackage, VariableSlot::VSF_SystemProtected);
/// setting the reason IN the package?
	tlFrame->varSlots[pkgIndex].setProtectionReason(protectionReason);
	ResolvePackageNamesVisitor rpnv(*this);
	typedef pair<string, intrusive_ptr<const AssignmentStatement> > Pair;
	vector<Pair> memberVars;
	for (intrusive_ptr<Statement> statement: this->statements->statements)
        {
          if (dynamic_pointer_cast<const AliasStatement>(statement))
            continue;
          if (const intrusive_ptr<const AssignmentStatement> assignment = dynamic_pointer_cast<const AssignmentStatement>(statement)) {
            intrusive_ptr<Identifier> idExp = intrusive_ptr_cast<Identifier>(assignment->leftExp);
            const string identifier = this->name+"."+idExp->identifier;
            memberVars.push_back(make_pair(identifier, assignment));
            const int index = runtimeEnv->setTopLevelVar(identifier, VoidValue::theInstance, VariableSlot::VSF_None);
            newPackage->ownedTopLevelIndices.push_back(index);
          } else if (const intrusive_ptr<const DefineStatement> define = dynamic_pointer_cast<const DefineStatement>(statement)) {
            define->funDecl->accept(&rpnv);
            const int index = runtimeEnv->setTopLevelVar(this->name+"."+define->name, new UserDefinedFunction(runtimeEnv, define->funDecl), VariableSlot::VSF_SystemProtected);
            newPackage->ownedTopLevelIndices.push_back(index);
            tlFrame->varSlots[index].setProtectionReason(protectionReason);
          } else
            assert(false);
	}
	for (const Token &tokName: this->exportedNames)
        {
          const string name(tokName.lexeme());
          int index = runtimeEnv->slotFor(name);
          if (index>=0) {
            if (tlFrame->varSlots[index].isProtected()) {
              runtimeEnv->interpreter->reportError("Cannot export the name \""+name+"\" because there is a top-level protected variable with the same name", tokName.getBegin(), tokName.getEnd());
              continue;
            }
            if (tlFrame->varSlots[index].isInPackage()) {
              intrusive_ptr<PackageValue> owner = FindOwnerPkg(runtimeEnv, name);
              assert(owner);
              runtimeEnv->interpreter->reportError("Cannot export the name \""+name+"\" because the same name was already exported by package "+owner->pkgName, tokName.getBegin(), tokName.getEnd());
              continue;
            }
          }
          const int indexFullyQualified = runtimeEnv->slotFor(this->name+"."+name);
          assert(indexFullyQualified>=0);
          index = runtimeEnv->setTopLevelVar(name, VoidValue::theInstance, VariableSlot::VSF_Package);
          newPackage->ownedTopLevelIndices.push_back(index);
          intrusive_ptr<Identifier> id(new Identifier(tokName));
          id->varData.depth = StaticEnv::TOP_LEVEL;
          id->varData.index = indexFullyQualified;
          tlFrame->varSlots[index].value = new VariableName(runtimeEnv, tlFrame, indexFullyQualified, id);
          tlFrame->varSlots[indexFullyQualified].setHasBeenExported();
	}
	for (const Pair &assignment: memberVars)
        {
          const intrusive_ptr<Expression> rightExp = assignment.second->rightExp;
          rightExp->accept(&rpnv);
          runtimeEnv->setTopLevelVar(assignment.first, rightExp->evalAs<RightValue>(runtimeEnv), VariableSlot::VSF_None);
	}
}

void DefineStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	const intrusive_ptr<UserDefinedFunction> fun(new UserDefinedFunction(runtimeEnv, this->funDecl));
	int slotIndex = runtimeEnv->slotFor(this->name);
	if (slotIndex<0)
		runtimeEnv->setTopLevelVar(this->name, fun);
	else {
		VariableSlot *vs = &(runtimeEnv->getTopLevelFrame()->varSlots[slotIndex]);
		CheckProtection(vs, this->expId, runtimeEnv, "define function ");
		vs->value = fun;
	}
}

intrusive_ptr<Value> FullyQualifiedIdentifier::implEval(RuntimeEnvironment *runtimeEnv) const {
	const int slot = runtimeEnv->slotFor(this->pkgName);
	if (slot<0)
		throw RuntimeException("Cannot find package \""+this->pkgName+"\"", this->tokPkgname);
	assert(runtimeEnv->getTopLevelFrame()->varSlots[slot].value);
	const intrusive_ptr<PackageValue> pv(intrusive_ptr_cast<PackageValue>(runtimeEnv->getTopLevelFrame()->varSlots[slot].value));
	return pv->toVariableName(runtimeEnv, this->id, this->tokId);
}

intrusive_ptr<Value> LambdaExpression::implEval(RuntimeEnvironment *runtimeEnv) const {
	return new UserDefinedFunction(runtimeEnv, this->funDecl);
}

void ReturnStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	throw Return(this->exp ? this->exp->evalAs<RightValue>(runtimeEnv) : VoidValue::theInstance);
}

void BreakStatement::implExecute(RuntimeEnvironment *) {
	throw Break(this->label);
}

void ContinueStatement::implExecute(RuntimeEnvironment *) {
	throw Continue(this->label);
}

void AssignmentStatement::implExecute(RuntimeEnvironment *runtimeEnv) {
	const intrusive_ptr<LeftValue> leftValue = this->leftExp->evalAs<LeftValue>(runtimeEnv);
	runtimeEnv->interpreter->checkForInterrupts(this->leftExp);
        leftValue->checkAssignable(runtimeEnv);  // see redmine 435
	intrusive_ptr<RightValue> rightValue = this->rightExp->evalAs<RightValue>(runtimeEnv);
	if (leftValue->assignmentNeedsOwnership())
		leftValue->obtainOwnership();
	leftValue->assign(rightValue, this->rightExp->getBegin(), this->rightExp->getEnd(), runtimeEnv);
}

} // namespace AST

namespace LexerNS {

using namespace InterpreterNS;

void ErrorReporter::reportError(const RuntimeException &exception) {
	this->reportError(exception.reason, exception.from, exception.to);
	const vector<SnapshotFrame>::size_type snapshotSize = exception.snapshot.size();
	if (snapshotSize) {
		this->printContext();
		int nesting=-1;
		bool dots=false;
		CharPointer from = exception.from;
		CharPointer to = exception.to;
		for (vector<SnapshotFrame>::size_type a=0; a<snapshotSize; ++a) {
			const SnapshotFrame &frame = exception.snapshot[a];
			++nesting;
			if ( (snapshotSize-a>10) && a>10 ) {
				if (dots) {
					--nesting;
					continue;
				}
				this->outputStream->print(string(nesting, ' '))->print("...\n");
				dots = true;
				continue;
			}
			if (nesting) {
				this->outputStream->print(string(nesting, ' '));
				this->printCalledBy();
			}
			const intrusive_ptr<const FunctionDeclaration> fnDecl(intrusive_ptr_cast<const FunctionDeclaration>(frame.block));
			if (fnDecl->fnName.length()) {
				this->outputStream->print("function ");
				this->printBold(fnDecl->fnName);
			} else
				this->printBold("anonymous-function");
			if (!this->reportLineNumberWhenMeaningful(from, to, false, false)) {
				if (frame.block)
					this->outputStream->print(" (previously defined at the prompt)");
				else
					this->outputStream->print(" at top-level");
			}
			this->outputStream->newline();
			from = frame.invocationExp->getBegin();
			to = frame.invocationExp->getEnd();
		}
		this->outputStream->print(string(++nesting, ' '))->print("called");
		if (!this->reportLineNumberWhenMeaningful(from, to, false, false))
			this->outputStream->print(" at top-level");
		this->outputStream->newline();
	}
}

  void ErrorReporter::reportInterrupt(const InterpreterNS::InterruptException &intr)
  {
    this->implReportInterrupt(intr.from, intr.to);
  }


} // namespace LexerNS

} // namespace CoCoA
