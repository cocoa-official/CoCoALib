//   Copyright (c)  2015  John Abbott,  Anna M Bigatti
//   Orig author:  2009  Giovanni Lagorio
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

#ifdef C5IDE

#include <QFileDialog>
#include <QFileInfo>
#include <QList>
// #include <QPlastiqueStyle>
#include <QFontDialog>
#include <QCoreApplication>
#include <QSettings>
#include <QStyleFactory>

#include "Banner.H"
#include "C5.H"
#include "C5Utils.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/interrupt.H"
#include "qce-config.h"
#include "qeditor.h"
#include "qlanguagefactory.h"
#include "qlinemarksinfocenter.h"
#include "qeditorinputbinding.h"
#include "widgets/qgotolinepanel.h"
#include "widgets/qfoldpanel.h"
#include "widgets/qstatuspanel.h"
#include "widgets/qlinemarkpanel.h"
#include "widgets/qsearchreplacepanel.h"
#include "widgets/qlinenumberpanel.h"
#include "widgets/qlinechangepanel.h"

using namespace std;
using namespace boost;
using namespace CoCoA;
using namespace CoCoA::LexerNS;
using namespace CoCoA::AST;
using namespace CoCoA::InterpreterNS;

void initQCodeEdit() {
	// the following are needed to statically (and successfully ;-) ) link QCodeEdit
	// note: this function must be outside any namespace (because of the macro Q_INIT_RESOURCE)
	QGotoLinePanel::_register();
	QFoldPanel::_register();
	QStatusPanel::_register();
	QLineMarkPanel::_register();
	QSearchReplacePanel::_register();
	QLineNumberPanel::_register();
	QLineChangePanel::_register();
	Q_INIT_RESOURCE(Edyuk);
}

namespace {
	const QString SETTING_FONT_FAMILY("font/family");
	const QString SETTING_FONT_POINTSIZE("font/pointSize");
	const QString SETTING_FONT_WEIGHT("font/weight");
	const QString SETTING_FONT_ITALIC("font/italic");
	const QString SETTING_COLORSCHEME_NAME("colorScheme/name");
	const QString SETTING_EMACS_LIKE_KEYBINDINGS("emacs/keyBindings");
}

namespace CoCoA {

namespace AST {
	intrusive_ptr<Value> Expression::eval(RuntimeEnvironment *re) const {
		Interpreter * const interpreter = re->interpreter;
		bool doStepOver = false;
		if (interpreter->singleStepExecution && !CoCoA::GetAndResetSignalReceived() && !this->skipDebugging())
			doStepOver = interpreter->pause(this);
		intrusive_ptr<Value> result;
		try	{
			result = this->implEval(re);
			if (doStepOver && interpreter->doStepOver)
				interpreter->singleStepExecution = true;
		} catch (InterruptException &) {
			throw;
		} catch (RuntimeException &) {
			if (doStepOver && interpreter->doStepOver)
				interpreter->singleStepExecution = true;
			throw;
		}
		return result;
	}

	bool ParsedObject::skipDebugging() const {
		return false;
	}

	bool EvalStatement::skipDebugging() const {
			return true;
	}

	bool Statements::skipDebugging() const {
			return true;
	}

	bool EmptyStatement::skipDebugging() const {
			return true;
	}

	bool InvocationExpression::skipDebugging() const {
			return this->targetExp->skipDebugging();
	}

	void Statement::execute(RuntimeEnvironment *re) {
		Interpreter * const interpreter = re->interpreter;
		bool doStepOver = false;
		if (interpreter->singleStepExecution && !CoCoA::GetAndResetSignalReceived() && !this->skipDebugging()) {
			doStepOver = interpreter->pause(this);
			if (doStepOver) {
				if (dynamic_cast<BreakStatement *>(this) || dynamic_cast<ContinueStatement *>(this))
					re->getCurrentFrame()->singleStepOnPop = true;
				if (dynamic_cast<ReturnStatement *>(this))
					re->getCurrentFunctionFrame()->singleStepOnPop = true;
			}
		}
		try	{
			this->implExecute(re);
			if (doStepOver && interpreter->doStepOver)
				interpreter->singleStepExecution = true;
		} catch (InterruptException &) {
			throw;
		} catch (RuntimeException &) {
			if (doStepOver && interpreter->doStepOver)
				interpreter->singleStepExecution = true;
			throw;
		}
	}
} // namespace AST

namespace {
	const int COL_NAME = 0;
	const int COL_KIND = 1;
	const int COL_VALUE = 2;
}

namespace InterpreterNS {
	bool Interpreter::pause(intrusive_ptr<const ParsedObject> po) {
		assert(std::this_thread::get_id()!=IDE::Console::guiThreadId);
		const Frame *tl = this->runtimeEnvironment.getTopLevelFrame();
		for(Frame *f=this->runtimeEnvironment.getCurrentFrame(); f>=tl;)
			(f--)->singleStepOnPop = false;
		this->doStepOver = false;
		this->pausedOn = po;
		InterpreterStatus status = this->status;
		assert(status==IS_RUNNING || status==IS_RUNNING_BUILTIN);
		this->status = IS_PAUSED;
		std::unique_lock<std::mutex> lock(this->mut);
		this->mustResume = false;
		do
			this->condVar.wait(lock);
		while (!this->mustResume);
		this->pausedOn = 0;
		this->status = status;
		this->runtimeEnvironment.getOutputStream()->flush(); // this is used to sync (that is, wait until the GUI notices the status has changed)
		this->checkForInterrupts(po);
		return this->doStepOver;
	}

	void Interpreter::resume() {
		assert(std::this_thread::get_id()==IDE::Console::guiThreadId);
		assert(this->status==IS_PAUSED);
		this->mustResume = true;
		this->condVar.notify_one();
	}

	void RightValue::initTreeWidget(QTreeWidgetItem *item) const {
		ostringstream os;
		os << this;
		item->setText(COL_VALUE, QString::fromStdString(os.str()));
	}

	void VoidValue::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, "Uninitialized");
	}

	void BuiltInFunction::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, "A built-in fn-proc");
	}

	void UserDefinedFunction::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, "A user-defined fn-proc");
	}

	void TYPE::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, "A type");
	}

	void OSTREAM::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, "An output-stream");
	}

	void ISTREAM::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, "An input-stream");
	}

	void LIST::initTreeWidget(QTreeWidgetItem *item) const {
		const ContainerType::size_type size = this->size();
		item->setText(COL_VALUE, QString::fromStdString("A "+lexical_cast<string>(size)+"-element list"));
		for(ContainerType::size_type a=0; a<size; ++a) {
			QTreeWidgetItem * const child = new QTreeWidgetItem();
			child->setText(COL_NAME, QString::fromStdString("["+lexical_cast<string>(a+1)+"]"));
			item->addChild(child);
			container[a]->initTreeWidget(child);
		}
	}

	void TaggedValue::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, QString::fromStdString("A tagged value"));
		QTreeWidgetItem * const child1 = new QTreeWidgetItem();
		child1->setText(COL_NAME, QString::fromStdString("Tag"));
		child1->setText(COL_VALUE, QString::fromStdString(this->tag));
		item->addChild(child1);
		QTreeWidgetItem * const child2 = new QTreeWidgetItem();
		child2->setText(COL_NAME, QString::fromStdString("Value"));
		item->addChild(child2);
		this->value->initTreeWidget(child2);
	}

	void RECORD::initTreeWidget(QTreeWidgetItem *item) const {
		const MapType::size_type size = this->numberOfFields();
		item->setText(COL_VALUE, QString::fromStdString("A "+lexical_cast<string>(size)+"-field record"));
		for(MapType::const_iterator pos = this->fields.begin(); pos!=this->fields.end(); ++pos) {
			QTreeWidgetItem * const child = new QTreeWidgetItem();
			child->setText(COL_NAME, QString::fromStdString("."+pos->first));
			item->addChild(child);
			pos->second->initTreeWidget(child);
		}
	}

	void IntMapValue::initTreeWidget(QTreeWidgetItem *item) const {
		item->setText(COL_VALUE, QString::fromStdString("A map"));
		for(MapType::const_iterator pos = this->map.begin(); pos!=this->map.end(); ++pos) {
			QTreeWidgetItem * const child = new QTreeWidgetItem();
			child->setText(COL_NAME, QString::fromStdString("["+lexical_cast<string>(pos->first)+"]"));
			item->addChild(child);
			pos->second->initTreeWidget(child);
		}
	}
} // namespace InterpreterNS

namespace IDE {

std::thread::id Console::guiThreadId;

void CocoaHighlighter::highlightBlock(const QString &text) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	if (!this->alwaysEnabled) {
		QTextBlockUserData *data = this->currentBlockUserData();
		if (data) {
			assert(dynamic_cast<HighlightingInfo *>(data));
			if (!static_cast<HighlightingInfo *>(data)->doHighlight)
				return;
		} else {
			this->setCurrentBlockUserData(new HighlightingInfo(this->enabled));
			if (!this->enabled)
				return;
		}
	}
	this->unclosedComment = this->unclosedStringLiteral = false;
	const int textLenght = text.length();
	setCurrentBlockState(0);
	int index=0;
	const int previous = previousBlockState();
	if (previous==IN_SINGLEQUOTE_STRING) {
		index = -1;
		goto findSingleQuoteClosing;
	}
	if (previous==IN_DOUBLEQUOTE_STRING) {
		index = -1;
		goto findDoubleQuoteClosing;
	}
	if (previous==IN_TRIPLEQUOTE_STRING) {
		index = -3;
		goto findTripleQuoteClosing;
	}
	if (previous==IN_COMMENT) {
		index = -2;
		goto findCommentClosing;
	}
	while(index<textLenght) {
		if (text.mid(index,2)==tr("--") || text.mid(index,2)==tr("//")) {
			int endPos = text.indexOf(this->singlelineCommentEndRE, index+2);
			int len;
			if (endPos==-1) {
				if (index<=0)
					index = 0;
				len = textLenght-index;
			} else {
				len = this->singlelineCommentEndRE.matchedLength();
				if (index>=0)
					len += 2;
				else
					index = 0;
			}
			setFormat(index, len, this->commentFormat);
			index += len;
			continue;
		}
		if (text.mid(index,2)==tr("/*")) {
findCommentClosing:
			setCurrentBlockState(IN_COMMENT);
			int endPos = text.indexOf(this->commentEndRE, index+2);
			int len;
			if (endPos==-1) {
				this->unclosedComment = true;
				if (index<=0)
					index = 0;
				len = textLenght-index;
			} else {
				len = this->commentEndRE.matchedLength();
				if (index>=0)
					len += 2;
				else
					index = 0;
				setCurrentBlockState(NONE);
			}
			setFormat(index, len, this->commentFormat);
			index += len;
			continue;
		}
		if (text.mid(index, 3)==tr("\"\"\"")) {
findTripleQuoteClosing:
			setCurrentBlockState(IN_TRIPLEQUOTE_STRING);
			int endPos = text.indexOf(this->tripleQuoteliteralStringEndRE, index+3);
			int len;
			if (endPos==-1) {
				this->unclosedStringLiteral = true;
				if (index<=0)
					index = 0;
				len = textLenght-index;
			} else {
				len = this->tripleQuoteliteralStringEndRE.matchedLength();
				if (index>=0)
					len += 3;
				else
					index = 0;
				setCurrentBlockState(NONE);
			}
			setFormat(index, len, this->literalStringFormat);
			index += len;
			continue;
		}
		if (text.mid(index, 1)==tr("\'")) {
findSingleQuoteClosing:
			setCurrentBlockState(IN_SINGLEQUOTE_STRING);
			int endPos = text.indexOf(this->singleQuoteliteralStringEndRE, index+1);
			int len;
			if (endPos==-1) {
				this->unclosedStringLiteral = true;
				if (index<=0)
					index = 0;
				len = textLenght-index;
			} else {
				len = this->singleQuoteliteralStringEndRE.matchedLength();
				if (index>=0)
					++len;
				else
					index = 0;
				setCurrentBlockState(NONE);
			}
			setFormat(index, len, this->literalStringFormat);
			index += len;
			continue;
		}
		if (text.mid(index, 1)==tr("\"")) {
findDoubleQuoteClosing:
			setCurrentBlockState(IN_DOUBLEQUOTE_STRING);
			int endPos = text.indexOf(this->doubleQuoteliteralStringEndRE, index+1);
			int len;
			if (endPos==-1) {
				this->unclosedStringLiteral = true;
				if (index<=0)
					index = 0;
				len = textLenght-index;
			} else {
				len = this->doubleQuoteliteralStringEndRE.matchedLength();
				if (index>=0)
					++len;
				else
					index = 0;
				setCurrentBlockState(NONE);
			}
			setFormat(index, len, this->literalStringFormat);
			index += len;
			continue;
		}
		const char firstChar = text.mid(index, 1).toStdString()[0];
		if (IsBlank(firstChar) && text.indexOf(this->blanksRE, index)==index) {
			index += this->blanksRE.matchedLength();
			continue;
		}
		if (firstChar=='_' || firstChar=='$')
			goto try_identifier;
		if (firstChar=='?')
			goto try_keyword;
		if (isalpha(firstChar)) {
			if (text.indexOf(this->typeRE, index)==index) {
				const int len = this->typeRE.matchedLength();
				setFormat(index, len, this->typeFormat);
				index += len;
				continue;
			}
			if (text.indexOf(this->constantRE, index)==index) {
				const int len = this->constantRE.matchedLength();
				setFormat(index, len, this->constantFormat);
				index += len;
				continue;
			}
try_keyword: if (text.indexOf(this->keywordRE, index)==index) { // note: RECORD is not a keyword (while all other casings of "record", are)
				const int len = this->keywordRE.matchedLength();
				setFormat(index, len, this->keywordFormat);
				index += len;
				continue;
			}
try_identifier: if (text.indexOf(this->identifierRE, index)==index) {
				const int len = this->identifierRE.matchedLength();
				setFormat(index, len, this->identifierFormat);
				index += len;
				continue;
			}
		}
		if (isdigit(firstChar) && text.indexOf(this->numericLiteralRE, index)==index) {
			const int len = this->numericLiteralRE.matchedLength();
			setFormat(index, len, this->numericLiteralFormat);
			index += len;
			continue;
		}
		if (text.indexOf(this->operatorRE, index)==index) {
			int len = this->operatorRE.matchedLength();
			if (text.indexOf(this->parenthesisRE, index)==index) {
				len = this->parenthesisRE.matchedLength();
				setFormat(index, len, this->operatorFormat);
			}
			index += len;
			continue;
		}
		if (this->highlightUnknownChars)
			setFormat(index, 1, this->unknownFormat);
		++index;
	}
}

CocoaHighlighter::CocoaHighlighter(QTextDocument *parent, bool highlightUnknownChars, bool alwaysEnabled) :
	QSyntaxHighlighter(parent),
	highlightUnknownChars(highlightUnknownChars),
	blanksRE("\\s+", Qt::CaseInsensitive, QRegExp::RegExp2),
	keywordRE(
			"[\\?]|\\b(?:Alias|And|Block|Break|Catch|Ciao|Clear|Continue|Define|"
			"DegLex|DegRevLex|Delete|Describe|Destroy|Do|Elif|Elim|"
			"Else|End|EndAlias|EndBlock|EndCatch|EndDefine|EndForeach|"
			"EndFor|EndFunc|EndIf|EndPackage|EndRepeat|EndTry|"
			"EndUsing|EndWhile|Exit|Export|False|Foreach|For|Func|If|"
			"ImportByRef|ImportByValue|In|IsDefined|IsIn|Lex|On|Opt|"
			"Ord|Or|Package|PosTo|PrintLn|Print|Protect|Quit|Record|"
			"Ref|Repeat|Return|Set|Skip|Source|Step|Then|Time|To|"
			"TopLevel|ToPos|True|Try|Unset|UponError|Until|Unprotect|"
			"Use|Using|Var|Weights|While|Xel)\\b", Qt::CaseInsensitive, QRegExp::RegExp2),
	constantRE("\\b(?:True|False|Lex|Xel|DegLex|DegRevLex|ToPos|PosTo)\\b", Qt::CaseInsensitive, QRegExp::RegExp2),
	singleQuoteliteralStringEndRE("(?:[^\"\\\\]|(?:\\\\(?:[nrta'\"\\\\]|(?:x[0-9a-fA-F]{2}))))*\'", Qt::CaseSensitive, QRegExp::RegExp2),
	doubleQuoteliteralStringEndRE("(?:[^\"\\\\]|(?:\\\\(?:[nrta'\"\\\\]|(?:x[0-9a-fA-F]{2}))))*\"", Qt::CaseSensitive, QRegExp::RegExp2),
	tripleQuoteliteralStringEndRE("(?:[^\"\\\\]|(?:\\\\(?:[nrta'\"\\\\]|(?:x[0-9a-fA-F]{2}))))*\"\"\"", Qt::CaseSensitive, QRegExp::RegExp2),
	identifierRE("[a-z$_][a-z0-9_]*", Qt::CaseInsensitive, QRegExp::RegExp2),
	commentEndRE("(?:[^\\*]|(?:\\*[^/]))*\\*/", Qt::CaseInsensitive, QRegExp::RegExp2),
	singlelineCommentEndRE("[^\\n]*", Qt::CaseInsensitive, QRegExp::RegExp2),
	numericLiteralRE("\\d+(?:\\.\\d+)?", Qt::CaseInsensitive, QRegExp::RegExp2),
	typeRE("\\b(?:BOOL|FUNCTION|LIST|INT|RAT|RECORD|TYPE|STRING|VOID|ERROR|OSTREAM"
			"RINGELEM|RATFUN|MODULEELEM|IDEAL|MODULE|MAT|RING|PACKAGE|RINGHOM|INTMAP)\\b", Qt::CaseSensitive, QRegExp::RegExp2),
	parenthesisRE("[()\\[\\]]|\\$\\{|\\}\\$", Qt::CaseInsensitive, QRegExp::RegExp2), // parenthesisRE is a quick-hack because QCodeEdit cannot highlight operators correctly
	operatorRE("[\\+\\-\\*/:<>=()\\[\\]|%^;,]|<=|>=|<>|:{1,2}=|<<|><|::|\\$\\{|\\}\\$|\\.{1,3}", Qt::CaseInsensitive, QRegExp::RegExp2),
	unclosedComment(false),
	unclosedStringLiteral(false),
	enabled(true),
	alwaysEnabled(alwaysEnabled)
{}

void CocoaHighlighter::setFormats(const ColorScheme &scheme) {
	this->literalStringFormat.setForeground(scheme.stringLiterals.foreground);
	if (scheme.stringLiterals.background.isValid())
		this->literalStringFormat.setBackground(scheme.stringLiterals.background);
	this->constantFormat.setForeground(scheme.constants.foreground);
	if (scheme.constants.background.isValid())
		this->constantFormat.setBackground(scheme.constants.background);
	this->commentFormat.setForeground(scheme.comments.foreground);
	if (scheme.comments.background.isValid())
		this->commentFormat.setBackground(scheme.comments.background);
	this->keywordFormat.setForeground(scheme.keywords.foreground);
	if (scheme.keywords.background.isValid())
		this->keywordFormat.setBackground(scheme.keywords.background);
	this->typeFormat.setForeground(scheme.types.foreground);
	if (scheme.types.background.isValid())
		this->typeFormat.setBackground(scheme.types.background);
	this->numericLiteralFormat.setForeground(scheme.numericalLiterals.foreground);
	if (scheme.numericalLiterals.background.isValid())
		this->numericLiteralFormat.setBackground(scheme.numericalLiterals.background);
	this->unknownFormat.setBackground(scheme.unknownChars.foreground);
	if (scheme.unknownChars.background.isValid())
		this->unknownFormat.setForeground(scheme.unknownChars.background);
}

IdeErrorReporter::IdeErrorReporter(LexerNS::WarningSeverity warningLevel, Console *console, const boost::intrusive_ptr<IdeOutputStream> outputStream) :
	ErrorReporter(warningLevel, outputStream),
	console(console)
{
	this->errorFormat.setFontWeight(QFont::Bold);
	this->errorFormat.setForeground(Qt::red);
	this->warningFormat.setFontWeight(QFont::Bold);
	this->warningFormat.setForeground(Qt::blue);
	this->calledByFormat.setFontWeight(QFont::Bold);
	this->whereFormat.setFontWeight(QFont::Bold);
	this->contextFormat.setFontWeight(QFont::Bold);
	this->boldFormat.setFontWeight(QFont::Bold);
	this->highlightErrorFormat.setUnderlineColor(Qt::red);
	this->highlightErrorFormat.setFontUnderline(true);
	this->highlightErrorFormat.setUnderlineStyle(QTextCharFormat::WaveUnderline);
}

void Console::closeEvent(QCloseEvent *event) {
	event->ignore();
	this->showMinimized();
}

void Debugger::closeEvent(QCloseEvent *event) {
	event->ignore();
	this->showMinimized();
}

IdeOutputStream::IdeOutputStream(Console *console) :
	OSTREAM(false, false),
	console(console)
{
	boldFormat.setFontWeight(QFont::Bold);
}

void FlushQC::execute(Console *) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	std::unique_lock<std::mutex> lock(this->mut);
	this->flushed = true;
	this->condVar.notify_one();
}

intrusive_ptr<FlushQC> FlushQC::theInstance(new FlushQC);

void IdeOutputStream::flush() {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	intrusive_ptr<FlushQC> flushQC = FlushQC::theInstance;
	std::unique_lock<std::mutex> lock(flushQC->mut);
	flushQC->flushed = false;
	this->console->postCommand(flushQC);
	do
		flushQC->condVar.wait(lock);
	while (!flushQC->flushed);
}

void IdeOutputStream::print(const string &s, bool highlight, QTextCharFormat format) {
	this->console->postCommand(new PrintQC(s, highlight, format));
}

intrusive_ptr<OSTREAM> IdeOutputStream::print(const string &s) {
	this->console->postCommand(new PrintQC(s, false, normalFormat));
	return this;
}

intrusive_ptr<OSTREAM> IdeOutputStream::print(intrusive_ptr<const RightValue> v) {
	ostringstream ss;
	ss << v;
	this->console->postCommand(new PrintQC(ss.str(), false, normalFormat));
	return this;
}

intrusive_ptr<PrintQC> PrintQC::newline(new PrintQC("\n", false, IdeOutputStream::normalFormat));
intrusive_ptr<PrintQC> PrintQC::newlineHL(new PrintQC("\n", true, IdeOutputStream::normalFormat));

intrusive_ptr<OSTREAM> IdeOutputStream::newline() {
	this->console->postCommand(PrintQC::newline);
	return this;
}

intrusive_ptr<OSTREAM> IdeOutputStream::newlineHL() {
	this->console->postCommand(PrintQC::newlineHL);
	return this;
}

void ClearReportedLocationsQC::execute(Console *console) {
	console->clearReportedLocations();
}

intrusive_ptr<ClearReportedLocationsQC> ClearReportedLocationsQC::theInstance(new ClearReportedLocationsQC);

bool IdeLineProvider::doReadNextLine(const LexerStatus &ls, const ParserNS::ParserStatus &ps, string &chars)
{
  assert(std::this_thread::get_id()==this->console->interpreterThreadId);
  this->console->outputStream->print(prompt(ls, ps), true, IdeOutputStream::boldFormat);
  std::unique_lock<std::mutex> lock(this->mut);
  for(;;)
  {
    const Interpreter::InterpreterStatus status = this->console->interpreter->getStatus();
    assert(status==Interpreter::IS_WAITING_FOR_COMMAND || status==Interpreter::IS_WAITING_FOR_COMMAND_COMPLETION);
    if (!this->enteredLines.empty())
    {
      chars = this->enteredLines.front();
      this->enteredLines.pop_front();
      // if (!IsWhiteSpace(chars)) ...
      if (chars.find_first_not_of(" \t") != string::npos)
        this->console->interpreter->UpdateStatusToWFCC();
      break;
    }
    if (status==Interpreter::IS_WAITING_FOR_COMMAND && !this->enteredFullCommands.empty())
    {
      chars = this->enteredFullCommands.front();
      this->enteredFullCommands.pop_front();
      this->console->postCommand(ClearReportedLocationsQC::theInstance);
      break;
    }
    this->condVar.wait(lock);
  }
  this->console->outputStream->print(chars, true, IdeOutputStream::normalFormat);
  this->console->outputStream->newline();
  chars += '\n';
  return true;
}

void IdeLineProvider::enterLine(const string &line, bool isFullCommand) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	std::unique_lock<std::mutex> lock(mut);
	(isFullCommand ? this->enteredFullCommands : this->enteredLines).push_back(line);
	this->condVar.notify_one();
}

void IdeErrorReporter::implReportWarning(const string &msg) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	this->printWarning();
	this->outputStream->print(msg)->newline();
}

void IdeErrorReporter::outputHighlightedChars(const CharPointer &from, const CharPointer &to) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	IdeOutputStream * const ios = static_pointer_cast<IdeOutputStream>(this->outputStream).get();
	if (!*from) {
		ios->print("\n// <End of file>", true, IdeOutputStream::normalFormat);
		ios->newline();
		return;
	}
	intrusive_ptr<const Line> fromLine = from.getLine();
	intrusive_ptr<const Line> toLine = to.getLine();
	if (fromLine==toLine)
		return outputHighlightedLine(fromLine, from.getPositionInLine(), to.getPositionInLine(), false);
	if (!*to) {
		ios->print("\n// ... <End of file>", true, IdeOutputStream::normalFormat);
		this->outputStream->newline();
		return;
	}
	this->outputHighlightedLine(fromLine, from.getPositionInLine(), fromLine->chars.length()-2, true); // the last char of every line is \n and we don't want to underline it
	int skippedLines=0;
	intrusive_ptr<const Line> currLine = fromLine;
	while ( (currLine = currLine->getNextLineInBuffer())!= toLine )
		++skippedLines;
	if (skippedLines) {
		static_pointer_cast<IdeOutputStream>(this->outputStream)->print("// ... (", true, IdeOutputStream::normalFormat);
		if (skippedLines==1)
			ios->print("another source line");
		else
			ios->print("other ")->print(lexical_cast<string>(skippedLines))->print(" source lines");
		ios->print(") ...", true, IdeOutputStream::normalFormat);
	}
	ios->newlineHL();
	this->outputHighlightedLine(toLine, 0, to.getPositionInLine(), false);
}

void IdeErrorReporter::outputHighlightedLine(intrusive_ptr<const Line> line, size_t from, size_t to, bool keepsHilighting) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	constexpr size_t MAX_PREFIX = 40; // these are from the tty version, should we change them?!?
	constexpr size_t MAX_SUFFIX = 30;
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	IdeOutputStream * const ios = static_pointer_cast<IdeOutputStream>(this->outputStream).get();
	string prefix, suffix("\n");
	const string &chars = line->chars;
    size_t printingStart=0, printingEnd=chars.length()==0 ? 0 : chars.length()-1;
    if (from==to) {
    	if (from>0)
    		--from;
    	if (to<printingEnd)
    		++to;
    }
	if (from>=MAX_PREFIX) {
		printingStart = from-MAX_PREFIX;
		prefix = "... ";
	}
	if ( (printingEnd-to)>=MAX_SUFFIX ) {
		printingEnd = to+MAX_SUFFIX;
		suffix = " ...\n";
	}
	assert(printingStart<=printingEnd);
	assert(printingEnd<chars.length());
	ios->print(prefix, keepsHilighting, IdeOutputStream::normalFormat);
	for(size_t i=printingStart; i<=printingEnd;++i) {
		char c = chars[i];
		if (c=='\n')
			break;
		ios->print(lexical_cast<string>(c=='\t' ? ' ':c), true, (from<=i && i<=to) ? this->highlightErrorFormat : IdeOutputStream::normalFormat);
	}
	ios->print(suffix, keepsHilighting, IdeOutputStream::normalFormat);
}

void IdeErrorReporter::implReportWarning(const string &msg, const CharPointer &from, const CharPointer &to) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	this->printWarning();
	this->outputStream->print(msg);
	this->reportLineNumberWhenMeaningful(from, to, true, true);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print("\n", true, IdeOutputStream::normalFormat);
	this->outputHighlightedChars(from, to);
}

void IdeErrorReporter::implReportInterrupt(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	this->outputStream->print("*** CoCoA-5 Interrupted ***")->newline();
	this->reportLineNumberWhenMeaningful(from, to, true, true);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print("\n", true, IdeOutputStream::normalFormat);
	this->outputHighlightedChars(from, to);
}

void IdeErrorReporter::implReportError(const string &msg) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	this->printError();
	this->outputStream->print(msg)->newline();
}

void IdeErrorReporter::implReportError(const string &msg, const CharPointer &from, const CharPointer &to) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	this->printError();
	this->outputStream->print(msg);
	this->reportLineNumberWhenMeaningful(from, to, true, true);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print("\n", true, IdeOutputStream::normalFormat);
	this->outputHighlightedChars(from, to);
}

void Console::print(const string &s, bool highlight, QTextCharFormat format) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->outputHL->enabled = highlight;
	QTextCursor cursor(this->outputTextEdit->document());
	cursor.movePosition(QTextCursor::End, QTextCursor::MoveAnchor);
	cursor.setCharFormat(format);
	this->outputTextEdit->setTextCursor(cursor);
	cursor.insertText(QString::fromStdString(s));
	this->outputTextEdit->ensureCursorVisible();
}

void Console::onEnterCommandClicked()
{
  assert(std::this_thread::get_id()==Console::guiThreadId);
  const string input(this->inputTextEdit->toPlainText().toStdString());
  if (input.length()==0)
    return;
  if (this->inputHL->thereIsUnclosedComment())
  {
    QMessageBox::critical(this, tr("Unclosed comment"), tr("The command cannot be entered because it contains an unclosed comment"), QMessageBox::Ok, QMessageBox::Ok);
    return;
  }
  if (this->inputHL->thereIsUnclosedStringLiteral())
  {
    QMessageBox::critical(this, tr("Unclosed string literal"), tr("The command cannot be entered because it contains an unclosed string literal"), QMessageBox::Ok, QMessageBox::Ok);
    return;
  }
  this->history.push_back(input);
  if (this->history.size()>static_cast<deque<string>::size_type>(MAX_HISTORY))
    this->history.pop_front();
  this->currentHistoryPosition = static_cast<int>(this->history.size());
  this->clearReportedLocations();
  this->inputTextEdit->clear();
  this->inputTextEdit->setFocus();
  switch(this->interpreter->getStatus())
  {
  case Interpreter::IS_WAITING_FOR_COMMAND:
  case Interpreter::IS_WAITING_FOR_COMMAND_COMPLETION:
    break;
  case Interpreter::IS_RUNNING:
  case Interpreter::IS_RUNNING_BUILTIN:
  case Interpreter::IS_PAUSED:
    QMessageBox::warning(this, "Deferred executions", "The interpreter cannot run this command right away; your request has been queued and will be executed ASAP", QMessageBox::Ok);
    break;
  case Interpreter::IS_ENDED:
    QMessageBox::warning(this, "Interpreter is not running anymore", "The interpreter is not running anymore, so it cannot execute anything", QMessageBox::Ok);
    break;
  }
  if (this->debuggerCheckbox->isChecked())
    this->interpreter->singleStepExecution = true;

  // Input may comprise several lines; pass them separately to lineProvider
  string::size_type BOL = 0;
  while (BOL < input.size())
  {
    const string::size_type EOL = input.find_first_of('\n', BOL);
    this->lineProvider->enterLine(input.substr(BOL, EOL-BOL), false);
    if (EOL == string::npos) break;
    BOL = EOL+1;
  }
}


int Console::reportLocation(const std::string &filename, int lineNumber, int columnNumber) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	for (const ReportedLocation &rl: this->reportedLocations)
		if (rl.filename==filename && rl.lineNumber==lineNumber && rl.columnNumber==columnNumber)
			return -1;
	this->reportedLocations.push_back(ReportedLocation(filename, lineNumber, columnNumber));
	return static_cast<int>(this->reportedLocations.size())-1;
}

bool IdeErrorReporter::reportLineNumberWhenMeaningful(const CharPointer &fromPos, const CharPointer & /*toPos*/, bool printColumns, bool includeHeader) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	bool result = this->ErrorReporter::reportLineNumberWhenMeaningful(fromPos, fromPos, printColumns, includeHeader);
	if (result) {
		const intrusive_ptr<const Line> line = fromPos.getLine();
		const intrusive_ptr<const FileLineProvider> p = intrusive_ptr_cast<const FileLineProvider>(line->provider);
		this->console->postCommand(new ReportErrorQC(p->myFileName(), line->number, fromPos.getPositionInLine()));
	}
	return result;
}

void ReportErrorQC::execute(Console *console) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	const int index = console->reportLocation(this->filename, this->lineNumber, this->columnNumber);
	if (index>=0) {
		console->locationComboBox->setEnabled(true);
		assert(this->lineNumber>=1);
		const string strippedFilename(QFileInfo(QString::fromStdString(this->filename)).baseName().toStdString());
		console->locationComboBox->addItem(QString::fromStdString("Line "+lexical_cast<string>(this->lineNumber)+" (col. "+ lexical_cast<string>(this->columnNumber+1) +") of "+strippedFilename), QVariant(index));
		console->openInEditorButton->setEnabled(true);
	}
}

void Console::clearReportedLocations() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->reportedLocations.clear();
	this->openInEditorButton->setEnabled(false);
	this->locationComboBox->clear();
	this->locationComboBox->setEnabled(false);
}

void Console::postCommand(intrusive_ptr<QueuedCommand> command) {
	int nCommands;
	{
		std::lock_guard<std::mutex> lock(this->mut);
		this->commands.push_back(command);
		nCommands = this->commands.size();
	}
	constexpr int tooManyQueuedCommands = 16;  // !!arbitrary!!
	if (nCommands>tooManyQueuedCommands && std::this_thread::get_id()==this->interpreterThreadId)
		std::this_thread::sleep_for(std::chrono::milliseconds(25));
}

bool Console::eventFilter(QObject *target, QEvent *event) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	if (target == this->inputTextEdit && event->type() == QEvent::KeyPress) {
		QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
		const int key = keyEvent->key();
		if ((keyEvent->modifiers() & Qt::ControlModifier)==Qt::ControlModifier) {
			if (key==Qt::Key_Up) {
				if (this->currentHistoryPosition>0) {
					this->inputTextEdit->setText(QString::fromStdString(this->history[--this->currentHistoryPosition]));
					goto cursorAtTheEnd;
				}
				return true;
			}
			if (key==Qt::Key_Down) {
				if (this->currentHistoryPosition < static_cast<int>((this->history.size()-1))) {
					this->inputTextEdit->setText(QString::fromStdString(this->history[++this->currentHistoryPosition]));
cursorAtTheEnd:
					QTextCursor cursor(this->inputTextEdit->document());
					cursor.movePosition(QTextCursor::End, QTextCursor::MoveAnchor);
					this->inputTextEdit->setTextCursor(cursor);
				}
				return true;
			}
		} else if (key == Qt::Key_Return) {
			this->onEnterCommandClicked();
			return true;
		}
	}
	return QWidget::eventFilter(target, event);
}

void Console::clearOutputWindow(bool AskUser) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	if (!AskUser || QMessageBox::question(this, tr("Deletion confirmation"), tr("Are you sure you want to delete the contents of the output window?"),
                                              QMessageBox::Yes|QMessageBox::No, QMessageBox::No)==QMessageBox::Yes) {
		this->outputTextEdit->clear();
	}
}

void PrintQC::execute(Console *console) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	console->print(this->s, this->highlight, this->format);
}

namespace {
	class RunTheInterpreter {
		intrusive_ptr<Interpreter> interpreter;
		Console *console;
	public:
		explicit RunTheInterpreter(Console *console, intrusive_ptr<Interpreter> interpreter) :
			interpreter(interpreter),
			console(console)
		{}
		void operator()() {
			this->console->interpreterThreadId = std::this_thread::get_id();
			this->interpreter->run(console);
		}
	};
}

Console::Console(QWidget *parent, MainWindow *mainWindow, WarningSeverity warningLevel, bool warnAboutCocoa5, const vector<string> & packageList, bool fullCoCoALibError) :
		QWidget(parent),
		currentHistoryPosition(-1),
		mainWindow(mainWindow),
		packageList(packageList),
		outputStream(new IdeOutputStream(this)),
		lineProvider(new IdeLineProvider(this)),
		errorReporter(new IdeErrorReporter(warningLevel, this, this->outputStream))
{
	guiThreadId = std::this_thread::get_id();
	this->setupUi(this);
	this->inputHL = new CocoaHighlighter(this->inputTextEdit->document(), true, true);
	this->outputHL = new CocoaHighlighter(this->outputTextEdit->document(), false, false);
	this->setHighlighterFormats();
	this->outputHL->enabled = false;
	connect(this->enterButton, SIGNAL(clicked()), this, SLOT(onEnterCommandClicked()));
	connect(this->clearOutputButton, SIGNAL(clicked()), this, SLOT(clearOutputWindow()));
	this->inputTextEdit->installEventFilter(this);
	this->interpreter = new Interpreter(warnAboutCocoa5, this->lineProvider, this->errorReporter, this->outputStream, fullCoCoALibError, 5000); // 5000 is "arbitrary" MaxStackSize
	std::thread intThread(RunTheInterpreter(this, this->interpreter));
        intThread.detach();
	QTimer *timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this, SLOT(processCommands()));
	timer->start(20);
	connect(this->openInEditorButton, SIGNAL(clicked()), this, SLOT(onOpenInEditorClicked()));
	connect(this->pauseButton, SIGNAL(clicked()), this, SLOT(onPauseClicked()));
	connect(this->interruptButton, SIGNAL(clicked()), this, SLOT(onInterruptClicked()));
	this->packageLoadingProgressBar->setMinimum(0);
	this->packageLoadingProgressBar->setMaximum(packageList.size()+1);
	this->packageLoadingProgressBar->setValue(0);
}

void Console::onOpenInEditorClicked() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	int currentIndex = this->locationComboBox->currentIndex();
	if (currentIndex<0)
		return;
	QVariant variant = this->locationComboBox->itemData(currentIndex);
	assert(variant.type()==QVariant::Int);
	int i = variant.toInt();
	assert(0<=i && i<static_cast<int>(this->reportedLocations.size()));
	QFileInfo file(QString::fromStdString(this->reportedLocations[i].filename));
	if (!file.exists()) {
		QMessageBox::critical(this, tr("File not found"), tr("Cannot find the selected file"), QMessageBox::Ok);
		return;
	}
	QString absPath(file.absoluteFilePath());
	SourceEditor *sourceEditor = this->mainWindow->editorFor(absPath, true);
	if (!sourceEditor) {
		sourceEditor = this->mainWindow->onFileNew();
		sourceEditor->load(absPath);
	}
	QEditor *const editor = sourceEditor->getEditor();
	const int lineNumber = this->reportedLocations[i].lineNumber;
	const int columnNumber = this->reportedLocations[i].columnNumber;
	assert(lineNumber>=1);
	editor->setCursor(QDocumentCursor(editor->document(), lineNumber - 1, columnNumber));
}

SourceEditor *MainWindow::editorFor(QString filename, bool activateWindow) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QString absPath = QFileInfo(filename).absoluteFilePath();
	QList<QMdiSubWindow *> wList = this->mdiArea->subWindowList();
	for (int i = 0; i < wList.size(); ++i) {
		QMdiSubWindow *const win = wList.at(i);
		SourceEditor *ed = dynamic_cast<SourceEditor *>(win->widget());
		if (ed && QFileInfo(ed->getEditor()->fileName()).absoluteFilePath()==absPath) {
			if (activateWindow) {
				win->raise();
				ed->getEditor()->setFocus();
			}
			return ed;
		}
	}
	return 0;
}

void MainWindow::closeEvent(QCloseEvent *e) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	if (this->interpreter->getStatus()!=Interpreter::IS_ENDED &&
			QMessageBox::question(
				this,
				tr("Exit confirmation"),
				tr("Are you sure you want to quit C5?"),
				QMessageBox::Yes|QMessageBox::No,
				QMessageBox::No)==QMessageBox::No) {
		e->ignore();
		return;
	}
        this->console->clearOutputWindow(false); // false -> do not ask for confirmation
	this->statusLabel->setText(" Exiting... "); // next cmd makes the update appear
        QMessageBox::information(this, "Ciao", "Ciao :-)");
	this->mdiArea->closeAllSubWindows();
	if (mdiArea->subWindowList().size()!=2 /* that is, the console and the debugger */)
		e->ignore();
	else {
		if (this->interpreter->getStatus()==Interpreter::IS_ENDED)
			QMessageBox::information(this, "Ciao", "Ciao :-)");
		e->accept();
	}
}

MainWindow *MainWindow::findMainWindow(QObject *o) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	while (o) {
		if (MainWindow *mw = dynamic_cast<MainWindow *>(o))
			return mw;
		o = o->parent();
		assert(o);
	}
	return 0; // to make the compiler happy
}

Console *MainWindow::findConsole(QObject *w) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	return findMainWindow(w)->console;
}

void Console::processCommands() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->mainWindow->updateStatusLabel(); // it's important to do this before executing the commands (resume() uses flush to sync)
	std::lock_guard<std::mutex> lock(this->mut);
	while (!this->commands.empty()) {
		this->commands.front()->execute(this);
		this->commands.pop_front();
	}
}

QTextCharFormat IdeOutputStream::normalFormat;
QTextCharFormat IdeOutputStream::boldFormat;

void IdeErrorReporter::printCalledBy() {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print(CalledbyPrefix, false, this->calledByFormat);
}

void IdeErrorReporter::printWhere() {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print(WherePrefix, false, this->whereFormat);
}

void IdeErrorReporter::printBold(const string &s) {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print(s, false, this->boldFormat);
}

void IdeErrorReporter::printWarning() {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print(WarningPrefix, false, this->warningFormat);
}

void IdeErrorReporter::printError() {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print(ErrorPrefix, false, this->errorFormat);
}

void IdeErrorReporter::printContext() {
	assert(std::this_thread::get_id()==this->console->interpreterThreadId);
	assert(dynamic_pointer_cast<IdeOutputStream>(this->outputStream));
	static_pointer_cast<IdeOutputStream>(this->outputStream)->print(ContextPrefix, false, this->contextFormat);
}

namespace {

	class PageUpDownCommand : public QEditorInputBinding::Command {
		const bool up;
		const QDocumentCursor::MoveMode mode;
		public:
		PageUpDownCommand(bool up, QDocumentCursor::MoveMode mode) :
			up(up),
			mode(mode)
		{}
		void exec(QEditor *e) {
			QDocumentCursor::MoveMode mode = this->mode;
			if ( e->flag(QEditor::LineWrap) && e->flag(QEditor::CursorJumpPastWrap) )
					mode |= QDocumentCursor::ThroughWrap;
			if (up)
				e->pageUp(mode);
			else
				e->pageDown(mode);
		}
	};

	class EmacsBinding : public QEditorInputBinding {
		QString id() const { return "emacs binding"; }
		QString name() const { return this->id(); }
	public:
		EmacsBinding() {
			this->setMapping(QKeySequence("Ctrl+A"), new QEditorInputBinding::MotionCommand(QDocumentCursor::StartOfLine, QDocumentCursor::MoveAnchor, 1));
			this->setMapping(QKeySequence("Ctrl+E"), new QEditorInputBinding::MotionCommand(QDocumentCursor::EndOfLine, QDocumentCursor::MoveAnchor, 1));
			this->setMapping(QKeySequence("Ctrl+V"), new PageUpDownCommand(false, QDocumentCursor::MoveAnchor));
			this->setMapping(QKeySequence("Ctrl+Shift+V"), new PageUpDownCommand(false, QDocumentCursor::KeepAnchor));
			this->setMapping(QKeySequence("Meta+V"), new PageUpDownCommand(true, QDocumentCursor::MoveAnchor));
			this->setMapping(QKeySequence("Meta+Shift+V"), new PageUpDownCommand(true, QDocumentCursor::KeepAnchor));
			// to add:
			// Meta+W copy
			// Ctrl+W cut
			// Ctrl+Y paste
			// Ctrl+_ undo
			// Ctrl+S search
			// Ctrl+X, Ctrl+S save
		}
	};

}

void ColorScheme::applyTo(QFormatScheme *fScheme) const {
	fScheme->setFormat("stringLiteral", this->stringLiterals);
	fScheme->setFormat("constant", this->constants);
	fScheme->setFormat("comment", this->comments);
	fScheme->setFormat("keyword", this->keywords);
	fScheme->setFormat("type", this->types);
	fScheme->setFormat("numericLiteral", this->numericalLiterals);
	fScheme->setFormat("unknownChar", this->unknownChars);
	fScheme->setFormat("braceMatch", this->braceMatch);
	fScheme->setFormat("braceMismatch", this->braceMismatch);
}

SourceEditor::SourceEditor(MainWindow *parent) :
		QWidget(parent),
		menuAction(new QAction(this))
{
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->setupUi(this);
	const bool emacsLike = parent->actionEmacsLike->isChecked();
	this->codeEdit = new QCodeEdit(!emacsLike, this);
	QEditor * const ed = this->getEditor();
	if (emacsLike) {
		ed->addInputBinding(new EmacsBinding());
	}
	//this->codeEdit->addPanel("Line Mark Panel", QCodeEdit::West, true)->setShortcut(QKeySequence("F6"));
	//this->codeEdit->addPanel("Line Number Panel", QCodeEdit::West, true)->setShortcut(QKeySequence("F11"));
	//this->codeEdit->addPanel("Fold Panel", QCodeEdit::West, true)->setShortcut(QKeySequence("F9"));
	this->codeEdit->addPanel("Line Change Panel", QCodeEdit::West, true);
	this->codeEdit->addPanel("Status Panel", QCodeEdit::South, true);
	this->codeEdit->addPanel("Goto Line Panel", QCodeEdit::South);
	this->codeEdit->addPanel("Search Replace Panel", QCodeEdit::South);
	this->outerVL->addWidget(ed);
	QFormatScheme *fScheme = new QFormatScheme(this);
	parent->getCurrentColorScheme()->applyTo(fScheme);
	QLanguageFactory *lFactory = new QLanguageFactory(fScheme, this);
	lFactory->addDefinitionPath(":/qxs");
	lFactory->setLanguage(this->getEditor(), "a.cocoa5");
	connect(this->saveButton, SIGNAL(clicked()), this, SLOT(save()));
	connect(this->saveAsButton, SIGNAL(clicked()), this, SLOT(saveAs()));
	connect(this->saveAndRunButton, SIGNAL(clicked()), this, SLOT(saveAndRun()));
	connect(ed, SIGNAL(titleChanged(const QString&)), this, SLOT(onEditorTitleChanged(const QString&)));
	connect(ed, SIGNAL(contentModified(bool)), this, SLOT(onEditorContentModified(bool)));
	ed->setTitle("unnamed [*]");
    this->menuAction->setCheckable(true);
    connect(this->menuAction, SIGNAL(triggered()), ed, SLOT(setFocus()));
}

SourceEditor *MainWindow::onFileNew() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	SourceEditor *sourceEditor = new SourceEditor(this);
	this->menuWindow->addAction(sourceEditor->menuAction);
	this->menuWindowActionGroup->addAction(sourceEditor->menuAction);
	QMdiSubWindow * const win = this->mdiArea->addSubWindow(sourceEditor);
	win->setAttribute(Qt::WA_DeleteOnClose);
	win->showNormal();
	return sourceEditor;
}

void SourceEditor::onEditorTitleChanged(const QString& title) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->setWindowTitle(title);
	QString fn = this->getEditor()->fileName();
	if (!fn.size())
		fn = "unnamed";
	this->menuAction->setText(fn);
}

void SourceEditor::onEditorContentModified(bool y) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->setWindowModified(y);
}

void SourceEditor::load(const QString& file) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QEditor *const ed = this->getEditor();
	QString filename = file.count() ? QFileInfo(file).absoluteFilePath() : file;
	if (filename.size() && QFile::exists(filename)) {
		ed->load(filename);
		// TODO: do something for recent files...
		// updateRecentFiles(filename);
		this->getEditor()->setTitle(QString("%1 [*]").arg(filename));
	} else {
		ed->setFileName("");
		ed->setText("");
		this->getEditor()->setTitle("unnamed [*]");
	}
}

void SourceEditor::save() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QEditor *const ed = this->getEditor();
	if (ed->fileName().size())
		ed->save();
	else
		this->saveAs();
}

void SourceEditor::saveAndRun() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->save();
	QEditor *const ed = this->getEditor();
	if (ed->isContentModified())
		return;
	ostringstream ss;
	intrusive_ptr<STRING> s(new STRING(ed->fileName().toStdString()));
	s->dumpAsString(ss);
	Console *console = MainWindow::findConsole(this);
	if (console->interpreter->getStatus()!=Interpreter::IS_WAITING_FOR_COMMAND)
		QMessageBox::warning(this, "Deferred executions", "The interpreter cannot run this file right away; your request has been queued and will be executed ASAP", QMessageBox::Ok);
	console->lineProvider->enterLine(this->debuggerCheckbox->isChecked() ? "debug("+ss.str()+") ;" : "Source "+ss.str()+";", true);
}

void MainWindow::onFileSaveOutputWindow() {
	QString fn = QFileDialog::getSaveFileName(this, tr("Save output-windows contents as..."));
	if (fn.isEmpty())
		return;
	QFile file(fn);
	if (!file.open(QIODevice::WriteOnly)) {
		QMessageBox::critical(this, "Cannot write", "Cannot open the file for writing", QMessageBox::Ok);
		return;
	}
	if (file.write(this->console->outputTextEdit->toPlainText().toUtf8()) < 0)
		QMessageBox::critical(this, "Cannot write", "Error while writing to the file", QMessageBox::Ok);
	file.close();
}

void SourceEditor::saveAs() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QEditor *const ed = this->getEditor();
	QString fn;
	for(;;) {
		fn = QFileDialog::getSaveFileName(this, tr("Save file as..."), ed->fileName(), tr("Sources (*.cocoa5 *.cpkg5)"));
		if (fn.isEmpty())
			return;
		if (!(fn.endsWith(tr(".cocoa5")) || fn.endsWith(tr(".cpkg5"))))
				fn.append(".cocoa5");
		SourceEditor *otherEditor = MainWindow::findMainWindow(this)->editorFor(fn, false);
		if (otherEditor==0 || otherEditor==this)
			break;
		QMessageBox::critical(this, "Name already in use", "This filename is already used by another edit-window; please choose another name", QMessageBox::Ok);
	}
	ed->save(fn);
	this->setWindowTitle(QString("[%1[*]]").arg(fn));
}

bool SourceEditor::maybeSave() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QEditor *const ed = this->getEditor();
	if (ed->isContentModified()) {
		int ret = QMessageBox::warning(
							this,
							tr("About to close"),
							tr("The file %1 contains unsaved modifications.\nWould you like to save it?").arg(ed->fileName()),
							QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel, QMessageBox::Yes);
		if (ret == QMessageBox::Cancel)
			return true;
		else if (ret == QMessageBox::Yes)
			ed->save();
	}
	return false;
}

void SourceEditor::closeEvent(QCloseEvent *e) {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	if ( this->maybeSave() ) {
		e->ignore();
		return;
	}
	this->QWidget::closeEvent(e);
}

void Debugger::setHighlighterFormats() {
	this->codeHL->setFormats(*this->mainWindow->getCurrentColorScheme());
	this->codeHL->rehighlight();
}

Debugger::Debugger(QWidget *parent, MainWindow *mainWindow) :
		QWidget(parent),
		mainWindow(mainWindow)
{
	this->setupUi(this);
	this->codeHL = new CocoaHighlighter(this->codeTextEdit->document(), false, true);
	this->setHighlighterFormats();
	connect(this->stepIntoButton, SIGNAL(clicked()), this, SLOT(onStepIntoClicked()));
	connect(this->stepOverButton, SIGNAL(clicked()), this, SLOT(onStepOverClicked()));
	connect(this->stepOutButton, SIGNAL(clicked()), this, SLOT(onStepOutClicked()));
	connect(this->stepOutFnProcButton, SIGNAL(clicked()), this, SLOT(onStepOutFnProcClicked()));
	connect(this->continueButton, SIGNAL(clicked()), this, SLOT(onContinueClicked()));
	connect(this->stopButton, SIGNAL(clicked()), this, SLOT(onStopClicked()));
	this->format.setBackground(QColor(255, 239, 0, 255));
	connect(this->callStackList, SIGNAL(itemSelectionChanged()), this, SLOT(onCallStackSelectionChanged()));
	connect(this->hideExportedNamesCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onHidingChoicesChanged()));
	connect(this->hideFunctionsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onHidingChoicesChanged()));
	connect(this->hideSystemNamesCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onHidingChoicesChanged()));
	connect(this->hideTypesCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onHidingChoicesChanged()));
	connect(this->autoFillLocalsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onAutoFillLocalsStateChanged()));
	connect(this->fillLocalsButton, SIGNAL(clicked()), this, SLOT(onFillLocalsClicked()));
}

void Debugger::onFillLocalsClicked() {
	this->fillLocalsButton->setEnabled(false);
	this->fillLocals();
}

void Debugger::onAutoFillLocalsStateChanged() {
	if (this->autoFillLocalsCheckBox->isChecked()) {
		this->fillLocals();
		this->fillLocalsButton->setEnabled(false);
	}
	else
		this->fillLocalsButton->setEnabled(true);
}

void Debugger::onCallStackSelectionChanged() {
	const int cRow = this->callStackList->currentRow();
	const bool isTopLevel = cRow==(this->callStackList->count()-1);
	const bool interpreterPaused = this->interpreter->getStatus() == Interpreter::IS_PAUSED;
	this->stepIntoButton->setEnabled(interpreterPaused && cRow==0);
	this->stepOverButton->setEnabled(interpreterPaused && cRow==0);
	this->stepOutFnProcButton->setEnabled(interpreterPaused && cRow==0 && !isTopLevel);
	this->stepOutButton->setEnabled(interpreterPaused && cRow==0);
	this->continueButton->setEnabled(interpreterPaused);
	this->stopButton->setEnabled(interpreterPaused);
	this->hideExportedNamesCheckBox->setEnabled(isTopLevel);
	this->hideSystemNamesCheckBox->setEnabled(isTopLevel);
	this->fillCode();
	if (this->autoFillLocalsCheckBox->isChecked()) {
		this->fillLocalsButton->setEnabled(false);
		this->fillLocals();
	} else {
		this->fillLocalsButton->setEnabled(true);
		this->localsTree->clear();
	}
}

void Debugger::onHidingChoicesChanged() {
	if (this->autoFillLocalsCheckBox->isChecked())
		this->fillLocals();
	else
		this->fillLocalsButton->setEnabled(true);
}

void Debugger::addLocal(const string &name, const StaticEnv::VarData &varData, const Frame *currentFrame) {
	if (varData.depth<0)
		return;
	const QString qname(QString::fromStdString(name));
	RuntimeEnvironment &re = this->interpreter->runtimeEnvironment;
	const Frame *f = varData.debuggerTryToFindFrame(re.getTopLevelFrame(), currentFrame);
	if (!f) {
		this->addLocalError(qname, "Its environment does not exist anymore");
		return;
	}
	if (varData.isCapturedValue) {
		assert(f->userdefinedFun);
		assert(static_cast<vector<intrusive_ptr<RightValue> >::size_type>(varData.index) < f->userdefinedFun->capturedValues.size());
		const intrusive_ptr<RightValue> v(f->userdefinedFun->capturedValues[varData.index]);
		if (dynamic_pointer_cast<FUNCTION>(v) && this->hideFunctionsCheckBox->isChecked())
			return;
		QTreeWidgetItem *item = new QTreeWidgetItem(this->localsTree);
		item->setText(COL_NAME, qname);
		item->setText(COL_KIND, "captured");
		v->initTreeWidget(item);
		this->localsTree->addTopLevelItem(item);
		return;
	}
	this->addLocal(name, f, varData.index);
}

void Debugger::addLocal(const std::string &name, const InterpreterNS::Frame *f, int index) {
	QString qname(QString::fromStdString(name));
	assert(index>=0 && index<static_cast<int>(f->varSlots.size()));
	const VariableSlot &vs = f->varSlots[index];
	if (vs.isSystemProtected() && !vs.isIterationVariable() && this->hideSystemNamesCheckBox->isEnabled() && this->hideSystemNamesCheckBox->isChecked())
		return;
	if (vs.isInPackage() && this->hideExportedNamesCheckBox->isEnabled() && this->hideExportedNamesCheckBox->isChecked())
		return;
	intrusive_ptr<Value> v = vs.value;
	if (!v) { // this can only happen when a package has been reloaded (and some of the previous version members are not defined anymore) or an indeterminate has been removed because of a Use statement 
		assert(f==this->interpreter->runtimeEnvironment.getTopLevelFrame());
		return;
	}
	intrusive_ptr<RightValue> rv;
	try {
		rv = v->asRightValue();
	} catch (const InterruptException &) {
		this->addLocalError(qname, "Evaluation has been interrupted");
		return;
	} catch (const RuntimeException &) {
		this->addLocalError(qname, "Reference cannot be evaluated");
		return;
	}
	if (dynamic_pointer_cast<TYPE>(rv) && this->hideTypesCheckBox->isChecked())
		return;
	if (dynamic_pointer_cast<FUNCTION>(rv) && this->hideFunctionsCheckBox->isChecked())
		return;
	QTreeWidgetItem *item = new QTreeWidgetItem(this->localsTree);
	item->setText(COL_NAME, qname);
	QString kind;
	if (vs.isIterationVariable())
		kind = "iteration";
	else if (vs.isInPackage())
		kind = "package export";
	else if (vs.isSystemProtected())
		kind = "system";
	else {
		if (dynamic_pointer_cast<LeftValue>(v))
			kind = vs.isProtected() ? "user protected, reference" : "reference";
		else if (vs.isProtected())
			kind = "user protected";
	}
	item->setText(COL_KIND, kind);
	rv->initTreeWidget(item);
	this->localsTree->addTopLevelItem(item);
}

void Debugger::addLocalError(const QString &name, const QString &message) {
	QTreeWidgetItem *item = new QTreeWidgetItem(this->localsTree);
	item->setText(COL_NAME, name);
	item->setText(COL_VALUE, message);
	this->localsTree->addTopLevelItem(item);
}

void Debugger::update() {
	switch (this->interpreter->getStatus()) {
	case Interpreter::IS_WAITING_FOR_COMMAND:
	case Interpreter::IS_WAITING_FOR_COMMAND_COMPLETION:
	case Interpreter::IS_PAUSED:
	case Interpreter::IS_ENDED:
		this->fillCallStack();
		this->setEnabled(true);
		break;
	case Interpreter::IS_RUNNING:
	case Interpreter::IS_RUNNING_BUILTIN:
		this->setEnabled(false);
		return;
	}
}

void Debugger::fillLocals() {
	this->localsTree->clear();
	if (this->listedFrames.size()==0)
		return;
	int selected = this->callStackList->currentRow();
	assert(selected>=0 && selected<static_cast<int>(this->listedFrames.size()));
	RuntimeEnvironment &re = this->interpreter->runtimeEnvironment;
	const Frame *const tlFrame = re.getTopLevelFrame();
	const Frame *f = selected==0 ? re.getCurrentFrame() : (this->listedFrames[selected-1]-1);
	for(;f!=tlFrame;--f) {
		assert(f>tlFrame);
		intrusive_ptr<StaticEnv> env = f->block->staticEnv;
		assert(env);
		assert(!f->userdefinedFun || f->userdefinedFun->fnDecl->staticEnv==f->block->staticEnv);
		for(StaticEnv::IdMap::const_iterator it = env->identifierMap.begin(); it!=env->identifierMap.end(); ++it)
			this->addLocal(it->first, it->second, f);
		if (f->userdefinedFun)
			break;
	}
	if (f==tlFrame) {
		for(map<string, int>::const_iterator it = re.topLevelIdentifiers.begin(); it!=re.topLevelIdentifiers.end(); ++it) {
			assert(it->first.length());
			assert(it->second>=0 && it->second<static_cast<int>(tlFrame->varSlots.size()));
			if (it->first[0]!='$')
				this->addLocal(it->first, f, it->second);
		}
	}
}

void Debugger::fillCode() {
	if (this->listedFrames.size()==0 || this->interpreter->getStatus()!=Interpreter::IS_PAUSED) {
		this->codeTextEdit->clear();
		this->codeGroupBox->setTitle("Code");
		return;
	}
	int selected = this->callStackList->currentRow();
	assert(selected>=0 && selected<static_cast<int>(this->listedFrames.size()));
	CharPointer begin(this->interpreter->pausedOn->getBegin());
	CharPointer end(this->interpreter->pausedOn->getEnd());
	if (selected) {
		const Frame *f = this->listedFrames[selected-1];
		assert(f->userdefinedFun);
		assert(f->invocationExp);
		begin = f->invocationExp->getBegin();
		end = f->invocationExp->getEnd();
	}
	const QTextCharFormat oldFormat = this->codeTextEdit->currentCharFormat();
        const intrusive_ptr<const LineProvider> provider = begin.getLine()->provider;
	if (provider->IamReadingFromFile()) {
          intrusive_ptr<const FileLineProvider> fileProvider = dynamic_pointer_cast<const FileLineProvider>(begin.getLine()->provider); // NASTY HACK
          this->codeGroupBox->setTitle(QString::fromStdString(provider->myFileName()));
		this->codeTextEdit->setText(QString::fromStdString(fileProvider->wholeFile));
		QTextCursor c(this->codeTextEdit->document());
		c.setPosition(0);
		const int beginLineNumber = begin.getLine()->number;
		const int endLineNumber = end.getLine()->number;
		assert(beginLineNumber>=1);
		assert(endLineNumber>=beginLineNumber);
		c.movePosition(QTextCursor::Down, QTextCursor::MoveAnchor, beginLineNumber-1);
		c.movePosition(QTextCursor::Right, QTextCursor::MoveAnchor, static_cast<int>(begin.getPositionInLine()));
		c.movePosition(QTextCursor::StartOfLine, QTextCursor::KeepAnchor);
		c.movePosition(QTextCursor::Down, QTextCursor::KeepAnchor, endLineNumber - beginLineNumber);
		c.movePosition(QTextCursor::Right, QTextCursor::KeepAnchor, static_cast<int>(end.getPositionInLine())+1);
		c.mergeCharFormat(this->format);
		c.clearSelection();
		c.setCharFormat(oldFormat);
		this->codeTextEdit->setTextCursor(c);
	} else {
		this->codeGroupBox->setTitle("Top-level code");
		this->codeTextEdit->clear();
		intrusive_ptr<const Line> line = begin.getLine();
		this->codeTextEdit->append(QString::fromStdString(line->chars));
		QTextCursor c(this->codeTextEdit->document());
		c.setPosition(0);
		c.movePosition(QTextCursor::Right, QTextCursor::MoveAnchor, static_cast<int>(begin.getPositionInLine()));
		c.movePosition(QTextCursor::StartOfLine, QTextCursor::KeepAnchor);
		intrusive_ptr<const Line> endLine = end.getLine();
		for(;;) {
			if (line==endLine) {
				c.movePosition(QTextCursor::Right, QTextCursor::KeepAnchor, static_cast<int>(end.getPositionInLine())+1);
				break;
			}
			line = line->getNextLineInBuffer();
			this->codeTextEdit->append(QString::fromStdString(line->chars));
			c.movePosition(QTextCursor::Down, QTextCursor::KeepAnchor);
		}
		c.mergeCharFormat(this->format);
		c.clearSelection();
		c.setCharFormat(oldFormat);
		this->codeTextEdit->setTextCursor(c);
	}
	this->codeTextEdit->ensureCursorVisible();
}

void Debugger::fillCallStack() {
	this->listedFrames.clear();
	this->callStackList->clear();
	RuntimeEnvironment &re = this->interpreter->runtimeEnvironment;
	const Frame *currentFrame = re.getCurrentFrame();
	const Frame *tlFrame = re.getTopLevelFrame();
	for(;currentFrame>tlFrame;--currentFrame) {
		if (currentFrame->userdefinedFun) {
			this->listedFrames.push_back(currentFrame);
			string name(currentFrame->userdefinedFun->fnDecl->fnName);
			if (name.empty())
				name = "<anonymous function>";
			this->callStackList->addItem(QString::fromStdString("Fn-proc "+name));
		}
	}
	this->callStackList->addItem("Top-level");
	assert(!tlFrame->userdefinedFun);
	this->listedFrames.push_back(tlFrame);
	this->callStackList->setCurrentRow(0);
}

void Debugger::clearAndDisable() {
	this->setEnabled(false);
	this->callStackList->clear();
	this->localsTree->clear();
	this->codeTextEdit->clear();
}

void Debugger::onStepIntoClicked() {
	this->clearAndDisable();
	this->interpreter->resume();
}

void Debugger::onStepOverClicked() {
	this->clearAndDisable();
	this->interpreter->singleStepExecution = false;
	this->interpreter->doStepOver = true;
	this->interpreter->resume();
}

void Debugger::onStepOutFnProcClicked() {
	this->clearAndDisable();
	Frame *f = this->interpreter->runtimeEnvironment.getCurrentFunctionFrameOrNull();
	assert(f);
	if (f) {
		f->singleStepOnPop = true;
		this->interpreter->singleStepExecution = false;
	}
	this->interpreter->resume();
}

void Debugger::onStepOutClicked() {
	this->clearAndDisable();
	this->interpreter->runtimeEnvironment.getCurrentFrame()->singleStepOnPop = true;
	this->interpreter->singleStepExecution = false;
	this->interpreter->resume();
}

void Debugger::onContinueClicked() {
	this->clearAndDisable();
	this->interpreter->singleStepExecution = false;
	this->interpreter->resume();
}

void Debugger::onStopClicked()
{
  this->clearAndDisable();
  this->interpreter->singleStepExecution = false;
  CoCoA::SetSignalReceived(2/*SIGINT*/); // actual signal does not matter (I believe)
  this->interpreter->resume();
}

void MainWindow::initColorSchemes() {
	static bool done;
	if (done)
		return;
	// Anna's color scheme
	csAnna.stringLiterals.foreground = QColor::fromRgb(0, 0x8b, 0);
	csAnna.constants.foreground = QColor::fromRgb(0xcd, 0x85, 0);
	csAnna.comments.foreground = QColor::fromRgb(0xcd, 0, 0);
	csAnna.keywords.foreground = QColor::fromRgb(0x36, 0x64, 0x8b);
	csAnna.types.foreground = QColor::fromRgb(0xcd, 0x69, 0xc9);
	csAnna.numericalLiterals.foreground = QColor::fromRgb(0, 0, 0);
	csAnna.unknownChars.foreground = QColor::fromRgb(0xff, 0xff, 0xff);
	csAnna.unknownChars.background = QColor::fromRgb(0xff, 0, 0);
	csAnna.braceMatch.foreground = QColor::fromRgb(0x8b, 0, 0x8b);
	csAnna.braceMatch.background = QColor::fromRgb(0xff, 0xff, 0);
	csAnna.braceMismatch.foreground = QColor::fromRgb(0xff, 0xff, 0xff);
	csAnna.braceMismatch.background = QColor::fromRgb(0xff, 0, 0);
	// Gio's color scheme
	csGio.stringLiterals.foreground = QColor::fromRgb(0xff, 0, 0x90);
	csGio.constants.foreground = QColor::fromRgb(0, 0, 0xff);
	csGio.comments.foreground = QColor::fromRgb(0, 0x70, 0);
	csGio.keywords.foreground = QColor::fromRgb(0, 0, 0xff);
	csGio.types.foreground = QColor::fromRgb(0x64, 0x95, 0xed);
	csGio.numericalLiterals.foreground = QColor::fromRgb(0, 0, 0x40);
	csGio.unknownChars.foreground = QColor::fromRgb(0xff, 0xff, 0xff);
	csGio.unknownChars.background = QColor::fromRgb(0xff, 0, 0);
	csGio.braceMatch.foreground = QColor::fromRgb(0x8b, 0, 0x8b);
	csGio.braceMatch.background = QColor::fromRgb(0xff, 0xff, 0);
	csGio.braceMismatch.foreground = QColor::fromRgb(0xff, 0xff, 0xff);
	csGio.braceMismatch.background = QColor::fromRgb(0xff, 0, 0);
	done = true;
}

void Console::setHighlighterFormats() {
	const ColorScheme *const cs = this->mainWindow->getCurrentColorScheme();
	this->inputHL->setFormats(*cs);
	this->inputHL->rehighlight();
	this->outputHL->setFormats(*cs);
	this->outputHL->rehighlight();
}

void MainWindow::setHighlighterFormats() {
	QList<QMdiSubWindow *> wList = this->mdiArea->subWindowList();
	for (int i = 0; i < wList.size(); ++i) {
		QMdiSubWindow *const win = wList.at(i);
		if (SourceEditor *se = dynamic_cast<SourceEditor *>(win->widget())) {
			QEditor * const ed = se->getEditor();
			this->currentColorScheme->applyTo(ed->document()->formatScheme());
			ed->highlight();
		}
	}
	this->console->setHighlighterFormats();
	this->debugger->setHighlighterFormats();
}

void MainWindow::onAnnaCS() {
	this->currentColorScheme = &csAnna;
	this->setHighlighterFormats();
	QSettings settings;
	settings.setValue(SETTING_COLORSCHEME_NAME, "Anna");
}

void MainWindow::onGioCS() {
	this->currentColorScheme = &csGio;
	this->setHighlighterFormats();
	QSettings settings;
	settings.setValue(SETTING_COLORSCHEME_NAME, "Gio");
}

ColorScheme MainWindow::csAnna;
ColorScheme MainWindow::csGio;

MainWindow::MainWindow(QWidget *parent, QApplication *application, WarningSeverity warningLevel, bool warnAboutCocoa5, const vector<string> & packageList, bool fullCoCoALibError) :
		QMainWindow(parent),
		application(application),
		mdiArea(new QMdiArea(this)),
		oldInterpreterStatus(static_cast<Interpreter::InterpreterStatus>(-1)),
		statusLabel(new QLabel(this)),
		menuWindowActionGroup(new QActionGroup(this)),
		colorSchemeActionGroup(new QActionGroup(this)),
		actionMenuConsole(new QAction("&Console",this)),
		actionMenuDebugger(new QAction("&Debugger", this)),
		currentColorScheme(&csAnna)
{
	initColorSchemes();
	this->setupUi(this);
	this->setCentralWidget(this->mdiArea);
	this->actionMenuConsole->setIcon(QIcon(":/images/utilities-terminal.png"));
	this->actionMenuDebugger->setIcon(QIcon(":/images/utilities-system-monitor.png"));
	this->toolBar->addAction(this->actionMenuConsole);
	this->toolBar->addAction(this->actionMenuDebugger);
	QSettings settings;
	if (settings.value(SETTING_COLORSCHEME_NAME, "Anna")=="Gio")
		this->currentColorScheme = &csGio;
	if (settings.value(SETTING_EMACS_LIKE_KEYBINDINGS, false).toBool())
		this->actionEmacsLike->setChecked(true);
	this->menuWindowActionGroup->addAction(this->actionMenuConsole);
	this->menuWindowActionGroup->addAction(this->actionMenuDebugger);
	this->colorSchemeActionGroup->addAction(this->actionAnna);
	this->actionAnna->setChecked(true);
	connect(this->actionAnna, SIGNAL(triggered()), this, SLOT(onAnnaCS()));
	this->colorSchemeActionGroup->addAction(this->actionGio);
	connect(this->actionGio, SIGNAL(triggered()), this, SLOT(onGioCS()));
	this->statusBar()->addWidget(this->statusLabel);
	this->menuWindow->addAction(this->actionMenuConsole);
	this->menuWindow->addAction(this->actionMenuDebugger);
	this->console = new Console(this, this, warningLevel, warnAboutCocoa5, packageList, fullCoCoALibError);
	this->debugger = new Debugger(this, this);
	this->debuggerMdiSubWin = this->mdiArea->addSubWindow(this->debugger);
	this->debuggerMdiSubWin->setWindowIcon(this->debugger->windowIcon());
	this->debuggerMdiSubWin->showMinimized();
	this->actionMenuDebugger->setCheckable(true);
	connect(this->actionMenuDebugger, SIGNAL(triggered()), this->debuggerMdiSubWin, SLOT(setFocus()));
	this->interpreter = this->debugger->interpreter = this->console->interpreter.get();
	this->consoleMdiSubWin = this->mdiArea->addSubWindow(this->console);
	this->consoleMdiSubWin->setWindowIcon(this->console->windowIcon());
	this->actionMenuConsole->setCheckable(true);
	connect(this->actionMenuConsole, SIGNAL(triggered()), this->console->inputTextEdit, SLOT(setFocus()));
	this->consoleMdiSubWin->show();
	connect(this->actionNew, SIGNAL(triggered()), this, SLOT(onFileNew()));
	connect(this->actionOpen, SIGNAL(triggered()), this, SLOT(onFileOpen()));
	connect(this->actionExit, SIGNAL(triggered()), this, SLOT(onFileExit()));
	connect(this->actionClose, SIGNAL(triggered()), this, SLOT(onWindowClose()));
	connect(this->actionClose_All, SIGNAL(triggered()), this, SLOT(onWindowCloseAll()));
	connect(this->actionTile, SIGNAL(triggered()), this, SLOT(onWindowTile()));
	connect(this->actionCascade, SIGNAL(triggered()), this, SLOT(onWindowCascade()));
	connect(this->actionNext, SIGNAL(triggered()), this, SLOT(onWindowNext()));
	connect(this->actionPrevious, SIGNAL(triggered()), this, SLOT(onWindowPrevious()));
	connect(this->actionAbout, SIGNAL(triggered()), this, SLOT(onHelpAbout()));
	connect(this->actionFont, SIGNAL(triggered()), this, SLOT(onOptionsFont()));
	connect(this->actionEmacsLike, SIGNAL(triggered()), this, SLOT(onOptionsEmacsLike()));
	connect(this->actionFileSaveOutputWindow, SIGNAL(triggered()), this, SLOT(onFileSaveOutputWindow()));
	connect(this->mdiArea, SIGNAL(subWindowActivated(QMdiSubWindow*)), this, SLOT(onSubWindowActivated(QMdiSubWindow*)));
}

void MainWindow::onOptionsFont() {
	 bool ok;
	 QFont font = QFontDialog::getFont(&ok, this);
	 if (ok) {
		 QDocument::setFont(font);
	     this->application->setFont(font);
	     QSettings settings;
	     settings.setValue(SETTING_FONT_FAMILY, font.family());
	     settings.setValue(SETTING_FONT_POINTSIZE, font.pointSize());
	     settings.setValue(SETTING_FONT_WEIGHT, font.weight());
	     settings.setValue(SETTING_FONT_ITALIC, font.italic());
	 }
}

void MainWindow::onOptionsEmacsLike() {
	 QSettings settings;
	 settings.setValue(SETTING_EMACS_LIKE_KEYBINDINGS, this->actionEmacsLike->isChecked());
	 if (this->mdiArea->subWindowList().size()!=2 /* that is, console and debugger */)
		 QMessageBox::information(this, "Warning", "Note: this choice affects only newly-opened editors, already opened editors keep their key-bindings.");
}

intrusive_ptr<IncrementProgressBarQC> IncrementProgressBarQC::theInstance(new IncrementProgressBarQC);

void IncrementProgressBarQC::execute(Console *console) {
	QProgressBar *bar = console->packageLoadingProgressBar;
	const int newValue = bar->value()+1;
	bar->setValue(newValue);
	if (newValue==bar->maximum()) {
		console->bottomLeftVL->removeWidget(console->packageLoadingProgressBar);
		console->inputTextEdit->setEnabled(true);
		console->enterButton->setEnabled(true);
		delete console->packageLoadingProgressBar;
		console->packageLoadingProgressBar = 0;
		console->inputTextEdit->setFocus();
	}
}

bool Console::myLoadPackages() {
	bool result = true;
	this->outputStream->print("", false, IdeOutputStream::normalFormat); // disable the HL for the following lines
        if (this->packageList.empty())
        {
          string LINE; LINE.resize(47, '-');  // 47 because it looks right on my screen
          this->outputStream->print(LINE)->newline();
          this->outputStream->print(PACKAGES_NOT_FOUND)->newline();
          this->outputStream->print(LINE)->newline();
        }
	for (const string &fullname: this->packageList)
        {
          this->postCommand(IncrementProgressBarQC::theInstance);
          if (result)
          {
// 2015-07-28 JAA			this->outputStream->print(PACKAGE_AUTOLOAD_LOADING, false, IdeOutputStream::normalFormat);
// 2015-07-28 JAA			this->outputStream->print(fullname)->newline();
            interpreter->readAndExecute(fullname, true, true);
            if (interpreter->errorReporter->getErrorCount())
            {
// JAA 2015-07-29				this->outputStream->print(PACKAGE_AUTOLOAD_FAILURE_MESSAGE)->flush();
              result = false;
            }
          }
          // 2015-07-28 suppress printing of "skipping package..."
          // else
          // {
          // 	this->outputStream->print(PACKAGE_AUTOLOAD_SKIPPING_PKG_DUE_TO_FAILURE, false, IdeOutputStream::normalFormat);
          // 	this->outputStream->print(fullname)->newline();
          // }
	}
  this->outputStream->print(CoCoA5BannerNonFixedWidthFonts())->newline();
  if (!result)
    this->outputStream->print(PACKAGE_AUTOLOAD_FAILURE_MESSAGE)->flush();
  this->postCommand(IncrementProgressBarQC::theInstance);
  return result;
}

void MainWindow::onSubWindowActivated(QMdiSubWindow *activatedWindow) {
	//cout << "activatedWindow " << static_cast<void *>(activatedWindow) << endl;
	const bool thereIsWin = activatedWindow!=0;
	this->actionClose->setEnabled(thereIsWin && activatedWindow!=this->consoleMdiSubWin && activatedWindow!=this->debuggerMdiSubWin);
	if (thereIsWin) {
		if (activatedWindow==this->consoleMdiSubWin) {
			this->actionMenuConsole->setChecked(true);
			this->console->inputTextEdit->setFocus();
		}
		else if (activatedWindow==this->debuggerMdiSubWin)
			this->actionMenuDebugger->setChecked(true);
		else if (SourceEditor *ed = dynamic_cast<SourceEditor *>(activatedWindow->widget())) {
			ed->menuAction->setChecked(true);
			ed->getEditor()->setFocus();
		} else {
			assert(false); // what's going on?!?
		}
	}
}

void MainWindow::updateStatusLabel() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	const Interpreter::InterpreterStatus status = this->console->interpreter->getStatus();
	if (status==this->oldInterpreterStatus)
		return;
	QString message;
	bool buttonsEnabled = false;
	switch ((this->oldInterpreterStatus = status)) {
	case Interpreter::IS_WAITING_FOR_COMMAND:
		this->interpreter->singleStepExecution = false;
		if (this->consoleMdiSubWin->isMinimized())
				this->consoleMdiSubWin->showNormal();
		this->console->inputTextEdit->setFocus();
		message = "The interpreter is <b>waiting</b> for a command";
		break;
	case Interpreter::IS_WAITING_FOR_COMMAND_COMPLETION:
		message = "The interpreter is <b>waiting</b> for the ending of the current command";
		break;
	case Interpreter::IS_RUNNING:
		message = "The interpreter is <b>running</b>";
		buttonsEnabled = true;
		break;
	case Interpreter::IS_RUNNING_BUILTIN:
		message = "The interpreter is <b>running</b> a built-in function";
		buttonsEnabled = true;
		break;
	case Interpreter::IS_PAUSED:
		message = "The interpreter is <b>paused</b>";
		if (this->debuggerMdiSubWin->isMinimized())
			this->debuggerMdiSubWin->showNormal();
		this->debuggerMdiSubWin->setFocus();
		break;
	case Interpreter::IS_ENDED:
		message = "The interpreter has <b>quit</b>";
		this->onFileExit();
		break;
	}
	this->console->pauseButton->setEnabled(buttonsEnabled);
	this->console->interruptButton->setEnabled(buttonsEnabled);
	this->statusLabel->setText(message);
	this->debugger->update();
}

void Console::onPauseClicked() {
	this->interpreter->singleStepExecution = true;
}

void Console::onInterruptClicked()
{
  if (QMessageBox::question(
        this,
        tr("Interrupt confirmation"),
        tr("Are you sure you want to abort the current computation?\n\nNote: some built-in fn-procs cannot be interrupted, so it might take a while before the computation actually halts"),
        QMessageBox::Yes|QMessageBox::No,
        QMessageBox::No)==QMessageBox::Yes)
    CoCoA::SetSignalReceived(2/*SIGINT*/); // actual signal does not matter (I believe)
}

void MainWindow::onFileOpen() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QString fn = QFileDialog::getOpenFileName(this, "Open CoCoA 5 source...", "", tr("Sources (*.cocoa5 *.cpkg5)"));
	if (!fn.size())
		return;
	SourceEditor *sourceEditor = this->editorFor(fn, true);
	if (!sourceEditor) {
		sourceEditor = this->onFileNew();
		sourceEditor->load(fn);
	}
}

void MainWindow::onFileExit() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	if (this->close())
		this->application->quit();
}

void MainWindow::onWindowClose() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QMdiSubWindow *w = this->mdiArea->activeSubWindow();
	if (w)
		w->close();
}

void MainWindow::onWindowCloseAll() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->mdiArea->closeAllSubWindows();
}

void MainWindow::onWindowTile() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->mdiArea->tileSubWindows();
}

void MainWindow::onWindowCascade() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->console->lower();
	this->mdiArea->cascadeSubWindows();
}

void MainWindow::onWindowNext() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->mdiArea->activateNextSubWindow();
}

void MainWindow::onWindowPrevious() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	this->mdiArea->activatePreviousSubWindow();
}

void MainWindow::onHelpAbout() {
	assert(std::this_thread::get_id()==Console::guiThreadId);
	QMessageBox::about(this, "About C5",
			QString::fromStdString("<H1 align=center>C5</H1>"
			"<H3 align=center>CoCoA 5 Integrated Development Environment</H3>"
			"<p><p align=center>Alpha-version, for testing purposes only."
			"<p align=center>Please, <font color=red><strong>DO NOT DISTRIBUTE</strong></font> this file."
			"<p><p><p>For more information about CoCoA, please visit the official <a href=\"http://cocoa.dima.unige.it/\">website</a>"));
}

namespace {
	void loadSettings(QApplication &app) {
		QSettings settings;
		QVariant fontFamily, fontPointSize, fontWeight, fontItalic;
		fontFamily = settings.value(SETTING_FONT_FAMILY);
		if (!fontFamily.isNull()) {
			fontPointSize = settings.value(SETTING_FONT_POINTSIZE);
			fontWeight = settings.value(SETTING_FONT_WEIGHT);
			fontItalic = settings.value(SETTING_FONT_ITALIC);
			//cout << "setting font " << fontFamily.toString().toStdString() << ", point-size=" << fontPointSize.toInt()
			//		<< ", weight=" << fontWeight.toInt() << ", italic=" << fontItalic.toBool() << endl;
			app.setFont(QFont(fontFamily.toString(), fontPointSize.toInt(), fontWeight.toInt(), fontItalic.toBool()));
		}
	}
}

int launchTheIDE(WarningSeverity warningLevel, bool warnAboutCocoa5, const vector<string> &packageList, bool fullCoCoALibError) {
	QCoreApplication::setOrganizationName("CoCoA Team");
	QCoreApplication::setOrganizationDomain("cocoa.dima.unige.it");
	QCoreApplication::setApplicationName("C5");
	char *argv[] = {const_cast<char *>("C5"), 0};
	int argc = 1;
	// QApplication::setStyle(new QPlastiqueStyle);
	QApplication::setStyle(QStyleFactory::create("Fusion"));
	QApplication app(argc, argv);
	initQCodeEdit();
	app.setWindowIcon(QIcon(":/images/CoCoALogo-icon.png"));
	loadSettings(app);
	QEditor::setDefaultFlags(QEditor::defaultFlags()|QEditor::LineWrap);
	QDocument::setFont(app.font());
	CoCoA::IDE::MainWindow main(0, &app, warningLevel, warnAboutCocoa5, packageList, fullCoCoALibError);
	main.showNormal();
	return app.exec();
}

bool loadPackages(Console *console) {
	return console->myLoadPackages();
}

} // namespace IDE

} // namespace CoCoA

#endif // #ifdef C5IDE
