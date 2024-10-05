#ifndef CoCoA5_C5IDE_H
#define CoCoA5_C5IDE_H

//   Copyright (c) 2010,2022   John Abbott,  Anna M. Bigatti
//   Original author: 2010 Giovanni Lagorio (lagorio@disi.unige.it)
//                    2020-2022  revitalized by Ulrich von der Ohe 
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

// Testing shows that the next 3 includes are not needed -- why???
#include <thread>
//using std::thread;
#include <mutex>
//using std::mutex;
#include <condition_variable>
//using std::condition_variable;

#include <deque>
#include <QApplication>
#include <QDialog>
#include <QTimer>
#include <QEvent>
#include <QKeyEvent>
#include <QMessageBox>
#include <QSyntaxHighlighter>
#include <QMdiSubWindow>
#include <QMdiArea>
#include <QFormLayout>
#include <QTreeWidget>
#include "qcodeedit.h"
#include "qformat.h"
#include "qformatscheme.h"

#include "Interpreter.H"
#include "ui_Console.h"
#include "ui_MainWindow.h"
#include "ui_SourceEditor.h"
#include "ui_Debugger.h"

namespace CoCoA {
namespace IDE {

struct ColorScheme {
	QFormat stringLiterals;
	QFormat constants;
	QFormat comments;
	QFormat keywords;
	QFormat types;
	QFormat numericalLiterals;
	QFormat unknownChars;
	QFormat braceMatch;
	QFormat braceMismatch;
	void applyTo(QFormatScheme *scheme) const;
};

class CocoaHighlighter : public QSyntaxHighlighter {
	struct HighlightingInfo : public QTextBlockUserData {
		bool doHighlight;
		int errorFrom, errorLen;
		HighlightingInfo(bool doHighlist) :
			doHighlight(doHighlist),
			errorFrom(-1),
			errorLen(-1)
		{}
	};
	Q_OBJECT
	bool highlightUnknownChars;
	QRegExp blanksRE;
	QTextCharFormat unknownFormat;
	QTextCharFormat keywordFormat;
	QRegExp keywordRE;
	QTextCharFormat constantFormat;
	QRegExp constantRE;
	QTextCharFormat literalStringFormat;
	QRegExp singleQuoteliteralStringEndRE;
	QRegExp doubleQuoteliteralStringEndRE;
	QRegExp tripleQuoteliteralStringEndRE;
	QTextCharFormat identifierFormat;
	QRegExp identifierRE;
	QTextCharFormat commentFormat;
	QRegExp commentEndRE;
	QRegExp singlelineCommentEndRE;
	QTextCharFormat numericLiteralFormat;
	QRegExp numericLiteralRE;
	QTextCharFormat typeFormat;
	QRegExp typeRE;
	QTextCharFormat operatorFormat;
	QRegExp parenthesisRE;
	QRegExp operatorRE;
	enum {
		NONE,
		IN_SINGLEQUOTE_STRING,
		IN_DOUBLEQUOTE_STRING,
		IN_TRIPLEQUOTE_STRING,
		IN_COMMENT
	};
	bool unclosedComment, unclosedStringLiteral;
	void highlightBlock(const QString &text);
public:
	bool enabled;
	const bool alwaysEnabled;
	inline bool thereIsUnclosedComment() const { return this->unclosedComment; }
	inline bool thereIsUnclosedStringLiteral() const { return this->unclosedStringLiteral; }
	CocoaHighlighter(QTextDocument *parent, bool highlightUnknownChars, bool alwaysEnabled);
	void setFormats(const ColorScheme &scheme);
};

class Console;

struct QueuedCommand : public LexerNS::ReferenceCountedObject {
	virtual void execute(Console *console) = 0;
	virtual ~QueuedCommand() throw () {}
};

class PrintQC : public QueuedCommand {
	std::string s;
	bool highlight;
	QTextCharFormat format;
public:
	PrintQC(const std::string &s, bool highlight, QTextCharFormat format) :
		s(s),
		highlight(highlight),
		format(format)
	{}
	void execute(Console *console);
	static boost::intrusive_ptr<PrintQC> newline, newlineHL;
	virtual ~PrintQC() throw () {}
};

class FlushQC : public QueuedCommand {
	FlushQC() {}
public:
	volatile bool flushed;
	static boost::intrusive_ptr<FlushQC> theInstance; // this is ok as long as we have a single interpreter per IDE (which is the case ;-) )
	void execute(Console *console);
	std::mutex mut;
	std::condition_variable condVar;
	virtual ~FlushQC() throw () {}
};

class IncrementProgressBarQC : public QueuedCommand {
	IncrementProgressBarQC() {}
public:
	static boost::intrusive_ptr<IncrementProgressBarQC> theInstance;
	void execute(Console *console);
	virtual ~IncrementProgressBarQC() throw () {}
};

class ClearReportedLocationsQC : public QueuedCommand {
	ClearReportedLocationsQC() {}
public:
	static boost::intrusive_ptr<ClearReportedLocationsQC> theInstance;
	void execute(Console *console);
	virtual ~ClearReportedLocationsQC() throw () {}
};

class ReportErrorQC : public QueuedCommand {
	std::string filename;
	int lineNumber;
	int columnNumber;
public:
	ReportErrorQC(const std::string &filename, int lineNumber, int columnNumber) :
		filename(filename),
		lineNumber(lineNumber),
		columnNumber(columnNumber)
	{}
	void execute(Console *console);
	virtual ~ReportErrorQC() throw () {}
};

class IdeOutputStream : public InterpreterNS::OSTREAM {
	Console *console;
public:
	IdeOutputStream(Console *console);
	static QTextCharFormat normalFormat;
	static QTextCharFormat boldFormat;
	boost::intrusive_ptr<InterpreterNS::OSTREAM> print(const std::string &s);
	boost::intrusive_ptr<InterpreterNS::OSTREAM> print(boost::intrusive_ptr<const InterpreterNS::RightValue> v);
	boost::intrusive_ptr<InterpreterNS::OSTREAM> newline();
	boost::intrusive_ptr<InterpreterNS::OSTREAM> newlineHL();
	void print(const std::string &s, bool highlight, QTextCharFormat format);
	void flush(); // NOTE: flush must *NOT* be called from the GUI thread
	virtual ~IdeOutputStream() throw () {}
};

class IdeLineProvider : public LexerNS::InteractiveLineProvider {
	Console *console;
	std::mutex mut;
	std::condition_variable condVar;
	std::deque<std::string> enteredLines;
	std::deque<std::string> enteredFullCommands;
protected:
	bool doReadNextLine(const LexerNS::LexerStatus &ls, const ParserNS::ParserStatus &ps, std::string &chars);
public:
	IdeLineProvider(Console *console) : console(console) {}
	void enterLine(const std::string &line, bool isFullCommand);
	~IdeLineProvider() throw () {}
};

class IdeErrorReporter : public LexerNS::ErrorReporter {
	Console *console;
	QTextCharFormat errorFormat;
	QTextCharFormat warningFormat;
	QTextCharFormat calledByFormat;
	QTextCharFormat whereFormat;
	QTextCharFormat contextFormat;
	QTextCharFormat boldFormat;
	QTextCharFormat highlightErrorFormat;
	void outputHighlightedChars(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void outputHighlightedLine(boost::intrusive_ptr<const LexerNS::Line> line, std::size_t from, std::size_t to, bool keepsHilighting);
protected:
	void implReportWarning(const std::string &msg);
	void implReportWarning(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void implReportInterrupt(const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
	void implReportError(const std::string &msg);
	void implReportError(const std::string &msg, const LexerNS::CharPointer &from, const LexerNS::CharPointer &to);
public:
	IdeErrorReporter(LexerNS::WarningSeverity warningLevel, Console *console, const boost::intrusive_ptr<IdeOutputStream> outputStream);
	bool reportLineNumberWhenMeaningful(const LexerNS::CharPointer &fromPos, const LexerNS::CharPointer &toPos, bool printColumns, bool includeHeader);
	void printWhere();
	void printError();
	void printWarning();
	void printCalledBy();
	void printContext();
	void printBold(const std::string &s);
	~IdeErrorReporter() throw () {}
};

class MainWindow;

class Debugger : public QWidget, public Ui::Debugger {
	Q_OBJECT
	MainWindow *mainWindow;
	QTextCharFormat format;
	std::vector<const InterpreterNS::Frame *> listedFrames;
	CocoaHighlighter *codeHL;
public:
	Debugger(QWidget *parent, MainWindow *mainWindow);
	InterpreterNS::Interpreter *interpreter;
	InterpreterNS::RuntimeEnvironment *re;
	void update();
	void setHighlighterFormats();
protected:
	void closeEvent(QCloseEvent *event);
private slots:
	void onStepIntoClicked();
	void onStepOverClicked();
	void onStepOutClicked();
	void onStepOutFnProcClicked();
	void onContinueClicked();
	void onStopClicked();
	void fillCallStack();
	void fillCode();
	void fillLocals();
	void clearAndDisable();
	void onCallStackSelectionChanged();
	void onAutoFillLocalsStateChanged();
	void onFillLocalsClicked();
	void onHidingChoicesChanged();
private:
	void addLocal(const std::string &name, const AST::StaticEnv::VarData &varData, const InterpreterNS::Frame *currentFrame);
	void addLocal(const std::string &name, const InterpreterNS::Frame *f, int index);
	void addLocalError(const QString &name, const QString &message);
};

class Console : public QWidget, public Ui::Console {
	Q_OBJECT
	std::deque<std::string> history;
	int currentHistoryPosition;
	CocoaHighlighter *inputHL;
	CocoaHighlighter *outputHL;
	static const int MAX_HISTORY = 128;
	std::mutex mut;
	std::deque<boost::intrusive_ptr<QueuedCommand> > commands;
	struct ReportedLocation {
		std::string filename;
		int lineNumber;
		int columnNumber;
		ReportedLocation(const std::string &filename, int lineNumber, int columnNumber) :
			filename(filename),
			lineNumber(lineNumber),
			columnNumber(columnNumber)
		{}
	};
	std::vector<ReportedLocation> reportedLocations;
	MainWindow *mainWindow;
	friend bool loadPackages(Console *);
	const std::vector<std::string> packageList;
	bool myLoadPackages();
public:
	Console(QWidget *parent, MainWindow *mainWindow, LexerNS::WarningSeverity warningLevel, bool warnAboutCocoa5, const std::vector<std::string> & packageList, bool fullCoCoALibError);
	const boost::intrusive_ptr<IdeOutputStream> outputStream;
	const boost::intrusive_ptr<IdeLineProvider> lineProvider;
	const boost::intrusive_ptr<IdeErrorReporter> errorReporter;
	void print(const std::string &s, bool highlight, QTextCharFormat format);
	void postCommand(boost::intrusive_ptr<QueuedCommand> command);
	int reportLocation(const std::string &filename, int lineNumber, int columnNumber);
	void clearReportedLocations();
	boost::intrusive_ptr<InterpreterNS::Interpreter> interpreter;
	std::thread::id interpreterThreadId;
	static std::thread::id guiThreadId;
	void setHighlighterFormats();
protected:
	void closeEvent(QCloseEvent *event);
public slots:
        void clearOutputWindow(bool AskUser=true); // 2015-07-27 JAA made this public so it can be called by MainWindow::closeEvent
private slots:
	void onEnterCommandClicked();
	void processCommands();
	void onOpenInEditorClicked();
	void onPauseClicked();
	void onInterruptClicked();
private:
	bool eventFilter(QObject *target, QEvent *event);
};

class MainWindow;

class SourceEditor : public QWidget, public Ui::SourceEditor {
	Q_OBJECT
	QCodeEdit *codeEdit;
public:
	SourceEditor(MainWindow *parent);
	QEditor *getEditor() { return this->codeEdit->editor(); }
	QAction * const menuAction;
public slots:
	void save();
	void saveAs();
	void saveAndRun();
	void load(const QString& file);
protected:
	void closeEvent(QCloseEvent *e);
private:
	bool maybeSave();
private slots:
	void onEditorTitleChanged(const QString& title);
	void onEditorContentModified(bool y);
};

class MainWindow : public QMainWindow, public Ui::MainWindow {
	Q_OBJECT
	QApplication * const application;
	Console *console;
	QMdiSubWindow *consoleMdiSubWin;
	Debugger *debugger;
	QMdiSubWindow *debuggerMdiSubWin;
	QMdiArea * const mdiArea;
	InterpreterNS::Interpreter::InterpreterStatus oldInterpreterStatus;
	QLabel * const statusLabel;
	QActionGroup * const menuWindowActionGroup;
	QActionGroup * const colorSchemeActionGroup;
	QAction *actionMenuConsole;
	QAction *actionMenuDebugger;
	static ColorScheme csAnna;
	static ColorScheme csGio;
	static void initColorSchemes();
	ColorScheme *currentColorScheme;
	void setHighlighterFormats();
public:
	MainWindow(QWidget *parent, QApplication *application, LexerNS::WarningSeverity warningLevel, bool warnAboutCocoa5, const std::vector<std::string> & packageList, bool fullCoCoALibError);
	static MainWindow *findMainWindow(QObject *);
	static Console *findConsole(QObject *);
	SourceEditor *editorFor(QString filename, bool activateWindow);
	InterpreterNS::Interpreter *interpreter;
	void updateStatusLabel();
	const ColorScheme *getCurrentColorScheme() { return this->currentColorScheme; }
public slots:
	SourceEditor *onFileNew();
	void onFileOpen();
	void onFileExit();
	void onWindowClose();
	void onWindowCloseAll();
	void onWindowTile();
	void onWindowCascade();
	void onWindowNext();
	void onWindowPrevious();
	void onHelpAbout();
protected:
	void closeEvent(QCloseEvent *e);
private slots:
	void onSubWindowActivated(QMdiSubWindow *activatedWindow);
	void onAnnaCS();
	void onGioCS();
	void onOptionsFont();
	void onOptionsEmacsLike();
	void onFileSaveOutputWindow();
};

} // namespace IDE
} // namespace CoCoA

#endif // #ifdef C5IDE

// Close include guard:
#endif
