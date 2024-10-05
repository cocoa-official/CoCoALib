# Makefile for CoCoALib/doc/ directory.

COCOA_ROOT=..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=doc/

SRC_DIR=txt
TEX_DIR=tex
TASK_DIR=CoCoALib-tasks
HTML_DIR=html
AUX_DIR=aux-txt2tags
DOC_SRCS = $(wildcard $(SRC_DIR)/*.txt)
DOC_HTML = $(patsubst $(SRC_DIR)/%.txt,$(HTML_DIR)/%.html,$(DOC_SRCS))
TASKS_OUT = CoCoALib-tasks.html  CoCoALib-CompletedTasks.html

.SUFFIXES:

# Rule for compiling *.txt files into *.html in HTML_DIR
$(HTML_DIR)/%.html:  $(SRC_DIR)/%.txt
	@echo "Generating $*.html"
	@txt2tags -t html -o $@ -i $<  >/dev/null


.PHONY : alldoc
alldoc: examples-index  htmldoc  texdoc  tasktables

$(TEX_DIR):
	@/bin/mkdir -p $(TEX_DIR)

$(TEX_DIR)/CoCoALib.tex:  $(DOC_SRCS)  $(TEX_DIR)  $(AUX_DIR)/DocTeX.t2t
	@echo "Generating CoCoALib.tex..."
	@sed -e "s/VERSION_NUMBER/${COCOALIB_VERSION}/" -e "s/DATE_OF_CREATION/`date +%Y-%m-%d`/"  $(AUX_DIR)/DocTeX.t2t  > $(AUX_DIR)/DocTeX-VERSION.t2t
	@txt2tags -o $(TEX_DIR)/CoCoALib.tex  $(AUX_DIR)/DocTeX-VERSION.t2t   >/dev/null

$(TEX_DIR)/cocoalib-doc.sty:  $(TEX_DIR)  $(AUX_DIR)/cocoalib-doc.sty
	@/bin/cp -p $(AUX_DIR)/cocoalib-doc.sty $(TEX_DIR)


CoCoALib.pdf:  $(COCOA_ROOT)/configuration/version  $(DOC_SRCS)  $(AUX_DIR)/DocTeX.t2t  $(AUX_DIR)/cocoalib-doc.sty
	@command -v txt2tags >/dev/null 2>&1 ; \
	 if [ $$? -ne 0 ]; \
	 then \
	   echo ">>>>>>>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<"; \
	   echo ">>>  Cannot build CoCoALib PDF documentation!     <<<"; \
	   echo ">>>  Install txt2tags from http://txt2tags.org/   <<<"; \
	   echo ">>>  (make sure you get version 3.0 or 2.7-DEV)   <<<"; \
	   echo ">>>>>>>>>>>>>>>>>>>> END WARNING <<<<<<<<<<<<<<<<<<<<"; \
	   echo; \
	   sleep 2; \
	   exit; \
	 fi; \
	 $(MAKE)  $(TEX_DIR)/CoCoALib.tex  $(TEX_DIR)/cocoalib-doc.sty; \
	 if [ \! -f "$(TEX_DIR)/CoCoALib.tex" ]; then echo; echo "ERROR!!!  FAILED TO MAKE $(TEX_DIR)/CoCoALib.tex"; echo; exit 1; fi; \
	 echo "Doing LaTeX-->PDF"; \
	 cd $(TEX_DIR); \
	 echo "[[CoCoALib:LaTeX-->PDF]]  Pass 1 of 3"; pdflatex -halt-on-error CoCoALib.tex < /dev/null > /dev/null  2>&1; \
	 if [ $$? -ne 0 ]; then echo "!!!LaTeX error!!! -- see file tex/CoCoALib.log"; exit 1; fi; \
	 echo "[[CoCoALib:LaTeX-->PDF]]  Pass 2 of 3"; pdflatex -halt-on-error CoCoALib.tex > /dev/null  2>&1; \
	 echo "[[CoCoALib:LaTeX-->PDF]]  Final pass";  pdflatex -halt-on-error CoCoALib.tex > /dev/null  2>&1; \
	 grep -F "LaTeX Warning:" CoCoALib.log; \
	 /bin/mv CoCoALib.pdf ..

# examples/index.html is really part of the documentation...
.PHONY: examples-index
examples-index:
	@cd $(COCOA_ROOT)/examples; $(MAKE) -f Makefile-indexhtml


.PHONY: texdoc
texdoc: CoCoALib.pdf


.PHONY: htmldoc
htmldoc:  $(COCOA_ROOT)/configuration/version  $(HTML_DIR)/cocoalib-doc.css  $(HTML_DIR)/index.html  $(HTML_DIR)/index-alpha.html
	@$(MAKE) -q $(DOC_HTML); \
	 if [ $$? -ne 0 ]; \
	 then \
	   command -v txt2tags >/dev/null 2>&1 ; \
	   if [ $$? -ne 0 ]; \
	   then \
	     echo ">>>>>>>>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<<<<<<<"; \
	     echo ">>>  Cannot build CoCoALib HTML documentation!    <<<"; \
	     echo ">>>  Install txt2tags from http://txt2tags.org/   <<<"; \
	     echo ">>>  (make sure you get version 3.0 or 2.7-DEV)   <<<"; \
	     echo ">>>>>>>>>>>>>>>>>>>> END WARNING <<<<<<<<<<<<<<<<<<<<"; \
	     echo; \
	     sleep 2; \
	     exit; \
	   fi; \
	   $(MAKE) -s $(DOC_HTML); \
	 fi


$(HTML_DIR):
	@/bin/mkdir -p $(HTML_DIR)

$(HTML_DIR)/cocoalib-doc.css:  $(HTML_DIR)  $(AUX_DIR)/cocoalib-doc.css
	@/bin/cp -p  $(AUX_DIR)/cocoalib-doc.css  $(HTML_DIR)

$(HTML_DIR)/index-alpha.html:  $(HTML_DIR)  $(AUX_DIR)/index-alpha.html
	@/bin/rm -f $(HTML_DIR)/index-alpha.html
	@sed -e "s/VERSION_NUMBER/${COCOALIB_VERSION}/"  $(AUX_DIR)/index-alpha.html  > $(HTML_DIR)/index-alpha.html

$(HTML_DIR)/index.html:  $(HTML_DIR)  $(AUX_DIR)/index.html
	@/bin/rm -f $(HTML_DIR)/index.html
	@sed -e "s/VERSION_NUMBER/${COCOALIB_VERSION}/"  $(AUX_DIR)/index.html  > $(HTML_DIR)/index.html

.PHONY: tasktables
tasktables:
	@cd $(TASK_DIR); $(MAKE) -s


.PHONY:  clean  clean-local  clean-subdirs
clean:  clean-local  clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f ./*~  ./.*~  ./.\#*
	@/bin/rm -rf $(TEX_DIR)

clean-subdirs:
	@cd $(TASK_DIR);  $(MAKE) -s clean
	@cd $(SRC_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*
	@cd $(AUX_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*  DocTeX-VERSION.t2t


.PHONY:  veryclean  veryclean-local  veryclean-subdirs
veryclean: veryclean-local  veryclean-subdirs
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local:  clean-local
	@/bin/rm -f $(HTML_DIR).tgz  CoCoALib.pdf  $(TASKS_OUT)  TIMESTAMP
	@/bin/rm -rf $(HTML_DIR)

veryclean-subdirs:
	@cd $(TASK_DIR);  $(MAKE) -s veryclean
	@cd $(SRC_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*
	@cd $(AUX_DIR);   /bin/rm -f  ./*~  ./.*~  ./.\#*  DocTeX-VERSION.t2t


# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/doc/Makefile,v 1.65 2024/03/18 20:54:56 abbott Exp $Log
# $Log: Makefile,v $
# Revision 1.65  2024/03/18 20:54:56  abbott
# Summary: INcluded version info for txt2tags (redmine 1662)
#
# Revision 1.64  2022/11/16 19:24:46  abbott
# Summary: Replaced fgrep (obsol) by grep -F
#
# Revision 1.63  2022/11/04 15:22:44  abbott
# Summary: Removed some cruft
#
# Revision 1.62  2022/04/10 16:14:08  abbott
# Summary: First run pdflatex takes input from /dev/null (o/w seems to hang)
#
# Revision 1.61  2022/04/06 16:34:53  abbott
# Summary: KISS: just try using builtin txt2tags; if fails, print warning. (redmine 1676)
#
# Revision 1.60  2022/04/04 13:31:57  abbott
# Summary: Removed debugging print -- first day back after holiday (sigh!)
#
# Revision 1.59  2022/04/04 13:29:58  abbott
# Summary: Oops: should have changed two places previously (hopefully now all OK).
#
# Revision 1.58  2022/04/04 13:26:47  abbott
# Summary: Fixed a coupld of bugs (see also redmine 1676)
#
# Revision 1.57  2022/03/25 16:27:11  bigatti
# Summary: new date procedure for LaTeX; added my personal path for t2t
#
# Revision 1.56  2022/03/20 21:34:40  abbott
# Summary: Moved location of script txt2tags
#
# Revision 1.55  2022/03/18 09:19:36  abbott
# Summary: Cleaned up generation of CoCoALib.pdf (copied from CoCoAManual/Makefil); redmine 1674
#
# Revision 1.54  2022/03/18 08:01:08  bigatti
# Summary: added  fgrep "LaTeX Warning:" CoCoALib.log;  as for cocoa5
#
# Revision 1.53  2022/03/17 18:26:58  abbott
# Summary: Neater progress mesgs
#
# Revision 1.52  2022/03/17 15:53:03  abbott
# Summary: Corrected syntax problem with test []
#
# Revision 1.51  2022/03/16 18:15:33  abbott
# Summary: Improved comment
#
# Revision 1.50  2022/03/11 14:32:05  bigatti
# Summary: added check of readability of local file txt2tags
#
# Revision 1.49  2022/03/09 15:33:34  abbott
# Summary: Introduced TXT2TAGS variable (since my built-in version does not work)
#
# Revision 1.48  2022/03/07 10:11:01  abbott
# Summary: Added space
#
# Revision 1.47  2020/06/19 14:45:25  abbott
# Summary: Improved progress mesgs about compiling latex to pdf
#
# Revision 1.46  2020/06/16 17:59:52  abbott
# Summary: Prints better progress mesgs
#
# Revision 1.45  2020/03/07 14:44:05  abbott
# Summary: clean now removes  DocTeX-VERSION.t2t
#
# Revision 1.44  2019/12/18 08:16:48  abbott
# Summary: Revised configure; now must explicitly request to build Qt GUI; changed VERSION into COCOALIB_VERSION; conseq changes to Makefiles
#
# Revision 1.43  2019/10/10 16:41:01  bigatti
# -- added dependency on configuration/version
#
# Revision 1.42  2018/08/02 14:55:42  abbott
# Summary: Replaced hack for examples/index.html by a clean impl
#
# Revision 1.41  2018/04/23 09:17:30  abbott
# Summary: Added dirty hack so that CopyInfo can be compiled even when libcocoa.a does not exist.
#
# Revision 1.40  2018/04/18 15:38:22  abbott
# Summary: Added dirty trick to avoid compilation problem with examples/CopyInfo.C
#
# Revision 1.39  2018/02/16 10:52:22  abbott
# Summary: Replaced /usr/bin/sed by sed (for better portability)
#
# Revision 1.38  2017/12/18 10:56:49  abbott
# Summary: Changed recursive calls of "make" into $(MAKE)
#
# Revision 1.37  2017/12/15 16:13:00  bigatti
# removed pointless backslash
#
# Revision 1.36  2017/12/15 15:26:14  bigatti
# added version in tex doc
#
# Revision 1.35  2017/11/26 08:32:06  abbott
# Summary: now redirect stderr into log file when running latex
#
# Revision 1.34  2017/11/16 14:18:46  abbott
# Summary: Added new feature: builting alldoc also creates examples/index.html
#
# Revision 1.33  2017/02/14 17:06:28  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.32  2014/07/28 14:42:33  abbott
# Summary: Improved the cleaning targets
# Author: JAA
#
# Revision 1.31  2014/03/14 10:59:55  abbott
# Summary: clean and veryclean targets now delete *~ and .*~ files too
# Author: JAA
#
# Revision 1.30  2014/01/29 17:57:05  abbott
# Summary: Rectified silly mistake
# Author: JAA
#
# Revision 1.29  2014/01/29 17:56:09  abbott
# Summary: Improved veryclean-local target
# Author: JAA
#
# Revision 1.28  2013/10/16 14:12:58  abbott
# Added check that pdflatex has run without error; prints mesg & bails if error occurs.
#
# Revision 1.27  2012/10/15 08:52:26  abbott
# Added clean-local to veryclean-local target.
#
# Revision 1.26  2012/05/30 11:00:54  bigatti
# -- added more dependencies
#
# Revision 1.25  2012/02/02 16:02:40  bigatti
# -- removed creation of index-new.html (had been renamed)
#
# Revision 1.24  2011/12/23 14:51:44  bigatti
# -- copying also the two new index pages
#
# Revision 1.23  2011/10/18 12:03:44  abbott
# Redirected stderr when calling "which" (twice)
#
# Revision 1.22  2011/09/23 13:28:21  abbott
# Corrected unified command script for making CoCoALib.pdf
# (removed some "cd" comands which are no longer needed).
#
# Revision 1.21  2011/09/22 16:22:29  abbott
# Minor correction: the CoCoALib.pdf rule is now a single shell script.
#
# Revision 1.20  2011/09/22 15:52:07  abbott
# Considerably changed, and cleaned.
# This Makefile now prints out a "warning" message if some documentation
# needs to be updated but txt2tags is absent.
#
# Revision 1.19  2011/08/31 09:37:18  bigatti
# -- undid last change about target tasktables
#
