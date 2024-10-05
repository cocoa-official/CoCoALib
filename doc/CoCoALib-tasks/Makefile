# Makefile for CoCoALib/doc/CoCoALib-tasks/ directory.

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk # sets CXX
CWD=doc/CoCoALib-tasks

EXEC = HTMLTasks
C_SRCS = HTMLTasks.C
XML_SRCS = tasks.xml
HTML_OUT = CoCoALib-tasks.html  CoCoALib-CompletedTasks.html


.PHONY: default
.NOTPARALLEL: default
default: $(HTML_OUT)

$(HTML_OUT): $(XML_SRCS)  $(C_SRCS)
	$(MAKE) $(EXEC)
	@echo "Generating HTML pages ..."
	@./$(EXEC)
	@/bin/cp $(HTML_OUT)  ../

$(EXEC): $(C_SRCS)
	@echo "Generating executable" $(EXEC) "..."
	@$(CXX) -g -Wextra  HTMLTasks.C  -o $(EXEC)
	/bin/rm -rf  $(EXEC).dSYM

.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f  ./*~  ./.*~  ./.\#*
	@/bin/rm -f  HTMLTasks  HTMLTasks.o	

.PHONY: veryclean  veryclean-local
veryclean: veryclean-local
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local: clean-local
	@/bin/rm -f  CoCoALib-tasks.html  CoCoALib-CompletedTasks.html


# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/doc/CoCoALib-tasks/Makefile,v 1.18 2022/11/04 15:23:57 abbott Exp $Log
# $Log: Makefile,v $
# Revision 1.18  2022/11/04 15:23:57  abbott
# Summary: Replaced hardwired "g++" by $CXX -- fixed by Ulrich
#
# Revision 1.17  2022/03/16 18:16:47  abbott
# Summary: Added some spaces
#
# Revision 1.16  2022/02/22 09:35:57  abbott
# Summary: Corrected first comment
#
# Revision 1.15  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.14  2014/07/28 14:42:51  abbott
# Summary: Improved the cleaning targets
# Author: JAA
#
# Revision 1.13  2014/03/14 11:00:27  abbott
# Summary: clean and veryclean targets now delete *~ and .*~ files too
# Author: JAA
#
# Revision 1.12  2011/08/31 09:38:06  bigatti
# -- made target "default" NOTPARALLEL
#
