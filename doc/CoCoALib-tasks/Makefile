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
