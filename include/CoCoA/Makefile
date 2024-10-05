# Makefile for directory CoCoALib/include/CoCoA/.

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=include/CoCoA/

# Default is to rebuild the unified header (in library.H)
.PHONY: default
default:
	@./MakeUnifiedHeader.sh --check


.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	/bin/rm -f  .mutex  ./*~  ./.*~  ./.\#*   PREPROCESSOR_DEFNS.H-old  library.H.old

.PHONY: veryclean veryclean-local
veryclean: veryclean-local
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local: clean-local
	/bin/rm -f  library.H  PREPROCESSOR_DEFNS.H
