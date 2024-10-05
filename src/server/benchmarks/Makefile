# Makefile for CoCoALib/src/server/benchmark directory

COCOA_ROOT=../../..
CWD=src/server/benchmarks/

.PHONY: benchmarks
benchmarks: ../CoCoAServer
	@if [ ! -e "inputs" -o "inputs.tgz" -nt "inputs" ] ; \
	 then \
	   tar xzf inputs.tgz; touch inputs; \
	 fi
	@$(MAKE) gb-charp
	@$(MAKE) gb-char0

.PHONY: gb-char0
gb-char0: ../CoCoAServer
	@./Benchmark0H.command;

.PHONY: gb-charp
gb-charp: ../CoCoAServer
	@./Benchmark32003H.command;

../CoCoAServer:
	@. $(COCOA_ROOT)/configuration/shell-fns.sh; echobox "Making CoCoAServer to run benchmarks"
	@(cd $(COCOA_ROOT); $(MAKE) -s server)

.PHONY: clean
clean:
	@/bin/rm -rf  inputs/
	@echo "Cleaned CoCoALib/$(CWD)"

.PHONY: veryclean
veryclean: clean

#############################################################################

# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/benchmarks/Makefile,v 1.5 2019/10/02 17:59:25 abbott Exp $
# $Log: Makefile,v $
# Revision 1.5  2019/10/02 17:59:25  abbott
# Summary: Replaced source command by more portable . (dot)
#
# Revision 1.4  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.3  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" bu calls to "$(MAKE)".
#
# Revision 1.2  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.1.1.1  2007/03/09 15:16:12  abbott
# Imported files
#
# Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
# Imported files
#
# Revision 1.5  2006/01/23 12:07:09  cocoa
# Further cleaning of Makefiles in response to a reported problem
# about spurious entries in src/Makefile_dependencies.
# Small change to doc for SmallFpImpl.
#
# Revision 1.4  2005/11/22 17:31:32  cocoa
# -- added: "tar xzf inputs.tgz" to target benchmarks
#
# Revision 1.3  2005/11/22 17:28:55  cocoa
# -- fixed : "/." path for calling scripts
#
# Revision 1.2  2005/11/22 16:20:46  cocoa
# -- added date info
#
# Revision 1.1  2005/11/22 15:56:27  cocoa
# -- first import
#
