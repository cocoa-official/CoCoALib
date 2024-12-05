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
