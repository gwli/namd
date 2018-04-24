
#TCLDIR=/Projects/namd2/tcl/tcl8.5.9-linux-x86_64
TCLDIR=/usr/include/tcl8.6
TCLINCL=-I$(TCLDIR)
#TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl
#TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl
TCLLIB=-ltcl8.6 -ldl -lpthread
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

