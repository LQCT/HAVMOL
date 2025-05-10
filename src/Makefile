#
# Makefile for HAVMOL for DOS
#
# The following variables will be set when the Makefile is invoked with
# a specific machine type.  The variables correspond to:
# CPP     => The C preprocessor on the machine
# DEBUG   => The flag used to invoke compiling in debug mode
# EXE     => The name of the executable that will be generated
# FC      => The name of the Fortran compiler
# FFLAGS  => A list of compile options
# LDFLAGS => A list of link or load options
# EXTRAS  => 'Special' routines used on this machine
#
#CPP	= cpp -P -traditional
DEBUG	= 
EXE	= havmol.exe
FC	= gfortran
FFLAGS  = -O3 -march=pentium -ffast-math -fbounds-check \
 -fomit-frame-pointer
#FFLAGS	= -O2 -march=athlon64 -ffast-math -fbounds-check \
# -fomit-frame-pointer 
#FFLAGS  = -O6 -march=pentium -malign-double \
# -fexpensive-optimizations -ffast-math -funroll-loops -fomit-frame-pointer \
# -fno-automatic
LDFLAGS = -static 
#LDFLAGS =
MAKE	= make
EXTRAS  =
#
# This variable is set to the list of "include" files which is just SIZES.
#
INCLUDES	= dimmm 
#
# OBJS is the list of all object files that need to be linked into the
# executable.  If a new source file is added, then the corresponding
# object file name must be added to the list below.
#
OBJS= main.o main1.o start1.o start2.o timit.o init.o hfscf.o integ.o integ1.o scfrun1.o scfrun2.o scfrun3.o jkder.o jkder1.o optx.o hfgrad.o hfgrad1.o forcx.o forcx1.o
#
# For Sun arcstation 2 - version 4.1.3 of OS
#all:
#	@echo "Compiling for SPARC..."
#	@$(MAKE) FC='$(SPARCFC)' FFLAGS='$(SPARCFFLAGS)' CPP='$(SPARCCPP)' \
#	        EXE='$(SPARCEXE)' EXTRAS='$(SPARCEXTRAS)' linkit
#
# For Linux
all:
	@echo "Compiling for LINUX..."
	@$(MAKE) FC='$(FC)' FFLAGS='$(FFLAGS)' \
		EXE='$(EXE)'  linkit

#
# The default method of how to make object files from Fortran files.
#
.f.o:
#	$(CPP) $< TEMP.f
#	$(FC) -c $(FFLAGS) TEMP.f
	$(FC) -c $(FFLAGS) $<
#	ren TEMP.o $*.o
#	del TEMP.[of]
#
#
linkit: $(OBJS)
	@echo " "
	@echo "Linking $(EXE)..."
	$(FC) $(LDFLAGS) -o $(EXE) $(OBJS)
	@echo " "
	@dir $(EXE)
#
# This removes all object files when "make clean" is entered by the user.
#
clean:
	del *.o
	del *.il
#
# Next comes all of the source file names and their dependencies (if any)
# and specific compile options (if any).  If no specific compile options
# are listed, then the default specified above is used.  If a new source file
# is added to HAVMOL, then an entry should be included below
main.o:		MAIN.f	$(INCLUDES) 
main1.o:	MAIN1.f	$(INCLUDES)  
start1.o:	START1.f	$(INCLUDES)  
start2.o:	START2.f	$(INCLUDES)  
timit.o:	TIMIT.f	$(INCLUDES)  
init.o:		INIT.f	$(INCLUDES)  
hfscf.o:	HFSCF.f	$(INCLUDES)  
integ.o:	INTEG.f	$(INCLUDES)  
integ1.o:	INTEG1.f	$(INCLUDES)  
scfrun1.o:	SCFRUN1.f	$(INCLUDES) 
scfrun2.o:	SCFRUN2.f	$(INCLUDES)  
scfrun3.o:	SCFRUN3.f	$(INCLUDES)  
jkder.o:	JKDER.f	$(INCLUDES)  
jkder1.o:	JKDER1.f	$(INCLUDES)  
optx.o:		OPTX.f	$(INCLUDES)  
hfgrad.o:	HFGRAD.f	$(INCLUDES)  
hfgrad1.o:	HFGRAD1.f	$(INCLUDES)  
forcx.o:	FORCX.f	$(INCLUDES)  
forcx1.o:	FORCX1.f	$(INCLUDES) 
