#### clear the list of suffices.
.SUFFIXES:
#### We want to specify suffix rules for .c and .o
.SUFFIXES: .c .o

VERSION = 20140305

# downloaded from CPAN
MATH_RANDOM_VERSION = 0.71

####  Installation related
PREFIX = $(HOME)
LIB_DIR = $(PREFIX)/lib/msbayes

#### compilation options
MAKE = make
CC = gcc

# for debugging
#CCOPT = -O0 -g  -Wunused
#-Wall -Werror
#LDFLAGS =  -g         # for debugging but not profiling

# for performance
CCOPT = -O3 -Wall -fomit-frame-pointer -funroll-loops
LDFLAGS =

INCLS = -I.
DEFS =

# With mac OS-X, static linking is discouraged.  So if you invoke make by
#  make STATIC_LINK=YES
# It statically link the gsl libs, but the rests are dynamically linked
# For linux, STATIC_LINK=YES will cause addtion of -static to linker.
STATIC_LINK := NO
GSL_LIB = -lgsl -lgslcblas

# with fink in OS-X, we may need to specify the location of GSL
Q_FINK := $(shell if [ -e /sw ]; then echo YES; else echo NO; fi)
ifeq (YES,$(Q_FINK))
LDFLAGS := $(LDFLAGS) -L/sw/lib
INCLS := $(INCLS) -I/sw/include
ifeq (YES,$(STATIC_LINK))
GSL_LIB := /sw/lib/libgsl.a /sw/lib/libgslcblas.a
CCOPT := $(CCOPT) -mmacosx-version-min=10.5
LDFLAGS := $(LDFLAGS) -mmacosx-version-min=10.5
endif
else
ifeq (YES,$(STATIC_LINK))
LDFLAGS := $(LDFLAGS) -static
endif
endif

CFLAGS = $(DEFS) $(INCLS) $(CCOPT)

#### A suffix rule to make .o
.c.o:; $(CC) -c -o $@ $(CFLAGS) $<

TARGET = sumstatsvector msDQH msprior msReject
PERL_SCRIPTS = acceptRej.pl convertIM.pl dataSummary.pl msbayes.pl msCombModels.pl obsSumStats.pl

LIB_EXEC_FILES = acceptRej.r make_pd2005.r loc2plot.r calmod.r
# LIBEXEC files are used by acceptRej.r, these has to be installed
# inside of ../lib/msbayes/ relative to binaries

MSPRIOR_O = msprior.o setup.o initvars.o whiteSpaces.o hashtab.o stringUtils.o
MSDQH_O = msQHarbpop.o QHarbpop-streec.o QHsubs.o QHarbpop-getpars.o rand2new.o
MSREJECT_O = rejectinC.o
SUMSTATSVECTOR = msStatsDQH.o sumStatsVector.o hashtab.o whiteSpaces.o stringUtils.o

#all: $(TARGET) Math-Random
all: $(TARGET)

sumstatsvector: $(SUMSTATSVECTOR) 
	$(CC) $(LDFLAGS) -o $@ $(SUMSTATSVECTOR) $(GSL_LIB) -lm

msDQH: $(MSDQH_O)
	$(CC) $(LDFLAGS) -o $@ $(MSDQH_O) $(GSL_LIB) -lm 

msReject: $(MSREJECT_O)
	$(CC) $(LDFLAGS) -o $@ $(MSREJECT_O) -lm

msprior: $(MSPRIOR_O)
	$(CC) $(LDFLAGS) -o $@ $(MSPRIOR_O) $(GSL_LIB) -lm

Math-Random:
	mkdir -p Math-Random
	cd Math-Random-$(MATH_RANDOM_VERSION); \
	perl Makefile.PL; make; \
	cp -r blib/lib/Math blib/arch/auto ../Math-Random
	rm -f Math-Random/Math/.exists Math-Random/auto/Math/Random/.exists

install: all
	mkdir -p $(PREFIX)/bin
	install -m 755 $(TARGET) $(PERL_SCRIPTS) $(PREFIX)/bin
	mkdir -p $(LIB_DIR)
	install -m 644 $(LIB_EXEC_FILES) $(LIB_DIR)
#	cp -rp Math-Random $(LIB_DIR)
	@echo
	@echo "##############################################################"
	@echo "Installation complete"
	@echo
	@echo "The executables ($(TARGET) $(PERL_SCRIPTS))"
	@echo "  are installed to $(PREFIX)/bin"
	@echo "R-scripts ($(LIB_EXEC_FILES))"
	@echo "  are installed to $(LIB_DIR)"
	@echo
	@echo "Make sure $(PREFIX)/bin is listed in your PATH"
	@echo

#	@echo "The binaries are copied to ../bin. Copy these to an appropriate"
#	@echo "location.  The appropriate location means some directory "
#	@echo "included in you PATH.  To check your PATH, type:"
#	@echo '   echo $$PATH '
#	@echo "If you have root access, you might want to copy them to"
#	@echo "/usr/local/bin.  If not, ~/bin may be the appropriate place."

help:
	@echo
	@echo " all             compile all binaries"
	@echo " install         install the binaries to $(PREFIX)/bin"
	@echo "                 and R-scripts to $(LIB_DIR)"
	@echo "                 To install in other location, use PREFIX, e.g.,"
	@echo "                   PREFIX=/usr/local make install"
	@echo "                 will install in /usr/local/bin and /usr/local/lib/msbayes"
	@echo " clean           Clean up everything"
	@echo " cleanObj        Clean up junks (*.o) but leave the executables"
	@echo
	@echo "Note that gsl (GNU Scientific Library) should be installed."
	@echo " make STATIC_LINK=YES will create staticly linked binary."
	@echo " In Mac OS-X, only GSL get linked statically."
	@echo

clean: cleanObj
#	cd Math-Random-$(MATH_RANDOM_VERSION); make clean
	rm -rf $(TARGET)
	rm -rf Math-Random

cleanObj:
	rm -f *.o core

distribution: clean
	cd ../../; ln -s msbayes-code msbayes-$(VERSION) && tar czvf msbayes-$(VERSION).tgz --exclude *~  --exclude .svn --dereference msbayes-$(VERSION) && rm msbayes-$(VERSION)

#### Generate dependencies (from GNU info on make)

%.d: %.c
	#$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) -I$(INC) $< \
	        #| sed '\''s/$*\\.o[ :]*/& $@/g'\'' > $@'

ROUND1 := $(SRCS:.c=.d)
ROUND2 := $(ROUND1:.l=.d)
ROUND3 := $(ROUND2:.y=.d)
#include $(ROUND3)
