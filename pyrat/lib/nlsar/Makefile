##
## Makefile: compilation script for th NL-SAR Toolbox
##
## This file is part of NL-SAR Toolbox version 0.6.
##
## Copyright Charles-Alban Deledalle (2013)
## Email charles-alban.deledalle@math.u-bordeaux1.fr
##
## This software is a computer program whose purpose is to provide a
## suite of tools to manipulate SAR images.
##
## This software is governed by the CeCILL license under French law and
## abiding by the rules of distribution of free software. You can use,
## modify and/ or redistribute the software under the terms of the CeCILL
## license as circulated by CEA, CNRS and INRIA at the following URL
## "http://www.cecill.info".
##
## As a counterpart to the access to the source code and rights to copy,
## modify and redistribute granted by the license, users are provided only
## with a limited warranty and the software's author, the holder of the
## economic rights, and the successive licensors have only limited
## liability.
##
## In this respect, the user's attention is drawn to the risks associated
## with loading, using, modifying and/or developing or reproducing the
## software by the user in light of its specific status of free software,
## that may mean that it is complicated to manipulate, and that also
## therefore means that it is reserved for developers and experienced
## professionals having in-depth computer knowledge. Users are therefore
## encouraged to load and test the software's suitability as regards their
## requirements in conditions enabling the security of their systems and/or
## data to be ensured and, more generally, to use and operate it in the
## same conditions as regards security.
##
## The fact that you are presently reading this means that you have had
## knowledge of the CeCILL license and that you accept its terms.
##
##
## Started on  Wed Jul 24 14:55:39 2013 Charles-Alban Deledalle
## Last update Wed Nov 20 17:55:50 2013 Charles-Alban Deledalle
##

MAKEFILE_RULES=Makefile.rules
MAKEFILE_DEPS=Makefile.deps

include $(MAKEFILE_RULES)

CFLAGS=$(CFLAGS_LIB) $(CFLAGS_COMMON)
LDFLAGS=$(LDFLAGS_LIB) $(LDFLAGS_COMMON)
# gprof -pg

SRC=	tools/sarerror.c tools/sarprintf.c tools/sarwaitbar.c \
	tools/matrixtools.c tools/inv_func.c \
	data/iorat.c data/iobin.c data/ioxima.c data/iosar.c \
	data/ionetpbm.c \
	data/sardata.c data/rgbdata.c data/fltdata.c \
	data/sar2rgb.c \
	algos/carfilter/sarboxcar.c \
	algos/carfilter/sardiskcar.c \
	algos/carfilter/sargausscar.c \
	algos/carfilter/fltboxcar.c \
	algos/nlsar/nlsar.c \
	algos/nlsar/sarsimstats.c \
	algos/nlsar/sarsim_glrwishart.c \
	algos/nlsar/phi.c \
	algos/noisegen/noisegen.c

DST=	$(SRC:.c=.o)

.SUFFIXES : .c .o .$(LDEXT)

TARGETS=lib/libnlsartoolbox.a lib/libnlsartoolbox.$(LDEXT)

all:	$(TARGETS)
	[ "$(DOCLI)" != "1" ]       || make -C interfaces/cli
	[ "$(DOIDL)" != "1" ]       || make -C interfaces/idl
	[ "$(DOMEX)" != "1" ]       || make -C interfaces/mex
	[ "$(DOPYTHON)" != "1" ]    || make -C interfaces/python
	[ "$(DOPOLSARPRO)" != "1" ] || make -C interfaces/polsarpro
	[ "$(DODOC)" != "1" ]       || make -C doc

-include $(MAKEFILE_DEPS)

lib/libnlsartoolbox.a: $(DST)
	mkdir -p lib
	ar cr $@ $(DST)
	ranlib $@

lib/libnlsartoolbox.$(LDEXT): $(DST)
	mkdir -p lib
	$(CC) -o $@ $(DST) -shared $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

install:
	install -m 0755 $(TARGETS) $(PREFIX)/lib/
	[ "$(DOCLI)" != "1" ]       || make -C interfaces/cli install
	[ "$(DOIDL)" != "1" ]       || make -C interfaces/idl install
	[ "$(DOMEX)" != "1" ]       || make -C interfaces/mex install
	[ "$(DOPYTHON)" != "1" ]    || make -C interfaces/python install
	[ "$(DOPOLSARPRO)" != "1" ] || make -C interfaces/polsarpro install
	make -C interfaces/lib install

depend:
	$(CC) -MM $(CFLAGS) $(SRC) > $(MAKEFILE_DEPS)
	[ "$(DOCLI)" != "1" ]       || make -C interfaces/cli depend
	[ "$(DOIDL)" != "1" ]       || make -C interfaces/idl depend
	[ "$(DOMEX)" != "1" ]       || make -C interfaces/mex depend
	[ "$(DOPYTHON)" != "1" ]    || make -C interfaces/python depend
	[ "$(DOPOLSARPRO)" != "1" ] || make -C interfaces/polsarpro depend

clean:
	rm -f  *.o *~ $(DST)
	[ "$(DOCLI)" != "1" ]       || make -C interfaces/cli clean
	[ "$(DOIDL)" != "1" ]       || make -C interfaces/idl clean
	[ "$(DOMEX)" != "1" ]       || make -C interfaces/mex clean
	[ "$(DOPYTHON)" != "1" ]    || make -C interfaces/python clean
	[ "$(DOPOLSARPRO)" != "1" ] || make -C interfaces/polsarpro clean
	[ "$(DODOC)" != "1" ]       || make -C doc clean

distclean: clean
	[ "$(DOCLI)" != "1" ]       || make -C interfaces/cli distclean
	[ "$(DOIDL)" != "1" ]       || make -C interfaces/idl distclean
	[ "$(DOMEX)" != "1" ]       || make -C interfaces/mex distclean
	[ "$(DOPYTHON)" != "1" ]    || make -C interfaces/python distclean
	[ "$(DOPOLSARPRO)" != "1" ] || make -C interfaces/polsarpro distclean
	[ "$(DODOC)" != "1" ]       || make -C doc distclean
	rm -f $(TARGETS) $(MAKEFILE_RULES) $(MAKEFILE_DEPS)
