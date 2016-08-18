# makefile for LGPL-licensed LAStools
#
#COPTS    = -g -Wall -Wno-deprecated -DDEBUG 
#COPTS     = -O3 -Wall -Wno-deprecated -DNDEBUG 
COPTS     = -O2 -w -DNDEBUG ${DEBUG_OPTS}
#COMPILER  = CC
#COMPILER  = g++
#LINKER  = g++
COMPILER ?= mpic++
LINKER ?= mpic++

#BITS     = -64

#LIBS     = -L/usr/lib64
#LIBS     = -L/usr/lib32
#INCLUDE  = -I/usr/include



LASLIBS     = -L../LASlib/lib
LASINCLUDE  = -I../LASzip/src -I../LASlib/inc 


all: p_lasmerge

p_lasmerge: lasmerge.o geoprojectionconverter.o
	${LINKER} ${BITS} ${COPTS} lasmerge.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin


.cpp.o: 
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $(LASINCLUDE) $< -o $@

.c.o: 
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $(LASINCLUDE) $< -o $@

clean:
	rm -rf *.o
	rm -rf laszip lasinfo lasprecision las2txt txt2las las2las lasdiff lasmerge lasindex p_lasmerge ../bin/p_lasmerge

clobber:
	rm -rf *.o
	rm -rf laszip lasinfo lasprecision las2txt txt2las las2las lasdiff lasmerge lasindex p_lasmerge ../bin/p_lasmerge
	rm -rf *~

