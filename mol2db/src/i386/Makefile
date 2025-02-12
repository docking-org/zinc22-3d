# mol2db linux makefile
# trivial change

ARCH = -m32
CC = gcc
CFLAGS = -g $(ARCH)

ALLOBJS =	add_atms.o color.o confhier.o fixcrgs.o getds.o \
		hiergen.o lookup_amb.o index.o other.o read_amber.o \
		readsolv.o rotateH.o skipmol.o term_n_s_o_h.o ui.o \
		vdwtype.o write_charges.o write_hierarchy.o

mol2db: $(ALLOBJS) mol2hier.o
	$(CC) $(CFLAGS) -o mol2db mol2hier.o $(ALLOBJS) -lm

clean:
	/bin/rm -f *.o mol2db

%.o: ../%.c ../mol2hier.h
	$(CC) -c $(CFLAGS) -o $*.o ../$*.c
