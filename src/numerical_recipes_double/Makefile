LIBOBJFILES = chsone.o  expdev.o  gammln.o  gcf.o    mrqcof.o  polint.o   qromb.o   qsimp2.o  rtbis.o   trapzd.o \
			  chstwo.c  fgauss.c  gammq.c   gser.o   mrqmin.o  polint1.o  qromb1.o  qtrap.o   sort.o    trapzd1.o \
			  covsrt.o  fit.o     gasdev.o  ksone.o  nrutil.o  polint2.o  qromb2.o  ran1.o    spline.o  trapzd2.o \
			  erfcc.o   gaussj.o  kstwo.o  piksrt.o  probks.o  qsimp.o    ratint.o  splint.o  zbrac.o	ran2.o

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

default: libnr.a

all: libnr.a

libnr.a: $(LIBOBJFILES)
	ar r $@ $^

clean:
	rm -f *.o libnr.a

tidy:
	rm -f *.o
