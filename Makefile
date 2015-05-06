default:
	make -C src

all:
	make -C src all

main:
	make -C src main

detfrac:
	make -C src detfrac

rfdetfrac:
	make -C src rfdetfrac

abdetfrac:
	make -C src abdetfrac

detfracflux:
	make -C src detfracflux

detfracboth:
	make -C src detfracboth

test:
	make -C src test

genpop:
	make -C src genpop

clean:
	make -C src clean

cleanall:
	make -C src cleanall

tidy:
	make -C src tidy

tidyall:
	make -C src tidyall
