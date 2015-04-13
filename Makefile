default: main

all: main test genpop

main:
	make -C src main

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
