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

tidy:
	make -C src tidy
