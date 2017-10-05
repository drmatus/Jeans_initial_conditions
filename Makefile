all:
	gcc inicial.c ran2.c ran1.c rank.c indexx.c gasdev.c odeint.c nrutil.c rkck.c rkqs.c -o inicial -lm -g

install:
	cp inicial ~/bin
