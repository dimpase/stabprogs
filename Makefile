CC=gcc
CFLAGS=-Os

all: stabil1 stabil2 stabcol

stabcol: germ1.c 
	$(CC) $(CFLAGS) -o stabcol germ1.c 

stabil1: rus1.c
	$(CC) $(CFLAGS) -o stabil1 rus1.c

wl: wl.c wl.h
	$(CC) $(CFLAGS) -fPIC -shared wl.c -o libwl.so

stabil2: wl rus2.c wl.h
	$(CC) $(CFLAGS) -o stabil2 rus2.c -L. -lwl

test: all
	-./stabcol input1; ./stabil1 input1; ./stabil2 input1

clean:
	-rm stabil1 stabil2 stabcol *.o *.out *.exe

