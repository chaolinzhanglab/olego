include ../common.mk
#XINC = -I$(MYSQLINC)

CC = gcc -m64
#CFLAGS=     -g -Wall

O = common.o sqlNum.o memalloc.o osunix.o errabort.o hash.o linefile.o dlist.o obscure.o wildcmp.o dystring.o localmem.o pipeline.o verbose.o

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${LOWELAB_DEFS} ${HG_DEFS} ${HG_WARN} -o $@ -c $<

lib:kent.a

kent.a: $(O)
	ar rcus $@ $(O)

test:
	cd tests && ${MAKE} test

cleanlocal:
	rm -rf ${O} kent.a

clean:cleanlocal

#cd tests && ${MAKE} clean
