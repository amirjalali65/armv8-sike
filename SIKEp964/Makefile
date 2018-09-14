####  Makefile for compilation on Linux  ####

OPT=-O3     # Optimization option by default

CC=gcc
ifeq "$(CC)" "gcc"
    COMPILER=gcc
else ifeq "$(CC)" "clang"
    COMPILER=clang
endif

ARCHITECTURE=_ARM64_
USE_OPT_LEVEL=_OPTIMIZED_FAST_

AR=ar rcs
RANLIB=ranlib

CFLAGS=$(OPT) -pie -std=gnu99 -D $(ARCHITECTURE) -D __LINUX__ -D $(USE_OPT_LEVEL)
LDFLAGS=-lm
EXTRA_OBJECTS_964=objs964/fp_arm64.o objs964/fp_arm64_asm.o
OBJECTS_964=objs964/P964.o $(EXTRA_OBJECTS_964) objs/random.o objs/fips202.o

all: lib964 tests KATS

objs964/%.o: %.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs964/fp_arm64.o: ARM64/fp_arm64.c
	$(CC) -c $(CFLAGS) ARM64/fp_arm64.c -o objs964/fp_arm64.o

objs964/fp_arm64_asm.o: ARM64/fp_arm64_asm.S
	$(CC) -c $(CFLAGS) ARM64/fp_arm64_asm.S -o objs964/fp_arm64_asm.o

objs/random.o: random/random.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) random/random.c -o objs/random.o

objs/fips202.o: sha3/fips202.c
	$(CC) -c $(CFLAGS) sha3/fips202.c -o objs/fips202.o

lib964: $(OBJECTS_964)
	rm -rf sike
	mkdir sike
	$(AR) sike/libsike.a $^
	$(RANLIB) sike/libsike.a

tests: lib964
	$(CC) -pie $(CFLAGS) -L./sike tests/test_SIKEp964.c tests/test_extras.c -lsike $(LDFLAGS) -o sike/test_KEM_964

# AES
AES_OBJS=objs/aes.o objs/aes_c.o

objs/%.o: tests/aes/%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

lib964_for_KATs: $(OBJECTS_964) $(AES_OBJS)
	$(AR) sike/libsike_for_testing.a $^
	$(RANLIB) sike/libsike_for_testing.a

KATS: lib964_for_KATs
	$(CC) $(CFLAGS) -L./sike tests/PQCtestKAT_kem.c tests/rng/rng.c -lsike_for_testing $(LDFLAGS) -o sike/PQCtestKAT_kem

check: tests

.PHONY: clean

clean:
	rm -rf *.req objs964 objs sike

