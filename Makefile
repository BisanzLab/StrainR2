CC = gcc
CFLAGS = -g -I./
CFLAGS += -Wall -Werror -Wno-unused-function -Wno-unused-parameter -Wcast-align
CFLAGS += -Wshadow -Wpointer-arith -Wwrite-strings -Wunreachable-code -pedantic
LDFLAGS = -lz
OBJS = hashcounter.o subcontig.o

all: subcontig hashcounter

release: CFLAGS += -O3 # release flags
release: clean all

debug: CFLAGS += -O0 -fsanitize=undefined # debug flags
debug: clean all

subcontig: subcontig.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

hashcounter: hashcounter.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c %.h kseq.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	@rm subcontig hashcounter $(OBJS) 2> /dev/null || true

test: subcontig hashcounter
	@../test/test.sh
