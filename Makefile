prog := qdmd
prog_objs := main.o md.o util.o tersoff1.o \
			 tersoff2.o tersoff2_forces.o tersoff2_params.o

CC := gcc
CFLAGS := -Wall -pedantic -g -O2 
#CFLAGS := -Wall -Wextra -pedantic -g -O2 
#CFLAGS := -Wall -Wtraditional-conversion -pedantic -g -O2 
#CFLAGS := -Wall -Wconversion -pedantic -g -O2 
LDFLAGS := -lm

.PHONY: all clean sim

all: $(prog) sim

sim: $(prog)
	cp -f $(prog) bin && chmod +x bin/$(prog)
	ctags *

$(prog): $(prog_objs)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

main.o: main.c
md.o: md.c md.h
tersoff1.o: tersoff1.c tersoff1.h
tersoff2.o: tersoff2.c tersoff2.h
tersoff2_forces.o: tersoff2_forces.c tersoff2_forces.h
tersoff2_params.o: tersoff2_params.c tersoff2_params.h
util.o: util.c util.h

clean:
	@rm -rf *.o $(prog)
