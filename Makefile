CC = gcc
CFLAGS = -O2 -Wall -I/usr/local/include/
LDFLAGS = -lgsl -lm

DEPS = mylib.h draine.h
OBJ = dust3d.o mylib.o draine.o

all: dust3d

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o $@ $<

dust3d: $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f *o