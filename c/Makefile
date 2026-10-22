# Build the C table lookup. Set HDF5 to the HDF5 installation prefix (1.8 or later),
# e.g.  make HDF5=/opt/homebrew/opt/hdf5   or   make HDF5=/usr
HDF5   ?= /usr/local
CC     ?= cc
CFLAGS ?= -O2 -std=c99 -Wall
CPPFLAGS += -I$(HDF5)/include
LDFLAGS  += -L$(HDF5)/lib -Wl,-rpath,$(HDF5)/lib
LDLIBS   += -lhdf5 -lm

all: champ_example

champ_example: champ_example.o champ_table.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c champ_table.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

clean:
	rm -f *.o champ_example
