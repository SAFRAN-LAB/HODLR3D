CC			=/usr/local/bin/g++-10
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I../
LDFLAGS		=-fopenmp -std=c++17 -I../
SOURCES		=./testHODLR3D.cpp
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testHODLR3D

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(KERNEL) $(HOMOG) $< -o $@

clean:
	rm a.out testHODLR3D *.o
