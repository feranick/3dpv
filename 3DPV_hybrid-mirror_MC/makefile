INCLUDE = -I/usr/include
LIBDIR  = -L/usr/lib

CC=g++
LIBRARIES = -static-libgcc -static-libstdc++

CFLAGS= -c -Wall $(INCLUDE)
LDFLAGS= $(LIBDIR)

SOURCES=src/main.cpp src/Cell_3.cpp src/chromosome.cpp src/crossover.cpp src/ga.cpp src/geometry.cpp src/globalSetup.cpp src/individual.cpp src/localsearch.cpp src/nsgapopulation.cpp src/population.cpp src/random.cpp src/selection.cpp src/spa.cpp src/AABBTriangleOverlapTest.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=solar3d

all: $(SOURCES) $(EXECUTABLE)

#$(EXECUTABLE): $(OBJECTS)
#	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBDIR) $< $(LIBRARIES)
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -O2 -o $@ $(LDFLAGS) $(LIBRARIES)

.cpp.o:
	$(CC) $(CFLAGS) $< -O2 -o $@ $(LDFLAGS)

clean:
	rm src/*.o
