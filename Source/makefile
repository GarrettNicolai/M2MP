GNUREC=-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time -pthread
GO=$(GNUREC)
CC=g++ $(GO) 

STLPORT=($HOME)/DTL/STLport-5.2.1  ## edit here

INCLUDES=-I./tclap-1.2.1/include/ -I/DTL/STLport-5.2.1/include/stlport
LIBS=-L/DTL/STLport-5.2.1/lib
CFLAGS=-c $(INCLUDES) -USESTLPORT
LDFLAGS=$(LIBS) 
INLIBS= -lgcc_s -lpthread -lc -lm

alignSrc=mmAligner.cpp mmEM.cpp
SOURCES=$(alignSrc)
alignObj=$(alignSrc:.cpp=.o)
OBJECTS=$(alignObj)

EXECUTABLE=M2M+

all: $($SOURCES) $(EXECUTABLE)

$(EXECUTABLE):	$(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -g -o $@ $(INLIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:	
	rm -f $(EXECUTABLE) $(OBJECTS)
