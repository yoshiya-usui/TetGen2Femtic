CXX           = icpc
CC            = icpc
CXXFLAGS      = -O2 \
                -DNDEBUG
DEST          = ./
LIBS          = 
OBJS          = TopographicData.o MeshData.o main.o
PROGRAM       = TetGen2Femtic

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;		rm -f *.o *~ $(PROGRAM)
