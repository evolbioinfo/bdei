CXXFLAGS =   -pthread -O3 -std=c++17 #-g  -fsanitize=address -O0
LIBS_NLOPT = -lnlopt_cxx
all: BDEI
	
BDEI:BDEI.o
	$(CXX) -o $@ $^ $(LIBS_NLOPT) $(CXXFLAGS)
clean:
	-rm -rf *.o BDEI
