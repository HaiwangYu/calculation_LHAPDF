
all: main.cc
	g++ -O3 -c main.cc -I $(shell lhapdf-config --incdir) -I $(ROOTSYS)/include
	g++ -O3 main.o -o main \
  -I $(shell lhapdf-config --incdir) \
  -L $(shell lhapdf-config --libdir) -lLHAPDF \
  -Wl,--no-as-needed $(shell root-config --glibs)

clean:
	rm -rf *.o main *~
