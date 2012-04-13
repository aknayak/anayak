
files = $(wildcard *.cc *.cpp)
OBJS = $(addsuffix .o,$(basename $(files)))

#OBJS = MyEvent.o MyMCParticle.o MyGlobalPoint.o MyMeasurement1D.o MyMET.o \
#       MyTrack.o MyVertex.o MyCaloTower.o MyJet.o MyImpactParameter.o \
#       eventdict.o etaPhiFunctions.o

OPT = -O -Wall -fPIC -D_REENTRANT
INC = -I${ROOTSYS}/include -I../../..

.cc.o:
	$(CXX) $(OPT) $(INC) -c $*.cc

all:
	@$(MAKE) eventdict.cc; $(MAKE) libMyEvent.so

eventdict.cc: ../interface/MyEvent.h ../interface/MyMCParticle.h ../interface/MyTrack.h ../interface/MyVertex.h ../interface/MyJet.h ../interface/MyMET.h ../interface/MyElectron.h ../interface/MyMuon.h ../interface/MyTau.h ../interface/SampleInfo.h ../interface/Reader.h ../interface/ObjectSelector.hh ../interface/HistogramPlotter.hh LinkDef.h
	rootcint -f eventdict.cc -c -I${PWD}/../../.. -p \
	MiniTree/Selection/interface/MyEvent.h \
	MiniTree/Selection/interface/MyMET.h \
	MiniTree/Selection/interface/MyMCParticle.h \
	MiniTree/Selection/interface/MyTrack.h \
	MiniTree/Selection/interface/MyVertex.h \
	MiniTree/Selection/interface/MyJet.h \
	MiniTree/Selection/interface/MyElectron.h \
	MiniTree/Selection/interface/MyMuon.h \
	MiniTree/Selection/interface/MyTau.h \
	MiniTree/Selection/interface/SampleInfo.h \
	MiniTree/Selection/interface/Reader.h \
	MiniTree/Selection/interface/ObjectSelector.hh \
	MiniTree/Selection/interface/HistogramPlotter.hh \
	MiniTree/Selection/src/LinkDef.h
	
compile: $(OBJS)

libMyEvent.so: $(OBJS)
	$(CXX) -shared -O *.o -o libMyEvent.so

clean:
	rm -f *.o; rm -f libMyEvent.so; rm -f eventdict.*
