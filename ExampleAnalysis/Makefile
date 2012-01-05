LIBS =  -L$(ROOTSYS)/lib  -lCint -lNew \
        -lGraf -lGraf3d -lHist -lHtml -lMatrix  -lMinuit  -lTree -lTreeViewer -lTreePlayer  \
        -lPostscript -lProof -lTree -lPhysics -lCore -lGX11 -lMathCore -lMathMore -lGenVector \
        -lGpad -lGui  \
        -L/usr/X11R6/lib -lm -ldl \
        -lm -lpthread -rdynamic \
	-L/soft/lip-sw/cmssw/users/anayak/CMS/NewFWK/myAnalysis/MiniTree/Selection/src -lMyEvent 
OPTCOMP =  -g  -Wno-deprecated
CXXFLAGS = --exceptions -I. -I/$(ROOTSYS)/include -I/usr/include -I/soft/lip-sw/cmssw/users/anayak/CMS/NewFWK/myAnalysis/MiniTree/Selection/interface -I/soft/lip-sw/cmssw/users/anayak/CMS/NewFWK/myAnalysis/MiniTree/Selection/src
NOM = myAnalysis
all:    $(NOM)
$(NOM): $(NOM).o
	g++ $(OPTCOMP) $(CXXFLAGS) -o $(NOM) -g $(NOM).o $(LIBS)
 $(NOM).o: $(NOM).C
	g++ $(OPTCOMP) $(CXXFLAGS)  -DNOT_INTER -DLINUX -Df2cFortran -c $(NOM).C
clean:
	rm -f *.o
	rm -f $(NOM)
	rm -f *~;
	rm -f core;
