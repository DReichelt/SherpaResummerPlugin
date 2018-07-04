CC = clang++

Color_SOURCES = Color.C Main.C CBasis.C 
Color_HEADERS = Tools.H Node.H Color.H CMatrix.H Auto_Basis.H Custom_Basis.H CBasis.H

CXXFLAGS = -O2 -DUSING__Color_only

EXTRA_DIST = COPYING colortest.sh

Color_DISTDIR = Color-1.0

Color: $(addsuffix .o,$(basename $(Color_SOURCES)))

clean: 
	-rm *.o Color

dist:
	mkdir $(Color_DISTDIR) && \
	cp $(Color_SOURCES) $(Color_HEADERS) \
	   $(EXTRA_DIST) Makefile $(Color_DISTDIR) && \
	tar -czf $(Color_DISTDIR).tar.gz $(Color_DISTDIR) && \
	rm -r $(Color_DISTDIR)

all: Color
