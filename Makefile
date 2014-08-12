
SOURCES = pancyclic_ccp.c shared/multicode_base.c shared/multicode_base.h\
          shared/multicode_input.c shared/multicode_input.h\
          shared/multicode_output.c shared/multicode_output.c\
          Makefile COPYRIGHT.txt LICENSE.txt README.md

MULTICODE_SHARED = shared/multicode_base.c shared/multicode_input.c\
                   shared/multicode_output.c

NAUTY_FILES = nauty/nauty.c nauty/nautil.c nauty/naugraph.c\
              nauty/schreier.c nauty/naurng.c

all: build/pancyclic_ccp

clean:
	rm -rf build
	rm -rf dist

build/pancyclic_ccp: pancyclic_ccp.c $(MULTICODE_SHARED) $(NAUTY_FILES)
	mkdir -p build
	cc -o $@ -O4 -Wall -DMAXN=1000 $^

sources: dist/pancyclic-sources.zip dist/pancyclic-sources.tar.gz

dist/pancyclic-sources.zip: $(SOURCES)
	mkdir -p dist
	zip dist/pancyclic-sources $(SOURCES)

dist/pancyclic-sources.tar.gz: $(SOURCES)
	mkdir -p dist
	tar czf $@ $(SOURCES)
