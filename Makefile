
SOURCES = pancyclic.c shared/multicode_base.c shared/multicode_base.h\
          shared/multicode_input.c shared/multicode_input.h\
          shared/multicode_output.c shared/multicode_output.c\
          Makefile COPYRIGHT.txt LICENSE.txt README.md

MULTICODE_SHARED = shared/multicode_base.c shared/multicode_input.c\
                   shared/multicode_output.c

all: build/pancyclic

clean:
	rm -rf build
	rm -rf dist

build/pancyclic: pancyclic.c $(MULTICODE_SHARED)
	mkdir -p build
	cc -o $@ -O4 -Wall $^

sources: dist/pancyclic-sources.zip dist/pancyclic-sources.tar.gz

dist/pancyclic-sources.zip: $(SOURCES)
	mkdir -p dist
	zip dist/pancyclic-sources $(SOURCES)

dist/pancyclic-sources.tar.gz: $(SOURCES)
	mkdir -p dist
	tar czf $@ $(SOURCES)
