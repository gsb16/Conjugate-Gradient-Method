CFLAGS = -Wall -O3 -mavx -march=native -DLIKWID_PERFMON
CPPFLAGS = -DDEBUG -g
LDLIBS = -lm -llikwid

all: cgSolver.c utils.c utils.h
	gcc utils.c cgSolver.c $(CFLAGS) $(LDLIBS) -o cgSolver

doc:
	doxygen config
	ln -s doc/html/index.html doc.html

debug: cgSolver.c utils.c utils.h
	gcc utils.c cgSolver.c $(CPPFLAGS) $(CFLAGS) $(LDLIBS) -o cgSolver

clean:
	rm -rf *.o *.bak doc/ doc.html cgSolver
	