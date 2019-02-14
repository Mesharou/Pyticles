


# compiler
     FC = f2py
    CXX = cpp

# compile flags
FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"

# source files and objects
   SRCSP = pyticles_3d_sig_sa.F
    RCSP = $(SRCSP:.F=.f)

# program name
PROGRAMP = pyticles_3d_sig_sa


all: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	$(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

pyticles: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	$(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

clean:
	rm -f *.f *.so  *.o *.mod
