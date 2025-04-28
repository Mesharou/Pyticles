


# compiler
     FC = f2py
    CXX = cpp

# compile flags
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --fcompiler=intelem --compiler=intelem -fPIC --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -DF2PY_REPORT_ON_ARRAY_COPY=1 -c --f90flags="-g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f90flags="-g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f77flags="-fallow-argument-mismatch" 
FCFLAGS = -I/usr/include/x86_64-linux-gnu -c 

# source files and objects
   SRCS = R_tools_fort.F
    RCS = $(SRCS:.F=_.f)
   SRCSP = pyticles_3d_sig_sa.F
    RCSP = $(SRCSP:.F=.f)

# program name
PROGRAM   = R_tools_fort
PROGRAMP  = pyticles_3d_sig_sa


all:
	$(CXX) $(SRCS) $(RCS)
	$(FC) $(FCFLAGS) -m $(PROGRAM) $(RCS)
	
pyticles: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	CFLAGS="--std=c99" $(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

clean:
	rm -f *.f *.so  *.o *.mod
