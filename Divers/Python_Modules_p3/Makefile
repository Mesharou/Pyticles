


# compiler
     FC = f2py
    CXX = cpp

# compile flags
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --fcompiler=intelem --compiler=intelem -fPIC --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -DF2PY_REPORT_ON_ARRAY_COPY=1 -c --f90flags="-g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f90flags="-g -check all -CA -CB -CS"
FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f77flags="-fallow-argument-mismatch" 


# source files and objects
   SRCS = R_tools_fort.F
    RCS = $(SRCS:.F=_.f)
   SRCSGU = R_tools_fort_gula.F
    RCSGU = $(SRCSGU:.F=_.f)
   SRCSGI = R_tools_fort_gigatl.F
    RCSGI = $(SRCSGI:.F=_.f)
   SRCSC = R_tools_fort_cuc.F
    RCSC = $(SRCSC:.F=_.f)
   SRCSP = pyticles_3d_sig_sa.F
    RCSP = $(SRCSP:.F=_.f)

# program name
PROGRAM   = R_tools_fort
PROGRAMGU = R_tools_fort_gula
PROGRAMGI = R_tools_fort_gigatl
PROGRAMC  = R_tools_fort_cuc
PROGRAMP  = pyticles_3d_sig_sa


all:
	$(CXX) $(SRCS) $(RCS)
	$(CXX) $(SRCSGU) $(RCSGU)
	$(CXX) $(SRCSGI) $(RCSGI)
	$(FC) $(FCFLAGS) -m $(PROGRAM) $(RCS)
	$(FC) $(FCFLAGS) -m $(PROGRAMGU) $(RCSGU)
	$(FC) $(FCFLAGS) -m $(PROGRAMGI) $(RCSGI)
#	$(CXX) $(SRCSC) $(RCSC)
#	$(FC) $(FCFLAGS) -m $(PROGRAMC) $(RCSC)
	f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m romstoolsfort_old romstoolsfort_old.F
	
gigatl:
	$(CXX) $(SRCSGI) $(RCSGI)
	$(FC) $(FCFLAGS) -m $(PROGRAMGI) $(RCSGI)
	
particles: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	$(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

clean:
	rm -f *.f *.so  *.o *.mod
