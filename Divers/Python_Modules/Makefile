


# compiler
     FC = f2py
    CXX = cpp

# compile flags
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -DF2PY_REPORT_ON_ARRAY_COPY=1 -c --fcompiler=intelem --compiler=intelem --f90flags="-fPIC -g -check all -CA -CB -CS"
#FCFLAGS = -c --fcompiler=intelem --compiler=intelem  --f90flags="-g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --fcompiler=intelem --compiler=intelem --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu --debug-capi -c --f90flags="-g -check all -CA -CB -CS"
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -DF2PY_REPORT_ON_ARRAY_COPY=1 -c --f90flags="-g -check all -CA -CB -CS"
FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f90flags="-fPIC -g -check all -CA -CB -CS"

#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --fcompiler=intelem --compiler=intelem --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"


# source files and objects
   SRCS = R_tools_fort.F
    RCS = $(SRCS:.F=.f)
   SRCSG = R_tools_fort_gula.F
    RCSG = $(SRCSG:.F=.f)
   SRCSC = R_tools_fort_cuc.F
    RCSC = $(SRCSC:.F=.f)
   SRCSP = pyticles_3d_sig_sa.F
    RCSP = $(SRCSP:.F=.f)

# program name
PROGRAM = R_tools_fort
PROGRAMG = R_tools_fort_gula
PROGRAMC = R_tools_fort_cuc
PROGRAMP = pyticles_3d_sig_sa


all: $(RCS)
	$(CXX) $(SRCS) $(RCS)
	$(CXX) $(SRCSG) $(RCSG)
	$(CXX) $(SRCSC) $(RCSC)
	$(FC) $(FCFLAGS) -m $(PROGRAM) $(RCS)
	$(FC) $(FCFLAGS) -m $(PROGRAMG) $(RCSG)
	$(FC) $(FCFLAGS) -m $(PROGRAMC) $(RCSC)

particles: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	$(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

clean:
	rm -f *.f *.so  *.o *.mod
