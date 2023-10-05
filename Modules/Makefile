# Pyticles and R_tools compilation receipes are within the same Makefile
# but **should not be compiled using the same Python environment.**
# Pyticles core is very sensitive to numpy version which is incompatible 
# with a changing postprocessing environment.
#
# compiler
     FC = f2py
    CXX = cpp

# compile flags
# for jupyter lemar (maybe lops too ?)
#FCFLAGS = -I/opt/linux/envs/envs/phyocean-2022.05/include -c --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
# usual linux path
#FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
FCFLAGS_DATARMOR = -I/appli/gcc/gcc-9.1.0__8.3.0/include/ -c --f90flags="-g -check all -CA -CB -CS"

# source files and objects
SRCSP = pyticles_3d_sig_sa.F
RCSP = $(SRCSP:.F=.f)
SRCS = R_tools_fort.F
RCS = $(SRCS:.F=.f)

# program name
PROGRAMP = pyticles_3d_sig_sa 

.PHONY: pass

pass:
	@echo "skip build"

# --
pyticles: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	CFLAGS="--std=c99" $(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

R_tools: $(RCS)
	$(CXX) $(SRCS) $(RCS)
	$(FC) $(FCFLAGS) -m $(PROGRAM) $(RCS)

# --
pyticles_datarmor: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	CFLAGS="--std=c99" $(FC) $(FCFLAGS_DATARMOR) -m $(PROGRAMP) $(RCSP)

R_tools_datarmor: $(RCS)
	$(CXX) $(SRCS) $(RCS)
	$(FC) $(FCFLAGS_DATARMOR) -m $(PROGRAM) $(RCS)

# --
clean:
	rm -f *.f *.so  *.o *.mod
