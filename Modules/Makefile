# compiler
     FC = f2py
    CXX = cpp

# compile flags
# for jupyter lemar (maybe lops too ?)
#FCFLAGS = -I/opt/linux/envs/envs/phyocean-2022.05/include -c --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
# usual linux path
FCFLAGS = -I/usr/include/x86_64-linux-gnu -c --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
FCFLAGS_DATARMOR = -I/appli/gcc/gcc-9.1.0__8.3.0/include/ -c --f90flags="-g -check all -CA -CB -CS"

# source files and objects
   SRCSP = pyticles_3d_sig_sa.F
    RCSP = $(SRCSP:.F=.f)

# program name
PROGRAMP = pyticles_3d_sig_sa 


all: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	CFLAGS="--std=c99" $(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

pyticles: $(RCSP)
	$(CXX) $(SRCSP) $(RCSP)
	CFLAGS="--std=c99" $(FC) $(FCFLAGS) -m $(PROGRAMP) $(RCSP)

clean:
	rm -f *.f *.so  *.o *.mod
