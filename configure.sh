#!/bin/bash

WHICH_PSI4=`which psi5`
if [[ "$WHICH_PSI4" = "" ]]; then
    echo "Error: psi4: command not found"
    exit 1
fi


TEMPDIR=`mktemp -d`
pushd . >/dev/null
cd $TEMPDIR
# run psi4 --new-plugin-makefile in a temporary directory
# and extract the configuration variables that we need to
# put into our Makefile
psi4 --new-plugin-makefile >/dev/null
make -pn Makefile > make.db.txt 2>/dev/null
while read var assign value; do
    if [[ ${assign} = '=' ]]; then
	if [[ ${var} = 'CXX' ]]; then
	    CXX="$value"
	elif [[ ${var} = 'CXXDEFS' ]]; then
	    CXXDEFS="$value"
	elif [[ ${var} = 'CXXFLAGS' ]]; then
	    CXXFLAGS="$value"
	elif [[ ${var} = 'LDFLAGS' ]]; then
 	    LDFLAGS="$value"
	elif [[ ${var} = 'INCLUDES' ]]; then
 	    INCLUDES="$value"
	elif [[ ${var} = 'OBJDIR' ]]; then
 	    OBJDIR="$value"
	fi
    fi
done < make.db.txt
popd >/dev/null
rm -rf $TEMPDIR

MAKEFILE_CONTENTS=$(cat <<'EOF'
# The name of your plugin. Taken from the directory name.
NAME := $(shell basename `pwd`)

# C++ source files for your plugin. By default we grab all *.cc files.
CXXSRC := $(wildcard *.cc) $(wildcard */*.cc)

# Used to determine linking flags.
UNAME = $(shell uname)

# Need to link against Psi4 plugin library
PSIPLUGIN = -L$(OBJDIR)/lib -lplugin

DEPENDINCLUDE = $(notdir $(wildcard *.h*))

PSITARGET = $(NAME).so

# Start the compilation rules
default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS += -shared
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS += -shared -undefined dynamic_lookup
    CXXFLAGS += -fno-common
endif

# The object files
BINOBJ := $(notdir $(CXXSRC:%.cc=%.o))

gitversion.hpp: .git/HEAD .git/index
	echo "const char *GIT_VERSION = \"$(shell git rev-parse --short HEAD)\";" > gitversion.hpp
mosignature.o: mosignature.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<
wavekernel.o: wavekernel.cc gitversion.hpp
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<
matrixutils.o: matrixutils.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<
fermilevel.o: fermilevel.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<
cnpy.o: cnpy/cnpy.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<
brent.o: brent/brent.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<

$(PSITARGET): $(BINOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(CXXDEFS) $(PSIPLUGIN)

clean:
	rm -f $(BINOBJ) $(PSITARGET) *.d *.pyc *.test output.dat psi.timer.dat

print-%:
	@echo '$*=$($*)'
EOF
)

rm Makefile
echo "CXX = ${CXX}" >> Makefile
echo "CXXDEFS = ${CXXDEFS}" >> Makefile
echo "CXXFLAGS = ${CXXFLAGS}" >> Makefile
echo "LDFLAGS = ${LDFLAGS}" >> Makefile
echo "INCLUDES = ${INCLUDES}" >> Makefile
echo "OBJDIR = ${OBJDIR}" >> Makefile
echo "" >> Makefile
echo "" >> Makefile
echo "$MAKEFILE_CONTENTS" >> Makefile

echo "Writing Makefile"
