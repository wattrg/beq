# aeolus
#Compiler and Linker
NVCC        := nvcc
CC          := g++
flavour     ?= debug
GPU         ?= 1

#The Target Binary Program
TARGET      := beq

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR     := src
INCDIR     := include
BUILDDIR   := build/obj
TARGETDIR  := build/target
INSTALLDIR := inst
LIBINSTDIR := $(LIBINSTDIR)/lib
LIBDIR     := build/lib
RESOURCES  := resources
RESDIR     := res
SRCEXT     := cu
DEPEXT     := d
OBJEXT     := o

#Flags, Libraries and Includes
CFLAGS      := -Wextra -Wall -fPIC -fopenmp -std=c++17
LIB         := -lm
INC         := -I$(INCDIR) -I/usr/local/include 
INCDEP      := -I$(INCDIR)


# optimisation
ifeq ($(flavour), debug)
	TARGET = beq_debug
	CFLAGS := $(CFLAGS) -g
else 
	CFLAGS := $(CFLAGS) -O3 -DNDEBUG
endif

# profiling
ifeq ($(profile), 1)
	ifneq ($(flavour), debug)
		CFLAGS := $(CFLAGS) -g
	endif
	CFLAGS := $(CFLAGS) -pg -lineinfo
endif


# clang needs some extra flags to make openmp work
ifeq ($(CC), clang++)
	CFLAGS "= $(CFLAGS) -fopenmp=libomp
endif


GIT_HASH := $(shell git describe --always --dirty)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
COMPILE_TIME := $(shell date +'%Y-%m-%d %H:%M:%S AEST')
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

SOURCES     := $(shell find $(SRCDIR) -path src/python_api -prune -o -type f -name *.$(SRCEXT) -print)
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
# LIBSOURCES  := $(shell find $(SRCDIR)/python_api -type f -name *.$(SRCEXT))
# PYBIND11    := $(shell python -m pybind11 --includes)

CUDA_CFLAGS := $(foreach option, $(CFLAGS), --compiler-options $(option))
CUDA_VERSION_FLAGS:= --compiler-options -DGIT_HASH="\"$(GIT_HASH)\"" --compiler-options -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" --compiler-options -DGIT_BRANCH="\"$(GIT_BRANCH)\""
# CFLAGS := $(CUDA_CFLAGS)
# VERSION_FLAGS := $(CUDA_VERSION_FLAGS)

#Default Make
all: directories $(TARGET) install
	@echo Finished

echo_sources:
	@echo $(SOURCES)

echo_objects:
	@echo $(OBJECTS)

echo_cflags:
	@echo $(CFLAGS)

echo_makeflags:
	@echo $(MAKEFLAGS)

echo_version_flags:
	@echo $(VERSION_FLAGS)

install: build $(TARGET) 
	@mkdir -p $(INSTALLDIR)/bin
	cp -r $(TARGETDIR)/* $(INSTALLDIR)/bin
	cp -p $(SRCDIR)/prep.py $(INSTALLDIR)/bin/beq_prep
	cp -p $(SRCDIR)/post.py $(INSTALLDIR)/bin/beq_post
	cp -r $(RESOURCES) $(INSTALLDIR)/resources
	@echo Finished installation

# Compile only
build: directories $(TARGET)
	@echo Finished build

#Remake
remake: cleaner all

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(LIBDIR)

#Clean only Objecst
clean:
	$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner:
	$(RM) -rf build
	$(RM) -rf inst

#Pull in dependency info for *existing* .o files
#-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(NVCC) $(CUDA_CFLAGS) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	$(NVCC) $(CUDA_CFLAGS) $(CUDA_VERSION_FLAGS) $(INC) -c -o $@ $< 
	$(NVCC) $(CUDA_CFLAGS) $(CUDA_VERSION_FLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

# make the dynamic libraries
# $(BUILDDIR)/lib/beq.so: $(OBJECTS)
# 	$(CC) $(CFLAGS) -shared $(VERSION_FLAGS) $(PYBIND11) $(SRCDIR)/python_api/lib.cpp -o $(LIBDIR)/beq.so $^ $(LIB)
#
# lib: directories $(BUILDDIR)/lib/beq.so
# 	@echo Finished building python library
#
# install_lib: lib
# 	mkdir -p $(LIBINSTDIR)
# 	cp -f $(LIBDIR)/* $(LIBINSTDIR)
# 	@echo Finished installing python library


#Non-File Targets
.PHONY: all remake clean cleaner resources
