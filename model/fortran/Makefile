FFLAGS := -fcheck=bounds -ffree-line-length-0
F_SOURCES := $(wildcard *.F95)
F_EXECUTABLE := $(F_SOURCES:.F95=.exe)

all:$(F_EXECUTABLE)

$(F_EXECUTABLE):$(F_SOURCES)
	gfortran $< $(FFLAGS) -o ../compiled/$@

