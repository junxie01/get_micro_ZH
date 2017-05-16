# call : make "name of program"
FC=gfortran
#FC=ifort
FFLAG1= -mcmodel=medium 
#-ffixed-line-length-none
objects=get_MICROZH.o clogc.o get_zh.o sacio.o zfour.o getper.o
executable=get_MICROZH
all: sacio.mod $(executable) 
%.o:%.f90
	$(FC) $^ -c
sacio.mod:sacio.f90
	$(FC) $^ -c
$(executable):$(objects)
	$(FC) $^ -o $@
install:
	cp $(executable) ../bin/$(executable)
uninstall:
	-rm ~/bin/$(excecutable)
clean:
	-rm $(objects) 
