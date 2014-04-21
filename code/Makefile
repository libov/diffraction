# Vladyslav Libov, vladyslav.libov@desy.de

FLAGS=`root-config --glibs --cflags` -g

FLAGS_OBJ=`root-config --cflags` -Wall -g

OBJ := tmp
SRC := src
INC := inc

vpath %.o   $(OBJ)
vpath %.cxx $(SRC)
vpath %.C   $(SRC)
vpath %.h   $(INC)

all: plots.exe

plots.o: plots.cxx
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

TVectorMeson.o: TVectorMeson.cxx TVectorMeson.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

plots.exe: plots.o TVectorMeson.o
	gcc -o plots.exe $(OBJ)/plots.o $(OBJ)/TVectorMeson.o  $(FLAGS) -lMinuit -lstdc++

.PHONY: clean
clean:
	-rm -rf $(OBJ)/*.o
	-rm -rf plots.exe