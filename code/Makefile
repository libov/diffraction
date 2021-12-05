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

all: plots.exe plot_R.exe

plots.o: plots.cxx
	gcc -o $(OBJ)/$@ -c $< -I$(TXMLPARSER)/inc -I inc/ $(FLAGS_OBJ)

plot_R.o: plot_R.cxx
	gcc -o $(OBJ)/$@ -c $< -I$(TXMLPARSER)/inc -I inc/ $(FLAGS_OBJ)

TVectorMeson.o: TVectorMeson.cxx TVectorMeson.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

TPlot.o: TPlot.cxx TPlot.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

TFitter.o: TFitter.cxx TFitter.h
	gcc -o $(OBJ)/$@ -c $< -I inc/ $(FLAGS_OBJ)

plots.exe: plots.o TVectorMeson.o TPlot.o TFitter.o
	gcc -o plots.exe $(OBJ)/plots.o $(OBJ)/TVectorMeson.o $(OBJ)/TPlot.o $(OBJ)/TFitter.o $(FLAGS) -lMinuit -lstdc++ -lMLP -lXMLIO -lTXMLParser -L$(TXMLPARSER)/lib

plot_R.exe: plot_R.o TVectorMeson.o TPlot.o TFitter.o
	gcc -o plot_R.exe $(OBJ)/plot_R.o $(OBJ)/TVectorMeson.o $(OBJ)/TPlot.o $(OBJ)/TFitter.o $(FLAGS) -lMinuit -lstdc++ -lMLP -lXMLIO -lTXMLParser -L$(TXMLPARSER)/lib

.PHONY: clean
clean:
	-rm -rf $(OBJ)/*.o
	-rm -rf *.exe
