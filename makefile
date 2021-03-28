OBJ = QBits.o QFT.o QComp.o intlib.o

main: main.o shor.o $(OBJ)
	g++ main.o shor.o $(OBJ)
fig1: fig1.o $(OBJ)
	g++ fig1.o $(OBJ)