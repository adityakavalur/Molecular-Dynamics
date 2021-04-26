CC = g++-10
CFLAGS = -Wall
DEPS = constants.h loadconfig.h binning.h neighcalculation.h initializenl.h potential_calc.h velverlet.h force_numerical.h
OBJ = neighborlist.o loadconfig.o binning.o neighcalculation.o initializenl.o potential_calc.o velverlet.o force_numerical.o
TARGET = md

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(OBJ) $(CFLAGS)

all: $(LIBS) $(TARGET)

clean:
	rm $(OBJ) $(TARGET)
