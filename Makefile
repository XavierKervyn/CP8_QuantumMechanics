CC = g++
CFLAGS = -std=c++0x -Wall -g
EXEC_NAME = Exercice8_solution
INCLUDES =
LIBS =
OBJ_FILES = Exercice8_solution.o

all : $(EXEC_NAME)

clean :
	rm $(EXEC_NAME) $(OBJ_FILES) *.out*

$(EXEC_NAME) : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

