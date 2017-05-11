# Use with GNU make.

CFLAGS = -W -Wall -g

OBJ = main.o hashtab.o

hashtab : $(OBJ)
	$(CC) $^ -o $@

main.o : main.c hashtab.h
hashtab.o : hashtab.c hashtab.h

.PHONY : clean
clean : 
	$(RM) $(OBJ) hashtab
