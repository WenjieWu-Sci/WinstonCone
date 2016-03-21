cc=g++
LIBS=$(shell root-config --libs)
INCLUDE=$(shell root-config --cflags)
exe=factor
obj=factor.o
$(exe):$(obj)
	$(cc) $(INCLUDE) -o $(exe) factor.C $(LIBS)
$(obj):factor.C
	$(cc) -c $(INCLUDE) factor.C
clean:
	rm -rf *.o $(exe)

