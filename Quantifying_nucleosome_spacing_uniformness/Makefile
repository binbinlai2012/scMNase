CXX = g++
CPPFLAGS = -O3
LDFLAGS = -O3

TARGET = ../bin/get_sc_nucl_nucl_distance
OBJECTS = operation.o get_sc_nucl_nucl_distance.o gene.o

.PHONY : all
all : $(TARGET) 

../bin/get_sc_nucl_nucl_distance : operation.o get_sc_nucl_nucl_distance.o gene.o
	$(CXX) $(LDFLAGS) -o ../bin/get_sc_nucl_nucl_distance operation.o get_sc_nucl_nucl_distance.o gene.o

operation.o :
gene.o : operation.h
get_sc_nucl_nucl_distance.o : operation.h gene.h
