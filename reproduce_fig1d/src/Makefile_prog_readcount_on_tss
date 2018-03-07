CXX = g++
CPPFLAGS = -O3
LDFLAGS = -O3

TARGET = ../extprog_readcount_on_tss
OBJECTS = enhancer.o gene.o histonmark_dynamic.o Histonmark.o main_tf_dynamics.o operation.o readfile.o ext_readcount_on_tss.o TF.o tfbinding_dynamic.o pet.o 


.PHONY : all
all : $(TARGET) 

$(TARGET) : $(OBJECTS)
	$(CXX) $(LDFLAGS) -o ../extprog_readcount_on_tss $(OBJECTS)

enhancer.o : operation.h
gene.o : 
pet.o : operation.h
histonmark_dynamic.o : databank.h operation.h
main_tf_dynamics.o : readfile.h tfbinding_dynamic.h histonmark_dynamic.h
operation.o : 
readfile.o : operation.h databank.h
ext_readcount_on_tss.o : enhancer.h gene.h histonmark_dynamic.h Histonmark.h main_tf_dynamics.h operation.h readfile.h TF.h tfbinding_dynamic.h
TF.o : operation.h
tfbinding_dynamic.o : databank.h operation.h

