CXX = g++
CPPFLAGS = -O3
LDFLAGS = -O3

TARGET = ../bin/cal_nucl_distribution_for_cell_HS ../bin/cal_center_region_overlap
OBJECTS = operation.o cal_nucl_distribution_for_cell_HS.o Nucl_space.o

.PHONY : all
all : $(TARGET) 


../bin/cal_nucl_distribution_for_cell_HS : operation.o cal_nucl_distribution_for_cell_HS.o Nucl_space.o
	$(CXX) $(LDFLAGS) -o ../bin/cal_nucl_distribution_for_cell_HS operation.o cal_nucl_distribution_for_cell_HS.o Nucl_space.o

../bin/cal_center_region_overlap : 	cal_center_region_overlap.cpp
	$(CXX) $(LDFLAGS) -o ../bin/cal_center_region_overlap cal_center_region_overlap.cpp


operation.o :
Nucl_space.o : operation.h
cal_nucl_distribution_for_cell_HS.o : operation.h Nucl_space.h 
