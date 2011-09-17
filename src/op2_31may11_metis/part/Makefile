
METISLIB 	= libmetis.a
METIS_PATH	= ../../Metis/metis-4.0.3
INC		:= -I $(METIS_PATH)/Lib
LIB		:= -L $(METIS_PATH)
LIBS 		= -lmetis -lm

#
# make all by default
#

all:		op_part.o op_part_metis.o op_part_lib.a


#
# C++ library used by sequential and OpenMP executables
#

op_part.o:	op_part.h op_part.cpp
		g++ -g -c op_part.cpp

op_part_metis.o:	op_part_metis.cpp op_part_metis.h op_part.h
		g++ -g -c  $(INC) $(LIB) op_part_metis.cpp $(LIBS)

op_part_rcb.o:op_part_rcb.cpp op_part_rcb.h op_part.h
		g++ -g -c op_part_rcb.cpp

op_part_lib.a:	op_part_metis.o op_part_rcb.o op_part.o
		ar rcs op_part_lib.a op_part_metis.o op_part_rcb.o op_part.o

#
# clean up
#

clean:
		rm -f op_part_lib.a op_part.o op_part_metis.o op_part_rcb.o
