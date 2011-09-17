METISLIB 	= libmetis.a
METIS_PATH	= $(HOME)/Dropbox/IC_DoC/3_IndividualProject/OP2/src/Metis/metis-4.0.3
INC_INPUT 	:= ../../airfoil/inputs
INC		:= -I$(METIS_PATH)/Lib -I$(INC_INPUT)
LIB		:= -L $(METIS_PATH)
LIBS 		= -lmetis -lm

metis:metis_main.cpp ../op_part_metis.cpp ../op_part_metis.h $(INC_INPUT)/input.cpp
	g++ -Wall -g $(INC) $(LIB) metis_main.cpp ../op_part_metis.cpp $(INC_INPUT)/input.cpp -o metis $(LIBS)

clean:
	rm -f *.o metis
