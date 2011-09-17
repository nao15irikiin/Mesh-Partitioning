INC_INPUT := ../../../airfoil/inputs
INC_PART  := ../..
INC := -I $(INC_INPUT)

all:rcb rcb_th rcb_mp

rcb:rcb_main.cpp $(INC_PART)/op_part_rcb.cpp $(INC_PART)/op_part_rcb.h $(INC_INPUT)/input.cpp $(INC_INPUT)/input.h
	g++ -Wall -g $(INC) rcb_main.cpp $(INC_PART)/op_part_rcb.cpp $(INC_INPUT)/input.cpp -o rcb

rcb_th:rcb_main.cpp $(INC_PART)/op_part_rcb_th.cpp $(INC_PART)/op_part_rcb.h $(INC_INPUT)/input.cpp $(INC_INPUT)/input.h
	g++ -Wall -g $(INC) rcb_main.cpp $(INC_PART)/op_part_rcb_th.cpp $(INC_INPUT)/input.cpp -o rcb_th -lpthread

rcb_mp:rcb_main.cpp $(INC_PART)/op_part_rcb_mp.cpp $(INC_PART)/op_part_rcb.h $(INC_INPUT)/input.cpp $(INC_INPUT)/input.h
	g++ -Wall -g $(INC) rcb_main.cpp $(INC_PART)/op_part_rcb_mp.cpp $(INC_INPUT)/input.cpp -o rcb_mp -fopenmp

clean:
	rm -f *.o rcb rcb_th rcb_mp
