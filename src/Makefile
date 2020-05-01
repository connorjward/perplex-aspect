PERPLEX_DIR=/home/connor/Documents/uni/project/code/perplex

test : meemum_wrapper.o test.o
	g++ -o test meemum_wrapper.o test.o \
		${PERPLEX_DIR}/resub.o ${PERPLEX_DIR}/*lib.o \
		-lm -lgfortran -g

meemum_wrapper.o : meemum_wrapper.f
	gfortran -c meemum_wrapper.f -I${PERPLEX_DIR} -g

test.o : test.cpp
	g++ -c test.cpp -g

clean :
	rm test meemum_wrapper.o test.o

