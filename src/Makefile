decompr-melt : decompr_melt.o meemum_wrapper.o
	g++ -o decompr-melt decompr_melt.o meemum_wrapper.o \
		${PERPLEX_DIR}/resub.o ${PERPLEX_DIR}/*lib.o \
		-lm -lgfortran -g

test : meemum_wrapper.o test.o
	g++ -o test meemum_wrapper.o test.o \
		${PERPLEX_DIR}/resub.o ${PERPLEX_DIR}/*lib.o \
		-lm -lgfortran -g

decompr_melt.o : decompr_melt.cpp
	g++ -c decompr_melt.cpp -g

meemum_wrapper.o : meemum_wrapper.f
	gfortran -c meemum_wrapper.f -I${PERPLEX_DIR} -g

test.o : test.cpp
	g++ -c test.cpp -g

clean :
	rm decompr-melt test decompr_melt.o meemum_wrapper.o test.o

