SHELL = /bin/sh

CC = /usr/local/bin/gcc-5 -std=c99 -fPIC
CFLAGS = -shared -dynamiclib -fopenmp -Wall 
INCLUDES = -I. -I${SEISMIC_PATH}/src/API/c -I/usr/local/include
LIBS =  -lm -lfftw3f -L/usr/local/lib

all: shotwem.dylib shotewem.dylib wavesep.dylib compute_angles.dylib shotillum.dylib shotwem_suboff.dylib shotewem_suboff.dylib

shotwem.dylib:shotwem.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
		${CC} ${INCLUDES} -o $@ ${CFLAGS} shotwem.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

shotewem.dylib:shotewem.c ewem.c ${SEISMIC_PATH}/src/API/c/seismic.c
		${CC} ${INCLUDES} -o $@ ${CFLAGS} shotewem.c ewem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

wavesep.dylib:wavesep.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
		${CC} ${INCLUDES} -o $@ ${CFLAGS} wavesep.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

compute_angles.dylib:compute_angles.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
		${CC} ${INCLUDES} -o $@ ${CFLAGS} compute_angles.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

shotillum.dylib:shotillum.c illum.c ${SEISMIC_PATH}/src/API/c/seismic.c
		${CC} ${INCLUDES} -o $@ ${CFLAGS} shotillum.c illum.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

shotwem_suboff.dylib:shotwem_suboff.c wem_suboff.c ${SEISMIC_PATH}/src/API/c/seismic.c
		${CC} ${INCLUDES} -o $@ ${CFLAGS} shotwem_suboff.c wem_suboff.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

shotewem_suboff.dylib:shotewem_suboff.c ewem_suboff.c ${SEISMIC_PATH}/src/API/c/seismic.c
		${CC} ${INCLUDES} -o $@ ${CFLAGS} shotewem_suboff.c ewem_suboff.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

clean:
		-rm -f *.o *.so *.dylib
