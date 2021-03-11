CXX      = clang++
CXXFLAGS = -std=c++11 -g3 -Wall -Wextra
LDFLAGS  = -g3
all: gibbs

gibbs: gibbs.o 
	${CXX} ${LDFLAGS} -o gibbs gibbs.o

gibbs.o: gibbs.cpp  

gibbsx: gibbsx.o
	${CXX} ${LDFLAGS} -o gibbsx gibbsx.o

gibbsx.o: gibbsbackup.cpp  	
