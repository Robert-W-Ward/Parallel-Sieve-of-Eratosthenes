all: sieveo3 sievempi
#all: sieve sieveo1 sieveo2 sieveo3 sieveomp


#sieve:
#g++ -O3 -o sieve sieve.cc

#sieveo1:
#g++ -O3 -o sieveo1 sieveo1.cc

#sieveo2:
#g++ -O3 -o sieveo2 sieveo2.cc

sieveo3:
#mpiCC -O0 -o sieveo3_mpi sieveo3_mpi.cc 
	mpic++ -O3 -o sieveo3 sieveo3.cc

#sieveomp:
#g++ -O3 -o sieveomp sieveomp.cc -fopenmp -lpthread

sievempi:
	mpic++ -O0 -o sievempi sievempi.cc 
clean:
#rm -f sieve sieveo1 sieveo2 sieveo3 sieveomp sievempi
	rm -f sieveo3 sievempi