#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std ;
// g++ -O3 -o sieveomp sieveomp.cc -fopenmp -lpthread
// OMP_NUM_THREADS=4 ./sieveomp

// Bit array to optimize the storage of booleans as bits
struct bitArray {
  vector<unsigned int> bit_array ;
  bitArray(size_t sz) {
    vector<unsigned int> tmp(sz,0) ;
    bit_array.swap(tmp) ;
  }
  bool operator[](size_t i) const {
    // i>>5 finds the index to the 32 bits that contain bit i
    // 1<<(i@0x1f) is the bit within the array
    return (bit_array[i>>5] & (1<<(i&0x1f)))!=0 ;
  }
  void setBit(size_t loc) {
    // works like above but sets the bit at location loc
    bit_array[loc>>5] |= (1<<(loc&0x1f)) ;
  }
  // zero out the bit array
  void zero() {
    for(size_t i=0;i<bit_array.size();++i)
      bit_array[i] = 0 ;
  }
} ;

//#define SAVE
int main() {
  int num_threads = 1 ;

  #pragma omp parallel
  {
    int thread_id = omp_get_thread_num() ;
    //    cout << "tid=" << thread_id << endl ;
    if(thread_id==0)
      num_threads = omp_get_num_threads() ;

  }
  cerr << "num_threads = " << num_threads << endl ;
#ifdef SAVE
  size_t N = 1024*1024*1024;
#else
  size_t N = 1024*1024*1024*16L ; // 30 bits
#endif
  // optimizations:
  // only need to odd primes (2 is only even prime)
  // Now our prime array is mapping odd numbers starting from 3:
  // primes[0] = 3
  // primes[1] = 5 ;
  // primes[2] = 7
  // primes[i] = (i+1)*2+1

  double start = omp_get_wtime() ;
  
  size_t blockSize = 1024*512 ;
  size_t hblockSize = blockSize>>1 ;
  // Now do a blocking version where we compute the sieve in blocks
  // of size blockSize.  Because we are only computing odd numbers it
  // is halve the size of blockSize.
  size_t hblockSizeb = (hblockSize+1)>>5 ;

  bitArray composite(hblockSizeb) ;
  
  vector<size_t> plist ;
  // estimate number of primes that will be found using prime counting function
  // This will improve performance of inserting primes on plist
  double num_est = 1.2*double(N)/log(double(N)) ;
  size_t est = size_t(num_est) ;
  cerr << "est = " << est << endl ;
  plist.reserve(est) ;

  // Compute the first block to get plist started
  for(size_t i=0;4*i*(i+3)+9<blockSize;++i) {
    // IF it is a prime, marke all odd factors as not prime
    if(!composite[i]) {
      // convert i to prime
      size_t p = (i+1)*2+1 ;
      // mark all odd prime factor
      for(size_t j=p+p+p;j<blockSize+2;j+= p+p) {
	size_t loc = (((j-1)>>1)-1) ;
	composite.setBit(loc) ;
      }
    }
  }
  // Now fill primes list (plist) with primes discovered from first block
  for(size_t i=0;i<hblockSize;++i)
    if(!composite[i]) {
      size_t p = (i+1)*2+1 ;
      plist.push_back(p) ;
    }
  // Note that since the first block computes all primes of size 1024*512
  // then this block is all that will be needed to compute primes to
  // 1024*1024*512*512 (or 1024*1024*1024*256).  So we can compute the all
  // of the blocks in parallel


  // Here we will process num_threads blocks at a time and will need some
  // bookkeeping arrays to process the information that comes from these blocks
  // counts gives the number of primes found by each block
  // starts gives where the blocks primes should be inserted into plist
  vector<size_t> counts(num_threads) ;
  vector<size_t> starts(num_threads) ;
  

  // Loop over blocks num_threads at a time
  for(size_t i=1;i*blockSize < N;i+=num_threads) {
    // parallelize over threads
#pragma omp parallel
    {
      // find thread id to index block
      int thread_id = omp_get_thread_num() ;
      // ii is the threads block
      size_t ii = i+thread_id ;
      // create prime marker for block (one for each thread)
      bitArray tcomposite(hblockSizeb) ;
      
      // compute the start and end of the block in the global index space
      size_t bstart = ii*blockSize+2 ;
      size_t bend = bstart+blockSize ;

      // process block if available
      if(ii*blockSize < N) {
	// process block
	for(size_t j=0;plist[j]*plist[j]<bend;++j) {
	  size_t p = plist[j] ;
	  size_t p2 = p*p ;
	  size_t k = (bstart)/p ;
	  // skip the the first odd prime within the block
	  size_t skip = p2>bstart?p2:(p*(k) + ((k&1)?0:p)) ;
	  while(skip < bstart)
	    skip += p+p ;
	  
	  // Mark factors in block
	  while(skip < bend) {
	    size_t indx = ((skip-1)>>1)-1 - ii*hblockSize ;
	    tcomposite.setBit(indx) ;
	    skip += p+p ;
	  }
	}
	// Now count the number of primes found in the block
	size_t cnt = 0  ;
	for(size_t k=0;k<hblockSize;++k)
	  if(!tcomposite[k])
	    cnt++ ;
	// set counts to allocate in plist
	counts[thread_id] = cnt ;
      } else
	counts[thread_id] = 0 ;

      // Wait for all blocks to be processed      
#pragma omp barrier
      if(thread_id == 0) {
	// Process the prefix sum operations in serial on thread 0
	size_t cntr = plist.size() ;
	// starts will be set to where each thread will be inserting its
	// primes
	for(int k=0;k<num_threads;++k) {
	  starts[k] = cntr ;
	  cntr +=counts[k] ;
	}
	// resize plist to make space to insert blocks into plist
	plist.resize(cntr) ;
      }
      // All threads wait for prefix sum step
#pragma omp barrier
      // now insert primes into plist at locations unique to each thread
      if(counts[thread_id] > 0) {
	// This thread starts inserting into plsit at cursor
	size_t cursor = starts[thread_id] ;
	for(size_t k=0;k<hblockSize;++k)
	  if(!tcomposite[k]) {
	    size_t p = (ii*hblockSize+k+1)*2+1 ;
	    plist[cursor] = p ;
	    cursor++ ;
	  }
      }
    }
  }

  // Now perform checksum
  size_t checksum = 2 ;
  size_t xcnt = 1+plist.size() ;
  for(size_t i=0;i<plist.size();++i) {
    checksum = checksum ^ plist[i]  ;
  }
  double end = omp_get_wtime() ;
  cerr << "checksum=" <<checksum << ", cnt=" << xcnt << endl ;
  cerr << "time to solve: " << (end-start) << endl ;
  
#ifdef SAVE
  cout << "2" << endl ;
  for(size_t i=0;i<plist.size();++i) {
    cout << plist[i] << endl ;
  }
#endif
}
  




  
  
