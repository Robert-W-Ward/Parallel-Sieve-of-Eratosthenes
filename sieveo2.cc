#include <iostream>
#include <vector>
using namespace std ;
//#define SAVE
int main() {

#ifdef SAVE
  size_t N = 1024*1024 ;
#else
  size_t N = 1024*1024*1024*16L ;
#endif
  // optimizations:
  // only need to odd primes (2 is only even prime)
  // Now our prime array is mapping odd numbers starting from 3:
  // primes[0] = 3
  // primes[1] = 5 ;
  // primes[2] = 7
  // primes[i] = (i+1)*2+1

  size_t blockSize = 1024*512 ;
  size_t hblockSize = blockSize>>1 ;

  // Now do a blocking version where we compute the sieve in blocks
  // of size blockSize.  Because we are only computing odd numbers it
  // is halve the size of blockSize.
  vector<bool> primes(hblockSize,true) ;

  // Now we are going to explicitly store primes (except 2)
  vector<size_t> plist ;

  // Compute the first block to get plist started
  for(size_t i=0;4*i*(i+3)+9<blockSize;++i) {
    // IF it is a prime, mark all odd factors as not prime
    if(primes[i]) {
      // convert i to prime
      size_t p = (i+1)*2+1 ;
      // mark all odd prime factor
      for(size_t j=p*p;j<blockSize+2;j+= p+p)
	primes[((j-1)/2)-1] = false ;
    }
  }
  // Now fill primes list (plist) with primes discovered from first block
  for(size_t i=0;i<hblockSize;++i)
    if(primes[i]) {
      size_t p = (i+1)*2+1 ;
      plist.push_back(p) ;
    }

  // now loop over remaining blocks to find remaining primes without
  // allocating a huge primes array.
  for(size_t i=1;i*blockSize < N;++i) {
    // reset primes
    for(size_t j=0;j<hblockSize;++j)
      primes[j] = true ;
    // compute the start and end of the block in the global index space
    size_t bstart = i*blockSize+2 ;
    size_t bend = bstart+blockSize ;
    // process block
    for(size_t j=0;plist[j]*plist[j]<=bend;++j) {
      size_t p = plist[j] ;
      size_t p2 = p*p ; 
      size_t k = bstart/p ;
      size_t skip = p2>bstart?p2:(p*k+((k&1)?0:p)) ;
      while(skip < bstart) {
	skip += p+p ;
      }

      // Mark factors in block
      while(skip < bend) {
	size_t indx = ((skip-1)/2)-1 - i*hblockSize ;
	primes[indx] = false ;
	skip += p+p ;
      }
    }
    
    // Enter primes found into the list of primes
    for(size_t k=0;k<hblockSize;++k)
      if(primes[k]) {
	size_t p = (i*hblockSize+k+1)*2+1 ;
	plist.push_back(p) ;
      }
  }
  // Compute checksum
  size_t checksum = 2 ;
  size_t cnt = 1+plist.size() ;
  for(size_t i=0;i<plist.size();++i) {
    checksum = checksum ^ plist[i]  ;
  }
  cerr << "checksum=" <<checksum << ", cnt=" << cnt << endl ;
  
#ifdef SAVE
  cout << "2" << endl ;
  for(size_t i=0;i<plist.size();++i) {
    cout << plist[i] << endl ;
  }
#endif
}
  
  
  
