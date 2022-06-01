#include <iostream>
#include <vector>
using namespace std ;
//#define SAVE
int main() {

#ifdef SAVE
  size_t N = 1024*1024*1024 ;
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
  // Size of our prime test array is halved in size
  size_t Nby2 = N/2 ;
  vector<bool> primes(Nby2,true) ;

  // Note transformation from i to candidate value is
  // p = (i+1)*2+1, no need to sieve for p*p > N
  for(size_t i=0;4*i*(i+3)+9<N;++i) {
    // If it is a prime, mark all odd factors as not prime
    if(primes[i]) {
      // covert i to prime
      size_t p = (i+1)*2+1 ;
      // mark all odd prime factor
      for(size_t j=p*p;j<N;j+= p+p)
	primes[((j-1)>>1)-1] = false ;
    }
  }
  // Form checksum of found primes.



  size_t checksum = 2 ;
  size_t cnt = 1 ;
  for(size_t i=0;(i+1)*2+1<N;++i) {
    if(primes[i]) {
      size_t p  = (i+1)*2+1  ;
      checksum = checksum ^ p ;
      cnt++ ;
    }
  }
  cerr << "checksum=" <<checksum << ", cnt=" << cnt << endl ;
#ifdef SAVE
  cout << "2" << endl ;
  for(size_t i=0;(i+1)*2+1<N;++i) {
    if(primes[i])
      cout << (i+1)*2+1 << endl ;
  }
#endif
}
  
  
  
