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
  vector<bool> primes(N,true) ;
  for(size_t i=2;i<N;++i) {
    if(primes[i]) {
      for(size_t j=i+i;j<N;j+= i)
	primes[j] = false ;
    }
  }
  size_t checksum =  0 ;
  size_t cnt = 0 ;
  for(size_t i=2;i<N;++i) {
    if(primes[i]) {
      checksum = checksum ^ i ;
      cnt++ ;
    }
  }
  cerr << "checksum=" <<checksum << ", cnt=" << cnt << endl ;
#ifdef SAVE
  for(size_t i=2;i<N;++i) {
    if(primes[i])
      cout << i << endl ;
  }
#endif
}
  
  
  
