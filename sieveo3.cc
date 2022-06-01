#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>
using namespace std;

// Bit array to optimize the storage of booleans as bits
struct bitArray
{
  vector<unsigned int> bit_array;
  bitArray(size_t sz)
  {
    vector<unsigned int> tmp(sz, 0);
    bit_array.swap(tmp);
  }
  bool operator[](size_t i) const
  {
    // i>>5 finds the index to the 32 bits that contain bit i
    // 1<<(i@0x1f) is the bit within the array
    return (bit_array[i >> 5] & (1 << (i & 0x1f))) != 0;
  }
  void setBit(size_t loc)
  {
    // works like above but sets the bit at location loc
    bit_array[loc >> 5] |= (1 << (loc & 0x1f));
  }
  // zero out the bit array
  void zero()
  {
    for (size_t i = 0; i < bit_array.size(); ++i)
      bit_array[i] = 0;
  }
};

// from Dr. Ed Luke untilies.h from previous project
double get_timer()
{
  static double to = 0;
  double tn, t;
  tn = MPI_Wtime();
  t = tn - to;
  to = tn;
  return t;
}

#define SAVE
int main()
{
get_timer();
#ifdef SAVE
  size_t N = 1024 * 1024 * 1024;
#else
  size_t N = 1024 * 1024 * 1024 * 16L;
#endif
  // optimizations:
  // only need to odd primes (2 is only even prime)
  // Now our prime array is mapping odd numbers starting from 3:
  // primes[0] = 3
  // primes[1] = 5 ;
  // primes[2] = 7
  // primes[i] = (i+1)*2+1

  size_t blockSize = 1024 * 512;
  size_t hblockSize = blockSize >> 1;
  // Now do a blocking version where we compute the sieve in blocks
  // of size blockSize.  Because we are only computing odd numbers it
  // is halve the size of blockSize.
  size_t hblockSizeb = (hblockSize + 1) >> 5;

  // Use bitarray to compress storage from vector of bools
  bitArray composite(hblockSizeb);

  vector<size_t> plist;
  // estimate number of primes that will be found using prime counting function
  // This will improve performance of inserting primes on plist
  double num_est = 1.2 * double(N) / log(double(N));
  size_t est = size_t(num_est);
  plist.reserve(est);

  // Compute the first block to get plist started
  for (size_t i = 0; 4 * i * (i + 3) + 9 < blockSize; ++i)
  {
    // IF it is a prime, marke all odd factors as not prime
    if (!composite[i])
    {
      // convert i to prime
      size_t p = (i + 1) * 2 + 1;
      // mark all odd prime factor
      for (size_t j = p + p + p; j < blockSize + 2; j += p + p)
      {
        size_t loc = (((j - 1) >> 1) - 1);
        composite.setBit(loc);
      }
    }
  }
  // Now fill primes list (plist) with primes discovered from first block
  for (size_t i = 0; i < hblockSize; ++i)
    if (!composite[i])
    {
      size_t p = (i + 1) * 2 + 1;
      plist.push_back(p);
    }
  // now loop over remaining blocks to find remaining primes without
  // allocating a huge primes array.
  for (size_t i = 1; i * blockSize < N; ++i)
  {
    // reset primes
    composite.zero();

    // compute the start and end of the block in the global index space
    size_t bstart = i * blockSize + 2;
    size_t bend = bstart + blockSize;

    // process block
    for (size_t j = 0; plist[j] * plist[j] < bend; ++j)
    {
      size_t p = plist[j];
      size_t p2 = p * p;
      size_t k = (bstart) / p;
      // skip the the first odd prime within the block
      size_t skip = p2 > bstart ? p2 : (p * (k) + ((k & 1) ? 0 : p));
      while (skip < bstart)
        skip += p + p;

      // Mark factors in block
      while (skip < bend)
      {
        size_t indx = ((skip - 1) >> 1) - 1 - i * hblockSize;
        composite.setBit(indx);
        skip += p + p;
      }
    }
    // Enter primes found into the list of primes
    for (size_t k = 0; k < hblockSize; ++k)
      if (!composite[k])
      {
        size_t p = (i * hblockSize + k + 1) * 2 + 1;
        plist.push_back(p);
      }
  }
  // Compute checksum
  size_t checksum = 2;
  size_t cnt = 1 + plist.size();
  for (size_t i = 0; i < plist.size(); ++i)
  {
    checksum = checksum ^ plist[i];
  }
  cerr << "checksum=" << checksum << ", cnt=" << cnt << " ";
  cerr << get_timer();
#ifdef SAVE
  //cout << "2" << endl;
  //for (size_t i = 0; i < plist.size(); ++i)
  //{
    //cout << plist[i] << endl;
  //}
#endif
}
