#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <stdio.h>
#include <mpi.h>
#include <fstream>
using namespace std;

const size_t N = (1024 * 1024 * 1024)*16L;
const size_t blockSize = (1024 * 512);
const int totalBlocks = (N / blockSize);
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

// A utility routine for measuring time
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

void processFirstBlock(size_t blockSize, bitArray &composite)
{
  // Compute the first block to get plist started
  for (size_t i = 0; 4 * i * (i + 3) + 9 < blockSize; ++i)
  {
    // IF it is a prime, mark all odd factors as not prime
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
}

// process block
void processBlock(vector<size_t> &plist, bitArray &composite, size_t blockSize, size_t hblockSize, size_t curBlock)
{
  // reset primes
  composite.zero();
  // compute the start and end of the block in the global index space
  size_t bstart = curBlock * blockSize + 2;
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
      size_t indx = ((skip - 1) >> 1) - 1 - curBlock * hblockSize;
      composite.setBit(indx);
      skip += p + p;
    }
  }
  // Enter primes found into the list of primes
  for (size_t k = 0; k < hblockSize; ++k)
    if (!composite[k])
    {
      size_t p = (curBlock * hblockSize + k + 1) * 2 + 1;
      plist.push_back(p);
    }
}

// partitions the blocks out so work is done evenly across processors
// weighted as higher numbers in the list of primes will be harder to compute
float partition(int processorCnt, int id)
{
  if (processorCnt == 2 && id == 0)
    return .75f;
  else if (processorCnt == 2 && id == 1)
    return .25f;

  if (id == processorCnt - 1)
    return partition(processorCnt, id - 1);
  if (id != 0)
    return (.5f) * partition(processorCnt, id - 1);
  else
    return .5f;
}

int getStartingBlockidx(int numBlocks, int id,int numProcessors)
{
  if (id == 0)
    return 1;
  else
    return (id * numBlocks)+(totalBlocks%numProcessors);
}
//#define SAVE
int main()
{
  //#ifdef SAVE
  // size_t N = 1024 * 1024 * 1024;
  // size_t blockSize = 1024*512;
  //#else
  // size_t N = 1024 * 1024 * 1024 * 16L;
  //#endif
  double start = get_timer();
  sleep(10);
  MPI_Init(NULL, NULL);

  int myId;
  int numProcessors;


  std::ofstream myfile;
  myfile.open("out.txt");
  MPI_Comm_size(MPI_COMM_WORLD, &numProcessors);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);

  // optimizations:
  // only need to odd primes (2 is only even prime)
  // Now our prime array is mapping odd numbers starting from 3:
  // primes[0] = 3
  // primes[1] = 5 ;
  // primes[2] = 7
  // primes[i] = (i+1)*2+1

  // size_t blockSize = 1024 * 512;
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

  processFirstBlock(blockSize, composite);
  // Now fill primes list (plist) with primes discovered from first block
  for (size_t i = 0; i < hblockSize; ++i)
    if (!composite[i])
    {
      size_t p = (i + 1) * 2 + 1;
      plist.push_back(p);
    }
  int firstBlockSize = plist.size();
  // now loop over remaining blocks to find remaining primes without
  // allocating a huge primes array.
  // measure here for each segment
  int numBlocks;
  int StartingBlockidx;

  if (totalBlocks % numProcessors == 0)
    numBlocks = (totalBlocks / numProcessors);
  else if (myId == 0)
  {
    int odd = totalBlocks % numProcessors;
    numBlocks = int(totalBlocks / numProcessors) + odd;
  }
  else
  {
    numBlocks = int(totalBlocks / numProcessors);
  }
  StartingBlockidx = getStartingBlockidx(numBlocks, myId,numProcessors);
  if (myId == 0)
    numBlocks -= 1;
  // for (size_t i = 1; i * blockSize <N; ++i)
  //{
  for (size_t i = StartingBlockidx; i < ( StartingBlockidx+numBlocks); ++i)
  {
    processBlock(plist, composite, blockSize, hblockSize, i);
  }
  std::cout<<*plist.end();
  // Compute checksum
  size_t localCheckSum = 2;
  size_t globalCheckSum = 0;
  size_t localCount = plist.size();
  size_t globalCount = 0;
  size_t partialCheckSums[numProcessors];
  for(int i =0 ;i<numProcessors;++i){
    partialCheckSums[i]=0;
  }
  if (myId == 0)
  {
    localCheckSum = 2;
    localCount += 1;
    for (size_t i = 0; i < plist.size(); ++i)
    {
      localCheckSum = localCheckSum ^ plist[i];
    }
    cerr<<"Local CheckSum: "<<localCheckSum<<endl;
  }
  else
  {
    localCheckSum = 0;
    localCount = plist.size() - firstBlockSize;
    for (size_t i = firstBlockSize; i < plist.size(); ++i)
    {
      localCheckSum = localCheckSum^plist[i];
    }
    cerr<<"Local CheckSum: "<<localCheckSum<<endl;

  }
  MPI_Allreduce(&localCount, &globalCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /*for (size_t i = 0; i < plist.size(); ++i)
  {
    checksum = checksum ^ plist[i];
  }*/
  MPI_Allgather(&localCheckSum,1,MPI_LONG,&partialCheckSums,1,MPI_LONG,MPI_COMM_WORLD);
  globalCheckSum = partialCheckSums[0];
  for(int i = 1;i<numProcessors;++i)
  {
    globalCheckSum = globalCheckSum ^ partialCheckSums[i];
  }
  myfile.close();

  //cerr << "checksum=" << localCheckSum << ", cnt=" << localCount << endl;
   cerr << "globalchecksum=" << globalCheckSum << ", globalCount=" << globalCount <<" ";
  //cout << "2" << endl;
  //for (size_t i = plist.size() - 10; i < plist.size() && myId == 1; ++i)
  //{
    //cout << plist[i] << endl;
  //}
  cerr <<"elapsed: "<<get_timer()<<endl;
  MPI_Finalize();
  return 0;
  // checksum for 1024*1024*1024*16 first half of blocks
  // correct checksum and prime count: checksum=11523766986, cnt=762939111 for 1024*1024*1024 * 16L
  // correct checksum and prime count: checksum=52560007, cnt=54,400,028 for 1024 *1024 *1024
}
