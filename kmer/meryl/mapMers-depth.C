#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bio++.H"
#include "seqCache.H"
#include "merStream.H"
#include "libmeryl.H"
#include "existDB.H"

//  Three outputs:
//
//  1) Number of kmers that span this position.  Count of each kmer is ignored.
//  2) Count of the kmer that begins at this position.
//  3) Stats of the counts of the kmers that span this position (e.g., ave, min, max, stddev).

int
main(int argc, char **argv) {
  uint32    merSize      = 0;
  char     *merylFile    = 0L;

  char     *fastaFile    = 0L;

  bool      outputCount  = false;
  bool      outputDepth  = false;
  bool      outputStats  = false;

  int arg=1;
  int err=0;
  while (arg < argc) {

    if      (strcmp(argv[arg], "-m") == 0)
      merSize = strtouint32(argv[++arg], 0L);

    else if (strcmp(argv[arg], "-mers") == 0)
      merylFile = argv[++arg];

    else if (strcmp(argv[arg], "-seq") == 0)
      fastaFile = argv[++arg];

    else if (strcmp(argv[arg], "-count") == 0)
      outputCount = true;
    else if (strcmp(argv[arg], "-depth") == 0)
      outputDepth = true;
    else if (strcmp(argv[arg], "-stats") == 0)
      outputStats = true;

    else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (merSize == 0)
    err++;
  if (fastaFile == 0L)
    err++;
  if (merylFile == 0L)
    err++;
  if (outputCount + outputDepth + outputStats != 1)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -mers MERYL -m MERSIZE -seq IN.FASTA [-count | -depth | -stats] > output\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "For sequence ordinal 's' and position in that sequence 'p':\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -count   - report the count (c) of the single kmer that starts at position (p).\n");
    fprintf(stderr, "             Format: 's p c'\n");
    fprintf(stderr, "  -depth   - report the number (n) of kmers that span position (p).  Format: 's p n'\n");
    fprintf(stderr, "  -stats   - report the min (m), max (M), ave (a) count of all mers that span\n");
    fprintf(stderr, "             position (p).  Format: 's p m M a t n'\n");
    fprintf(stderr, "             (also reports total count (t) and number of kmers (n))\n");
    fprintf(stderr, "\n");

    if (merSize == 0)
      fprintf(stderr, "ERROR:  No mer size (-m) suppled.\n");
    if (fastaFile == 0L)
      fprintf(stderr, "ERROR:  No fasta input (-seq) suppled.\n");
    if (merylFile == 0L)
      fprintf(stderr, "ERROR:  No meryl database (-mers) suppled.\n");
    if (outputCount + outputDepth + outputStats != 1)
      fprintf(stderr, "ERROR:  Exactly one of -count, -depth and -stats may be supplied.\n");

    exit(1);
  }

  //  Open the input sequences

  seqCache      *F = new seqCache(fastaFile);

  //  Load kmer counts from a meryl database.  existDBcompressBuckets is broken.

  existDB       *E = new existDB(merylFile, merSize, existDBcounts | existDBcompressCounts, 0, UINT32_MAX);

  //  For each sequence...

  for (uint32 Sid=0; Sid < F->getNumberOfSequences(); Sid++) {
    seqInCore  *S  = F->getSequenceInCore(Sid);
    merStream  *MS = new merStream(new kMerBuilder(merSize),
                                   new seqStream(S->sequence(), S->sequenceLength()),
                                   true, true);

    //  Build a lists of the min, max and total count at each position.

    uint32  *mincount = new uint32 [S->sequenceLength() + 1];
    uint32  *maxcount = new uint32 [S->sequenceLength() + 1];
    uint32  *totcount = new uint32 [S->sequenceLength() + 1];
    uint32  *numcount = new uint32 [S->sequenceLength() + 1];

    for (uint32 xx=0; xx<S->sequenceLength() + 1; xx++) {
      mincount[xx] = UINT32_MAX;
      maxcount[xx] = 0;
      totcount[xx] = 0;
      numcount[xx] = 0;
    }

    //  Scan the sequence, find the count.

    while (MS->nextMer()) {
      uint32   pos = MS->thePositionInSequence();
      uint32   cnt = E->count(MS->theFMer()) + E->count(MS->theRMer());

      if (cnt == 0)
        //  Mer doesn't exist in the database.
        continue;

      if (outputCount)
        totcount[pos] = cnt;

      if (outputDepth)
        for (uint32 xx=pos; xx<pos+merSize; xx++)
          numcount[xx]++;

      if (outputStats)
        for (uint32 xx=pos; xx<pos+merSize; xx++) {
          totcount[xx] += cnt;
          numcount[xx]++;

          if (cnt < mincount[xx])
            mincount[xx] = cnt;

          if (maxcount[xx] < cnt)
            maxcount[xx] = cnt;
        }
    }

    //  If there is no coverage, the min is still set to UINT32_MAX.

    for (uint32 x=0; x < S->sequenceLength(); x++)
      if (numcount[x] == 0)
        mincount[x] = 0;

    //  Report the single kmer count?

    if (outputCount) {
      for (uint32 x=0; x < S->sequenceLength(); x++)
        fprintf(stdout, "%u\t%u\t%u\n", Sid, x, totcount[x]);
    }

    //  Report the depth?

    if (outputDepth) {
      for (uint32 x=0; x < S->sequenceLength(); x++)
        fprintf(stdout, "%u\t%u\t%u\n", Sid, x, numcount[x]);
    }

    //  Report the min/max/ave count?

    if (outputStats) {
      for (uint32 x=0; x < S->sequenceLength(); x++)
        fprintf(stdout, "%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
                Sid, x,
                mincount[x],
                maxcount[x],
                (numcount[x] > 0) ? totcount[x] / numcount[x] : 0,
                totcount[x],
                numcount[x]);
    }

    delete [] mincount;
    delete [] maxcount;
    delete [] totcount;
    delete [] numcount;

    delete MS;
    delete S;
  }


  delete F;
  delete E;
}
