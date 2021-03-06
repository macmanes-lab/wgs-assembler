#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bio.h"

#define MAX_ESTS    (16 * 1024 * 1024)
#define MAX_HITS    (18474961)  //  for 20-03-5000-0.4

//#define SHOW_ONE

////////////////////////////////////////

struct aHit {
  uint32  _forward;
  uint32  _qsIdx;
  uint32  _dsIdx;
  uint32  _dsLo;
  uint32  _dsHi;
  uint32  _covered;
  uint32  _matched;
  uint32  _numMers;
  uint32  _yesno;
  uint32  _identity;
  uint32  _coverage;

  float    scoreCov;
  float    scoreMult;
};


void   ahit_writeBinary(aHit *a, FILE *F) {
  fwrite(a, sizeof(aHit), 1, F);
}

void   ahit_readBinary(aHit *a, FILE *F) {
  fread(a, sizeof(aHit), 1, F);
}

void   ahit_printASCII(aHit *a, FILE *F) {
#ifdef TRUE64BIT
  fprintf(F, "-%c -e %u -D %u %u %u -M %u %u %u %s %u %u\n",
          a->_forward ? 'f' : 'r',
          a->_qsIdx,
          a->_dsIdx,
          a->_dsLo,
          a->_dsHi,
          a->_covered,
          a->_matched,
          a->_numMers,
          a->_yesno ? "-Y" : "-N",
          a->_identity,
          a->_coverage);
#else
  fprintf(F, "-%c -e %lu -D %lu %lu %lu -M %lu %lu %lu %s %lu %lu\n",
          a->_forward ? 'f' : 'r',
          a->_qsIdx,
          a->_dsIdx,
          a->_dsLo,
          a->_dsHi,
          a->_covered,
          a->_matched,
          a->_numMers,
          a->_yesno ? "-Y" : "-N",
          a->_identity,
          a->_coverage);
#endif
}


void   ahit_parseString(aHit *a, char *b) {
  char *c = b+1;

  a->_forward = (*c == 'f');
  c += 1;

  if (c[2] != 'e')  fprintf(stderr, "'%s' didn't get -e\n", b);

  c += 4;
  a->_qsIdx     = (uint32)strtoul(c, &c, 10);

  if (c[2] != 'D')  fprintf(stderr, "'%s' didn't get -D\n", b);

  c += 4;
  a->_dsIdx     = (uint32)strtoul(c, &c, 10);
  a->_dsLo      = (uint32)strtoul(c, &c, 10);
  a->_dsHi      = (uint32)strtoul(c, &c, 10);

  if (c[2] == 'M') {
    c += 4;
    a->_covered   = (uint32)strtoul(c, &c, 10);
    a->_matched   = (uint32)strtoul(c, &c, 10);
    a->_numMers   = (uint32)strtoul(c, &c, 10);
  } else {
    //fprintf(stderr, "'%s' didn't get -M\n", b);
    a->_covered   = 0;
    a->_matched   = 0;
    a->_numMers   = 0;
  }

  a->_yesno    = 0;
  a->_identity = 0;
  a->_coverage = 0;

  if (c[2] == 'Y') {
    c += 4;
    a->_yesno     = 1;
    a->_identity  = (uint32)strtoul(c, &c, 10);
    a->_coverage  = (uint32)strtoul(c, &c, 10);
  }

#if 0
  if (c[2] == 'N') {
    c += 4;
    a->_yesno     = 0;
    a->_identity  = (uint32)strtoul(c, &c, 10);
    a->_coverage  = (uint32)strtoul(c, &c, 10);
  }
#endif
}

////////////////////////////////////////

int
hitCompare(const void *a, const void *b) {
  const aHit  *A = (const aHit *)a;
  const aHit  *B = (const aHit *)b;

  if (A->scoreCov > B->scoreCov)
    return(-1);
  else
    return(A->scoreCov < B->scoreCov);
}

int
hitCompareID(const void *a, const void *b) {
  const aHit  *A = (const aHit *)a;
  const aHit  *B = (const aHit *)b;

  if (A->_qsIdx < B->_qsIdx)
    return(-1);
  if (A->_qsIdx > B->_qsIdx)
    return(1);
  return(0);
}



int
main(int argc, char **argv) {
  aHit        *hits         = new aHit   [MAX_HITS];
  uint32       hitsLen      = 0;


  //  read all the hits from stdin -- assumes ascii format
  //
  char   hitLine[1025];

  while (!feof(stdin)) {
    fgets(hitLine, 1024, stdin);
    if (!feof(stdin)) {
      ahit_parseString(hits + hitsLen, hitLine);

      //  These are the scores used by the filter
      //
      hits[hitsLen].scoreCov  = (float)hits[hitsLen]._covered / (float)hits[hitsLen]._numMers;
      hits[hitsLen].scoreMult = (float)hits[hitsLen]._matched / (float)hits[hitsLen]._covered;

      //  aHit->_covered is in bases, but aHit->_numMers is the
      //  number of mers.  Possible for coverage to be > 1.0.
      //
      if (hits[hitsLen].scoreCov > 1.0)
        hits[hitsLen].scoreCov = 1.0;

      hitsLen++;

      if ((hitsLen & 0xff) == 0) {
        fprintf(stderr, "reading hits "uint32FMT"\r", hitsLen);
        fflush(stderr);
      }
    }
  }

  fprintf(stderr, "reading hits "uint32FMT"\n", hitsLen);


  //  Sort the hits by estid
  //
  fprintf(stderr, "sorting hits by cDNA\n");
  qsort(hits, hitsLen, sizeof(aHit), hitCompareID);


  //  Sort the hits by score (scoreCov), in decreasing order.
  //
  fprintf(stderr, "sorting hits by score\n");
  for (uint32 currentHit = 0; currentHit < hitsLen; ) {
    uint32 estOfInterest = hits[currentHit]._qsIdx;
    uint32 numHits       = 0;
    for (uint32 t=currentHit; (t < hitsLen) && (hits[t]._qsIdx == estOfInterest); t++)
      numHits++;

    qsort(hits + currentHit, numHits, sizeof(aHit), hitCompare);

    currentHit += numHits;
  }

  fprintf(stderr, "filtering hits\n");

  double L  = 0.0;
  double H  = 0.0;
  double V  = 0.1;
  double M  = 1.0;
  double MC = 0.0;
  uint32 ML = 0;

  double minIdentity = 98.0;
  double minCoverage = 96.0;

  for (uint32 Hcnt = 10; Hcnt <= 100; Hcnt += 10) {
    for (uint32 Lcnt = 10; Lcnt < Hcnt && Lcnt < 60; Lcnt += 10) {
      for (uint32 Vcnt = 10; Vcnt < 100; Vcnt += 10) {
#ifdef SHOW_ONE
        Lcnt = 30;
        Hcnt = 40;
        Vcnt = 100;
#endif
        L = Lcnt / 100.0;
        H = Hcnt / 100.0;
        V = Vcnt / 100.0;

        uint32 truepositive  = 0;
        uint32 falsepositive = 0;
        uint32 truenegative  = 0;
        uint32 falsenegative = 0;

        for (uint32 currentHit = 0; currentHit < hitsLen; ) {

          //  Find the number of hits for this ESTid
          //
          uint32 estOfInterest = hits[currentHit]._qsIdx;
          uint32 numHits       = 0;
          for (uint32 t=currentHit; (t < hitsLen) && (hits[t]._qsIdx == estOfInterest); t++)
            numHits++;

          double h = hits[currentHit].scoreCov - hits[currentHit + numHits - 1].scoreCov;
          double p = 0.0;

          if (h <= L)    p = 1.0;
          if (h >= H)    p = V;
          if (p == 0.0)  p = 1.0 - (1.0 - V) * (h - L) / (H - L);

          //  check p; it should be between V and 1.0
          if (p > 1.0) {
            fprintf(stderr, "error in p; p=%f h=%f (%f %f %f)\n", p, h, L, H, V);
            p = 1.0;
          }

          if (p < V) {
            fprintf(stderr, "error in p; p=%f h=%f (%f %f %f)\n", p, h, L, H, V);
            p = V;
          }

          //  Output the top p% hits, by score.
          //
          double cutL = hits[currentHit].scoreCov - p * h;

          if (cutL > M)
            cutL = M;

#ifdef SHOW_ONE
          fprintf(stdout, "LHV = %f %f %f  p=%f h=%f  cutL=%f\n", L, H, V, p, h, cutL);
#endif


          for (uint32 i=currentHit; i < currentHit + numHits; i++) {
            if ((cutL <= hits[i].scoreCov) &&
                ((MC <= hits[i].scoreCov) || (ML <= hits[i]._covered))) {
#ifdef SHOW_ONE
              fprintf(stdout, "POS: (%f)", hits[i].scoreCov);
              ahit_printASCII(hits+i, stdout);
#endif
              if ((hits[i]._yesno == 1) && (hits[i]._identity >= minIdentity) && (hits[i]._coverage >= minCoverage))
                truepositive++;
              else
                falsepositive++;
            } else {
#ifdef SHOW_ONE
              fprintf(stdout, "NEG: (%f)", hits[i].scoreCov);
              ahit_printASCII(hits+i, stdout);
#endif
              if ((hits[i]._yesno == 1) && (hits[i]._identity >= minIdentity) && (hits[i]._coverage >= minCoverage))
                falsenegative++;
              else
                truenegative++;
            }
          }

#ifdef SHOW_ONE
          fprintf(stdout, "----\n");
#endif
          currentHit += numHits;
        }

        //  Print L, H, V, sensitivity, specificity
        //
        fprintf(stdout, "%f %f %f  %6.4f %6.4f  "uint32FMT" "uint32FMT" "uint32FMT" "uint32FMT"\n",
                L, H, V,
                (double)truepositive / (truepositive + falsenegative),
                (double)truenegative / (truenegative + falsepositive),
                truepositive, falsepositive,
                truenegative, falsenegative);
        fflush(stdout);

#ifdef SHOW_ONE
        exit(0);
#endif
      }
    }
  }

  return(0);
}

