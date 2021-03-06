#include "snapper2.H"


#define MATCH            0
#define GAPA             1
#define GAPB             2
#define STOP             3

#define MATCHSCORE       2
#define GAPSCORE        -3
#define MISMATCHSCORE   -1


void
reverse(char *a, char *b, int len) {
  char   c=0;
  char  *s=a,  *S=a+len-1;
  char  *q=b,  *Q=b+len-1;

  while (s < S) {
    c    = *s;
    *s++ =  *S;
    *S-- =  c;

    c    = *q;
    *q++ = *Q;
    *Q-- =  c;
  }
}



class dpMatch {
public:
  dpMatch() {
    matches  = 0;
    alignLen = 0;

    begI = begJ = 0;
    endI = endJ = 0;
    lenA = lenB = 0;
  };

  int      matches;
  int      alignLen;

  int      begI, begJ;
  int      endI, endJ;
  int      lenA, lenB;

  char    *alignA;
  char    *alignB;
};


class dpMatrix {
private:
  typedef struct {
    unsigned int  score  : 30;
    unsigned int  action : 2;
  } dpCell;

public:
  dpMatrix() {
    aMax = 0;
    bMax = 0;

    alignA = 0L;
    alignB = 0L;
    matrix = 0L;
  };

  ~dpMatrix() {
    delete [] alignA;
    delete [] alignB;
    delete [] matrix;
  };

  void   dpMatrixInit(int lenA, int lenB) {

    if ((aMax <= lenA) || (bMax <= lenB)) {
      delete [] alignA;
      delete [] alignB;
      delete [] matrix;

      aMax = MAX(aMax, lenA) + 1000;
      bMax = MAX(bMax, lenB) + 1000;

      fprintf(stderr, "dpMatrix-- reallocate to "uint32FMT" x "uint32FMT"\n", aMax, bMax);

      alignA = new char   [aMax + bMax + 1];
      alignB = new char   [bMax + bMax + 1];
      matrix = new dpCell [aMax * bMax];
    }

    int i, j, p = 0;

    for (i=0; i<lenA+1; i++) {
      matrix[p].score = 1 << 29;
      matrix[p].action = STOP;
      p += bMax;
    }

    p = 0;
    for (j=0; j<lenB+1; j++) {
      matrix[p].score = 1 << 29;
      matrix[p].action = STOP;
      p++;
    }
  };

  int    dpMatrixCellGetScore(int a, int b) {
    return(matrix[a * bMax + b].score);
  };

  int    dpMatrixCellGetAction(int a, int b) {
    return(matrix[a * bMax + b].action);
  };

  void   dpMatrixCellSet(int a, int b, int score, int action) {
    dpCell  x;
    x.score  = score;
    x.action = action;
    matrix[a * bMax + b] = x;
  };

  dpMatch *dpAlign(char *stringA, int lenA,
                   char *stringB, int lenB,
                   dpMatch *match);
private:
  int      aMax;
  int      bMax;

  char    *alignA;
  char    *alignB;
  dpCell  *matrix;
};


dpMatch *
dpMatrix::dpAlign(char     *stringA,  int lenA,
                  char     *stringB,  int lenB,
                  dpMatch  *match) {

  int   i, j;

  dpMatrixInit(lenA, lenB);

  int   scoreMax  = 0;

  int   begI=0, endI=0, curI=0;
  int   begJ=0, endJ=0, curJ=0;

  for (i=1; i<=lenA; i++){
    for (j=1; j<=lenB; j++){

      //  Pick the max of these

      int ul = dpMatrixCellGetScore(i-1, j-1) + ((stringA[i-1] == stringB[j-1]) ? MATCHSCORE : MISMATCHSCORE);
      int lf = dpMatrixCellGetScore(i-1, j) + GAPSCORE;
      int up = dpMatrixCellGetScore(i, j-1) + GAPSCORE;

      // (i,j) is the beginning of a subsequence, our default behavior
      int sc = 1 << 29;
      int ac = STOP;

      if (sc < ul) {
        sc = ul;
        ac = MATCH;
      }

      if (sc < lf) {
        sc = lf;
        ac = GAPB;
      }

      if (sc < up) {
        sc = up;
        ac = GAPA;
      }

      dpMatrixCellSet(i, j, sc, ac);

      if (scoreMax < sc) {
        scoreMax  = sc;
        endI = curI = i;
        endJ = curJ = j;
      }
    }
  }

  //fprintf(stdout, "SCORE %d at %d,%d\n", scoreMax - (1 << 29), endI, endJ);

  int  alignLen  = 0;
  int  matches   = 0;
  int  terminate = 0;

  while (terminate == 0) {
    switch (dpMatrixCellGetAction(curI, curJ)) {
      case STOP:
        terminate = 1;
        break;
      case MATCH:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = stringB[curJ-1];

        if (alignA[alignLen] == alignB[alignLen]) {
          alignA[alignLen] = tolower(alignA[alignLen]);
          alignB[alignLen] = tolower(alignB[alignLen]);
          matches++;
        } else {
          //fprintf(stdout, "MIS  at %d\n", alignLen);
          alignA[alignLen] = toupper(alignA[alignLen]);
          alignB[alignLen] = toupper(alignB[alignLen]);
        }

        curI--;
        curJ--;
        alignLen++;
        break;
      case GAPA:
        //fprintf(stdout, "GAPA at %d\n", alignLen);
        alignA[alignLen] = '-';
        alignB[alignLen] = stringB[curJ-1];
        curJ--;
        alignLen++;
        break;
      case GAPB:
        //fprintf(stdout, "GAPB at %d\n", alignLen);
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = '-';
        curI--;
        alignLen++;
        break;
    }
  }

  begI = curI;
  begJ = curJ;

  alignA[alignLen] = 0;
  alignB[alignLen] = 0;

  reverse(alignA, alignB, alignLen);

  match->matches  = matches;
  match->alignLen = alignLen;
  match->begI     = begI;
  match->begJ     = begJ;
  match->endI     = endI;
  match->endJ     = endJ;
  match->lenA     = lenA;
  match->lenB     = lenB;

//warning alignA and alignB aliases to dpMatrix
  match->alignA   = alignA;
  match->alignB   = alignB;

  return(match);
}








void
doPolishDP(searcherState       *state,
           query               *qry) {

  //  For the autofilter
  uint64   successes    = uint64ZERO;
  uint64   successMask  = uint64MASK(config._afLength);
  uint32   attempts     = 0;

  if (qry->theHitsLen == 0)
    return;

  qry->theOutputLen = 0;
  qry->theOutputMax = 2 * 1024 * qry->theHitsLen;
  qry->theOutput    = new char [qry->theOutputMax];

  qry->theOutput[0] = 0;

  //  Move these to searcherState!

  if (state->DP == 0L)
    state->DP = new dpMatrix;

  dpMatch   match;
  dpMatrix *matrix = (dpMatrix *)state->DP;

  for (uint32 h=0; h<qry->theHitsLen; h++) {

    //  If the hit was discarded, move along.
    //
    if (qry->theHits[h]._status & AHIT_DISCARDED)
      continue;

    //  If the hit was filtered out, move along.
    //
    if ((config._doValidation == false) &&
        ((qry->theHits[h]._status & AHIT_POLISHABLE) == 0) &&
        ((qry->theHits[h]._status & AHIT_HAS_UNIQUE) == 0))
      continue;

    //  If our recent success rate is pretty terrible, continue.
    //
    if (config._afEnabled) {
      if (attempts > config._afInit) {
        double  rat = countNumberOfSetBits64(successes) / (double)((attempts < config._afLength) ? attempts : config._afLength);

        //  If we've hit the end of the good polishes, give up.  But
        //  still do all the stuff with unique mers in them.
        //
        if (((qry->theHits[h]._status & AHIT_HAS_UNIQUE) == 0) &&
            (rat < config._afThreshold))
          continue;
      }

      attempts++;
    }

    //
    //  Polish it up!
    //

    seqInCore            *QRYseq = qry->seq;
    seqInCore            *GENseq = genome->getSequenceInCore(qry->theHits[h]._dsIdx);
    uint32                GENlo  = qry->theHits[h]._dsLo;
    uint32                GENhi  = qry->theHits[h]._dsHi;

    char  *q = QRYseq->sequence();
    char  *g = GENseq->sequence() + GENlo;

    if (GENhi > GENseq->sequenceLength())
      GENhi = GENseq->sequenceLength();

    uint32   qlen = qry->seq->sequenceLength();
    uint32   glen = GENhi - GENlo;

    bool     doForward =  qry->theHits[h]._status & AHIT_DIRECTION_MASK;
    bool     doReverse = !doForward;

    if (doReverse) {
      reverseComplementSequence(q, qlen);
    }

#if 0
    fprintf(stderr, "align QRYlen="uint32FMT" GEN="uint32FMT"-"uint32FMT" GENlen="uint32FMT"\n",
            qlen, GENlo, GENhi, glen);
#endif

    //if ((qlen * 3 > glen) && ((qlen / 1024) * (glen / 1024) < 4 * 1024))

    matrix->dpAlign(q, qlen, g, glen, &match);

    if (doReverse) {
      reverseComplementSequence(q, qlen);

      uint32 x = match.begI;
      match.begI = qlen - match.endI;
      match.endI = qlen - x;
    }


    //  Build the proper match if it's even remotely good
    //
    if (match.matches > 0) {
      sim4polish      p;
      sim4polishExon  e;

      qry->theHits[h]._status |= AHIT_VERIFIED;

      p._estID    = QRYseq->getIID();
      p._estLen   = QRYseq->sequenceLength();
      p._estPolyA = 0;
      p._estPolyT = 0;

      p._genID            = GENseq->getIID();
      p._genRegionOffset  = GENlo;
      p._genRegionLength  = GENhi - GENlo;

      p._numMatches   = match.matches;
      p._numMatchesN  = 0;
      p._numCovered   = match.endI - match.begI;
      
      p._percentIdentity  = 0;
      p._querySeqIdentity = 0;

      p._matchOrientation  = (doReverse) ? SIM4_MATCH_COMPLEMENT : SIM4_MATCH_FORWARD;
      p._strandOrientation = SIM4_STRAND_UNKNOWN;

      p._comment    = NULL;
      p._estDefLine = QRYseq->header();
      p._genDefLine = GENseq->header();

      p._numExons = 1;
      p._exons    = &e;

      e._estFrom           = match.begI + 1;
      e._estTo             = match.endI;
      e._genFrom           = match.begJ + GENlo + 1;
      e._genTo             = match.endJ + GENlo;
      e._numMatches        = match.matches;
      e._numMatchesN       = 0;
      e._percentIdentity   = 0;
      e._intronOrientation = SIM4_INTRON_NONE;

      //  The alignments are needed for updateAlignmentScores().

      e._estAlignment      = match.alignA;  //  'e' DOES NOT own this, must reset the pointer later.
      e._genAlignment      = match.alignB;

      p.s4p_updateAlignmentScores();

      //  Since we're not using sim4, the normal method of ignoring aligns doesn't work.
      //  Do it explicitly.

      if (config._doAlignments == false) {
        e._estAlignment      = NULL;
        e._genAlignment      = NULL;
      }


      //  Save it if it is truely good.
      if ((p._percentIdentity  >= config._minMatchIdentity) &&
          (p._querySeqIdentity >= config._minMatchCoverage)) {
        char *pstr = p.s4p_polishToString(sim4polishStyleDefault);

        uint32 l = (uint32)strlen(pstr);

        if (qry->theOutputLen + l + 1 >= qry->theOutputMax) {
          qry->theOutputMax = qry->theOutputMax + qry->theOutputMax + l;
          char *o = 0L;
          try {
            o = new char [qry->theOutputMax];
          } catch (...) {
            fprintf(stderr, "doPolish()-- Can't reallocate space for the output string ("uint32FMT" bytes) in thread "uint64FMT"\n", qry->theOutputMax, state->threadID);
            abort();
          }
          memcpy(o, qry->theOutput, sizeof(char) * qry->theOutputLen);
          delete [] qry->theOutput;
          qry->theOutput = o;
        }

        memcpy(qry->theOutput + qry->theOutputLen, pstr, sizeof(char) * l);
        qry->theOutputLen += l;

        qry->theOutput[qry->theOutputLen] = 0;

        delete [] pstr;

        //  Save the best scores
        //
        uint32 pi = p._percentIdentity;
        uint32 pc = p._querySeqIdentity;

        qry->theHits[h]._status |= pi << 16;
        qry->theHits[h]._status |= pc << 24;

        successes <<= 1;
        if ((pi  >= config._minMatchIdentity) &&
            (pc >= config._minMatchCoverage)) {
          //fprintf(stderr, "GOOD "uint32FMT" "uint32FMT"\n", pi, pc);
          successes |= uint64ONE;
        } else {
          //fprintf(stderr, "BAD  "uint32FMT" "uint32FMT"\n", pi, pc);
          successes |= uint64ZERO;
        }
        successes  &= successMask;
      }

      //  Before sim4polish and sim4polishExon go out of scope, reset the pointers.  Ugly, but needs
      //  to be done else the destructors try to delete things it doesn't own (alignments) or
      //  allocated on the stack (sim4polishExon).
      p._estDefLine   = 0L;
      p._genDefLine   = 0L;
      p._exons        = 0L;
      e._estAlignment = 0L;
      e._genAlignment = 0L;
    }

    delete GENseq;
  }  //  over all hits
}

