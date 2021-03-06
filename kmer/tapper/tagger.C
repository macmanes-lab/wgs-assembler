#include "tapperTag.H"
#include "tapperResult.H"
#include "tapperAlignment.H"
#include "tapperHit.H"

#include "seqCache.H"

//  Convert reads from ASCI to tapper binary.
//
//  ASSUMPTIONS
//
//  1) User is smart enough to give the correct set of mated files.
//  Code doesn't check that an F tag goes with an R tag, just that the
//  tag coordinates agree.  It is possible to mate an F to an F if the
//  wrong inputs are given.
//
//  2) Tag coords are 16-bit integers.  File UIDs are 16-bit integers.
//


//  Define this to test the encode/decode functionality.
//#define TEST_ENCODING


int
tapperTagCompare(const void *a, const void *b) {
  tapperTag const *A = (tapperTag const *)a;
  tapperTag const *B = (tapperTag const *)b;
  if (A->tagID()  < B->tagID()) return(-1);
  return(A->tagID() != B->tagID());
}


bool
readTag(uint32 fileUID, FILE *seq, FILE *qlt, tapperTag *T) {
  static uint16  id[4];
  static char    seqhdr[1024];
  static char    seqseq[1024];
  static char    qlthdr[1024];
  static char    qltseq[1024];
  static uint64  qltnum[1024];
  static splitToWords  S;

  seqhdr[0] = 0;
  seqseq[0] = 0;
  qlthdr[0] = 0;
  qltseq[0] = 0;

  if (feof(seq) || feof(qlt))
    return(false);

  fgets(seqhdr, 1024, seq);
  while (seqhdr[0] == '#')
    fgets(seqhdr, 1024, seq);
  fgets(seqseq, 1024, seq);

  fgets(qlthdr, 1024, qlt);
  while (qlthdr[0] == '#')
    fgets(qlthdr, 1024, qlt);
  fgets(qltseq, 1024, qlt);

  if ((seqhdr[0] == 0) || (qlthdr[0] == 0))
    return(false);

  chomp(seqhdr);
  chomp(seqseq);
  chomp(qlthdr);
  chomp(qltseq);

  if (strcmp(seqhdr, qlthdr) != 0)
    fprintf(stderr, "WARNING:  Got unpaired seq '%s' and qlt '%s'\n", seqhdr, qlthdr);

  //  Assumes the header is >461_28_1918_F3
  //  -- copies it to the left by one to remove the >
  //  -- the loop below doesn't move the zero-terminator
  //  -- resulting string is "461 28 1918 F33"
  //
  for (uint32 i=1; seqhdr[i]; i++) {
    if (seqhdr[i] == '_')
      seqhdr[i] = ' ';
    seqhdr[i-1] = seqhdr[i];
  }

  S.split(seqhdr);

  id[0] = fileUID;
  id[1] = strtouint32(S[0], 0L);
  id[2] = strtouint32(S[1], 0L);
  id[3] = strtouint32(S[2], 0L);

  S.split(qltseq);

  //  Not sure why there are negative numbers here, but there are.
  //
  for (uint32 i=0; i<S.numWords(); i++) {
    qltnum[i] = (S[i][0] == '-') ? 0 : strtouint64(S[i], 0L);

#ifdef TEST_ENCODING
    //  We need to fudge the QV's here, so our tests pass.
    if (qltnum[i] > 31)
      qltnum[i] = 31;
#endif
  }

  T->encode(id, seqseq, qltnum);

#ifdef TEST_ENCODING
  {
    uint16  it[4];
    char    seqtst[1024];
    uint64  qlttst[1024];

    T->decode(it, seqtst, qlttst);

    uint32  len = strlen(seqtst);
    uint32  fail = 0;
    uint64  qltsum=0, tstsum=0;

    for (uint32 l=0; l<len; l++) {
      qltsum += qltnum[l];
      tstsum += qlttst[l];
      if ((seqseq[l] != seqtst[l]) || (qltnum[l] != qlttst[l]))
        fail++;
    }

    if ((id[0] != it[0]) ||
        (id[1] != it[1]) ||
        (id[2] != it[2]) ||
        (id[3] != it[3]) ||
        (fail)) {
      fprintf(stderr, "FAIL:  ("uint32FMT"_"uint32FMT"_"uint32FMT"_"uint32FMT",%s,"uint64FMT") != ("uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT",%s,"uint64FMT")\n",
              id[0], id[1], id[2], id[3], seqseq, qltsum,
              it[0], it[1], it[2], it[3], seqtst, tstsum);
      for (uint32 l=0; l<len; l++)
        fprintf(stderr, "  %2d -- "uint64FMT" "uint64FMT"\n", l, qltnum[l], qlttst[l]);
    }
  }
#endif

  return(true);
}



void
dumpTagFileStats(char *tagfile) {
  tapperTagFile  *TF = new tapperTagFile(tagfile, 'r');

  if (TF->metaData()->isPairedTagFile()) {
    fprintf(stdout, "%s\ttype\tmated tags\n", tagfile);
    fprintf(stdout, "%s\tlength\t"uint32FMT"\n", tagfile, TF->metaData()->tagSize());
    fprintf(stdout, "%s\tnumMates\t"uint64FMT"\n", tagfile, TF->numberOfMatePairs());
    fprintf(stdout, "%s\tmean\t"uint32FMT"\n", tagfile, TF->metaData()->mean());
    fprintf(stdout, "%s\tstddev\t"uint32FMT"\n", tagfile, TF->metaData()->stddev());
  } else {
    fprintf(stdout, "%s\ttype\tfragment tags\n", tagfile);
    fprintf(stdout, "%s\tlength\t"uint32FMT"\n", tagfile, TF->metaData()->tagSize());
    fprintf(stdout, "%s\tnumTags\t"uint64FMT"\n", tagfile, TF->numberOfFragmentTags());
 
  }
}


void
dumpTagFile(char *tagfile) {
  tapperTagFile  *TF = new tapperTagFile(tagfile, 'r');
  tapperTag       a, b;
  uint16          ida[4],    idb[4];
  char            seqa[265], seqb[256];
  char            quaa[256], quab[256];
  uint64          qvsa[256], qvsb[256];
  uint32          i;

  if (TF->metaData()->isPairedTagFile()) {
    while (TF->get(&a, &b)) {
      a.decode(ida, seqa, qvsa);
      b.decode(idb, seqb, qvsb);
      for (i=0; seqa[i+1]; i++)
        quaa[i] = qvsa[i] + '0';
      for (i=0; seqb[i+1]; i++)
        quab[i] = qvsb[i] + '0';
      fprintf(stdout, ">"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t%s/%s\t>"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t%s/%s\n",
              ida[0], ida[1], ida[2], ida[3], seqa, quaa,
              idb[0], idb[1], idb[2], idb[3], seqb, quab);
    }
  } else {
    while (TF->get(&a)) {
      a.decode(ida, seqa, qvsa);
      for (i=0; seqa[i+1]; i++)
        quaa[i] = qvsa[i] + '0';
      fprintf(stdout, ">"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t%s/%s\n",
              ida[0], ida[1], ida[2], ida[3], seqa, quaa);
    }
  }

  delete TF;
}



int
main(int argc, char **argv) {
  char  *prefix  = 0L;

  uint32  sampleSize    = 0;
  char   *sampleFile    = 0L;
  uint32  sampleErrors  = 3;
  uint32  sampleTagSize = 25;

  uint32  tagfuid = 0,   tagruid = 0;
  char   *tagfseq = 0L, *tagrseq  = 0L;
  char   *tagfqlt = 0L, *tagrqlt  = 0L;

  uint32  mean=0, stddev=0;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-tagout", 5) == 0) {
      prefix   = argv[++arg];

    } else if (strncmp(argv[arg], "-tags", 5) == 0) {
      tagfuid  = strtouint32(argv[++arg], 0L);
      tagfseq  = argv[++arg];
      tagfqlt  = argv[++arg];

    } else if (strncmp(argv[arg], "-ftags", 2) == 0) {
      tagfuid  = strtouint32(argv[++arg], 0L);
      tagfseq  = argv[++arg];
      tagfqlt  = argv[++arg];
    } else if (strncmp(argv[arg], "-rtags", 2) == 0) {
      tagruid  = strtouint32(argv[++arg], 0L);
      tagrseq  = argv[++arg];
      tagrqlt  = argv[++arg];

    } else if (strncmp(argv[arg], "-insertsize", 2) == 0) {
      mean   = strtouint32(argv[++arg], 0L);
      stddev = strtouint32(argv[++arg], 0L);

      if (mean > MAX_INSERT_SIZE)
        fprintf(stderr, "%s: insert size limited to at most %dbp.\n", argv[0], MAX_INSERT_SIZE), exit(1);
      if (stddev > MAX_INSERT_DEVIATION)
        fprintf(stderr, "%s: insert size limited to at most +- %dbp.\n", argv[0], MAX_INSERT_DEVIATION), exit(1);

    } else if (strcmp(argv[arg], "-sample") == 0) {
      sampleSize = strtouint32(argv[++arg], 0L);
      sampleFile = argv[++arg];
    } else if (strcmp(argv[arg], "-sampleerrors") == 0) {
      sampleErrors = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-sampletagsize") == 0) {
      sampleTagSize = strtouint32(argv[++arg], 0L);

    } else if (strncmp(argv[arg], "-stats", 3) == 0) {
      dumpTagFileStats(argv[++arg]);
      exit(0);

    } else if (strncmp(argv[arg], "-dump", 2) == 0) {
      dumpTagFile(argv[++arg]);
      exit(0);

    } else {
      err++;
    }
    arg++;
  }
  if (sampleFile == 0L) {
    if ((tagfseq == 0L) || (tagfqlt == 0L))  err++;
    if ((tagfseq != 0L) && (tagfqlt == 0L))  err++;
    if ((tagfseq == 0L) && (tagfqlt != 0L))  err++;
  }
  if ((err) || (prefix == 0L)) {
    fprintf(stderr, "usage: %s -tagout prefix  -tags fileUID xx.csfasta xx.qual\n", argv[0]);
    fprintf(stderr, "usage: %s -tagout prefix -ftags fileUID ff.csfasta ff.qual -rtags fileUID rr.csfasta rr.qual\n", argv[0]);
    fprintf(stderr, "usage: %s -dump file.tapperTags\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "unmated tags will be placed in 'prefix.frag.tapperTags'\n");
    fprintf(stderr, "  mated tags will be placed in 'prefix.mate.tapperTags'\n");
    exit(1);
  }

  uint64  numTagsF = 0, maxTagsF = 0;
  uint64  numTagsR = 0, maxTagsR = 0;
  uint64  numTagsM = 0;

  tapperTag      *TF = 0L;
  tapperTag      *TR = 0L;

  //  If given a sampleFile, generate some tags from there.
  if (sampleFile) {
    seqCache   *F = new seqCache(sampleFile);
    seqInCore  *s = F->getSequenceInCore();

    uint32  pos = 0;
    uint32  len = s->sequenceLength();

    uint16  id[4];
    char    cor[64] = {0};
    char    seq[64] = {0};
    uint64  qlt[64] = {0};

    char    acgt[4] = {'A', 'C', 'G', 'T'};

    mt_s   *mtctx = mtInit(time(0));

    maxTagsF = sampleSize;
    TF       = new tapperTag [maxTagsF];

    maxTagsR = sampleSize;
    TR       = new tapperTag [maxTagsR];

    for (uint32 i=0; i<sampleSize; i++) {
      pos = mtRandom32(mtctx) % (len - sampleTagSize);

      char  n = acgt[mtRandom32(mtctx) % 4];
      char  l = n;

      cor[0] = n;
      seq[0] = n;

      bool   doForward = (mtRandom32(mtctx) & 0x1000) == 0x1000;
      //doForward = false;

      if (doForward) {
        uint32 sp = pos;
        for (uint32 x=1; x<=sampleTagSize; x++) {
          n = s->sequence()[sp++];
          cor[x] = n;
          seq[x] = baseToColor[l][n];
          l = n;
        }
      } else {
        uint32 sp = pos + sampleTagSize - 1;
        for (uint32 x=1; x<=sampleTagSize; x++) {
          n = complementSymbol[s->sequence()[sp--]];
          cor[x] = n;
          seq[x] = baseToColor[l][n];
          l = n;
        }
      }

      //  Insert errors.

      char     errors[256] = {0};
      char     errort[256] = {0};
      uint32   nerrs = mtRandom32(mtctx) % (sampleErrors + 1);

      for (uint32 xx=0; xx<nerrs; xx++) {
        uint32 e = mtRandom32(mtctx) % (sampleTagSize-1) + 1;
        char   o = seq[e];
        seq[e] = seq[e] + 1;
        if (seq[e] > '3')
          seq[e] = '0';
        sprintf(errort, "\t%c->%c@%02d", o, seq[e], e);
        strcat(errors, errort);
      }

      id[0] = i;
      id[1] = 0;
      id[2] = 0;
      id[3] = 0;

      fprintf(stdout, "F\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t0\t"uint32FMT"\t%c\t%s%s\t%s\n",
              id[0], id[1], id[2], id[3],
              pos,
              (doForward) ? 'f' : 'r',
              cor+1,
              errors,
              seq);

      //  TF is NOT just storing the 'forward' reads, it's all the
      //  reads from the first half of the mate.  Since we're not
      //  mated, this is just all reads.

      TF[numTagsF++].encode(id, seq, qlt);
    }
  }

  //
  //  Suck in all the F tags.
  //
  if (tagfseq) {
    FILE *fseq = fopen(tagfseq, "r");
    FILE *fqlt = fopen(tagfqlt, "r");

    speedCounter *CT = new speedCounter(" reading F tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

    maxTagsF = sizeOfFile(tagfseq) / 44 + 1000000;
    TF       = new tapperTag [maxTagsF];

    while (readTag(tagfuid, fseq, fqlt, TF + numTagsF)) {
      numTagsF++;
      if (numTagsF >= maxTagsF)
        fprintf(stderr, "Too many F tags.  Boom.\n"), exit(1);
      CT->tick();
    }
    delete CT;

    fclose(fseq);
    fclose(fqlt);
  }

  //
  //  Suck in all the R tags.
  //
  if (tagrseq) {
    FILE *rseq = fopen(tagrseq, "r");
    FILE *rqlt = fopen(tagrqlt, "r");

    speedCounter *CT = new speedCounter(" reading R tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

    maxTagsR = sizeOfFile(tagrseq) / 44 + 1000000;
    TR       = new tapperTag [maxTagsR];;

    while (readTag(tagruid, rseq, rqlt, TR + numTagsR)) {
      numTagsR++;
      if (numTagsR >= maxTagsR)
        fprintf(stderr, "Too many R tags.  Boom.\n"), exit(1);
      CT->tick();
    }
    delete CT;

    fclose(rseq);
    fclose(rqlt);
  }

    maxTagsF = numTagsF;
    numTagsF = 0;

    maxTagsR = numTagsR;
    numTagsR = 0;

  //
  //  Sort them.
  //
  qsort_mt(TF, maxTagsF, sizeof(tapperTag), tapperTagCompare, 4, 4 * 1024 * 1024);
  qsort_mt(TR, maxTagsR, sizeof(tapperTag), tapperTagCompare, 4, 4 * 1024 * 1024);

  //
  //  Merge to find pairs, output.
  //
  char            fragout[FILENAME_MAX];
  char            mateout[FILENAME_MAX];

  sprintf(fragout, "%s.frag.tapperTags", prefix);
  sprintf(mateout, "%s.mate.tapperTags", prefix);

  tapperTagFile  *TOfrag = 0L;
  tapperTagFile  *TOmate = 0L;

  speedCounter *CF = new speedCounter(" writing frag tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);
  speedCounter *CM = new speedCounter(" writing mate tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

  while ((numTagsF < maxTagsF) && (numTagsR < maxTagsR)) {
    uint64   fID = TF[numTagsF].tagID() & uint64MASK(48);
    uint64   rID = TR[numTagsR].tagID() & uint64MASK(48);

    if (fID == rID) {
      if (TOmate == 0L)
        TOmate = new tapperTagFile(mateout, 'w');
      TOmate->put(TF + numTagsF, TR + numTagsR);
      numTagsF++;
      numTagsR++;
      numTagsM++;
      CM->tick();
    } else if (fID < rID) {
      if (TOfrag == 0L)
        TOfrag = new tapperTagFile(fragout, 'w');
      TOfrag->put(TF + numTagsF);
      numTagsF++;
      CF->tick();
    } else {
      if (TOfrag == 0L)
        TOfrag = new tapperTagFile(fragout, 'w');
      TOfrag->put(TR + numTagsR);
      numTagsR++;
      CF->tick();
    }
  }
  while (numTagsF < maxTagsF) {
    if (TOfrag == 0L)
      TOfrag = new tapperTagFile(fragout, 'w');
    TOfrag->put(TF + numTagsF);
    numTagsF++;
    CF->tick();
  }
  while (numTagsR < maxTagsR) {
    if (TOfrag == 0L)
      TOfrag = new tapperTagFile(fragout, 'w');
    TOfrag->put(TR + numTagsR);
    numTagsR++;
    CF->tick();
  }

  delete CF;
  delete CM;

  if (TOmate)
    TOmate->metaData()->setMeanStdDev(mean, stddev);

  delete TOmate;
  delete TOfrag;

  delete [] TR;
  delete [] TF;
}
