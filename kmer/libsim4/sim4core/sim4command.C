#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sim4command.H"

#include <algorithm>

using namespace std;


//  Run a single EST against a genomic range
//
//  XXX: We should pull out the EST and GEN from the seqCache,
//  and store them as the "two char*" method.
//
sim4command::sim4command(uint32      ESTid,
                         seqCache   *ESTs,
                         uint32      GENid,
                         uint32      GENlo,
                         uint32      GENhi,
                         seqCache   *GENs,
                         bool        doFor,
                         bool        doRev) {

  _estIdx = ESTid;

  _ESTs              = ESTs;
  _ESTloaded         = 0L;
  _ESTsequence       = 0L;
  _ESTsequenceLength = 0;

  _genIdx = GENid;
  _genLo  = GENlo;
  _genHi  = GENhi;

  _GENs              = GENs;
  _GENloaded         = 0L;
  _GENsequence       = 0L;
  _GENsequenceLength = 0;

  _doForward = doFor;
  _doReverse = doRev;

  _externalSeedsLen   = 0;
  _externalSeedsMax   = 0;
  _externalSeeds      = 0L;
}


sim4command::sim4command(seqInCore  *EST,
                         seqInCore  *GEN,
                         uint32      GENlo,
                         uint32      GENhi,
                         bool        doFor,
                         bool        doRev) {

  _estIdx = EST->getIID();

  _ESTs              = 0L;
  _ESTloaded         = EST;
  _ESTsequence       = 0L;
  _ESTsequenceLength = 0;

  _genIdx = GEN->getIID();
  _genLo  = GENlo;
  _genHi  = GENhi;

  _GENs              = 0L;
  _GENloaded         = GEN;
  _GENsequence       = 0L;
  _GENsequenceLength = 0;

  _doForward = doFor;
  _doReverse = doRev;

  _externalSeedsLen   = 0;
  _externalSeedsMax   = 0;
  _externalSeeds      = 0L;
}


//  Use two char*'s for sequence sources
//
sim4command::sim4command(char             *EST,
                         uint32            ESTlen,
                         char             *GEN,
                         uint32            GENlen,
                         uint32            GENlo,
                         uint32            GENhi,
                         bool              doFor,
                         bool              doRev) {
  _estIdx = 0;

  _ESTs              = 0L;
  _ESTloaded         = 0L;
  _ESTsequence       = EST;
  _ESTsequenceLength = ESTlen;

  _genIdx = 0;
  _genLo  = GENlo;
  _genHi  = GENhi;

  _GENs              = 0L;
  _GENloaded         = 0L;
  _GENsequence       = GEN;
  _GENsequenceLength = GENlen;

  _doForward = doFor;
  _doReverse = doRev;

  _externalSeedsLen   = 0;
  _externalSeedsMax   = 0;
  _externalSeeds      = 0L;
}


sim4command::~sim4command() {
  if (_ESTs)
    delete _ESTloaded;
  if (_GENs)
    delete _GENloaded;

  delete [] _externalSeeds;
}


//  Make absolutely sure that the genomic sequence start and end
//  positions are within the actual sequence.  Ideally, this should
//  be checked by whatever generates the input, but it probably
//  isn't.
//
//  If the end position is too big, make it the same as the sequence
//  length.
//
//  If the start position is bigger than the (corrected) end
//  position, make it 100K less than the end position.
//
//  This has the side-effect of loading the genomic sequence.
//
void
sim4command::finalize(void) {

  if (_genHi > getGENlength())
    _genHi = getGENlength();

  if (_genLo > _genHi)
    if (_genHi > 100000)
      _genLo = _genHi - 100000;
    else
      _genLo = 0;
}



//  get() routines have multple cases
//
//  if no fastaBase, they can quickly return
//  otherwise
//  if nothing loaded or the thing loaded isn't right:
//    delete the current
//    load the correct
//

void
sim4command::loadEST(void) {
  if ((_ESTloaded == 0L) ||
      (_ESTloaded->getIID() != _estIdx)) {
    delete _ESTloaded;
    _ESTloaded = _ESTs->getSequenceInCore(_estIdx);
  }
}


uint32
sim4command::getESTidx(void) {
  if (_ESTsequence)
    return(0);
  return(_estIdx);
}

char*
sim4command::getESTheader(void) {
  static char *xxx = "anonymous cDNA sequence";
  if (_ESTsequence)
    return(xxx);
  loadEST();
  return(_ESTloaded->header());
}

char*
sim4command::getESTsequence(void) {
  if (_ESTsequence)
    return(_ESTsequence);
  loadEST();
  return(_ESTloaded->sequence());
}

uint32
sim4command::getESTlength(void) {
  if (_ESTsequence)
    return(_ESTsequenceLength);
  loadEST();
  return(_ESTloaded->sequenceLength());
}





void
sim4command::loadGEN(void) {
  if ((_GENloaded == 0L) ||
      (_GENloaded->getIID() != _genIdx)) {
    delete _GENloaded;
    _GENloaded = _GENs->getSequenceInCore(_genIdx);
  }
}

char*
sim4command::getGENheader(void) {
  char *xxx = "anonymous genomic sequence";
  if (_GENsequence)
    return(xxx);
  loadGEN();
  return(_GENloaded->header());
}

char*
sim4command::getGENsequence(void) {
  if (_GENsequence)
    return(_GENsequence);
  loadGEN();
  return(_GENloaded->sequence());
}

uint32
sim4command::getGENlength(void) {
  if (_GENsequence)
    return(_GENsequenceLength);
  loadGEN();
  return(_GENloaded->sequenceLength());
}



////////////////////////////////////////
//
//  This expects base-based seeds.
//  This expects that the position of the seed is the base in the seed.
//  This expects that GENpos is relative to the genomic subsequence.
//
//  If reverse-complement match, the EST is reversed, the GEN is forward.
//
void
sim4command::addSeed(uint32 GENpos, uint32 ESTpos, uint32 length) {

  if (_externalSeedsLen >= _externalSeedsMax) {
    if (_externalSeedsMax == 0)
      _externalSeedsMax = 256;
    _externalSeedsMax *= 2;
    externalSeed *n = new externalSeed [_externalSeedsMax];
    memcpy(n, _externalSeeds, sizeof(externalSeed) * _externalSeedsLen);
    delete [] _externalSeeds;
    _externalSeeds = n;
  }

  _externalSeeds[_externalSeedsLen]._GENposition = GENpos;
  _externalSeeds[_externalSeedsLen]._ESTposition = ESTpos;
  _externalSeeds[_externalSeedsLen]._length      = length;

  //  fprintf(stderr, "sim4command::addSeed()-- GEN="uint32FMT" EST="uint32FMT" of length "uint32FMT"\n", GENpos, ESTpos, length);

  _externalSeedsLen++;
}



void
sim4command::sortExternalSeeds(void) {
  sort(_externalSeeds, _externalSeeds + _externalSeedsLen);
}
