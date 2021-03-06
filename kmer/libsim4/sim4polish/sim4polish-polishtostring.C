#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

//#define DEBUG_CIGAR

const char *mOriFWD = "forward";
const char *mOriCMP = "complement";
const char *mOriERR = "error";
const char *mOriDEF = "UNKNOWN";

const char *sOriFWD = "forward";
const char *sOriREV = "reverse";
const char *sOriUNK = "unknown";
const char *sOriINT = "intractable";
const char *sOriABT = "aborted";
const char *sOriERR = "error";
const char *sOriDEF = "UNKNOWN";

const char *iOriPOS = " ->";
const char *iOriNEG = " <-";
const char *iOriAMB = " --";
const char *iOriGAP = " ==";
const char *iOriERR = " ??";
const char *iOriNOO = "";


bool            sim4polishStyleSet     = false;
sim4polishStyle sim4polishStyleDefault = sim4polishS4DB;
uint32          sim4polishPolishID     = 0;


char *
encodeGap(char *ref, char *tgt) {

  if ((ref == 0L) || (tgt == 0L))
    return(0L);

  uint32 lenref = strlen(ref);
  uint32 lentgt = strlen(tgt);

  assert(lenref == lentgt);

  char   *gap    = new char [3 * lenref];
  char   *gpp    = gap;

  char    gaptyp = 0;
  uint32  gapcnt = 0;

  for (uint32 i=0; i<lenref; i++) {
    if        ((ref[i] == '-') && (tgt[i] != '-')) {
      if (gaptyp != 'I') {
        if (gaptyp != 0) {
          sprintf(gpp, "%c"uint32FMT" ", gaptyp, gapcnt);
          while (*gpp) gpp++;
        }
        gaptyp = 'I';
        gapcnt = 0;
      }
      gapcnt++;
    } else if ((ref[i] != '-') && (tgt[i] == '-')) {
      if (gaptyp != 'D') {
        if (gaptyp != 0) {
          sprintf(gpp, "%c"uint32FMT" ", gaptyp, gapcnt);
          while (*gpp) gpp++;
        }
        gaptyp = 'D';
        gapcnt = 0;
      }
      gapcnt++;
    } else if ((ref[i] == '-') && (tgt[i] == '-')) {
      assert(0);
    } else {
      if (gaptyp != 'M') {
        if (gaptyp != 0) {
          sprintf(gpp, "%c"uint32FMT" ", gaptyp, gapcnt);
          while (*gpp) gpp++;
        }
        gaptyp = 'M';
        gapcnt = 0;
      }
      gapcnt++;
    }
  }

  if (gaptyp != 0) {
    sprintf(gpp, "%c"uint32FMT"", gaptyp, gapcnt);
    while (*gpp) gpp++;
  }

#ifdef DEBUG_CIGAR
  fprintf(stderr, "REF=%s\n", ref);
  fprintf(stderr, "TGT=%s\n", tgt);
  fprintf(stderr, "GAP=%s\n", gap);
  fprintf(stderr, "---\n");
#endif

  return(gap);
}





char *
sim4polish::s4p_polishToString(sim4polishStyle style) {
  char *ret = NULL;

  if (_numExons == 0)
    return(ret);

  switch (style) {
    case sim4polishS4DB:
      ret = s4p_polishToStringS4DB();
      break;
    case sim4polishGFF3:
      ret = s4p_polishToStringGFF3();
      break;
    case sim4polishATAC:
      ret = s4p_polishToStringATAC();
      break;
    default:
      fprintf(stderr, "s4p_polishToString()-- unknown style='%d'\n",
              style);
      exit(1);
  }

  return(ret);
}




char *
sim4polish::s4p_polishToStringS4DB(void) {
  const char   *mOri = mOriDEF;
  const char   *sOri = sOriDEF;
  const char   *iOri = iOriERR;

  //  Make a decent estimate of how much space we'll need to store the string
  //
  uint32 spaceNeeded = (1024 + 128 * _numExons +
                        ((_comment)    ? strlen(_comment)    : 0) +
                        ((_estDefLine) ? strlen(_estDefLine) : 0) +
                        ((_genDefLine) ? strlen(_genDefLine) : 0));

  for (uint32 i=0; i<_numExons; i++)
    if (_exons[i]._estAlignment)
      spaceNeeded += 2 * strlen(_exons[i]._estAlignment);

  char *outs = new char [spaceNeeded];
  char *outc = outs;

  switch (_matchOrientation) {
    case SIM4_MATCH_FORWARD:     mOri = mOriFWD;  break;
    case SIM4_MATCH_COMPLEMENT:  mOri = mOriCMP;  break;
    case SIM4_MATCH_ERROR:       mOri = mOriERR;  break;
    default:
      fprintf(stderr, "sim4reader: Unknown matchOrientation '"uint32FMT"' in printPolish()\n", _matchOrientation);
      mOri = mOriDEF;
      break;
  }

  switch (_strandOrientation) {
    case SIM4_STRAND_POSITIVE:    sOri = sOriFWD;  break;
    case SIM4_STRAND_NEGATIVE:    sOri = sOriREV;  break;
    case SIM4_STRAND_UNKNOWN:     sOri = sOriUNK;  break;
    case SIM4_STRAND_INTRACTABLE: sOri = sOriINT;  break;
    case SIM4_STRAND_FAILED:      sOri = sOriABT;  break;
    case SIM4_STRAND_ERROR:       sOri = sOriERR;  break;
    default:
      fprintf(stderr, "sim4reader: Unknown strandOrientation '"uint32FMT"' in printPolish()\n", _matchOrientation);
      sOri = sOriDEF;
      break;
  }

  sprintf(outc, "sim4begin\n"uint32FMT"["uint32FMT"-"uint32FMT"-"uint32FMT"] "uint32FMT"["uint32FMT"-"uint32FMT"] <"uint32FMT"-"uint32FMT"-"uint32FMT"-%s-%s>\n",
          _estID, _estLen, _estPolyA, _estPolyT,
          _genID, _genRegionOffset, _genRegionLength,
          _numMatches, _numMatchesN, _percentIdentity, mOri, sOri);
  while (*outc)  outc++;

  if (_comment) {
    sprintf(outc, "comment=%s\n", _comment);
    while (*outc)  outc++;
  }

  if (_estDefLine) {
    sprintf(outc, "edef=%s\n", _estDefLine);
    while (*outc)  outc++;
  }

  if (_genDefLine) {
    sprintf(outc, "ddef=%s\n", _genDefLine);
    while (*outc)  outc++;
  }

  for (uint32 i=0; i<_numExons; i++) {
    switch (_exons[i]._intronOrientation) {
      case SIM4_INTRON_POSITIVE:    iOri = iOriPOS;  break;
      case SIM4_INTRON_NEGATIVE:    iOri = iOriNEG;  break;
      case SIM4_INTRON_AMBIGUOUS:   iOri = iOriAMB;  break;
      case SIM4_INTRON_GAP:         iOri = iOriGAP;  break;
      case SIM4_INTRON_ERROR:       iOri = iOriERR;  break;
      default:                      iOri = iOriNOO;  break;
    }

    sprintf(outc, ""uint32FMT"-"uint32FMT" ("uint32FMT"-"uint32FMT") <"uint32FMT"-"uint32FMT"-"uint32FMT">%s\n",
            _exons[i]._estFrom, _exons[i]._estTo,
            _exons[i]._genFrom, _exons[i]._genTo,
            _exons[i]._numMatches, _exons[i]._numMatchesN, _exons[i]._percentIdentity, iOri);

    while (*outc)  outc++;
  }

  for (uint32 i=0; i<_numExons; i++) {
    if (_exons[i]._estAlignment) {
      strcpy(outc, _exons[i]._estAlignment);
      while (*outc)  outc++;
      *outc++ = '\n';
    }
    if (_exons[i]._genAlignment) {
      strcpy(outc, _exons[i]._genAlignment);
      while (*outc)  outc++;
      *outc++ = '\n';
    }
  }

  strcpy(outc, "sim4end\n");

  return(outs);
}

char *
sim4polish::s4p_polishToStringGFF3(void) {

  //  9 columns, tab separated
  //  tab, newline, cr and control MUST be escaped
  //  reserved letters:  ; = % & ,
  //  spaces ARE ALLOWED in fields
  //  undefined values should use '.'
  //
  //  1 seqid, genome name (a-zA-Z0-9.:^*$@!+_?-|), no whitespace (??) and not begin with >
  //  2 source ("sim4db")
  //  3 type ("mRNA" or "exon")
  //  4 begin, 1-based
  //  5 end, zero-length start=end, to the right of this base
  //  6 score (percent identity)
  //  7 strand
  //  8 phase
  //  9 attributes
  //      ID        (unique within scope of file)
  //      Name      (display name)
  //      Parent    ()
  //      Target
  //      Gap
  //      Derives_from
  //      Note
  //      Dbxref
  //      Ontology_term
  //      Is_circular
  //      others, user-defined (lowercase first letter; see below)
  //
  //   Example:
  //     0:arm_2L        sim4db  mRNA    2372455 2373234 98      -       .       ID=sim4db0;Name=61728:gb|CA807305;Target=61728:gb|CA807305 22 685 +;targetLen=685;pA=0;pT=21;genRegion=2370482-2375223
  //     0:arm_2L        sim4db  exon    2372455 2372770 99      -       .       Parent=sim4db0;Target=61728:gb|CA807305 22 337 +;Gap=M316;nMatches=313;intron=<-
  //     0:arm_2L        sim4db  exon    2372830 2373076 96      -       .       Parent=sim4db0;Target=61728:gb|CA807305 338 584 +;Gap=M74 D1 M2 I1 M170;nMatches=238;intron=<-
  //     0:arm_2L        sim4db  exon    2373134 2373234 99      -       .       Parent=sim4db0;Target=61728:gb|CA807305 585 685 +;Gap=M101;nMatches=100
  //

  //  Make a decent estimate of how much space we'll need to store the string
  //
  uint32 spaceNeeded = (1024 + 128 * _numExons +
                        ((_comment)    ? strlen(_comment)    : 0) +
                        ((_estDefLine) ? strlen(_estDefLine) : 0) +
                        ((_genDefLine) ? strlen(_genDefLine) : 0));

  for (uint32 i=0; i<_numExons; i++)
    if (_exons[i]._estAlignment)
      spaceNeeded += 2 * strlen(_exons[i]._estAlignment);

  char   *outs = new char [spaceNeeded];
  char   *outc = outs;

  //  Find extents of this match.
  uint32  estbgn = _exons[0]._estFrom;
  uint32  estend = _exons[_numExons-1]._estTo;
  uint32  genbgn = _exons[0]._genFrom;
  uint32  genend = _exons[_numExons-1]._genTo;

  for (uint32 i=0; i<_numExons; i++) {
    if (_exons[i]._genFrom < genbgn)  genbgn = _exons[i]._genFrom;
    if (_exons[i]._genTo   < genbgn)  genbgn = _exons[i]._genTo;
    if (genend < _exons[i]._genFrom)  genend = _exons[i]._genFrom;
    if (genend < _exons[i]._genTo)    genend = _exons[i]._genTo;

    if (_exons[i]._estFrom < estbgn)  estbgn = _exons[i]._estFrom;
    if (_exons[i]._estTo   < estbgn)  estbgn = _exons[i]._estTo;
    if (estend < _exons[i]._estFrom)  estend = _exons[i]._estFrom;
    if (estend < _exons[i]._estTo)    estend = _exons[i]._estTo;
  }

  //  Find the orientation
  char    mOri = '?';

  if (_matchOrientation == SIM4_MATCH_FORWARD)     mOri = '+';
  if (_matchOrientation == SIM4_MATCH_COMPLEMENT)  mOri = '-';

  // Find the strand
  char    sOri = '?';
  switch (_strandOrientation) {
    case SIM4_STRAND_POSITIVE:    sOri = '+';  break;
    case SIM4_STRAND_NEGATIVE:    sOri = '-';  break;
    case SIM4_STRAND_UNKNOWN:
    case SIM4_STRAND_INTRACTABLE:
    case SIM4_STRAND_FAILED:
    case SIM4_STRAND_ERROR:       sOri = '.';  break;
    default:
      fprintf(stderr, "sim4reader: Unknown strandOrientation '"uint32FMT"' in printPolishGFF3()\n", _matchOrientation);
      sOri = '.';
      break;
  }

  //  Get rid of spaces in the names (and do it non-destructively).

  uint32  estDefSpace = 0;
  uint32  genDefSpace = 0;

  while ((_estDefLine[estDefSpace]) && (isspace(_estDefLine[estDefSpace]) == 0))
    estDefSpace++;
  while ((_genDefLine[genDefSpace]) && (isspace(_genDefLine[genDefSpace]) == 0))
    genDefSpace++;

  char estDefChar = _estDefLine[estDefSpace];
  char genDefChar = _genDefLine[genDefSpace];

  _estDefLine[estDefSpace] = 0;
  _genDefLine[genDefSpace] = 0;

  //  The main mRNA match line.

  sprintf(outc, uint32FMT":%s\tsim4db\tmRNA\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t%c\t.\t",
          _genID, _genDefLine, genbgn, genend, _percentIdentity, sOri);
  while (*outc)  outc++;

  sprintf(outc, "ID=sim4db"uint32FMT";Name="uint32FMT":%s;Target="uint32FMT":%s "uint32FMT" "uint32FMT" %c;",
          sim4polishPolishID, _estID, _estDefLine, _estID, _estDefLine, estbgn, estend, mOri);
  while (*outc)  outc++;

  sprintf(outc, "targetLen="uint32FMT";pA="uint32FMT";pT="uint32FMT";genRegion="uint32FMT"-"uint32FMT"\n",
          _estLen, _estPolyA, _estPolyT, _genRegionOffset, _genRegionOffset + _genRegionLength -1);
  while (*outc)  outc++;

  //  Exons.

  for (uint32 i=0; i<_numExons; i++) {
    char *gap = encodeGap(_exons[i]._genAlignment, _exons[i]._estAlignment);

    sprintf(outc, uint32FMT":%s\tsim4db\texon\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t%c\t.\t",
            _genID, _genDefLine, _exons[i]._genFrom, _exons[i]._genTo, _exons[i]._percentIdentity, sOri);
    while (*outc)  outc++;

    if (gap)
      sprintf(outc, "Parent=sim4db"uint32FMT";Target="uint32FMT":%s "uint32FMT" "uint32FMT" %c;Gap=%s;nMatches="uint32FMT"",
              sim4polishPolishID, _estID, _estDefLine, _exons[i]._estFrom, _exons[i]._estTo, mOri, gap, _exons[i]._numMatches);
    else
      sprintf(outc, "Parent=sim4db"uint32FMT";Target="uint32FMT":%s "uint32FMT" "uint32FMT" %c;nMatches="uint32FMT"",
              sim4polishPolishID, _estID, _estDefLine, _exons[i]._estFrom, _exons[i]._estTo, mOri, _exons[i]._numMatches);
    while (*outc)  outc++;

    delete [] gap;

    switch (_exons[i]._intronOrientation) {
      // +1 to exclude the front blank space
      case SIM4_INTRON_POSITIVE:    sprintf(outc, ";intron=%s\n", iOriPOS +1); break;
      case SIM4_INTRON_NEGATIVE:    sprintf(outc, ";intron=%s\n", iOriNEG +1); break;
      case SIM4_INTRON_AMBIGUOUS:   sprintf(outc, ";intron=%s\n", iOriAMB +1); break;
      case SIM4_INTRON_GAP:         sprintf(outc, ";intron=%s\n", iOriGAP +1); break;
      case SIM4_INTRON_ERROR:       sprintf(outc, ";intron=%s\n", iOriERR +1); break;
      default:                      sprintf(outc, "\n"); break;
    }

    while (*outc)  outc++;
  }

  sim4polishPolishID++;

  _estDefLine[estDefSpace] = estDefChar;
  _genDefLine[genDefSpace] = genDefChar;

  return(outs);
}

char *
sim4polish::s4p_polishToStringATAC(void) {
  return(0L);
}


