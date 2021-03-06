#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>


void
sim4polish::s4p_readPolishS4DB(readBuffer *rb) {

  //  Clear this polish.

  _numExons = 0;

  delete [] _comment;     _comment    = 0L;
  delete [] _estDefLine;  _estDefLine = 0L;
  delete [] _genDefLine;  _genDefLine = 0L;
  delete [] _exons;       _exons      = 0L;

  //  Decide the type of record we're reading.

  //  Read it.

  uint64    startPosition = rb->tell();

  uint64    thisLineMax  = 1048576;
  uint64    thisLineLen  = 0;
  char     *thisLine     = new char [thisLineMax];

  uint32    numLines = 10240;
  uint32    curLine  = 0;

  char    **lines   = new char * [numLines + 1];
  uint32   *lengths = new uint32 [numLines + 1];

  memset(lines,   0, sizeof(char *) * numLines);
  memset(lengths, 0, sizeof(uint32) * numLines);

  thisLineLen = rb->read(thisLine, thisLineMax, '\n');
  chompL(thisLine, thisLineLen);

  while (!rb->eof() && strcmp(thisLine, "sim4begin")) {
    fprintf(stderr, "sim4reader: Got '%s', expecting 'sim4begin' at byte "uint64FMT"\n",
            thisLine, startPosition);
    thisLineLen = rb->read(thisLine, thisLineMax, '\n');
    chompL(thisLine, thisLineLen);
  }

  //  Stash the 'sim4begin' line into the lines array.
  lines[curLine]   = new char [thisLineLen + 1];
  lengths[curLine] = thisLineLen;
  memcpy(lines[curLine++], thisLine, sizeof(char) * (thisLineLen + 1));

  //  Until we hit 'sim4end' stash lines into lines.  Yes, we test the previous line, then read the
  //  next.  At the end of the loop, we'll read 'sim4end', stash it in lines[], then test.

  while (!rb->eof() && strcmp(thisLine, "sim4end")) {
    thisLineLen = rb->read(thisLine, thisLineMax, '\n');
    chompL(thisLine, thisLineLen);

    if (curLine >= numLines) {
#warning LAZY PROGRAMMER did not extend an array
      fprintf(stderr, "ERROR: too many lines, lazy programmer.\n");
      exit(1);
    }

    //  Stash the line in the lines array.
    lines[curLine]   = new char [thisLineLen + 1];
    lengths[curLine] = thisLineLen;
    memcpy(lines[curLine++], thisLine, sizeof(char) * (thisLineLen + 1));
  }

  delete [] thisLine;

  if (numLines > 0)
    s4p_linesToPolishS4DB(startPosition, numLines, lines, lengths);

  for (uint32 i=0; i<curLine; i++)
    delete [] lines[i];

  delete [] lines;
  delete [] lengths;
}



void
sim4polish::s4p_readPolishGFF3(readBuffer *rb) {
  //  Clear this polish.

  _numExons = 0;

  delete [] _comment;     _comment    = 0L;
  delete [] _estDefLine;  _estDefLine = 0L;
  delete [] _genDefLine;  _genDefLine = 0L;
  delete [] _exons;       _exons      = 0L;

  //  Decide the type of record we're reading.

  //  Read it.
  uint64    startPosition = rb->tell();

  uint64    thisLineMax  = 1048576;
  uint64    thisLineLen  = 0;
  char     *thisLine     = new char [thisLineMax];

  uint32    numLines = 10240;
  uint32    curLine  = 0;

  bool      firstLine = true;

  char    **lines   = new char * [numLines + 1];
  uint32   *lengths = new uint32 [numLines + 1];

  memset(lines,   0, sizeof(char *) * numLines);
  memset(lengths, 0, sizeof(uint32) * numLines);

  thisLineLen = rb->read(thisLine, thisLineMax, '\n');
  chompL(thisLine, thisLineLen);

  while (!rb->eof() && (!strstr(thisLine, "\tsim4db\tmRNA") || (thisLine[0]=='#'))) {
    if (thisLine[0]!='#')
      fprintf(stderr, "sim4reader: Got '%s', expecting GFF3 mRNA line at byte "uint64FMT"\n",
              thisLine, startPosition);
    thisLineLen = rb->read(thisLine, thisLineMax, '\n');
    chompL(thisLine, thisLineLen);
  }

  //  Check the mRNA line (!), then stash into the lines array.
  lines[curLine]   = new char [thisLineLen + 1];
  lengths[curLine] = thisLineLen;
  memcpy(lines[curLine++], thisLine, sizeof(char) * (thisLineLen + 1));

  //  Read the GFF3 record, till the next mRNA line.
  //  We expect 'intron' on each exon line but the last; until we hit an intron-less line,
  //  stash lines into lines.  Yes, we test the previous line, then read the next.  
  //  At the end of the loop, we'll read the intron-less line, stash it in lines[], then test.

  while (!rb->eof() && (firstLine || strstr(thisLine, "\tsim4db\texon\t"))) {

    if ((firstLine == false) && !strstr(thisLine, "intron=")) break;

    thisLineLen = rb->read(thisLine, thisLineMax, '\n');
    chompL(thisLine, thisLineLen);

    if (curLine >= numLines) {
#warning LAZY PROGRAMMER did not extend an array
      fprintf(stderr, "ERROR: too many lines, lazy programmer.\n");
      exit(1);
    }

    //  If not a comment, stash the line in the lines array
    if (thisLine[0] == '#') continue;

    lines[curLine]   = new char [thisLineLen + 1];
    lengths[curLine] = thisLineLen;
    memcpy(lines[curLine++], thisLine, sizeof(char) * (thisLineLen + 1));

    firstLine = false;
  }

  delete [] thisLine;

  if (curLine > 0)
    s4p_linesToPolishGFF3(startPosition, numLines, lines, lengths);

  for (uint32 i=0; i<curLine; i++)
    delete [] lines[i];

  delete [] lines;
  delete [] lengths;
}



void
sim4polish::s4p_readPolishATAC(readBuffer *rb) {
}
