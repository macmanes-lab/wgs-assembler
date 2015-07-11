// This file is part of A2Amapper.
// Copyright (c) 2005, 2006 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bio++.H"
#include "atac.H"


atacMatchList::atacMatchList() {
  _matchesLen = 0;
  _matchesMax = 256;
  _matches    = new atacMatch [_matchesMax];
}

atacMatchList::~atacMatchList() {
  delete [] _matches;
}

void
atacMatchList::add(atacMatch &m) {

  if (_matchesLen >= _matchesMax) {
    _matchesMax <<= 2;
    atacMatch  *A = new atacMatch [_matchesMax];
    memcpy(A, _matches, sizeof(atacMatch) * _matchesLen);
    delete [] _matches;
    _matches = A;
  }

  memcpy(&_matches[_matchesLen], &m, sizeof(atacMatch));

  _matches[_matchesLen].matchiid = _matchesLen++;
}
