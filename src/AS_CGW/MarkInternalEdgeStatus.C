
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
static char *rcsid = "$Id: MarkInternalEdgeStatus.C 4371 2013-08-01 17:19:47Z brianwalenz $";

#include "AS_global.H"
#include "AS_UTL_Var.H"
#include "AS_UTL_interval.H"
#include "AS_CGW_dataTypes.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"
#include "ScaffoldGraphIterator_CGW.H"
#include "ChiSquareTest_CGW.H"


#undef VERBOSE_MARKING


static
void
findOrientAndDistance(PairOrient     &edgeOrient,
                      NodeCGW_T      *CI, SequenceOrient &CIorient,
                      NodeCGW_T      *MI, SequenceOrient &MIorient,
                      LengthT        &gapDistance) {

  assert(edgeOrient.isUnknown() == false);

  if (edgeOrient.isAB_BA()) {
    //      CI                                        MI
    //  A --------------------- B               B --------------------- A
    //    5'----->                                           <------5'
    CIorient.setIsForward();
    MIorient.setIsReverse();
    gapDistance.mean     = MI->offsetBEnd.mean     - CI->offsetBEnd.mean;
    gapDistance.variance = MI->offsetBEnd.variance - CI->offsetBEnd.variance;
  }

  if (edgeOrient.isAB_AB()) {
    //      CI                                        MI
    //  A --------------------- B               A --------------------- B
    //    5'----->                                           <------5'
    CIorient.setIsForward();
    MIorient.setIsForward();
    gapDistance.mean     = MI->offsetAEnd.mean     - CI->offsetBEnd.mean;
    gapDistance.variance = MI->offsetAEnd.variance - CI->offsetBEnd.variance;
  }

  if (edgeOrient.isBA_BA()) {
    //      CI                                        MI
    //  B --------------------- A               B --------------------- A
    //    5'----->                                           <------5'
    CIorient.setIsReverse();
    MIorient.setIsReverse();
    gapDistance.mean     = MI->offsetBEnd.mean     - CI->offsetAEnd.mean;
    gapDistance.variance = MI->offsetBEnd.variance - CI->offsetAEnd.variance;
  }

  if (edgeOrient.isBA_AB()) {
    //      CI                                        MI
    //  B --------------------- A               A --------------------- B
    //    5'----->                                           <------5'
    CIorient.setIsReverse();
    MIorient.setIsForward();
    gapDistance.mean     = MI->offsetAEnd.mean     - CI->offsetAEnd.mean;
    gapDistance.variance = MI->offsetAEnd.variance - CI->offsetAEnd.variance;
  }

  //  This condition should not occur but when it does it causes an assert in the PairwiseChiSquare
  //  subroutine so we will just mark the edge as untrusted so the program can keep going but this
  //  needs to be investigated and fixed!!!!  The PairwiseChiSquare test is conditional on variance
  //  being positive.

  if (gapDistance.variance <= 0)
    fprintf(stderr, "WARNING:  Negative variance (%.2f +- %.2f) from edgeOrient %c CI %d %c (%.2f +- %.2f -- %.2f +- %.2f) to MI %d %c (%.2f +- %.2f -- %.2f +- %.2f)\n",
            gapDistance.mean, gapDistance.variance,
            edgeOrient.toLetter(),
            CI->id,
            CIorient.toLetter(),
            CI->offsetAEnd.mean, CI->offsetAEnd.variance,
            CI->offsetBEnd.mean, CI->offsetBEnd.variance,
            MI->id,
            MIorient.toLetter(),
            MI->offsetAEnd.mean, MI->offsetAEnd.variance,
            MI->offsetBEnd.mean, MI->offsetBEnd.variance);
  //assert(gapDistance.variance > 0);
}




void
MarkInternalEdgeStatus(ScaffoldGraphT  *graph,
                       CIScaffoldT     *scaffold,
                       int              beLoose,
                       int              operateOnMerged,
                       double           pairwiseChiSquaredThreshold,
                       double           maxVariance) {

  CIScaffoldTIterator   CIs;
  NodeCGW_T            *CI;
  NodeCGW_T            *MI;
  int32                 numCIs;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  for (numCIs=0; (CI = NextCIScaffoldTIterator(&CIs)) != NULL; numCIs++)
    CI->indexInScaffold = numCIs;

  assert(numCIs == scaffold->info.Scaffold.numElements);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    GraphEdgeIterator     edges(ScaffoldGraph->ContigGraph, CI->id, ALL_END, ALL_EDGES);
    EdgeCGW_T            *edge;

    while ((edge = (operateOnMerged) ? edges.nextMerged() : edges.nextRaw()) != NULL) {
      MI = GetGraphNode(ScaffoldGraph->ContigGraph, (edge->idA == CI->id) ? edge->idB : edge->idA);

      //  Only label same-scaffold edges.
      if (MI->scaffoldID != CI->scaffoldID) {
#ifdef VERBOSE_MARKING
        //PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "Interscaffold", edge, -1);
#endif
        SetEdgeStatus(graph->ContigGraph, edge, INTER_SCAFFOLD_EDGE_STATUS);
        continue;
      }

      //  Only label canonical edges.  Skipping this marks all edges as bad orient (because we see
      //  the same edge twice, and one of those must be mis-oriented).
      if (MI->indexInScaffold <= CI->indexInScaffold) {
#ifdef VERBOSE_MARKING
        //PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "NonCanonical ", edge, -1);
#endif
        continue;
      }

#ifdef VERBOSE_MARKING
      PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "Examining", edge, -1);
#endif

      //  Check that the edge orientation is consistent with the CI positions
      //  and orientations within the scaffold.
      PairOrient     edgeOrient = GetEdgeOrientationWRT(edge, CI->id);
      SequenceOrient CIorient;
      SequenceOrient MIorient;
      LengthT        gapDistance = {0,0};

      findOrientAndDistance(edgeOrient, CI, CIorient, MI, MIorient, gapDistance);

      //  Throw out edges across bad variance gaps.
      //if (gapDistance.variance <= 0.0) {
      //  SetEdgeStatus(graph->ContigGraph, edge, UNTRUSTED_EDGE_STATUS);
      //  continue;
      //}

      //  Mark as untrusted an edge whose orientation does not agree
      //  with the orientation of the CIs in the scaffold.
      if ((beLoose < 3) &&
          ((GetNodeOrient(CI) != CIorient) ||
           (GetNodeOrient(MI) != MIorient))) {
#ifdef VERBOSE_MARKING
        //fprintf(stderr, "["F_CID","F_CID"]Bad orientation (%c,%c) (%c,%c)\n",
        //        CI->id, MI->id,
        //        GetNodeOrient(CI).toLetter(), CIorient.toLetter(), GetNodeOrient(MI).toLetter(), MIorient.toLetter());
        PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "BadOrient    ", edge, -1);
#endif
        SetEdgeStatus(graph->ContigGraph, edge, UNTRUSTED_EDGE_STATUS);
        continue;
      }

      //  Mark this edge as untrusted if the distance of the edge is not consistent with the
      //  estimated gap distance as judged by the Chi Squared Test.
      //
      //  If the gapDistance.variance is negative, this fails.  Previous to version 1.2 we would
      //  mark all edges across this gap as untrusted -- but doesn't that automatically disconnect
      //  the scaffold?!

      if ((beLoose < 2) &&
          (gapDistance.variance > 0)) {
        double  chiSquareResult;

        if (FALSE == PairwiseChiSquare(gapDistance.mean,
                                       gapDistance.variance,
                                       edge->distance.mean,
                                       edge->distance.variance,
                                       NULL,
                                       &chiSquareResult,
                                       pairwiseChiSquaredThreshold)) {
#ifdef VERBOSE_MARKING
          //fprintf(stderr, "              X=%f gap=(%f,%f) edge=(%f,%f)\n",
          //      chiSquareResult,
          //        gapDistance.mean, gapDistance.variance,
          //        edge->distance.mean, edge->distance.variance);
          //PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "BadChiSq     ", edge, -1);
#endif
          SetEdgeStatus(graph->ContigGraph, edge, UNTRUSTED_EDGE_STATUS);
          continue;
        }
      }

      //  If not raw, check for large variance.
      if ((beLoose < 1) &&
          (edge->edgesContributing > 1) &&
          (edge->distance.variance > maxVariance)) {
#ifdef VERBOSE_MARKING
        PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "LargeVariance", edge, -1);
#endif
        SetEdgeStatus(graph->ContigGraph,edge, LARGE_VARIANCE_EDGE_STATUS);
        continue;
      }

#ifdef VERBOSE_MARKING
      PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "TRUSTED      ", edge, -1);
#endif
      SetEdgeStatus(graph->ContigGraph,edge, TRUSTED_EDGE_STATUS);
    }
  }

  uint32  nInvalid  = 0;
  uint32  nUnknown  = 0;
  uint32  nUntrust  = 0;
  uint32  nUntrustT = 0;
  uint32  nTrustT   = 0;
  uint32  nTrust    = 0;
  uint32  nLargeVar = 0;
  uint32  nInterScf = 0;
  uint32  nError    = 0;
  uint32  nNonCanon = 0;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
    GraphEdgeIterator     edges(ScaffoldGraph->ContigGraph, CI->id, ALL_END, ALL_EDGES);
    EdgeCGW_T            *edge;

    while ((edge = (operateOnMerged) ? edges.nextMerged() : edges.nextRaw()) != NULL) {
      MI  = GetGraphNode(ScaffoldGraph->ContigGraph, (edge->idA == CI->id) ? edge->idB : edge->idA);

      if (MI->indexInScaffold <= CI->indexInScaffold) {
        nNonCanon++;
        continue;
      }

      switch (GetEdgeStatus(edge)) {
        case INVALID_EDGE_STATUS:
          nInvalid++;
          break;
        case UNKNOWN_EDGE_STATUS:
          nUnknown++;
          break;
        case UNTRUSTED_EDGE_STATUS:
          nUntrust++;
          break;
        case TENTATIVE_UNTRUSTED_EDGE_STATUS:
          nUntrustT++;
          break;
        case TENTATIVE_TRUSTED_EDGE_STATUS:
          nTrustT++;
          break;
        case TRUSTED_EDGE_STATUS :
          nTrust++;
          break;
        case LARGE_VARIANCE_EDGE_STATUS:
          nLargeVar++;
          break;
        case INTER_SCAFFOLD_EDGE_STATUS:
          nInterScf++;
          break;
        default:
          nError++;
          break;
      }
    }
  }

  //fprintf(stderr, "#Edges: invalid %u unknown %u untrusted %u (tentative %u) trusted %u (tentative %u) largevar %u interscf %u noncanon %u error %u\n",
  //        nInvalid, nUnknown, nUntrust, nUntrustT, nTrust, nTrustT, nLargeVar, nInterScf, nNonCanon, nError);

#warning is internalEdges and confirmedInternalEdges used
  scaffold->info.Scaffold.internalEdges          = nTrust + nTrustT + nUntrustT;
  scaffold->info.Scaffold.confirmedInternalEdges = nTrust + nTrustT;
}
