
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
/*************************************************
* Module:  OlapFromSeeds.h
* Description:
*   Declarations for olap-from-seeds program
*
*    Programmer:  A. Delcher
*       Started:   15 February 2007
*
* Assumptions:
*
* Notes:
*
*************************************************/

/* RCS info
 * $Id: OlapFromSeedsOVL.H 4600 2015-02-06 05:56:41Z brianwalenz $
 * $Revision: 1.27 $
*/


#ifndef  __OLAPFROMSEEDS_H_INCLUDED
#define  __OLAPFROMSEEDS_H_INCLUDED

static const char *rcsid_OLAPFROMSEEDS_H_INCLUDED = "$Id: OlapFromSeedsOVL.H 4600 2015-02-06 05:56:41Z brianwalenz $";

//**ALD determine if use new code to analyze true multialignments
#define  USE_NEW_STUFF  1

#include  "AS_global.H"
#include  "AS_UTL_reverseComplement.H"
#include  "AS_PER_gkpStore.H"
#include  "AS_OVS_overlapStore.H"
#include  "SharedOVL.H"

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>
#include  <pthread.h>

//  Constants

#define  DEFAULT_CORRECTION_FILENAME  "frag.cor"
  //  Default name of file where corrections are sent
#define  DEFAULT_DEGREE_THRESHOLD    2
  //  Default value for  Degree_Threshold
#define  DEFAULT_END_EXCLUDE_LEN     3
  //  Default value for  End_Exclude_Len
#define  DEFAULT_KMER_LEN            9
  //  Default value for  Kmer_Len
#define  DEFAULT_NUM_PTHREADS        2
  //  Default number of pthreads to use
#define  DEFAULT_VOTE_QUALIFY_LEN    9
  //  Default value for bases surrounding SNP to vote for change
#define  DISPLAY_WIDTH               60
  //  Number of characters per line in alignments
#define  EDIT_DIST_PROB_BOUND        1e-4
  //  Probability limit to "band" edit-distance calculation
  //  Determines  NORMAL_DISTRIB_THOLD
#define  ERATE_BITS                  16
  //  Number of bits to store integer versions of error rates
#define  ERRORS_FOR_FREE             1
  //  The number of errors that are ignored in setting probability
  //  bound for terminating alignment extensions in edit distance
  //  calculations
#define  EXPANSION_FACTOR            1.4
  //  Factor by which to grow memory in olap array when reading it
#define  FRAGS_PER_BATCH             100000
  //  Number of old fragments to read into memory-based fragment
  //  store at a time for processing
#define  HOMOPOLY_LEN_VARIATION      3
  //  Homopolymer run-length differences greater than this
  //  can be errors for correlated difference detection
#define  MAX_DEGREE                  32767
  //  Highest number of votes before overflow
#define  MAX_LINE                    1000
  //  Longest possible input line
#define  MAX_VOTE                    255
  //  Highest number of votes before overflow
#define  MIN_HAPLO_OCCURS            3
  //  This many or more votes at the same base indicate
  //  a separate haplotype
#define  NORMAL_DISTRIB_THOLD        3.62
  //  Determined by  EDIT_DIST_PROB_BOUND
#define  THREAD_STACKSIZE        (128 * 512 * 512)
  //  The amount of memory to allocate for the stack of each thread

#define USE_STORE_DIRECTLY_READ
  //  Use the store directly during the initial load -- we aren't doing
  //  random access, just streaming through and loading.  This lets us
  //  load about 2x the frags.

#define USE_STORE_DIRECTLY_STREAM
  //  Use the store directly during the stream -- good if you don't have
  //  lots of frags loaded.

//#define USE_STREAM_FOR_EXTRACT
  //  When loading frags during the stream, use a fragStream instead
  //  of random access to the store.  Useless unless USE_STORE_DIRECTLY_STREAM
  //  is enabled.


//  Type definitions

typedef struct
  {
   short int  sub;    // subscript of aligned sequence
   short int  diff_ct;
   unsigned  disregard : 1;
   char * seq;
  }  Difference_Signature_t;

typedef enum
  {HOMOPOLY, STANDARD}  Vote_Kind_t;

typedef struct
  {
   unsigned char  hp_char_ct [5];  // a,c,g,t,gap respectively for homopoly reads
   unsigned char  nonhp_char_ct [5];  // a,c,g,t,gap respectively for standard reads
   unsigned int  homopoly_ct : 16;
   unsigned int  non_hp_ct : 16;
   unsigned int  mutab : 1;  //  This was called 'mutable', but c++ sucked that up from us
   double  homopoly_sum;
   double  non_hp_sum;
  }  New_Vote_t;

typedef  struct
  {
   unsigned int  confirmed : 8;
   unsigned int  deletes : 8;
   unsigned int  a_subst : 8;
   unsigned int  c_subst : 8;
   unsigned int  g_subst : 8;
   unsigned int  t_subst : 8;
   unsigned int  no_insert : 8;
   unsigned int  a_insert : 8;
   unsigned int  c_insert : 8;
   unsigned int  g_insert : 8;
   unsigned int  t_insert : 8;
   // homopoly values are associated with the first base of a
   // consecutive run of the same base (i.e., a homopolymer run)
   // homopoly_ct is the number of reads with a matching run
   // homopoly_sum is the sum of the lengths of those runs
   unsigned int  homopoly_ct : 8;
   uint32  homopoly_sum;
  }  Vote_Tally_t;

typedef  struct
  {
   int  frag_sub;
   int  align_sub;
   Vote_Value_t  vote_val;
  }  Vote_t;

typedef  struct
  {
   char  * sequence;
   Vote_Tally_t  * vote;
#if USE_NEW_STUFF
   Sequence_Diff_t  * diff_list;
   uint32  num_diffs;
#endif
   unsigned  clear_len        : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  len              : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  trim_5p          : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  trim_3p          : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  left_degree      : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  shredded         : 1;    // True if shredded read
   unsigned  right_degree     : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  is_homopoly_type : 1;
  }  Frag_Info_t;

const int  INNIE = 0;
const int  NORMAL = 1;

typedef  struct
  {
   AS_IID a_iid, b_iid;
   signed int  a_hang    : AS_OVS_HNGBITS;
   signed int  b_hang    : AS_OVS_HNGBITS;
   signed int  orient    : 2;
   unsigned int k_count  : 8;
  }  Olap_Info_t;

typedef  struct
  {
   AS_IID id;
   unsigned  trim_5p          : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  trim_3p          : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  len              : AS_READ_MAX_NORMAL_LEN_BITS;
   unsigned  shredded         : 1;
   unsigned  is_homopoly_type : 1;
   int  start;              // position of beginning of sequence in  buffer
  }  Frag_List_Entry_t;

typedef  struct
  {
   Frag_List_Entry_t  * entry;
   char  * buffer;
   int  size, ct, buffer_size;
  }  Frag_List_t;

typedef  struct
  {
   int  thread_id;
   AS_IID  lo_frag, hi_frag;
   int  next_olap;
   int  failed_olaps;
   gkStream  * frag_stream;
   gkFragment     frag_read;
   Frag_List_t  * frag_list;
   char  rev_seq [AS_READ_MAX_NORMAL_LEN + 1];
   int  rev_id;
   int  ** edit_array;
   Homopoly_Match_Entry_t  ** homopoly_edit_array;
   HP_LV_Cell_t  ** hp_cell;
   unsigned  * can_look;
   int  * edit_space;
#if USE_NEW_STUFF
   Diff_Entry_t  diff_list [AS_READ_MAX_NORMAL_LEN+1];  //  only MAX_ERRORS needed
   Alignment_Cell_t  * banded_space;
#endif
  }  Thread_Work_Area_t;

typedef enum
  {TEXT_FILE, BINARY_FILE}  OVL_Output_t;



//  Static Globals

static int  Asymmetric_Olaps = FALSE;
  // If set true by the -a option, then only output overlaps
  // where a_iid < b_iid
static BinaryOverlapFile  * Binary_OVL_Output_fp = NULL;
  // Pointer for binary overlap outputs
static double  Char_Match_Value = 0.0;
  // Positive score for matching characters in computing alignments
  // Score for mismatching characters is (this value) - 1.0
  // This is set at runtime.
static int  Check_Correlated_Differences = FALSE;
  // If set true by the -w option then check for confirmed differences
  // in columns to eliminate overlaps
static char  * Correction_Filename = DEFAULT_CORRECTION_FILENAME;
  // Name of file to which correction information is sent
static int  Degree_Threshold = DEFAULT_DEGREE_THRESHOLD;
  // Set keep flag on end of fragment if number of olaps < this value
static int  Doing_Corrections = TRUE;
  // Determines whether error-corrections of reads will be
  // computed and output
static int  Doing_Partial_Overlaps = FALSE;
  // If set true by the G option (G for Granger)
  // then allow overlaps that do not extend to the end
  // of either read.
static int  * Edit_Array [AS_READ_MAX_NORMAL_LEN+1];  //  only MAX_ERRORS needed
  // Use for alignment calculation.  Points into  Edit_Space .
static int  Edit_Match_Limit [AS_READ_MAX_NORMAL_LEN+1] = {0};  //  only MAX_ERRORS needed
  // This array [e] is the minimum value of  Edit_Array [e] [d]
  // to be worth pursuing in edit-distance computations between guides
static int *Edit_Space = NULL;
  // Memory used by alignment calculation
static int  End_Exclude_Len = DEFAULT_END_EXCLUDE_LEN;
  // Length of ends of exact-match regions not used in preventing
  // sequence correction
static int  Error_Bound [AS_READ_MAX_NORMAL_LEN + 1];
  // This array [i]  is the maximum number of errors allowed
  // in a match between sequences of length  i , which is
  // i * MAXERROR_RATE .
static double  Error_Rate = 0.0;
  // Highest allowed error rate used in alignments
  // Can be set at run time.
static int  Extend_Fragments = FALSE;
  // If true, try to extend clear range of fragments.
  // Set by  -e  option
static int  Failed_Olaps = 0;
  // Counts overlaps that didn't make the error bound
static Frag_Info_t  * Frag;
  // Sequence and vote information for current range of fragments
  // being corrected
static Frag_List_t  Frag_List;
  // List of ids and sequences of fragments with overlaps to fragments
  // in  Frag .  Allows simultaneous access by threads.
static gkStore  *gkpStore;
  // Fragment store from which fragments are loaded
static gkStream  *Frag_Stream;
  // Stream to extract fragments from internal store
static char  * gkpStore_Path;
  // Name of directory containing fragment store from which to get fragments
static unsigned  X_Arg = 0;
  // Set to argument of special X option
static AS_IID Hi_Frag_IID;
  // Internal ID of last fragment in frag store to process
static gkStore  * Internal_gkpStore;
  // Holds partial frag store to be processed simultanously by
  // multiple threads
static int  Kmer_Len = DEFAULT_KMER_LEN;
  // Length of minimum exact match in overlap to confirm base pairs
static AS_IID Lo_Frag_IID;
  // Internal ID of first fragment in frag store to process
static int  Min_Olap_Len = AS_OVERLAP_MIN_LEN;
  // The minimum number of bp in each read to report an overlap
static time_t  Now;
  // Used to get current time
static uint32 Num_Frags;
  // Number of fragments being corrected
static uint64 Num_Olaps;
  // Number of overlaps being used
static int  Num_PThreads = DEFAULT_NUM_PTHREADS;
  // Number of pthreads to process overlaps/corrections;
static int  Offsets_WRT_Raw = TRUE;
  // Indicates if offset information of seeds is relative to the
  // entire raw read, or just the clear range
static Olap_Info_t  * Olap = NULL;
  // Array of overlaps being used
static FILE  * OVL_Output_fp = NULL;
  // File where output overlaps should be sent
static char  * OVL_Output_Path = NULL;
  // Name of file or store to which overlaps should be sent
static OVL_Output_t  OVL_Output_Type = TEXT_FILE;
  // Type of output for overlaps
static int  Seeds_From_Store = FALSE;
  // Indicates if overlap info comes from  get-olaps  or from
  // a binary overlap store
static char  * Olap_Path;
  // Name of file containing a sorted list of overlaps
static pthread_mutex_t  Print_Mutex;
  // To make debugging printout come out together
static int  Use_Haplo_Ct = TRUE;
  // Set false by  -h  option to ignore haplotype counts
  // when correcting
static int  Vote_Qualify_Len = DEFAULT_VOTE_QUALIFY_LEN;
  // Number of bases surrounding a SNP to vote for change



//  Static Functions

static void  Add_Homopoly_Vote
  (Vote_Tally_t * vote, int hp_len);
static void  Adjust_For_Compressed_Seeds
  (int * a_offset, int * b_offset, char * * a_part, char * * b_part);
static void  Adjust_For_Eliminated
  (char * ref, int * ref_len, Sequence_Diff_t * dp, int dp_ct);
static void  Analyze_Alignment
  (int delta [], int delta_len, char * a_part, char * b_part,
   int a_len, int b_len, int a_offset, int sub);
static void  Analyze_Frag
    (FILE * fp, int sub);
static void  Analyze_Diffs
    (void);
static int  Binomial_Bound
  (int, double, int, double);
static int  By_A_Lo
  (const void * a, const void * b);
static int  By_B_IID
  (const void * a, const void * b);
static int  By_Diffs_And_Alpha
  (const void * a, const void * b);
static void  Cast_Confirmation_Vote
  (Vote_Tally_t * vote);
static void  Cast_Delete_Vote
  (Vote_Tally_t * vote);
static void  Cast_Insert_Vote
  (Vote_Tally_t * vote, unsigned ch);
static void  Cast_New_Vote_Char
  (New_Vote_t * vp, char ch, Vote_Kind_t type);
static void  Cast_New_Vote_Code
  (New_Vote_t * vp, unsigned code, Vote_Kind_t type);
static void  Cast_No_Insert_Vote
  (Vote_Tally_t * vote);
static void  Cast_Substitution_Vote
  (Vote_Tally_t * vote, unsigned ch);
static void  Cast_Vote
  (Vote_Value_t val, int p, int sub);
static int  Char_To_Code
  (char ch);
static int  Char_Matches
  (char ch, unsigned code);
static void  Compute_Delta
  (int delta [], int * delta_len, int ** edit_array,
   int e, int d, int row);
static void  Convert_Delta_To_Diff
  (int delta [], int delta_len, char * a_part, char * b_part,
   int  a_len, int b_len, Sequence_Diff_t * diff, int errors,
   Thread_Work_Area_t * wa);
static void  Count_From_Diff
  (short int count [AS_READ_MAX_NORMAL_LEN] [5], const char * seq, int seq_len,
   int seq_is_homopoly, const Sequence_Diff_t * dp);
static void  Determine_Homopoly_Corrections
    (FILE * fp, int sub, New_Vote_t * vote, char * seq,
     char * correct, int len);
static void  Determine_Standard_Corrections
    (FILE * fp, int sub, New_Vote_t * vote, char * seq,
     char * correct, int len);
static void  Display_Alignment
  (char * a, int a_len, char * b, int b_len, int delta [], int delta_ct,
   int capitalize_start);
static void  Display_Diffs
  (FILE * fp, const Sequence_Diff_t * dp, const char * a, int a_len);
static void  Display_Frags
  (void);
static void  Display_Multialignment
  (FILE * fp, int sub, const char * ref, const char * anno, int ref_len,
   const Sequence_Diff_t * dp, int dp_ct, int skip_disregard);
static void  Display_Partial_Diff
  (FILE * fp, const Sequence_Diff_t * dp, int lo, int hi, const char * a);
static int  Eliminate_Correlated_Diff_Olaps
  (int sub, const char * ref, int ref_len, Sequence_Diff_t * dp, int dp_ct,
   int ref_is_homopoly);
static void  Extract_Needed_Frags
  (gkStore *store, AS_IID lo_frag, AS_IID hi_frag,
   Frag_List_t * list, int * next_olap);
static char  Filter
  (char ch);
static void  Get_Seeds_From_Store
  (char * path, int32 lo_id, int32 hi_id, Olap_Info_t ** olap, int * num);
static char  Homopoly_Should_Be
  (char curr, New_Vote_t * vp, int * ch_ct, int * tot);
static void  Init_Frag_List
  (Frag_List_t * list);
static void  Initialize_Globals
  (void);
static void  Init_Thread_Work_Area
  (Thread_Work_Area_t * wa, int id);
static char *  Insert_Gaps
  (char * seq, const short unsigned insert_size [], int n);
static int  Is_Homopoly_Type
  (gkFragment * fr, gkStore * gkp);
static Vote_Value_t  Matching_Vote
  (char ch);
static void  Modify_For_Inserts
  (Sequence_Diff_t * mod_dp, Sequence_Diff_t * dp, const short unsigned insert_size []);
static void  Output_Correction_Header
  (FILE * fp, int sub);
static void  Output_Corrections
  (FILE * fp);
static void  Output_Olap
  (Olap_Info_t * olap, int a_lo, int a_hi, int a_len,
   int b_lo, int b_hi, int b_len, int errors);
static void  Output_Olap_From_Diff
  (const Sequence_Diff_t * dp, int32 a_iid);
static void  Parse_Command_Line
  (int argc, char * argv []);
static void  Process_Seed
  (Olap_Info_t * olap, char * b_seq, unsigned b_len, char * rev_seq, int * rev_id,
   int shredded, int is_homopoly, Thread_Work_Area_t * wa);
static void  Read_Frags
  (void);
static void  Read_Seeds
  (void);
static void  Set_Diff_Entry
  (Diff_Entry_t ** de, int * pos, int * size, unsigned len, unsigned action,
   unsigned ch);
static void  Set_Homopoly_Votes_From_Diffs
  (int sub, Sequence_Diff_t * dp);
static void  Set_Insert_Sizes
  (short unsigned * insert_size, const Sequence_Diff_t * dp, const char * seq,
   int seq_len);
static void  Set_New_Homopoly_Votes
  (New_Vote_t * vote, const char * seq, int len, const Sequence_Diff_t * dp);
static void  Set_New_Self_Homopoly_Votes
  (New_Vote_t * vote, const char * seq, int len);
static void  Set_New_Self_Votes
  (New_Vote_t * vote, const char * seq, int len);
static void  Set_New_Standard_Votes
  (New_Vote_t * vote, const char * seq, int len, const Sequence_Diff_t * dp);
static void  Set_Self_Homopoly_Votes
  (int sub, int frag_len);
static void  Set_Self_Votes
  (int sub, int frag_len);
static void  Set_Signature
  (Difference_Signature_t * signature, short int diff_col [], int num_diff_cols,
   const Sequence_Diff_t * dp, const char * seq, int seq_len, int seq_is_homopoly);
static void  Set_Votes_From_Diffs
  (int sub, Sequence_Diff_t * dp);
static void  Show_Corrections
  (FILE * fp, const char * seq, const char * corr, int len);
static void  Show_Edit_Array
  (FILE * fp, int ** ea, int errs);
static void  Show_Frag_Votes
  (FILE * fp, int sub);
static void  Show_New_Votes
  (FILE * fp, const char * seq, const New_Vote_t * vp, int n);
static void  Show_Votes
  (FILE * fp);
static char  Standard_Should_Be
  (char curr, New_Vote_t * vp, int * ch_ct, int * tot);
static void  Stream_Old_Frags
  (void);
void *  Threaded_Process_Stream
  (void * ptr);
static void  Threaded_Stream_Old_Frags
  (void);
static void  Tidy_Up
  (void);
static void  Usage
  (char * command);
static int  Votes_For
  (char ch, const New_Vote_t * vp);


#endif
