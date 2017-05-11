#ifndef SUM_STATS_VECTOR_H
#define SUM_STATS_VECTOR_H

typedef struct
{
  int nsam;      /* number of total samples, from msDQH line */
  int segsites;  /* number of segregating sites, from line after // */
  char ** seqDat;  /* sequence data */
  int nsub;      /* nadv? default 0, -a can be used to set this value */
  int npops;        /* # of subpopulations, from -m or -D */
  int * n;          /* array of # seqs in each sub pops, npops elements, from -m or -D */
  double theta;     /* 4 Ne mu used for the sim, arg of -t from msDQH line,
		       mu is per locus, not per site, with -s it is
		       calculated by watterson's estimater from seg sites */
  int isNumSegSitesFixed;  /* 1: fixed numSegSites set by -s
		              O: theta is specified in the simulation
			         by -t (msbayes) */
  int Qbool;        /* 1: -Q used to specify ts/tv and base freq, 0: no -Q */
  int Fst_bool;     /* 1: multiple population sampled (npops > 0) */
  int replicateID;        /* simulation number (1 to howmany), was count */
  int numReplicates;       /* total number of repetitionreplication, (was howmany) */
  int taxonID;  /* serial number to identify taxon pair */
  int locusID;
  int taxonLocusID; /* sequential ID for each taxon:locus (1-NumPairs), was TAXAcount, from msDQH line */
  int NumTaxonLocusPairs; /* Total # of taxon:locus pairs, was NumTaxa, from msDQH line */
  int BasePairs;    /* Sequence length, from msDQH line */
} msOutput;

typedef struct
{
  msOutput *dat;
  int numElements;
  int allocatedSize;
  int numTaxonPairs;
  int numLoci;
  int numTaxonLocusPairs;
} msOutputArray;

struct SumStat
{  
	// 21 summary statistics + 2 imtermiate summary statistics + 2 spots for locsusID and taxonID
	double SS[25];
};

/* in msStatsDQH.c */
struct SumStat * CalcSumStats (msOutput *msOut);
void PrintSumStatsArray (struct SumStat **SumStat_list, int numTaxonLocusPairs, int numLoci, int numTaxonPairs);



/*
void printstats (struct SumStat **SumStat_list, int numTauClassesHyper, double *tau, int NumTaxa);
*/

#endif /* SUM_STATS_VECTOR_H */


