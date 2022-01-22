#ifndef MS_PRIOR_H
#define MS_PRIOR_H

/* Default values used for interactive setup */
#define DEFAULT_LOWER_THETA 0.00000000004 /* mu = 10^(-11), N = 1 */
#define DEFAULT_UPPER_THETA 0.04 /* mu = 10^(-8), N=10^6 */
#define DEFAULT_UPPER_TAU  1.0
#define DEFAULT_UPPER_MIG  0.0	/* when set 0, it is a constant */
#define DEFAULT_UPPER_REC  0.0
#define DEFAULT_UPPER_ANC_POPSIZE 0.25
#define DEFAULT_REPS  1000000
#define DEFAULT_MUT_FILE "SampleSize_MuModel_Vector"
#define DEFAULT_SUBPARAMCONSTRAIN "000000000"

#define MAX_FILENAME_LEN 1024
#define MAX_NAME_CHAR_LEN 1024
#define NUMBER_OF_CONPARAM 9

#include "hashtab.h"
#include "array_utils.h"
#include "partition_combinatorics_random.h"

typedef struct
{
  double upperTheta;		/* upper limit of prior dist'n for theta */
  double lowerTheta;		/* lower limit of prior dist'n for theta */
  double upperTau;		/* upper limit of prior dist'n for time of separation */
  double bufferTauClasses;	/* Whether or not to buffer distance between tau values */
  double concentrationShape;  /* Concentration parameter for dpp for numTauClasses */
  double concentrationScale;  /* Concentration parameter for dpp for numTauClasses */
  double upperMig;		/* upper limit of prior dist'n for migration rate */
  double upperRec;		/* upper limit of prior dist'n for recombination rate */
  double upperAncPopSize;	/* upper limit of prior dist'n for ancestral pop size */
  unsigned long long reps;
  unsigned int numTaxonLocusPairs; /* total number of taxon:locus pairs */
  unsigned int numTaxonPairs;     /* number of unique taxon pairs */ 
  unsigned int numLoci;         /* number of unique loci */
  unsigned int numTauClasses;
  long prngSeed;
  char configFile[MAX_FILENAME_LEN];
  char scratchFile[MAX_FILENAME_LEN];
  unsigned int constrain;
  int printConf;                /* 1: -i opt print the config and exit (default: 0) */
  char subParamConstrain[NUMBER_OF_CONPARAM + 1]; /* +1 for string terminating null */
} runParameters;

extern runParameters gParam;

typedef struct
{
  char taxonName[MAX_NAME_CHAR_LEN];
  char locusName[MAX_NAME_CHAR_LEN];
  unsigned int taxonID;  /* integer (0-) representation of taxon pair name */
  unsigned int locusID;  /* integer (0-) representation of locus name */
  /* double thetaScaler; */
  double NScaler;     /* scaler of pop size, e.g. 0.25 for mtDNA of dioecy */
  double mutScaler;   /* mutation rate multiplier, relative to nuclear */
  unsigned int numPerTaxa;
  unsigned int sample[2];
  double tstv[2];
  double gamma;
  unsigned int seqLen;
  double freqA;
  double freqC;
  double freqG;
  double freqT;
  char filename[MAX_FILENAME_LEN];
} mutParameter;

typedef struct {
  int **tbl;
  unsigned int numTaxon;  /* number of rows (not 0-offset index) */
  unsigned int numLoci;   /* number of columns */
} lociTbl;

/* keeps track of sample size mutation model settings */
typedef struct
{
  mutParameter *data;
  int numElements;
  int numAllocated;
  int numSpecies;
  int numLoci;
  hashtab_t *taxonIDTbl;  /* taxon name string to -> unique ID hash table */
  hashtab_t *locusIDTbl;  /* locus name string to -> unique ID hash table */
  lociTbl *locTbl;
} mutParameterArray;

extern mutParameterArray gMutParam;

typedef struct
{
  double conTau;
  double conBottPop1;
  double conBottPop2;
  double conBottleTime;
  double conMig;		// migration rate
  double conTheta;
  double conN1;			//current population size1
  double conNanc;		// ancestral pupulation size
  double conRec;		// recombination rate
} constrainedParameter;

typedef struct
{
  constrainedParameter *conData;
  int conNumElements;
  int conNumAllocated;
} constrainedParameterArray;

extern constrainedParameterArray gConParam;

#endif
