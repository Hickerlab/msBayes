/*
 * sumStatsVector.c
 *
 * Copyright (C) 2006  Michael Hickerson, Naoki Takebayashi
 *
 * This file is a part of sumstatsvector, distributed with msBayes.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

/*
  Change Log
  * Fri May 12 2006 Naoki Takebayashi <ffnt@uaf.edu>
  - Added ParseCommandLine().
  - tauequalizer of msprior and this program has to be synchronized.  So
    now this program takes upper bound of theta as an option (-T).
    This value enables us to calculate the tauequalizer.
    All of this calc is not needed any more.
  - msbayes.pl gets this value correctly from msprior, and pass it with -T.
  - nadv used to be a commandline argument.  But I moved it to be an option.
    To pass nadv, we use -a N option now.
  - It will print the usage message with -h option.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>		/* for getopt_long() */
#include <ctype.h>		/* for isspace() */
#include "msprior.h"		/* MAX_FILENAME_LEN */
#include "whiteSpaces.h"
#include "stringUtils.h"        /* cmatrix */
#include "sumStatsVector.h"

#define MAX_LNSZ 1000

/* posit[] was originally double, but I thougt int makes more sense, Naoki
 * Change the typedef to double if I'm wrong */
typedef unsigned int tPositionOfSegSites;


runParameters gParam;
int gNadv = 0;
int gPrintHeader = 0;

int sizeMemAllocated = 0;


/* Function prototypes of external functions */
int multiplepopssampledfrom (int nsam, int npops, int *config);
double thetaW (int, int), mu, N;
void PrintSumStatNames (void);

/* function prototypes for the functions in this file */
static int biggerlist (int nsam, unsigned nmax, char **list);

static int FindNumPopsAndSubpopSampleSizes (const char line[],
					    int **subPopSampleSize);
static void ReadInPositionOfSegSites (const char *line,
				      tPositionOfSegSites * positionArray,
				      int numSegSites);

static void ParseCommandLine (int argc, char *argv[]);
static void ReadInMSoutputAndCalculateSummaryStatistics (FILE * pfin);

int gSortPattern = 7;


int
main (int argc, char *argv[])
{
  /* check the command line argument */
  ParseCommandLine (argc, argv);

  /* read the data from STDIN */
  ReadInMSoutputAndCalculateSummaryStatistics(stdin);

  return (0);
}

static void
ReadInMSoutputAndCalculateSummaryStatistics (FILE * pfin){
  msOutputArray *msOutArr;
  msOutput *outPtr;
  struct SumStat **sumStatArr;
  struct SumStat **sumStatArrTemp;

  int initialArrSize;
  int msbayesFormat, taxonID, locusID;

  int maxsites = 1000; /* max number of seg sites, used for size of data mat */
  int nsam, i, howmany, npops, *config;
  char **list, line[MAX_LNSZ], longline[32000], *mutscanline;
  char dum[100];
  tPositionOfSegSites *posit;
  double theta;
  int segsites, count, numTaxonLocusPairs, BasePairs, taxonLocusID;
  int Fst_bool = 0, Qbool = 0;
  int isNumSegSitesConst = 0;	/* 1 with -s, the number of segregating sites
				 *    will be constant in each sample
				 * 0 with -t, varies between samples
				 */

  /* read in first line of output (command line options of msDQH) */
  fgets (line, MAX_LNSZ, pfin);

  if (! (msOutArr = malloc(sizeof(msOutputArray)))) {
    perror ("ERROR: not enough memory in ReadInMSoutput\n");
    exit(EXIT_FAILURE);
  }
  
  /* processing the header information */
  if (strcmp(line, "# BEGIN MSBAYES\n") == 0) { /* process info from msbayes.pl */
    fgets(line, MAX_LNSZ,pfin); /* get next line */
    sscanf(line, "# numTaxonLocusPairs %d numTaxonPairs %d numLoci %d",
	   &(msOutArr->numTaxonLocusPairs), &(msOutArr->numTaxonPairs),
	   &(msOutArr->numLoci));
    /* if taxon:locus matrix is required, we can process here  */

    fgets(line, MAX_LNSZ,pfin); /* get next line */
    msbayesFormat = 1;
    initialArrSize = 500;
  } else {
    msOutArr->numTaxonPairs=msOutArr->numTaxonLocusPairs = 1;
    msOutArr->numLoci = 1;
    msbayesFormat = 0;
    initialArrSize = 1;
  }

  /* allocate memory to msOutArr */
  if (! (msOutArr->dat = (msOutput *)calloc(initialArrSize, sizeof(msOutput)))) {
    perror ("ERROR: not enough memory 2 in ReadInMSoutput\n");
    exit(EXIT_FAILURE);
   }
  msOutArr->numElements = 0;
  msOutArr->allocatedSize = initialArrSize;

   /* go through the array of each msDQH run and get sum stats  */
  if ((sumStatArr = calloc(initialArrSize,
			    sizeof(struct SumStat *))) == NULL) {
    perror("ERROR: No mem in main ");
    exit(EXIT_FAILURE);
  }
				
   /* go through the array of each msDQH run and get sum stats  */
  if ((sumStatArrTemp = calloc(msOutArr->numTaxonLocusPairs,
			    sizeof(struct SumStat *))) == NULL) {
    perror("ERROR: No mem in main ");
    exit(EXIT_FAILURE);
  }
					
  int sumStatCounter = 0;					
  do {
    
    int endOfFile = 0;
    if (msbayesFormat) {
      while (strncmp(line, "# taxonID ", 10) != 0) {
	if (! fgets(line, MAX_LNSZ,pfin)) { /* basically skipping empty line */
	  endOfFile = 1;
	  break;
	}
      }
      if(endOfFile)
	break;
      /* Before msDQH output, numerial IDs for taxon and locus are inserted */
      sscanf(line, "# taxonID %d locusID %d\n", &taxonID, &locusID);
      fgets(line, MAX_LNSZ,pfin);
    } else {
      /* CHECK THIS WELL, NAOKI, SWRS */
      while (BlankCharStringQ(line)) {
	if (! fgets(line,MAX_LNSZ,pfin)) {
	  endOfFile=1;
	  break;
	}
      }
      if (endOfFile)
	break;
      taxonID = locusID = 1;
    }

    /*
     * Get the following variables from the command line options
     * NOTE: this line has to match with system() line of msbayes.pl
     *
     * nsam:      number of total samples
     * howmany:   how many simulations were run
     * THETA:     4 Ne mu used for the simulation
     *            This is removed, and getting this in more prper way
     * BasePairs: sequence length
     * taxonLocusID: sequential ID for each taxon:locus pair
                     (1 to # of taxon:locus pairs)
     * numTaxonLocusPairs: total number of taxon:locus pairs per 1 set of sims.
     */
    numTaxonLocusPairs = taxonLocusID = BasePairs = -1;
    sscanf (line,
	    " %s %s %d %d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %u %s %s %s %d %s %s %s %u ",
	    dum, dum, &nsam, &howmany, dum, dum, dum, dum, dum, dum, dum,
	    dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum,
	    dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum,
	    dum, dum, dum, dum, dum, dum, dum, dum, &BasePairs, dum, dum,
	    dum, &taxonLocusID, dum, dum, dum, &numTaxonLocusPairs);
	
    if(!msbayesFormat) {
      msOutArr->numTaxonPairs=msOutArr->numTaxonLocusPairs = 1;
      taxonLocusID=numTaxonLocusPairs=1;
    }
    /*
     * of course, I have to put a prior generator in the actual sample
     * generator for theta and tau down below for each count
     */

    /* Find the theta or number of segregating sites from -t or -s */
    mutscanline = strstr (line, "-s");
    if (mutscanline != NULL)
      {
	/* number of segregating sites is constant */
	sscanf (mutscanline, " %d", &segsites);
	isNumSegSitesConst = 1;
	theta = thetaW (nsam, segsites);
      }
    else
      {
	mutscanline = strstr (line, "-t");
	if (mutscanline != NULL)
	  sscanf (mutscanline, "-t %s", dum);
	else
	  {
	    fprintf (stderr, "\nmutscanline problem -s or -t not found \n");
	    exit (1);
	  }
	theta = atof (dum);

	/* -Q will tell transition transversion rate ratio and base freqs */
	if ((mutscanline = strstr (line, "-Q")) != NULL) {
	  Qbool = 1;
	}
      }

    mutscanline = strstr(line, "-r");
    if(mutscanline != NULL) {
      sscanf(mutscanline, "-r %s %d", dum, &BasePairs);
      // dum contains recomb rate
    }

    /*
     * config become an array with npops elements,
     * it contains subpop sample sizes
     */
    npops = FindNumPopsAndSubpopSampleSizes (line, &config);

    if (npops == 1) {
      config[0] = nsam;
    }

    /* Checking if 0 < config[i] < nsam for all i */
    if ((npops > 1) && (multiplepopssampledfrom (nsam, npops, config)))
      Fst_bool = 1;

    /* prepare the storage for segregating sites data */
    if (isNumSegSitesConst)
      maxsites = segsites;

    list = cmatrix (nsam, maxsites + 1);

    posit =
      (tPositionOfSegSites *) calloc (maxsites, sizeof (tPositionOfSegSites));

    if (list == NULL)
      {
	    fprintf (stderr, "No mem for segregating sites data, couldn't allocated memory for list\n");
	    exit (EXIT_FAILURE);
      }

	if (posit == NULL)
      {
	    fprintf (stderr, "No mem for segregating sites data, couldn't allocated memory for posit\n");
	    exit (EXIT_FAILURE);
      }

    /* Start to process the data */
    count = 0;
    while (howmany - count++)
      {
	  
	/* The line after "//" is the beginning of simulation data */
	while (strcmp (line, "//\n") != 0)
	  fgets (line, MAX_LNSZ, pfin);

	/* Number of segregating sites line */
	fgets (line, MAX_LNSZ, pfin);
	if (!isNumSegSitesConst)
	  {
	    sscanf (line, "segsites: %d\n", &segsites);

	    if (segsites >= maxsites)	/* readjust the size of data matrix */
	      {
		maxsites = segsites + 10;	/* extra 10 elements */
		posit = (tPositionOfSegSites *)
		  realloc (posit, maxsites * sizeof (tPositionOfSegSites));
		/*printf("PRE %d %d %d\n", segsites, maxsites, nsam); */
		if (posit == NULL || biggerlist (nsam, maxsites, list) != 0)
		  {
		    fprintf (stderr,
			     "Not enough memory for reallocating char matrix\n");
		    exit (EXIT_FAILURE);
		  }
	      }
	  }

	/* get rid of base frequency line */
	if (Qbool)
	  {
	    fgets (line, MAX_LNSZ, pfin);
	    sscanf (line, "freqACGT: %s %s %s %s", dum, dum, dum, dum);
	  }

	if (segsites > 0)
	  {
	    /* read in position of segregating sites */
	    fgets (longline, 32000, pfin);

	    /* posit array initialized */
	    ReadInPositionOfSegSites (longline, posit, segsites);

	    /* list[][] get initialized with character states */
	    for (i = 0; i < nsam; i++)
	      fscanf (pfin, " %s", list[i]);

	  }
	/* what do we do if segsites = 0?, Naoki */

	/* insert the data into the array */
	msOutArr->numElements ++;
	if (msOutArr->numElements > msOutArr->allocatedSize) {
	  /* reallocate the memory */
	  msOutArr->allocatedSize += 1000;
	  msOutArr->dat = realloc(msOutArr->dat,
				  sizeof(msOutput) * (msOutArr->allocatedSize));
	  if(msOutArr->dat ==NULL) {
	    perror("Realloc of msOutArr->dat failed\n");
	    exit(EXIT_FAILURE);
	  }
	  
	  sumStatArr = realloc(sumStatArr,
				  msOutArr->allocatedSize * sizeof(struct SumStat *));
	  if(sumStatArr ==NULL) {
	    perror("Realloc of sumStatArr failed\n");
	    exit(EXIT_FAILURE);
	  }
	}
	outPtr = & msOutArr->dat[msOutArr->numElements-1];

	outPtr->nsam = nsam;
	outPtr->segsites = segsites;
	outPtr->seqDat = list;
	outPtr->nsub = gNadv;
	outPtr->npops = npops;
	outPtr->n = config;
	outPtr->theta = theta;
	outPtr->isNumSegSitesFixed = isNumSegSitesConst;
	outPtr->Qbool = Qbool;
	outPtr->Fst_bool = Fst_bool;
	outPtr->replicateID = count;
	outPtr->numReplicates = howmany;
	outPtr->taxonID = taxonID;
	outPtr->locusID = locusID;
	outPtr->taxonLocusID = taxonLocusID;
	outPtr->NumTaxonLocusPairs = numTaxonLocusPairs;
	outPtr->BasePairs = BasePairs;
	
	sumStatArr[sumStatCounter++] = CalcSumStats (outPtr);

	}
	 	 
    freeCMatrix (nsam, list);
	free(posit);
	
  } while (fgets(line, MAX_LNSZ,pfin));

 
  int k, j;
  for (k = 0; k < msOutArr->numElements; k++) {
    j = k % msOutArr->numTaxonLocusPairs;
    sumStatArrTemp[j] = sumStatArr[k];

    if (j == msOutArr->numTaxonLocusPairs - 1) {
      // we got sumStats for 1 set 
      PrintSumStatsArray(sumStatArrTemp, msOutArr->numTaxonLocusPairs, msOutArr->numLoci, msOutArr->numTaxonPairs);
    }
  }
  
  free(msOutArr->dat); 
  free(msOutArr);
  free(sumStatArr);
 
}


/*
 * Take a character string line[] as the argument.  Find -m (migration
 * rate) or -D (demographic history?) from this line[], and extract
 * the information about sub-population sample sizes.
 *
 * The line[] is something like this:
 *
./msDQH 15 1 -t 22.052669 -Q 4.765400 0.207500 0.231800 0.215000 0.345700 -H 999.000000 -r 137.043900 639 -D 5 2 2 13 0 I 0.000000 1.723360 1.723360 0.276640 0.276640 2.522721 2 1 0 0 1 0 I 0.000000 Nc 0.352105 0.673049 4.170096 1 Nc 0.215405 1.1 1 Nc 0.215405 639 1 Nc 0.215405 2
 *
 *
 * Returned value: numSubPops
 *     The number of subpopulations (1st number after these options)
 *
 * Side effect:
 *   Mem is allocated to the pointer *subpopSampleSize (numSubPops elements),
 *   and this array contains the sample sizes of sub-populations.
 *
 */
static int
FindNumPopsAndSubpopSampleSizes (const char line[], int **subPopSampleSize)
{
  char *mscanline, *Dscanline, *charPtr;
  int numSubPops, i;
  int *arrayPtr;

  mscanline = strstr (line, "-m");
  Dscanline = strstr (line, "-D");

  if ((mscanline == NULL) && (Dscanline == NULL))
    {
      *subPopSampleSize = (int *) malloc (sizeof (int));
      if (*subPopSampleSize == NULL)
	{
	  fprintf (stderr,
		   "No memory in FindNumPopsAndSubpopSampleSizes()\n");
	  exit (EXIT_FAILURE);
	}

      arrayPtr = (int *) calloc (1, sizeof (int));
      arrayPtr[0] = -1;
      *subPopSampleSize = arrayPtr;
      return 1;			/* no -D nor -m, so only 1 population */
    }

  /* set the char pointer to the beginning of arguments for the option */
  if (mscanline != NULL)
    {
      if (Dscanline != NULL)  /* Both -m & -D specified, can be a problem, Naoki */
	fprintf (stderr, "WARN: both -m and -D are specified, ignoring -D\n");

      charPtr = mscanline;

    }
  else if (Dscanline != NULL)
    {
      charPtr = Dscanline;
    }

  /* the 2nd argument after -D contains the number of sub populations. */
  for (i = 0; i < 2; i++)
    {
      charPtr = FindFirstSpace (charPtr);
      if (!charPtr)
	{
	  fprintf (stderr, "ERROR: no arguments are given after -m or -D");
	  exit (EXIT_FAILURE);
	}
      charPtr = RmLeadingSpaces (charPtr);
    }
  numSubPops = (int) strtol (charPtr, &charPtr, 10);

  /* subPopSampleSize is an array storing the sub-population sample sizes */
  arrayPtr = (int *) calloc (numSubPops, sizeof (int));
  /* This sucks up the subpop sample sizes into arrayPtr[] */
  for (i = 0; i < numSubPops; i++)
    {
      arrayPtr[i] = (int) strtol (charPtr, &charPtr, 10);

      if (errno == ERANGE)
	{
	  fprintf (stderr, "WARN: out of range in population sizes. "
		   "Continuing... But the results are probably wrong\n");
	}
      if (errno == EINVAL)
	{
	  fprintf (stderr, "ERROR: The arguments to -m options is weird:\n"
		   "%s\n", line);
	  exit (EXIT_FAILURE);
	}
    }

  *subPopSampleSize = arrayPtr;	/* returning the address of the array */
  return numSubPops;
}


/*
 * Process a line containing positions of segregating sites, and
 * assign the values to an array
 *
 * Arguments:
 *   line: A character string which contains the position of segregating sites
 *         This line starts with "positions:" and contains integers
 *         delimited by spaces.  There should be numSegSites integers.
 *
 *   positionArray: The values extracted from line will be assigned to this
 *                  array.  Assumes correct size mem is allocated
 *   numSegSites: number of segregating sites.  The size of positionArray
 *                should be equal to (or greater than) this value.
 *
 * Returns: nothing
 *
 * Side effects: values will be assigned to positionArray
 *
 */
static void
ReadInPositionOfSegSites (const char *line,
			  tPositionOfSegSites * positionArray,
			  int numSegSites)
{
  int i;
  char *charPtr, *tempPtr;

  charPtr = strstr (line, "positions:");
  if (charPtr == NULL)
    {
      fprintf (stderr, "ERROR: positions line not found, ignoring\n");
      exit (EXIT_FAILURE);
    }

  charPtr = FindFirstSpace (charPtr);
  if (charPtr == NULL)
    {
      fprintf (stderr, "ERROR: positions line doesn't contain a space\n");
      exit (EXIT_FAILURE);
    }

  for (i = 0; i < numSegSites; i++)
    {
      positionArray[i] = (int) strtol (charPtr, &tempPtr, 10);
      if (charPtr == tempPtr)
	{
	  fprintf (stderr, "ERROR: reached to the end, while processing "
		   "positions line\n");
	  exit (EXIT_FAILURE);
	}
      charPtr = tempPtr;
    }
  return;
}

/*****************  Character matrix functions **********************/

/*
 * Arguments:
 *   list[][]: a matrix storing segregating sites data, allocated by cmatrix()
 *   nsam:     number of rows in list[][]
 *   nmax:     number of columns
 *
 * Side effect:
 *   grow the size (character string length) to nmax.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int
biggerlist (int nsam, unsigned nmax, char **list)
{
  int i;

  for (i = 0; i < nsam; i++)
    {
      list[i] = (char *) realloc (list[i], nmax * sizeof (char));
      if (list[i] == NULL)
	{
	  perror ("realloc error. in biggerlist()");
	  return -1;
	}
    }
  return 0;
}

/* Macro to move string from s to d */
#define strMove(d,s) memmove(d,s,strlen(s)+1)

/*
 * Print out the usage, and exit
 */
static void
PrintUsage (char *progname)
{
  char *p;

  /*
   * Strip the path from the program name for when we
   * print the usage.
   */
  p = strrchr (progname, '/');
  if (p)
    p++;
  else
    p = progname;

  fprintf (stderr,
	   "\nUsage: %s [--help] [--header] [--name] [--sort N] [--nadv N] [< output_line_of_msDQH]\n\n"
	   "        help: Print this usage function (-h)\n"
	   "      header: Print column header (-H)\n"
	   "        name: Print names of available summary statistics (-n)\n"
	   "        sort: specify the sorting pattern (-s), the choices of N are:\n"
	   "               0: no sort\n"
	   "               1: simple sort (default)\n"
	   "               2: group by locus then sort by pi.b within each locus(?)\n"
	   "               3: group by locus then sort by the average of pi.b, can use up to 4 moments\n"
	   "               4: group by taxon then sort by the average of pi.b, use first 4 moments\n"
	   "               5: group by taxon then sort by the average of pi.b, use first 3 moments\n"
	   "               6: group by taxon then sort by the average of pi.b, use first 2 moments\n"
	   "               7: group by taxon then sort by the average of pi.b, use the first moment\n"
	   "        nadv: Specify nadv (-a)\n"
	   "stdin is used to read in a single line of msDQH output "
	   "(output_line_of_msDQH)" "\n\n", p);
  exit (EXIT_FAILURE);
}


static struct option sim_opts[] = {
  {"help", 0, NULL, 'h'},	/* list options */
  {"header", 0, NULL, 'H'},
  {"name", 0, NULL, 'n'},
  {"sort", 7, NULL, 's'},
  {"nadv", 1, NULL, 'a'},
  {NULL, 0, NULL, 0}
};

static void
ParseCommandLine (int argc, char *argv[])
{
  while (1)
    {
      int opt;
      int opt_index;

      opt = getopt_long (argc, argv, "hHns:a:", sim_opts, &opt_index);
      if (opt < 0)
	break;

      switch (opt)
	{
	case 'h':		/* Print usage and exit */
	  PrintUsage (argv[0]);	/* This function will exit */
	  break;
	case 'n':
	  PrintSumStatNames ();
	  break;
	case 'H':		/* Print header */
	  gPrintHeader = 1;
	  break;
	case 's':
	  if (!optarg)
	    {
	      fprintf (stderr, "Must select sort Pattern with --sort option\n");
	      PrintUsage (argv[0]);
	    }
	  gSortPattern = strtol(optarg, NULL, 10);
	  if (errno || (gSortPattern < 0) || (gSortPattern > 7))
	    {
	      fprintf (stderr, "Invalid sort: %s\n", optarg);
	      PrintUsage (argv[0]);
	    }
	  break;
	case 'a':		/* specify nadv */
	  if (!optarg)
	    {
	      fprintf (stderr, "Must select nadv value " "of theta\n");
	      PrintUsage (argv[0]);
	    }
	  gNadv = strtol (optarg, NULL, 10);
	  if (errno || (gNadv < 0))
	    {
	      fprintf (stderr, "Invalid nadv: %s\n", optarg);
	      PrintUsage (argv[0]);
	    }
	  break;
	default:
	  PrintUsage (argv[0]);	/* This function will exit */
	  break;
	}
    }
}
