/*
 * msStatsDQH.c
 *
 * Copyright (C) 2006  Richard Hudson, Eli Stahl,
 *                     Michael Hickerson, Naoki Takebayashi
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

#include "hashtab.h"
#include "msprior.h"  /* I'm not sure why this is included 03/31/09 */
#include "sumStatsVector.h"
#include <string.h>

/* #include <float.h> */
/* #if defined (__SVR4) && defined (__sun) */
/*   int isinf(double x) { return !finite(x) && x==x; } */
/* #endif */

#define MAX_LEN_COLUMN_NAME 128	/* Used for header. This is the maximum char length
				   of names for each column */
#define HIGHEST_MOMENT 3

extern int gSortPattern;

double nucdiv (int, int, char **);
double nucdiv_w (int, int, char **, int, int *, int);
double nucdiv_bw (int, int, char **, int, int *);
double tajddenominator (int, double);

double thetaW (int, int);
double thetah (int, int, char **);
void FuLi (double *D, double *F, int, int, char **, double pi);
void FrequencyDistrInfSites (int *freqdist, int nsam, int segsites,
			     char **list);
void FrequencyDistrQ (int *freqdist, int nsam, int segsites, char **list);
int segsub (int nsub, int segsites, char **list);
void segsubs (int *segwithin, int segsites, char **list, int npops, int *n);
int multiplepopssampledfrom (int nsam, int npops, int *n);	/*zzz n[] is called config[] in samples.c zzz */
static int SS_comp (const void *, const void *);
#if 0
static int compare_doubles (const void *a, const void *b);
#endif
int pwDiffAtOneSite (int, int, char **);
int cntBase (char, int, int, char **);


void shannonIndex (char **list, int *config, double **shannonIndexArray);
int charCount (char *arr);
double S_w (int ss, int nsam, char **list, int np, int *n, int pop, double D_w, int basePairs);
double S_xy (int ss, int nsam, char **list, int np, int *n, int pop1,int pop2, double D_xy, int basePairs);
double JWakely_Psi(int* n, double Sx, double Sy, double Sxy, double dx, double dy, double k);
void getUniqueItems(int *Array, int *uniqueItemArray, int totalSize, int uniqueSize);
void getSumStat(struct SumStat** base, double *data, int ssIndex, int size);
double moment(double *data, int size, int order);


extern int gPrintHeader;	/* boolean 1 means print header (column names), 0 = no header
				   -H option invoke printing of the header */
static void PrintHeader (char priorNames[][MAX_LEN_COLUMN_NAME],
			 int numPriors,
			 char sumStatNames[][MAX_LEN_COLUMN_NAME],
			 int numSumStats, int numTaxonPairs);
static void PrintHeader_fourMoments2 (char priorNames[][MAX_LEN_COLUMN_NAME], int numPriors,
	     char sumStatNames[][MAX_LEN_COLUMN_NAME], int numSumStats,
	     int numLoci);
		 
static void PrintHeader_fourMoments (char priorNames[][MAX_LEN_COLUMN_NAME], int numPriors,
	     char sumStatNames[][MAX_LEN_COLUMN_NAME], int numSumStats,
	     int numTaxonPairs);

/***** MAKE SURE the following two lines are up-to-date *****/
int numSumStats = 21;
char ssNameVect[][MAX_LEN_COLUMN_NAME] =
  { "pi.b", "pi.w", "pi", "wattTheta", "pi.net", "tajD", "tajD.denom",
  "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1",
  "tajD.denomPop2", "tajD.denomPop1", "ShannonsIndex.Between",
  "ShannonsIndex.Net", "ShannonsIndex.Pop1", "ShannonsIndex.Pop2",
  "S_x", "S_y", "S_xy","JWakeley_Psi"
};

enum SS //summmary statistics enum, the order is the same as in ssNameVect
{
  taxonID = 23,
  locusID = 24,
  PI_b = 0,
  PI_w,	// 1
  PI,	// 2
  TW,	// 3
  PI_Net,	// 4
  TD,	// 5
  TDD,	//6
  PI_w2,	// 7
  PI_w1,	// 8
  TW2,	// 9
  TW1,	// 10
  TDD2,	// 11
  TDD1,	// 12
  si1,	// 13 ?
  si2,	// 14
  si3,	// 15
  si4,	// 16
  Sx,	// 17
  Sy,	// 18
  Sxy,	// 19
  JW_Psi,	// 20
  // we are using summary statistics from pi_b to JW_PSI
  TD1,	// 21
  TD2	// 22
};


/* Print out the available summary stats and Exit */
void
PrintSumStatNames (void)
{
  int i;
  for (i = 0; i < numSumStats; i++)
    {
      printf ("%s\n", ssNameVect[i]);
    }
  exit (0);
}

/*
  Takes a pointer to structure which contains the results of SINGLE run of
  msDQH and calculates the summary statistics.

  Side effect: memory is allocated to store all SumStats,
  Returns: a pointer to the SumStats which contains the results

  nsam: number of samples
  segsites: number of segregating sites
  list: character matrix (A,T,G,orC) containing nsam rows and segsites columns
  nsub: gNadv, default 0, but can be specified by --nadv option
  npops: number of sub-populations
  n: a vector of sub-population sizes, there are npops elements
  theta: 4 Ne mu used for the simulations, it comes from the command line
         option (-t) to msDQH
*/
struct SumStat *
CalcSumStats (msOutput *msOut)
{
  int nsam, npops, *n, BasePairs, segsites;
  int *segwithin;
  double  FuLiD, FuLiF, Fst;
  // int i, tW_w_npops;
  // double tW_w;
  // int *freqdist;
  // double Nm;

  double  Dx,Dy, Dxy, K;
  double sx, sy, sxy, jw_psi;

  /* double h, th, ObsvVARD, ObsvVARPi_Net, ObsvEPi_Net, ObsvCV_pi_net_tW; */
  char **list;

  /* struct SumStat SumStat_list[NumTaxa]; */
  struct SumStat * resultSS;

  /* copying frequently used values for easy access */
  nsam = msOut->nsam;
  npops = msOut->npops;
  n = msOut->n;
  BasePairs = msOut->BasePairs;
  segsites = msOut->segsites;
  list = msOut->seqDat;

  /* SIDE EFFECT: allocating memory to return results */
  if (!(resultSS = malloc (sizeof(struct SumStat)))) {
    perror("ERROR: No Mem in CalcSumStats\n");
    exit(EXIT_FAILURE);
  }

  //initialize PI_b, PI_Net, PI_w2, PI_w1 and PI_w with -1 
  (resultSS->SS)[PI_b] = (resultSS->SS)[PI_Net] = (resultSS->SS)[PI_w2] = (resultSS->SS)[PI_w1] = (resultSS->SS)[PI_w] = -1;
 
  //initialize Sx, Sy, Sxy and JW_Psi with -1
  (resultSS->SS)[Sx] = (resultSS->SS)[Sy] = (resultSS->SS)[Sxy] = (resultSS->SS)[JW_Psi] = -1;

#if 0
  /* Not used */
  if (! (freqdist = (int *) malloc (nsam * sizeof (int)))) {
	perror("ERROR: No Mem 2 in CalcSumStats\n");
	exit(EXIT_FAILURE);
  }
  if (msOut->Qbool)
    FrequencyDistrQ (freqdist, nsam, segsites, list);
  else
    FrequencyDistrInfSites (freqdist, nsam, segsites, list);
#endif

  /* this is actually already taken care of, so this is redundant */
  if (msOut->isNumSegSitesFixed)
    (resultSS->SS)[TW] = msOut->theta; 
  else
    (resultSS->SS)[TW] = thetaW (nsam, segsites); 

  (resultSS->SS)[PI] = nucdiv (nsam, segsites, list); // pi per gene

  /* Tajima's D denominator PER SITE */
  (resultSS->SS)[TDD] = tajddenominator (nsam, (double) segsites/BasePairs); 
  /* Tajima's D PER GENE */
  (resultSS->SS)[TD] = ((resultSS->SS)[PI] - (resultSS->SS)[TW]) / tajddenominator (nsam, (double) segsites);   
  /* converting Pi to PER SITE, This should occur after TajD calc */
  (resultSS->SS)[PI] = (resultSS->SS)[PI] / BasePairs;  
  
  if (msOut->Fst_bool) /* stats for the simulation with more than 1 subpop */
    {
      if (!(segwithin = (int *) malloc (npops * sizeof (int)))) {
	perror("ERROR: No Mem 3 in CalcSumStats\n");
	exit(EXIT_FAILURE);
      }

      segsubs (segwithin, segsites, list, npops, n);

#if 0
      /* not used */
      tW_w = 0.;
      tW_w_npops = 0;
      for (i = 0; i < npops; i++)
	if (n[i] > 1)
	  {
	    tW_w += thetaW (n[i], segwithin[i]);
	    tW_w_npops++;
	  }
      tW_w /= tW_w_npops;  /* watterson's theta averaged over subpops */
#endif

      (resultSS->SS)[TW1] =  thetaW (n[0], segwithin[0]); 
	  (resultSS->SS)[TW2] = thetaW (n[1], segwithin[1]); 
      (resultSS->SS)[PI_w1] =  nucdiv_w (nsam, segsites, list, npops, n, 0); 
      (resultSS->SS)[PI_w2] =  nucdiv_w (nsam, segsites, list, npops, n, 1); 
#if 0
      (resultSS->SS)[PI_w] = 
	nucdiv_w (nsam, segsites, list, npops, n, -1) ; 
                         /* -1 signals Average of pi's within subpop */
#endif
      /* This is equivalent to the commented out PI_w above, but faster.
       * This basically pi_within for each sub pop, and taking the
       * weighted averages.  I confirmed that this gives the identical
       * results as above.
       */
      (resultSS->SS)[PI_w] = ((resultSS->SS)[PI_w1] * n[0] * (n[0] - 1)  +
			(resultSS->SS)[PI_w2] * n[1] * (n[1] - 1))/
	( n[0] * (n[0] - 1) +  n[1] * (n[1] - 1)); 

      (resultSS->SS)[PI_b] = nucdiv_bw (nsam, segsites, list, npops, n);

      /* Tajima's D denominator PER SITE for each sub pop */
      (resultSS->SS)[TDD1] = tajddenominator (n[0], (double) segwithin[0]/BasePairs);
      (resultSS->SS)[TDD2] = tajddenominator (n[1], (double) segwithin[1]/BasePairs);
      
      /* Tajima's D PER GENE for each sub pop*/
      (resultSS->SS)[TD1] = ((resultSS->SS)[PI_w1] - (resultSS->SS)[TW1]) / 
	tajddenominator (n[0], (double) segwithin[0]);
      (resultSS->SS)[TD2] = ((resultSS->SS)[PI_w2] - (resultSS->SS)[TW2]) / 
	tajddenominator (n[1], (double) segwithin[1]);

      /* converting Pi to PER SITE , This should occur after TajD, but
	 before jw_psi calc */
      (resultSS->SS)[PI_w] = (resultSS->SS)[PI_w] / BasePairs;
      (resultSS->SS)[PI_w1] = (resultSS->SS)[PI_w1] / BasePairs;
      (resultSS->SS)[PI_w2] = (resultSS->SS)[PI_w2] / BasePairs;
      (resultSS->SS)[PI_b] = (resultSS->SS)[PI_b] / BasePairs;

      /* expected values of pairwise differences */
      /* All of PI are PER SITE below */
      Dx = (resultSS->SS)[PI_w1];  //dx  (Pop1)
      Dy = (resultSS->SS)[PI_w2];  //dy  (Pop2)
      Dxy = (resultSS->SS)[PI_b];  //dxy (between)
      K = (resultSS->SS)[PI];      //k   (overall)

      /* These are coming from Wakely 1996a,b Thoretical Population
	 Biology 49:39-57 and 49:369-386 */
      sx = S_w(segsites, nsam, list, npops, n, 0, Dx, BasePairs);
      sy = S_w(segsites, nsam, list, npops, n, 1, Dy, BasePairs);
      sxy = S_xy(segsites, nsam, list, npops, n, 0,1,Dxy, BasePairs);
      /* expression (11) of Wakely 1996b */
      jw_psi = JWakely_Psi(n, sx, sy, sxy, Dx, Dy, K); 

      (resultSS->SS)[Sx] = sx;
	  (resultSS->SS)[Sy] = sy;
	  (resultSS->SS)[Sxy] = sxy;
      (resultSS->SS)[JW_Psi] = jw_psi;

      (resultSS->SS)[PI_Net] = (resultSS->SS)[PI_b] - (resultSS->SS)[PI_w]; /* PER SITE */

#if 0
      /* not used */
      Fst = 1. - (resultSS->SS)[PI_w] / (resultSS->SS)[PI_b];
      if (Fst < 0)
	  {
	    Fst = 0.;
	    // Nm = -1.;  /* Nm doesn't seems to be used any where */
	  }
      else
	  {
	    // Nm = (1. / Fst - 1.) / 4.;
	  }
#endif
    }
  /*   th = thetah(nsam, segsites, list) ; */
  
  
/*  FuLi(&FuLiD,&FuLiF,nsam,segsites,list,resultSS->PI);*/
/*   h = resultSS->PI-th ; */
#if 0
  /* not used NT */
  if (msOut->nsub > 0)
    nsegsub = segsub (msOut->nsub, segsites, list);
#endif

  /* Convert watterson's theta to PER SITE */
  (resultSS->SS)[TW] = (resultSS->SS)[TW] / BasePairs;
  if (msOut->Fst_bool) {
    (resultSS->SS)[TW1] = (resultSS->SS)[TW1] / BasePairs;
    (resultSS->SS)[TW2] = (resultSS->SS)[TW2] / BasePairs;
  }

  if (segsites < 1) {
    (resultSS->SS)[PI_b] =  (resultSS->SS)[TD] = (resultSS->SS)[PI] = (resultSS->SS)[TW] = 
      (resultSS->SS)[TDD] = FuLiD = FuLiF = 0;
    if (msOut->Fst_bool) {
      (resultSS->SS)[PI_Net] = (resultSS->SS)[TD1] = (resultSS->SS)[TD2] = 
	(resultSS->SS)[PI_w] = (resultSS->SS)[PI_w1] = (resultSS->SS)[PI_w2] = 
	(resultSS->SS)[TW1] = (resultSS->SS)[TW2] = (resultSS->SS)[TDD1] = (resultSS->SS)[TDD2] = 
	Fst = 0;
    }
  }
  if (segwithin[0] < 1)
    (resultSS->SS)[TD1] = 0;
  if (segwithin[1] < 1)
    (resultSS->SS)[TD2] = 0;
  if ((resultSS->SS)[PI_Net] < 0)
    (resultSS->SS)[PI_Net] = 0;
  if ((resultSS->SS)[PI_b] < 0)
    (resultSS->SS)[PI_b] = 0;

  /* the following Shannon Index probably should be done after testing
   Fst_bool Naoki Nov 2, 2009 */
  double *shannonIndexArray;

  shannonIndexArray = (double *) malloc (4 * sizeof (double));
  shannonIndex (list, n, &shannonIndexArray);
  (resultSS->SS)[si1] = shannonIndexArray[0];
  (resultSS->SS)[si2] = shannonIndexArray[1];
  (resultSS->SS)[si3] = shannonIndexArray[2];
  (resultSS->SS)[si4] = shannonIndexArray[3];

  (resultSS->SS)[taxonID] = msOut->taxonID;
  (resultSS->SS)[locusID] = msOut->locusID;
    
  free (shannonIndexArray);
  free(segwithin);
#if 0
  free(freqdist);
#endif

  return resultSS;
}

#if 0				/* commented out since it is not used */
/*
 * used for qsort to compare two doubles.
 * Takes two pointers to double as the arguments.
 *
 * Returns: 1  if a > b
 *          0  if a == b
 *         -1  if a < b
 */

static int
compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}
#endif

/*
 * Used for qsort to sort the array of SumStat
 *   Uses PI_b as the criteria of sorting
 * Returns: 1  if p1 > p2
 *          0  if p1 == p2
 *         -1  if p1 < p2
 */
static int
SS_comp (const void *p1, const void *p2)
{
  struct SumStat **sp1 = (struct SumStat **) p1;
  struct SumStat **sp2 = (struct SumStat **) p2;

  return ((((*sp1)->SS)[PI_b]) > (((*sp2)->SS)[PI_b])) - ((((*sp1)->SS)[PI_b]) < (((*sp2)->SS)[PI_b]));
}

/*
 * Used for qsort to sort the array of SumStat_moments
 *   Uses PI_b_mean as the criteria of sorting
 * Returns: 1  if p1 > p2
 *          0  if p1 == p2
 *         -1  if p1 < p2
 */
static int
SS_comp_moments (const void *p1, const void *p2)
{
  double **sp1 = (double **) p1;
  double **sp2 = (double **) p2;
  
  // *sp1[0][PI_b]
  return (((*sp1)[PI_b]) > ((*sp2)[PI_b])) - (((*sp1)[PI_b]) < ((*sp2)[PI_b]));
}

/*
  This function gets unique taxonID from all taxaLoci by looking at each taxonLocus's taxonID
*/
void getUniqueItems(int *Array, int *uniqueItemArray, int totalSize, int UniqueSize)
{

   int a, b, c = 0, notFound = 1;
   // fill unqiueLoci with dummy value
   for (b= 0; b < UniqueSize; b++)
   {
      uniqueItemArray[b] = -1;
   }


   for (a = 0; a < totalSize; a++ )
   {
      for(b = 0; b < UniqueSize; b++)
	  {
	     if(Array[a] == uniqueItemArray[b])
		 {
		    notFound = 0;
			break;
		 }
	  }// for (uniqueItem)

	  if(notFound == 1)
	  {
	     uniqueItemArray[c++] = Array[a];
	  }

	  // reinitiallize notFound for each item
	  notFound = 1;
   }//for (all items)

}

void getSumStat(struct SumStat **base, double *data, int ssIndex, int size)//numTaxa can also mean numLoci, depending on the structure of base)
{
   int i;

   if((ssIndex < PI_b) || (ssIndex > JW_Psi))
   {
		perror("ERROR: ssName index out of bound\n");
        exit(EXIT_FAILURE);
   }

   for(i = 0; i < size; i ++)
	{	data[i] = (base[i]->SS)[ssIndex]; }
	
}

double moment(double *data, int size, int order)
{
    if ((size <= 0)  || (order < 1))
    {
	perror("ERROR: correct range expected for number of taxa and order to calculate moment\n");
        exit(EXIT_FAILURE);
    }

    int a, b;
    double sum = 0.0, val = 1.0;

    for(a = 0; a < size; a++)
    {
	for(b = 1; b <= order; b++)
	{
	       val *= data[a];
	       //fprintf(stderr, "val : %lf\n", val);
	}
       sum += val;
       val = 1.0;
    }

    //fprintf(stderr, "sum : %lf \n", sum);
    sum /= ((double)size);
    
    //fprintf(stderr, "answer : %lf \n", sum);
    //fprintf(stderr, "%lf %lf %d\n", data[0], sum, size);
    return sum;
}

/*
*/
void
PrintSumStatsArray (struct SumStat **SumStat_list, int numTaxonLocusPairs, int numLoci, int numTaxonPairs)
{
  int a, b;

  /****** NOTE ******
   *
   * (A) If new summary stat is added or the print order is
   *     changed, please modify the global: numStats and
   *     ssNameVect (top of this file).  numSumStats should be
   *     the number of summary statistics used for each taxon
   *     pair.
   *
   * (B) If new prior is added or the print order is changed,
   *     please modify numPriorColumns and priorNameVect.  For
   *     prior names, start with "PRI."
   *
   * ORDER of names is important!
   */

  
  /* WORK HERE: Is this sorting appropriate? */
  if (gSortPattern == 0) {
    ; /* no sorting */
  }
  else if(gSortPattern == 1)
  {
    /* Simple sorting of summary statistics without paying attention to
     * taxonPairs:genes.
     * The taxonPairs:gene with the largerest pi.b (divergence between the
     * pair) becomes the 1st column (left most).
     */	 
   qsort(SumStat_list, numTaxonLocusPairs, sizeof(SumStat_list[0]), SS_comp);
  }
  else if(gSortPattern == 2)
  {  
     /*
     * SumStats vectors of each locus will be sorted seperately by PI_b
     */
    int *uniqueLociArray =  calloc (numLoci, sizeof (int));
    int *lociArray = calloc (numTaxonLocusPairs, sizeof(int));

    struct SumStat **tempSumStat = calloc(numTaxonLocusPairs, sizeof(struct SumStat*));

    if ((uniqueLociArray == NULL) || (lociArray == NULL) || (tempSumStat == NULL))
	{
	    fprintf(stderr, "Not enough memory for sorting summary statistics array\n");
        exit(0);
	}

    for (a = 0; a <numTaxonLocusPairs; a++)
    {
       lociArray[a] = (int)((SumStat_list[a]->SS)[locusID]);
    }

			
    getUniqueItems(lociArray, uniqueLociArray, numTaxonLocusPairs, numLoci);
	
	int headCursor = -1, tailCursor = -1, c = 0;

	for(b = 0; b < numLoci; b++)
	{
	   for (a = 0; a <numTaxonLocusPairs; a++)
	     {
		    if((int)((SumStat_list[a]->SS)[locusID]) == uniqueLociArray[b])
			{
               if(headCursor == -1)
               {
				   headCursor = tailCursor = c;
			   }
               tempSumStat[c++] = SumStat_list[a];

               tailCursor ++;

			}// if
		 }// for a

		//head = &tempSumStat[headCursor];		
		qsort(&tempSumStat[headCursor], tailCursor - headCursor, sizeof(tempSumStat[headCursor]), SS_comp);
			
		headCursor = tailCursor = -1;
		//head = NULL;
	}// for b

		
	for (a = 0; a < numTaxonLocusPairs; a++)
	{
		SumStat_list[a] = tempSumStat[a];
	}
		
    free(uniqueLociArray);
    free(lociArray);
	free(tempSumStat);
  }
  else if(gSortPattern == 3) /* moment across taxa (by locus)*/
  {
	int *uniqueLociArray =  calloc (numLoci, sizeof (int));
    int *lociArray = calloc (numTaxonLocusPairs, sizeof(int));

    //struct SumStat **head;
    struct SumStat **tempSumStat = calloc(numTaxonLocusPairs, sizeof(struct SumStat*));

	if ((uniqueLociArray == NULL) || (lociArray == NULL) || (tempSumStat == NULL))
	{
	    fprintf(stderr, "Not enough memory for sorting summary statistics array\n");
        exit(0);
	}

	double **SumStat_mean = calloc(numLoci, sizeof(double *)); //1st central moment
	double **SumStat_var = calloc(numLoci, sizeof(double *)); //2nd central moment
	double **SumStat_skewness = calloc(numLoci, sizeof(double *)); //3rd central moment
	double **SumStat_kurtosis = calloc(numLoci, sizeof(double *)); //4th central moment

	if((SumStat_mean == NULL) || (SumStat_var == NULL) || (SumStat_skewness == NULL) || (SumStat_kurtosis == NULL))
	{
	   perror("ERROR:Not enough memory for 4-moment summary statistics vector across taxa\n");
	   exit(EXIT_FAILURE);
	}

	// allocate memory for SumStat_mean, SumStat_var, SumStat_skewness and SumStat_kurtosis
    for (a = 0; a < numLoci; a++)
	{
	   SumStat_mean[a] = calloc(numSumStats, sizeof(double));
	   SumStat_var[a] = calloc(numSumStats, sizeof(double));
	   SumStat_skewness[a] = calloc(numSumStats, sizeof(double));
	   SumStat_kurtosis[a] = calloc(numSumStats, sizeof(double));
	   if ((SumStat_mean[a] == NULL) || (SumStat_var[a] == NULL) || (SumStat_skewness[a] == NULL) || (SumStat_kurtosis[a] == NULL))
	   {
	      perror("ERROR: No Mem for 4 moments SS across taxa\n");
          exit(EXIT_FAILURE);
	   }
	}

    for (a = 0; a <numTaxonLocusPairs; a++)
    {
		lociArray[a] = (SumStat_list[a]->SS)[locusID];
    }

    getUniqueItems(lociArray, uniqueLociArray, numTaxonLocusPairs, numLoci);

	int headCursor = -1, tailCursor = -1, count = 0;
	enum SS currentSS;
	
	for(b = 0; b < numLoci; b++)
	{
	   for (a = 0; a <numTaxonLocusPairs; a++)
	     {
		    if(((SumStat_list[a]->SS)[locusID]) == uniqueLociArray[b])
			{
               if(headCursor == -1)
               {
				   headCursor = tailCursor = count;
			   }
               tempSumStat[count++] = SumStat_list[a];

               tailCursor ++;

			}// if
		 }// for a

        //head = &tempSumStat[headCursor];

        double *data = calloc(tailCursor - headCursor, sizeof(double));

		for (currentSS = PI_b; currentSS <= JW_Psi; currentSS ++ )
		{
			getSumStat(&tempSumStat[headCursor], data, currentSS, tailCursor - headCursor);
			SumStat_mean[b][currentSS] = moment(data,tailCursor - headCursor, 1 );
			SumStat_var[b][currentSS] = moment(data, tailCursor - headCursor, 2);
			SumStat_skewness[b][currentSS] = moment(data, tailCursor - headCursor, 3);
			SumStat_kurtosis[b][currentSS] = moment(data, tailCursor - headCursor, 4);
		}

		headCursor = tailCursor = -1;
		//head = NULL;

		free(data);
	}// for b

    if (gPrintHeader)
    {
      int numPriorColumns = 0; /* Prior is not printed anymore */
      char priorNameVect2[][MAX_LEN_COLUMN_NAME] = { "PRI.Psi", "PRI.var.t", "PRI.E.t", "PRI.omega" };
      PrintHeader_fourMoments2 (priorNameVect2, numPriorColumns, ssNameVect,
		   numSumStats, numLoci);

      gPrintHeader = 0; /* lower the flag to print only once */

    }

	for(currentSS = PI_b; currentSS <= JW_Psi; currentSS ++)
	{
		for(a = 0; a < numTaxonPairs; a++)
		{	
			printf("%lf\t", SumStat_mean[a][currentSS]);
			printf("%lf\t", SumStat_var[a][currentSS]);
			printf("%lf\t", SumStat_skewness[a][currentSS]);
			printf("%lf\t", SumStat_kurtosis[a][currentSS]);
			
		}
	}

	printf ("\n");

    free(uniqueLociArray);
    free(lociArray);
	free(tempSumStat);

    // free memory for SumStat_mean
    for (a = 0; a < numLoci; a++)
	{
	   free(SumStat_mean[a]);
	   free(SumStat_var[a]);
	   free(SumStat_skewness[a]);
	   free(SumStat_kurtosis[a]);
	}

    free(SumStat_mean);
	free(SumStat_var);
	free(SumStat_skewness);
	free(SumStat_kurtosis);

  } 
  else if((gSortPattern >= 4) && (gSortPattern <= 7)) // moment across loci (by taxon), default gSortPattern == 1
  {
    
    int *uniqueTaxaArray =  calloc (numTaxonPairs, sizeof (int));
    int *taxaArray = calloc (numTaxonLocusPairs, sizeof(int));

        struct SumStat **tempSumStat = calloc(numTaxonLocusPairs, sizeof(struct SumStat*));
		
	if ((uniqueTaxaArray == NULL) || (taxaArray == NULL) || (tempSumStat == NULL))
	{
	    fprintf(stderr, "Not enough memory for sorting summary statistics array\n");
        exit(0);
	}
	
	double **SumStat_mean = calloc(numTaxonPairs, sizeof(double *)); //1st central moment
	double **SumStat_var = calloc(numTaxonPairs, sizeof(double *)); //2nd central moment
	double **SumStat_skewness = calloc(numTaxonPairs, sizeof(double *)); //3rd central moment
	double **SumStat_kurtosis = calloc(numTaxonPairs, sizeof(double *)); //4th central moment
	
	if((SumStat_mean == NULL) || (SumStat_var == NULL) || (SumStat_skewness == NULL) || (SumStat_kurtosis == NULL))
	{
	   perror("ERROR:Not enough memory for 4 moment summary statistics vectr\n");
	   exit(EXIT_FAILURE);
	}

	// allocate memory for SumStat_mean
	for (a = 0; a < numTaxonPairs; a++)
	{
	   SumStat_mean[a] = calloc(numSumStats, sizeof(double));
	   SumStat_var[a] = calloc(numSumStats, sizeof(double));
	   SumStat_skewness[a] = calloc(numSumStats, sizeof(double));
	   SumStat_kurtosis[a] = calloc(numSumStats, sizeof(double));
	   if ((SumStat_mean[a] == NULL) || (SumStat_var[a] == NULL) || (SumStat_skewness[a] == NULL) || (SumStat_kurtosis[a] == NULL))
	   {
	      perror("ERROR: No Mem for moments SumStat vectors\n");
	      exit(EXIT_FAILURE);
	   }
	   
	}
	
    for (a = 0; a < numTaxonLocusPairs; a ++)
    {
		taxaArray[a] = (int)((SumStat_list[a]->SS)[taxonID]);
    }

   getUniqueItems(taxaArray, uniqueTaxaArray, numTaxonLocusPairs, numTaxonPairs);

	int headCursor = -1, tailCursor = -1, count = 0;

	for(b = 0; b < numTaxonPairs; b++)
	{
	   for (a = 0; a <numTaxonLocusPairs; a++)
	     {
		    if((int)((SumStat_list[a]->SS)[taxonID]) == uniqueTaxaArray[b])
			{
               if(headCursor == -1)
               {
				   headCursor = tailCursor = count;
			   }
               tempSumStat[count++] = SumStat_list[a];

               tailCursor ++;

			}// if
		 }// for a


		double *data = calloc(tailCursor - headCursor, sizeof(double));

		// calculate moment for each summary statistics
		enum SS currentSS;
		for (currentSS = PI_b; currentSS <= JW_Psi; currentSS ++ )
		{
			getSumStat(&tempSumStat[headCursor], data, currentSS, tailCursor - headCursor);
			SumStat_mean[b][currentSS] = moment(data,tailCursor - headCursor, 1 );
			
			//fprintf(stderr, "%d %d %lf\n", b, currentSS,SumStat_mean[b][currentSS]);
		
			SumStat_var[b][currentSS] = moment(data, tailCursor - headCursor, 2);
			SumStat_skewness[b][currentSS] = moment(data, tailCursor - headCursor, 3);
			SumStat_kurtosis[b][currentSS] = moment(data, tailCursor - headCursor, 4);
			
		}

		headCursor = tailCursor = -1;
		//head = NULL;

		free(data);
	}// for b
	
	// fill entries of matrixToPrint
	double **matrixToPrint = calloc(numTaxonPairs, sizeof(double*));
	int rowSize = (8 - gSortPattern)*numSumStats;
	int ii; 
	enum SS currentSS;
	for(ii = 0; ii < numTaxonPairs; ii++)
			{	matrixToPrint[ii] = calloc(rowSize, sizeof(double));	}
	
	// scale is:
	// 4 if gSortPattern = 4, we use mean, variance, skewness and kurtosis
	// 3 if gSortPattern = 5, we use mean, variance, skewness
	// 2 if gSortPattern = 6, we use mean, variance
	// 1 if gSortPattern = 7, we use mean
	int scale = 8 - gSortPattern;
	for(ii = 0; ii < numTaxonPairs; ii++)
	{		
		for (currentSS = PI_b; currentSS <= JW_Psi; currentSS ++)
		{	
			matrixToPrint[ii][currentSS*scale + 0] = SumStat_mean[ii][currentSS];
			//fprintf(stderr, "%d %d %lf\n", ii, currentSS, matrixToPrint[ii][currentSS*scale + 0]);
			if((gSortPattern == 4) || (gSortPattern == 5) || (gSortPattern == 6))
				{	matrixToPrint[ii][currentSS*scale + 1] = SumStat_var[ii][currentSS]; }
			if((gSortPattern == 4) || (gSortPattern == 5))
				{	matrixToPrint[ii][currentSS*scale + 2] = SumStat_skewness[ii][currentSS]; }
			if(gSortPattern == 4)
				{	matrixToPrint[ii][currentSS*scale + 3] = SumStat_kurtosis[ii][currentSS]; }
		}
	}
	
	
	// free allocated arrays that we don't need anymore
	free(uniqueTaxaArray);
	free(taxaArray);
	free(tempSumStat);

    // free memory for four moment SumStat
    for (a = 0; a < numTaxonPairs; a++)
	{
	   free(SumStat_mean[a]);
	   free(SumStat_var[a]);
	   free(SumStat_skewness[a]);
	   free(SumStat_kurtosis[a]);
	}

    free(SumStat_mean);
	free(SumStat_var);
	free(SumStat_skewness);
	free(SumStat_kurtosis);

	// sort rows of ss_moments by PI_b_mean(across loci)
	qsort(matrixToPrint, numTaxonPairs, sizeof(matrixToPrint[0]), SS_comp_moments);
    
	if (gPrintHeader)
    {
      int numPriorColumns = 0; // Prior is not printed anymore 
      char priorNameVect2[][MAX_LEN_COLUMN_NAME] = { "PRI.Psi", "PRI.var.t", "PRI.E.t", "PRI.omega" };
      PrintHeader_fourMoments (priorNameVect2, numPriorColumns, ssNameVect,
		   numSumStats, numTaxonPairs);

      gPrintHeader = 0; // lower the flag to print only once

    }

	for(currentSS = PI_b; currentSS <= JW_Psi; currentSS ++)
	{
		for(a = 0; a < numTaxonPairs; a++)
		{	
			int aa;
			for(aa = 0; aa < scale; aa ++)
				{	printf("%lf\t", matrixToPrint[a][currentSS*scale + aa]);	}
		}
	}


	printf ("\n");
	
	    // free memory for four moment SumStat
    for (a = 0; a < numTaxonPairs; a++)
	{
	   free(matrixToPrint[a]);
	}

    free(matrixToPrint);

  }

  if (gPrintHeader)
	{
		int numPriorColumns = 0; /* Prior is not printed anymore */
		char priorNameVect[][MAX_LEN_COLUMN_NAME] =
			{ "PRI.Psi", "PRI.var.t", "PRI.E.t", "PRI.omega" };
		PrintHeader (priorNameVect, numPriorColumns, ssNameVect,
		   numSumStats, numTaxonLocusPairs);
		gPrintHeader = 0; /* lower the flag to print only once */
	}

	if(gSortPattern < 3)
	{
		/* start to print sum stats */
		enum SS currentSS;
		for(currentSS = PI_b; currentSS <= JW_Psi; currentSS ++)
		{
			for(a = 0; a < numTaxonLocusPairs; a ++)
			{
				printf("%lf\t", (SumStat_list[a]->SS)[currentSS]);
			}
		}
	
		printf("\n");
		
	}// if (gSortPattern < 3)

}

static void
PrintHeader (char priorNames[][MAX_LEN_COLUMN_NAME], int numPriors,
	     char sumStatNames[][MAX_LEN_COLUMN_NAME], int numSumStats,
	     int numTaxonPairs)
{
  int i, a;
  for (i = 0; i < numPriors; i++)
    {
      printf ("%s\t", priorNames[i]);
    }
  for (i = 0; i < numSumStats; i++)
    {
      for (a = 0; a < numTaxonPairs; a++)
	{
	  printf ("%s.%d\t", sumStatNames[i], a + 1);
	}
    }
  printf ("\n");
  return;
}


static void
PrintHeader_fourMoments2 (char priorNames[][MAX_LEN_COLUMN_NAME], int numPriors,
	     char sumStatNames[][MAX_LEN_COLUMN_NAME], int numSumStats,
	     int numLoci)
{
  int i, a, count;
  for (i = 0; i < numPriors; i++)
    {
      printf ("%s\t", priorNames[i]);
    }
  for (i = 0; i < numSumStats; i++)
    {
	  count = 1;
      for (a = 0; a < numLoci; a++)
	  {
	     printf ("%s.%d\t", sumStatNames[i], count++);
	     printf ("%s.%d\t", sumStatNames[i], count++);
	     printf ("%s.%d\t", sumStatNames[i], count++);
		 printf ("%s.%d\t", sumStatNames[i], count++);
	  }
    }
  printf ("\n");
  return;
}

static void
PrintHeader_fourMoments (char priorNames[][MAX_LEN_COLUMN_NAME], int numPriors,
	     char sumStatNames[][MAX_LEN_COLUMN_NAME], int numSumStats,
	     int numTaxonPairs)
{
  int i, a, count;
  for (i = 0; i < numPriors; i++)
    {
      printf ("%s\t", priorNames[i]);
    }
  for (i = 0; i < numSumStats; i++)
    {
	  count = 1;
      for (a = 0; a < numTaxonPairs; a++)
	  {
	     printf ("%s.%d\t", sumStatNames[i], count++);
		 if((gSortPattern == 4) || (gSortPattern == 5) || (gSortPattern == 6))
			printf ("%s.%d\t", sumStatNames[i], count++);
		 if((gSortPattern == 4) || (gSortPattern == 5))
			printf ("%s.%d\t", sumStatNames[i], count++);
		 if(gSortPattern == 4)
			printf ("%s.%d\t", sumStatNames[i], count++);
	  }
    }
  printf ("\n");
  return;
}


/*
 * Checks that sub population sample sizes n[] are reasonable.
 * Arguments:
 *   nsam:     number of total samples in the simulation
 *   npops:    number of subpopulations
 *   n[npops]: sub-population sample sizes
 *
 * Returns: 1 if all subpop sample sizes are bet. 0 and nsam (ends exclusive)
 *            and elements of n[] add up to nsam
 *          0 otherwise
 */
int
multiplepopssampledfrom (int nsam, int npops, int *n)
{
  int i, sum = 0;
  for (i = 0; i < npops; i++)
    {
      sum += n[i];

      if ((n[i] <= 0) || (n[i] >= nsam))
	return (0);
    }
  /* This function was checking only the first n[i], and returning 1
     if at least 1 element is 0 < n[i] < nsam.
     I don't think this is the intention, so I corrected it. Naoki
   */

  /* the additional check is added by Naoki */
  if (sum != nsam)
    return 0;

  return (1);
}

int
pairwisediffc_w (int ss, int nsam, char **list, int np, int *n, int pop)
{
  int n1, n2, diffc = 0;
  int popi, startn = 0;
  int s;
  if (n[pop] > 1)
    {
      for (popi = 0; popi < pop; popi++)
	startn += n[popi];
      for (n1 = startn; n1 < (startn + n[pop] - 1); n1++)
	for (n2 = n1 + 1; n2 < (startn + n[pop]); n2++)
	  {
	    for (s = 0; s < ss; s++)
	      if (list[n1][s] != list[n2][s])
		diffc++;
	  }
    }
  /*printf("piW: %d\n", diffc);  test print */
  return (diffc);
}


int
pairwisediffc_b (int ss, int nsam, char **list, int np, int *n, int pop1,
		 int pop2)
{
  int n1, n2, diffc = 0;
  int popi, startn1, startn2;
  int s;
  if ((n[pop1] > 0) && (n[pop2] > 0))
    {
      startn1 = startn2 = 0;
      for (popi = 0; popi < pop1; popi++)
	startn1 += n[popi];
      for (popi = 0; popi < pop2; popi++)
	startn2 += n[popi];
      for (n1 = startn1; n1 < (startn1 + n[pop1]); n1++)
	for (n2 = startn2; n2 < (startn2 + n[pop2]); n2++)
	  {
	    for (s = 0; s < ss; s++)
	      if (list[n1][s] != list[n2][s])
		diffc++;
	    /*printf("diffc: %d", diffc);  test print */
	  }
    }
  /*printf("piB: %d\n", diffc);  test print */
  return (diffc);

}

//Square root of variance of pairwise differences for one population
double
S_w (int ss, int nsam, char **list, int np, int *n, int pop, double D_w, int basePairs)
{
  int n1, n2, diffc;
  int popi, startn = 0;
  int s;
  double diffcPerSite, diff_square = 0.0;
  double var_D_w = 0.0, Sw = 0.0;
  
  if (n[pop] > 1)
    {
      for (popi = 0; popi < pop; popi++)
	{
	  startn += n[popi];
	}
      
      for (n1 = startn; n1 < (startn + n[pop] - 1); n1++)
        for (n2 = n1 + 1; n2 < (startn + n[pop]); n2++)
    	  {
    	    diffc = 0;
    	    for (s = 0; s < ss; s++)
    	      if (list[n1][s] != list[n2][s])
		diffc++;
	    
	    diffcPerSite = (double) diffc / basePairs;
    	    diff_square = (diffcPerSite - D_w) * (diffcPerSite - D_w);
    	    var_D_w += diff_square;
    	  }//for n2
    }// if
    else
      return -1.0;

  /* denominator is number of pairwise comparison n[pop] choose 2 */
  var_D_w = var_D_w/ (n[pop] * (n[pop]-1)/2);
  Sw = sqrt(var_D_w);
  //fprintf(stderr, "Sw, ss = %d, nsam = %d, np = %d, pop = %d, Sw = %lf\n", ss, nsam, np, pop, Sw);

  return (Sw);
}



//variance of pairwise difference between
double
S_xy (int ss, int nsam, char **list, int np, int *n, int pop1,
      int pop2, double D_xy, int basePairs)
{
  int n1, n2, diffc;
  int popi, startn1, startn2;
  int s;

  double diffcPerSite, diff_square = 0.0;
  double var_D_b = 0.0, Sxy = 0.0;
  if ((n[pop1] > 0) && (n[pop2] > 0))
    {
      startn1 = startn2 = 0;
      for (popi = 0; popi < pop1; popi++)
	startn1 += n[popi];
      for (popi = 0; popi < pop2; popi++)
	startn2 += n[popi];
      
      for (n1 = startn1; n1 < (startn1 + n[pop1]); n1++)
	for (n2 = startn2; n2 < (startn2 + n[pop2]); n2++)
	  {
	    diffc = 0;
	    for (s = 0; s < ss; s++)
	      if (list[n1][s] != list[n2][s])
		diffc++;
	    
	    diffcPerSite = (double) diffc / basePairs;
	    diff_square = (diffcPerSite - D_xy)*(diffcPerSite - D_xy);
	    var_D_b += diff_square;
	  }
    }

  var_D_b = var_D_b / (n[pop1] * n[pop2]);
  Sxy = sqrt(var_D_b);
  //fprintf(stderr, "S_xy: ss = %d, nsam = %d, np = %d, pop1 = %d, pop2 = %d, Sxy = %lf\n", ss, nsam, np, pop1, pop2, Sxy);

  return (Sxy);

}

//John Wakely's Psi
double JWakely_Psi(int* n, double S_x, double S_y, double S_xy, double d_x, double d_y, double d_k)
{
       //fprintf(stderr, "S_x = %lf, S_y = %lf, S_xy = %lf,d_x = %lf, d_y = %lf, d_k = %lf\n", S_x, S_y, S_xy, d_x, d_y, d_k);
	double JW_Psi = 0.0;
	double item0, item1, item2, item3, item4;
	item0 = item1 = item2 = item3 = item4 = 0.0;
	int nx, ny, nTotal;
	nx = ny = nTotal = 0;
        nx = n[0], ny = n[1];
        nTotal = nx + ny;

        item0 = nTotal * (nTotal - 1);
	item1 = 1 / item0;
      	item2 = (d_x > 0)? nx*(nx - 1)*S_x/d_x : 0;
	item3 = (d_y > 0)? ny*(ny - 1)*S_y/d_y : 0;
	item4 = (d_k > 0)? 2*nx*ny*S_xy/d_k : 0;

	JW_Psi = item1 * (item2 + item3 + item4);

    //fprintf(stderr, "nx = %d, ny = %d, item1 = %lf, item2 = %lf, item3 = %lf, item4 = %lf, JW_Psi = %lf\n", nx, ny, item1, item2, item3, item4, JW_Psi);
	return JW_Psi;
}


/*yyy BELOW  Eli 05/15/06 yyy*/
/*
 * Arguments:
 *   segwReturnArr: results get returned here
 *   segsites: number of segSite in the entire data = # columns in list matrix
 *   list: sequence data array
 *   numSubpops: number of sub populations
 *   n: array with numSubpops elements.  i-th element is # seqs in i-th subpops
 *
 * Assumes that memory is allocated to segwReturnArr (numSubpops * sizeof(int)).
 *
 * Side effect: segReturnArr get filled up by # seg sites within each subpops.
 */
void
segsubs (int *segwReturnArr, int segsites, char **list, int numSubpops, int *n)
{
  int i, count, npi, n0 = 0, ni, npops_gt1 = 0;

  /* initialize the result counter */
  memset(segwReturnArr, 0, numSubpops *  sizeof(int));

  for (npi = 0; npi < numSubpops; npi++)
    {
      if (n[npi] > 1)
	{
	  ++npops_gt1;
	  count = 0;
	  for (i = 0; i < segsites; i++)
	    {
	      for (ni = n0 + 1; ni < n0 + n[npi]; ni++)
		{
		  if (list[ni][i] != list[n0][i])
		    {
		      segwReturnArr[npi]++;
		      break;
		    }
		}
	    }
	}
      else
	{
	  segwReturnArr[npi] = 0;
	}
      n0 += n[npi];
    }
}

/*yyy ABOVE  Eli 05/15/06 yyy*/


/*
 * Calculate the pi (per gene and not per site) within sub populations.
 * Specify the index (0-offset) of sub population from which pi are calculated
 * in targetPop.
 *
 * If a negative value of targetPop is given, it calculate the pi within subpop
 * for each sub population and average is taken:  Let's say two subpops with
 * n0 and n1 samples, and pi's within each subpops are pi0 and pi1.
 * Weighted average is (C(n0,2) * n0 + C(n1,2) * n1) /(C(n0,2)+C(n1,2)),
 * where C(x,y) denote combinatorial: Choose y from x.
 */
double
nucdiv_w (int nsam, int segsites, char **list, int np, int *n, int targetPop)
{
  int pop, pairwisediffc_w (int, int, char **, int, int *, int);
  int beginPop, endPop;
  double pi, nd;
  double num_comps;

  pi = 0.0;

  nd = nsam;
  num_comps = 0.;

  if (targetPop < 0)
    {
      beginPop = 0;
      endPop = np;
    }
  else
    {
      beginPop = targetPop;
      endPop = targetPop + 1;
    }

  for (pop = beginPop; pop < endPop; pop++)
    {
      pi += pairwisediffc_w (segsites, nsam, list, np, n, pop);
      num_comps += (double) n[pop] * ((double) n[pop] - 1.) / 2.;
    }

  pi /= num_comps;
  return (pi);
}


double
nucdiv_bw (int nsam, int segsites, char **list, int np, int *n)
{
  int pop1, pop2, pairwisediffc_b (int, int, char **, int, int *, int, int);
  double pi, nd;
  double num_comps;

  pi = 0.0;

  nd = nsam;
  num_comps = 0;
  for (pop1 = 0; pop1 < (np - 1); pop1++)
    for (pop2 = (pop1 + 1); pop2 < np; pop2++)
      {
	pi += pairwisediffc_b (segsites, nsam, list, np, n, pop1, pop2);
	/*printf("piB: %lf\n", pi);  test print */
	num_comps += (double) n[pop1] * (double) n[pop2];
      }
  pi /= num_comps;
  /*printf("piB-FINAL: %lf\n", pi);  test print */
  return (pi);
}


void
FrequencyDistrInfSites (int *freqdist, int n, int S, char **list)
{
  int i;
  for (i = 0; i < n; i++)
    freqdist[i] = 0;
  for (i = 0; i < S; i++)
    {
      freqdist[cntBase ('1', i, n, list)]++;	/* probably bogus for ACGT */
    }
}


void
FrequencyDistrQ (int *freqdist, int n, int S, char **list)
{
  int i, f, cntBase (char, int, int, char **);
  for (i = 0; i < n; i++)
    freqdist[i] = 0;
  for (i = 0; i < S; i++)
    {
      f = cntBase (list[0][i], i, n, list);
      freqdist[f < n / 2 + 0.0001 ? f : n - f]++;	/* probably bogus for ACGT */
    }
}


void
FuLi (D, F, n, S, list, pi)
     int n, S;
     char **list;
     double *D, *F, pi;
{
  int k, s, etae, cntBase (char, int, int, char **);
  double n1, S1, vD, uD, uF, vF, an, bn, cn;
  n1 = (double) n;
  S1 = (double) S;
  for (s = etae = 0; s < S; s++)
    if (cntBase ('1', s, n, list) == 1)
      etae++;
  for (k = 1, an = bn = 0.; k < n; k++)
    {
      an += 1. / (double) k;
      bn += 1. / (double) (k * k);
    }
  if (n == 2)
    cn = 1.;
  else
    cn = 2. * (n1 * an - 2. * (n1 - 1)) / ((n1 - 1) * (n1 - 2));
/*   printf("an:\t%9.7f\tbn:\t%9.7f\tcn:\t%9.7f\t",an,bn,cn); */
  vD = 1. + an * an / (bn + an * an) * (cn - (n1 + 1) / (n1 - 1));
  uD = an - 1. - vD;
  *D = S1 - an * etae;
  *D /= sqrt (uD * S1 + vD * S1 * S1);
  vF = cn + 2. * (n1 * n1 + n1 + 3) / (9. * n1 * (n1 - 1)) - 2. / (n1 - 1);
  vF /= an * an + bn;
  uF = 1. + (n1 + 1) / (3. * (n1 - 1)) - 4. * (n1 + 1) /
    ((n1 - 1) * (n1 - 1)) * (an + 1. / n1 - 2. * n1 / (n1 + 1));
  uF /= an;
  uF -= vF;
/*   printf("vF:\t%9.7f\tuF:\t%9.7f\t",vF,uF); */
  *F = pi - etae;
  *F /= sqrt (uF * S1 + vF * S1 * S1);
}

/* 
 * Argument:
 *  n: number of samples
 *  S: number of segregating sites
 * S can be floating points, e.g. S/seqLen to calculate "PER SITE" value
 */
double
tajddenominator (int n, double S)
{
  int i;
  double n1, a1, a2, b1, b2, e1, e2, denom;
  n1 = (double) n;

  a1 = a2 = 0.;
  for (i = n-1; i > 0; i--)
    {
      a1 += 1. / i;
      a2 += 1. / (i * i);
    }
  b1 = (n1 + 1) / (3. * (n1 - 1));
  b2 = 2. * ((n1 * n1) + n1 + 3) / (9. * n1 * (n1 - 1));
  e1 = (b1 - 1. / a1) / a1;
  e2 = (b2 - (n1 + 2.) / (a1 * n1) + a2 / (a1 * a1)) / (a1 * a1 + a2);
  denom = sqrt (e1 * S + e2 * S * (S - 1));
  return (denom);
}

double
thetaW (int n, int S)
{
  int i;
  double a1=0.;

  for (i = n-1; i > 0; i--)
    a1 += 1. / i;

  return ((double) S / a1);
}


/* This calculate overall average pairwise differences (pi per gene)
 * ignoring sub population designation  */
double
nucdiv (int nsam, int segsites, char **list)
{
  int s;
  double pi = 0.0, denom;

  for (s = 0; s < segsites; s++)
    {
      pi += pwDiffAtOneSite (s, nsam, list);
      /* pwDiffAtOneSite() returns the number of pair wise differences at site s
       * from all pairwise comparison */
    }

  /* denom is  # of ways to choose a pair: nsam choose 2 */
  denom = nsam * (nsam - 1) / 2;
  return (pi / denom);
}


/*
 * count the occurence of "allele" at the "site"
 * nsam: number of samples
 * For example, 
 *   cntBase('A', 0, 10, list)
 * goes through the 1st polymorphic site (10 samples), and 
 * returns number of 'A' at this site
 */
int
cntBase (char base, int site, int nsam, char **list)
{
  int i, count = 0;
  for (i = 0; i < nsam; i++)
    count += (list[i][site] == base);
  return (count);
}



/*
 * Count the number of pairwise differences at the site.
 * nsam * (nsam - 1) / 2 pairs are compared.
 *
 * Arguments:
 *   site: i-th segregating sites
 *   nsam: total number of samples
 *   list: character matrix of segregating sites
 *
 * Returns:  the number of pairwise differences at the site
 */
int
pwDiffAtOneSite (int site, int nsam, char **list)
{
  char allele1;			/*7/27/04; Hickerson */
  int i, n, denom, count = 0;

  denom = 0;

  for (n = 0; n < nsam; n++)
    {
      allele1 = list[n][site];

      for (i = n + 1; i < nsam; i++)
	{
	  if (list[i][site] != allele1)
	    count = count + 1;
	}

    }
  return (count);
}

int
frequencySING (char base, int site, int nsam, char **list)	/* in progress Hickerson 7/29/04 */
{
  char allele1;			/*7/27/04; Hickerson */
  int i, n, denom, singleton, count = 0;

  denom = 0;
  singleton = 0;
  for (n = 0; n < nsam; n++)
    {				/*7/27/04; Hickerson */
      allele1 = list[n][site];	/*7/27/04; Hickerson */
      count = 0;


      for (i = n; i < nsam; i++)
	{			/*7/27/04; Hickerson */
	  if (list[i][site] == allele1)
	    count = count;	/*7/27/04; Hickerson */
	  else
	    count = count + 1;	/*7/27/04; Hickerson */

	}
      if ((count = nsam - 1))
	{
	  singleton++;
	}
    }

  return (singleton);
}

/*  thetah - pi   */
/* 	double */
/* hfay( int nsam, int segsites, char **list) */
/* { */
/* 	int s, frequency( char, int, int, char**); */
/* 	double pi, p1, nd, nnm1  ; */

/* 	pi = 0.0 ; */

/* 	nd = nsam; */
/* 	nnm1 = nd/(nd-1.0) ; */
/*    	for( s = 0; s <segsites; s++){ */
/* 		p1 = frequency('1', s,nsam,list)/nd ; */
/* 		pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ; */
/* 		} */
/* 	return( pi ) ; */
/* } */

/* Fay's theta_H  */
double
thetah (int nsam, int segsites, char **list)
{
  int s;
  double pi, p1, nd, nnm1;

  pi = 0.0;

  nd = nsam;
  nnm1 = nd / (nd - 1.0);
  for (s = 0; s < segsites; s++)
    {
      p1 = cntBase ('1', s, nsam, list);

      pi += p1 * p1;
    }
  return (pi * 2.0 / (nd * (nd - 1.0)));
}




int
segsub (int nsub, int segsites, char **list)
{
  int i, count = 0, c1;

  for (i = 0; i < segsites; i++)
    {
      c1 = cntBase ('1', i, nsub, list);
      if ((c1 > 0) && (c1 < nsub))
	count++;
    }
  return (count);
}

void
shannonIndex (char **list, int *config, double **shannonIndexArray)
{
  int i, sizeOfSp1, sizeOfSp2, sizeAll, unit=1;
  double sHa1 = 0, sHa2 = 0, sHu = 0, sHua = 0, temp;
  hashtab_iter_t iHash;
  hashtab_t *subPop1, *subPop2, *pool;

  sizeOfSp1 = config[0];
  sizeOfSp2 = config[1];

  subPop1 = ht_init (sizeOfSp1, NULL);
  subPop2 = ht_init (sizeOfSp2, NULL);
  pool = ht_init (sizeOfSp1 + sizeOfSp2, NULL);

  if (subPop1 == NULL || subPop2 == NULL || pool == NULL) {
    fprintf(stderr, "ERROR: no memory in shannonIndex\n");
    exit(EXIT_FAILURE);
  }

  sizeOfSp1 = config[0];
  sizeOfSp2 = config[1];
  sizeAll = sizeOfSp1 + sizeOfSp2;

  /* insert allele-count pair into hashtables for subPop1
     key: allele as string, value: number of allele as int) */
  for (i = 0; i < sizeOfSp1; i++)
    {
      int *thisCnt = (int *) ht_search (subPop1, list[i], charCount (list[i]));
      if ( thisCnt == NULL) { /* new key */
	ht_insert (subPop1, list[i], charCount (list[i]), &unit, sizeof (int));
	ht_insert(pool,list[i], charCount (list[i]), &unit, sizeof (int));
      } else {/* this key have seen before */
	(*thisCnt)++;
	(*((int *) ht_search (pool, list[i], charCount (list[i])))) ++;
      }
    }

  /* For the subPop2 */
  for (i = 0; i < sizeOfSp2; i++)
    {
      int *thisCnt = (int*) ht_search (subPop2, list[sizeOfSp1 + i],
				       charCount (list[sizeOfSp1 + i]));

      if (thisCnt == NULL)
	ht_insert (subPop2, list[sizeOfSp1 + i], charCount(list[sizeOfSp1 + i]),
		   &unit, sizeof (int));
      else
	(*thisCnt)++;

      /* deal with the pooled hashTbl */
      thisCnt = (int*) ht_search (pool, list[sizeOfSp1 + i],
				  charCount (list[sizeOfSp1 + i]));

      if (thisCnt == NULL)
	ht_insert (pool, list[sizeOfSp1 + i], charCount(list[sizeOfSp1 + i]),
		   &unit, sizeof (int));
      else
	(*thisCnt)++;

    }

  // initialize hash table iterator, and calculate sHua1
  for (ht_iter_init (subPop1, &iHash); iHash.key != NULL; ht_iter_inc (&iHash))
    {
      temp = (double) (*((int *) (iHash.value)));
      temp = temp / (double) sizeOfSp1;
      sHa1 += (-(log (temp) / log ((double) 2) * temp));
    }

  // calculate sHua2
  for (ht_iter_init (subPop2, &iHash); iHash.key != NULL; ht_iter_inc (&iHash))
    {
      temp = (double) (*((int *) (iHash.value)));
      temp = temp / sizeOfSp2;
      sHa2 += (-(log (temp) / log ((double) 2) * temp));
    }

  // calculate sHu
  for (ht_iter_init (pool, &iHash); iHash.key != NULL; ht_iter_inc (&iHash))
    {
      temp = (double) (*((int *) (iHash.value)));
      temp = temp / sizeAll;
      sHu += (-((log (temp) / log ((double) 2)) * temp));
    }

  // calculate sHua
  sHua = sHu - (sizeOfSp1 * sHa1 - sizeOfSp2 * sHa2)/sizeAll;

  // throw values into double array shannonIndexArray
  *(*shannonIndexArray + 0) = sHu, *(*shannonIndexArray + 1) =
    sHua, *(*shannonIndexArray + 2) = sHa1, *(*shannonIndexArray + 3) = sHa2;

  ht_destroy (pool);
  ht_destroy (subPop1);
  ht_destroy (subPop2);
}

/*
 * Count the number of characters in a cstring (size of the char*)
 *
 * Argument:
 *   arr: the cstring whose number of characters to be counted
 *
 * Returns: the size of the string
 *
 */
int
charCount (char *arr)
{
  int k = 0;
  while (arr[k] != '\0')
    ++k;
  return k;
}				//int charCount(char*)
