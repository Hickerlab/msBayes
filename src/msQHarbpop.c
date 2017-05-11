/***** msarbpop.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number 
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.   
        usage: ms seed nsam howmany  -t  theta  [options]
		or
	       ms seed nsam howmany -s segsites  [options] 

	   nsam is the number of gametes per sample.
	   howmany is the number of samples to produce.
	   With -t the numbers of segregating sites will randomly vary 
		from one sample to the next.
	   with -s segsites,  the number of segregating sites will be
		segsites in each sample.

           Other options:

              -r rho nsites     (recombination)
              -c f track_len    (gene conversion:  To use this option,
				 the -r option must also be used, and rho must be > 0. )   
              -m mig_rate npop n1 n2 ...  (Island model with symmetric migration.)
              -d  nintn  nr1 np1 t1  [nr2 np2 t2 ... ]  (population size changes) 

*	  Arguments of the options are explained here:

           npop:  Number of subpopulations which make up the total population
           ni:  the sample size from the i th subpopulation (all must be 
		specified.) The output will have the gametes in order such that
		the first n1 gametes are from the first island, the next n2 are
		from the second island, etc.
           nsites: number of sites between which recombination can occur.
           theta: 4No times the neutral mutation rate 
           rho: recombination rate between ends of segment times 4No
	   f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
	   track_len:  mean length of conversion track in units of sites.  The 
		       total number of sites is nsites, specified with the -r option.
           mig_rate: migration rate: the fraction of each subpop made up of
                 migrants times 4No. 
           howmany: howmany samples to generate.
	   nintn:  number of intervals with specified demographic parameters.
		(Following nintn there must be 3*nintn arguments specifying the
		values of the following three parameters for each interval.
		starting with the most recent interval and proceeding into the past.
		The population growth rate implied by the values of the last
		set of three values will apply into the indefinite past.)
	   nri:  the population size (relative to No) at the beginning (most 
		recent point) of the i-th  interval.
	   npi:  the population size (relative to No) at the end (farthest in 
		the past point) of the i-th interval.
	   ti:  The duration (in units of 4No) of the i-th interval.

	Note:  In the above definition, No is the total diploid population if
		npop is one, otherwise, No is the diploid population size of each
		subpopulation. However, if the -d option is used the population sizes
		are No times various factors specified by the -d arguments.
	The following behavior of seedms is removed for msbayes.  seed is
                fed from the command line argument
	A seed file called "seedms" will be created  if it doesn't exist. The
		seed(s) in this file will be modified by the program. 
		So subsequent runs
		will produce new output.  The initial contents of seedms will be
		printed on the second line of the output.
        Output consists of one line with the command line arguments and one
	 	line with the seed(s).
		The samples appear sequentially following that line.
		Each sample begins with "//", then the number of segregating sites, the positions
		of the segregating sites (on a scale of 0.0 - 1.0). On the following
		lines are the sampled gametes, with mutants alleles represented as
		ones and ancestral alleles as zeros.
	To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
		or:  cc -o ms ms.c streec.c rand2.c -lm
	 (Of course, gcc would be used instead of cc on some machines.  And -O3 or 
		some other optimization switches might be usefully employed with some 
		compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

example runs:  ms 7 3 -t 2.0 -r 1.0 1000 -d 1 100. 1. .75  -m .5 2 4 3 >sample3.out
 	   with no seedms file present   
	produces the results shown in the file sample3.out .
	These are samples from a two-island (symmetric migration) with 4Nom = 0.5.
	Sample size 4 from one subpopulation, sample size 3 from the other.  There 
	1000 sites between which recombination can occur.  rho  = 1.0.  4Nou (the 
	mutation parameter) = 2.0.  The population exponentially shrinks as one
	goes back in time, such that the population size .75*4No generations ago
	it is No, while the present population size is 100No. The exponential
	 shrinking continues in this case indefinitely into the past. 

	An additional example of the use of -d:
       ms 10 1000 -t 2.0 -d 3 .5 .5 1.0  2.0 2.0 2.0  1.0 1.0 1.0 
	will produce 1000 samples of size 10, for a demographic history in which
	 the population size is
	No/2 for 4No generations, then is  2No for 2*4No generations, then is No
	indefinitely into the past beyond that. 
*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
	Modifications to add gene conversion 5 Mar 2001
	Added demographic options -d  13 Mar 2001
	Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
	Changed cleftr() to check for zero rand() .13 June 2001
	Move seed stuff to subroutine seedit()  11 July 2001
	Modified streec.c to handle zero length demographic intervals 9 Aug 2001
	Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
***************************************************************************/
/* 1 Nov 01 Eli changes */
/* new options:
   -Q [Rtstvratio [Ytstvratio] [Afreq Cfreq Gfreq Tfreq]] (0,1,2,4,5 or 6 args)
   -H gammaalpha
compile with QHsubs.c in addition to the standard files. */
/* end Eli */


/*
 * msQHarbpop.c
 *
 * Copyright (C) 2006  Richard Hudson and Eli Stahl
 *
 * This file is a part of msDQH, distributed with msBayes.
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
/* modified by Michael Hickerson and Naoki Takebayashi */

#include <math.h>
#include <string.h>
#include <time.h>

#include "msQHarbpop.h"
#include "QHsubs.h"
#include "msprior.h"  /* I'm not sure why this is included, Naoki 3/31/09 */

#define SITESINC 10



/* 1 Nov 01 Eli changes */
/*   added *QH, to hold most new vars (and other stuff) with minimal changes to function calls */
/*   modified main(), gensam( + &seQmut) and getpars( + &seQmut) */
/*   added functions getQpars(), Qparamscheck(), makeQHgametes(), assign_ndes(), newttime() */
void getQpars (int argc, char *argv[], int *parg, struct QHparameters *QH);
void Qparamscheck (int argc, struct QHparameters *QH);
int makeQHgametes (int n, int st, int end, struct node *tre, int ns,
		   struct QHparameters *QH, char **);
void assign_ndes (int nsam, struct node *ptree);
double newttime (struct node *ptree, int nsam, int nodenum);
/*   also added functions in QHsubs.c */
/*   added function prototypes below */
int biggerlist (int n, char **list);
/* int gensam(struct parameters , char ** , double * , double * ); */
int gensam (struct parameters, char **);
/*   fprintf(stderr,"HERE "); */
/*   fprintf(stderr,"HERE2\n"); */
/* end Eli */

int order (int n, double pbuf[]);
int ranvec (int n, double pbuf[]);
int tdesn (struct node *ptree, int tip, int node);
int pickb (int nsam, struct node *ptree, double tt);
int ordran (int n, double pbuf[]);
int mnmial (int n, int nclass, double p[], int rv[]);
int locate (int n, double beg, double len, double *ptr);
int make_gametes (int nsam, struct node *ptree, double tt, int newsites,
		  int ns, char **list);
int poisso (double u);

int seedit (unsigned long seed);
int cleanPRNG (void);

extern void assign_QHsiterate (double theta, int nsites,
			       struct QHparameters *QH);

unsigned maxsites = SITESINC;
double *posit;

int
main (argc, argv)
     int argc;
     char *argv[];
{
  int i;
  long count;
  unsigned long seed;
  char **list, **cmatrix ();
  FILE *pf, *fopen ();
  int segsites;
  struct parameters P;
  int getpars (int argc, char *argv[], struct parameters *);
  int printpars (char *, struct parameters);
/* 	unsigned int seedit( const char * ) ; */
  /* double t, ttot; */
/* Aug 11, 2004 Eli CHANGES */
/* 	int seQmut; -PLACED IN struct parameters */
/*	struct QHparameters *QH; -PLACED IN struct parameters */
/* 	clock_t starttime=clock(); */

  for (i = 0; i < argc; i++)
    printf ("%s ", argv[i]);
  printf ("\n");
  if (getpars (argc, argv, &P))
    {
      printpars (argv[0], P);
      exit (0);
    }
  /* seedit( "t"); */
  seed = strtoul (argv[1], NULL, 10);
  seedit (seed);
  printf ("\n%lu\n", seed);
  pf = stdout;

/* Aug 11, 2004 Eli changes */
  if (P.seQmut)
    {
      assign_QHsiterate (P.theta, P.nsites, P.QH);
      P.QH->MRCAseq =
	(char *) malloc ((unsigned) ((P.nsites) * sizeof (char)));
      P.QH->treeseq =
	(char *) malloc ((unsigned) ((2 * P.nsam - 1) * sizeof (char)));
      posit = (double *) malloc ((unsigned) (maxsites * sizeof (double)));
      if (P.QH->output)
	list = cmatrix (P.nsam, P.nsites + 1);
      else
	list = cmatrix (P.nsam, maxsites + 1);
    }
  else
    {
/* end Eli */
      if (P.theta > DBL_EPSILON)
	{
	  list = cmatrix (P.nsam, maxsites + 1);
	  posit = (double *) malloc ((unsigned) (maxsites * sizeof (double)));
	}
      else
	{
	  list = cmatrix (P.nsam, P.segsites + 1);
	  posit =
	    (double *) malloc ((unsigned) (P.segsites * sizeof (double)));
	}
    }				//Aug 11, 2004 Eli changes

  count = 0;
  while (P.howmany - count++)
    {
/* 	    segsites = gensam(P, list, &t,&ttot); */
      segsites = gensam (P, list);
      fprintf (pf, "\n//");
      fprintf (pf, "\nsegsites: %d", segsites);
/* 	    fprintf(pf," tMRCA: %9.7f",t); */
/* 	    fprintf(pf," ttot: %9.7f",ttot); */
/* Aug 11, 2004 Eli changes */
      if ((P.seQmut))
	{
	  /*QHoption- output nt frequencies */
	  fprintf (pf, "\nfreqACGT:");
	  for (i = 0; i < 4; i++)
	    fprintf (pf, " %6.4f", P.QH->obsfreq[i]);
	  /*QHoption- output polymorphism positions (format like integers) */
	  if ((segsites > 0))
	    fprintf (pf, "\npositions: ");
	  for (i = 0; i < segsites; i++)
	    fprintf (pf, "%6.0f ", posit[i]);
	}
      else if (!P.seQmut)
	{
/* end Eli */
	  if (segsites > 0)
	    fprintf (pf, "\npositions: ");
	  for (i = 0; i < segsites; i++)
	    {
	      // if( i%10 == 0 ) fprintf(pf,"\n");
	      fprintf (pf, "%6.4lf ", posit[i]);
	    }
	}			// Aug 11 2004 Eli
      fprintf (pf, "\n");
      if (segsites > 0)
	for (i = 0; i < P.nsam; i++)
	  {
	    fprintf (pf, "%s\n", list[i]);
	  }
    }
  fprintf (pf, "\n");
/* 	seedit( "end" ); */
  cleanPRNG ();
/* 	fprintf(pf,"\nCPUtime: %4.2f\n",(double)(clock()-starttime)/CLOCKS_PER_SEC); */
  exit (EXIT_SUCCESS);
}


/* int gensam(struct parameters P, char **list, double*pt,double *pttot) */
int
gensam (struct parameters P, char **list)
{
  int nsegs, i, k, seg, ns, start, end, len, segsit, *ss;
  struct segl *seglst, *segtre_mig (struct parameters, int *);	/* used to be: [MAXSEG];  */
  double nsinv, tseg, tt, *pk, ttime ();

  nsinv = 1. / P.nsites;
  seglst = segtre_mig (P, &nsegs);

/* 1 Nov 01 Eli changes */
  if (P.seQmut)
    {
      ns = 0;
      for (i = 0; i < 4; i++)
	P.QH->obsfreq[i] = 0.;
      createrandomseQ (P.QH->MRCAseq, P.nsites, P.QH);
      for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++)
	{
	  start = seglst[seg].beg;
	  end =
	    (k < nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : P.nsites - 1);
	  assign_ndes (P.nsam, seglst[seg].ptree);
	  segsit =
	    makeQHgametes (P.nsam, start, end, seglst[seg].ptree, ns, P.QH,
			   list);
	  free (seglst[seg].ptree);
	  ns += segsit;
	}
      for (i = 0; i < 4; i++)
	P.QH->obsfreq[i] /= (P.nsam * P.nsites);
    }
  else
    {
/* end Eli */
      if (P.theta > DBL_EPSILON)
	{
	  ns = 0;
	  for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++)
	    {
/* 	      if (seg==0) *pt = (seglst[0].ptree + 2* (P.nsam) -2)->time; */
	      end =
		(k <
		 nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : P.nsites - 1);
	      start = seglst[seg].beg;
	      len = end - start + 1;
	      tseg = (double) len *(P.theta / (double) P.nsites);
	      tt = ttime (seglst[seg].ptree, P.nsam);
/* 	      if (seg==0) *pttot = tt; */
	      segsit = poisso (tseg * tt);
	      if ((segsit + ns) >= maxsites)
		{
		  maxsites = segsit + ns + SITESINC;
		  posit =
		    (double *) realloc (posit, maxsites * sizeof (double));
		  biggerlist (P.nsam, list);
		}
	      make_gametes (P.nsam, seglst[seg].ptree, tt, segsit, ns, list);
	      free (seglst[seg].ptree);
	      locate (segsit, start * nsinv, len * nsinv, posit + ns);
	      ns += segsit;
	    }
	}
      else
	{
	  pk = (double *) malloc ((unsigned) (nsegs * sizeof (double)));
	  ss = (int *) malloc ((unsigned) (nsegs * sizeof (int)));
	  if ((pk == NULL) || (ss == NULL))
	    perror ("malloc error. gensam.2");
	  tt = 0.0;
	  for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++)
	    {
	      end =
		(k <
		 nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : P.nsites - 1);
	      start = seglst[seg].beg;
	      len = end - start + 1;
	      tseg = (double) len / (double) P.nsites;
	      pk[k] = ttime (seglst[seg].ptree, P.nsam) * tseg;
	      tt += pk[k];
	    }
	  for (k = 0; k < nsegs; k++)
	    pk[k] /= tt;
	  mnmial (P.segsites, nsegs, pk, ss);
	  ns = 0;
	  for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++)
	    {
	      end =
		(k <
		 nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : P.nsites - 1);
	      start = seglst[seg].beg;
	      len = end - start + 1;
	      tseg = (double) len / (double) P.nsites;
	      make_gametes (P.nsam, seglst[seg].ptree, tt * pk[k] / tseg,
			    ss[k], ns, list);
	      free (seglst[seg].ptree);
	      locate (ss[k], start * nsinv, len * nsinv, posit + ns);
	      ns += ss[k];
	    }
	  free (pk);
	  free (ss);
	}
    }				// Aug 11 2004  Eli
  for (i = 0; i < P.nsam; i++)
    list[i][(P.seQmut) && (P.QH->output) ? P.nsites : ns] = '\0';	// Aug 11 2004  Eli
  return (ns);
}


int
biggerlist (nsam, list)
     int nsam;
     char **list;
{
  int i;
/*  fprintf(stderr,"maxsites: %d\n",maxsites);  */
  for (i = 0; i < nsam; i++)
    {
      list[i] = (char *) realloc (list[i], maxsites * sizeof (char));
      if (list[i] == NULL)
	perror ("realloc error. bigger");
    }
  return 0;
}


/* allocates space for gametes (character strings) */
char **
cmatrix (nsam, len)
     int nsam, len;
{
  int i;
  char **m;
  if (!(m = (char **) malloc ((unsigned) nsam * sizeof (char *))))
    perror ("alloc error in cmatrix");
  for (i = 0; i < nsam; i++)
    {
      if (!(m[i] = (char *) malloc ((unsigned) len * sizeof (char))))
	perror ("alloc error in cmatric. 2");
    }
  return (m);
}


int
locate (n, beg, len, ptr)
     int n;
     double beg, len, *ptr;
{
  int i;
  ordran (n, ptr);
  for (i = 0; i < n; i++)
    ptr[i] = beg + ptr[i] * len;

  return 0;
}




/************ make_gametes.c  *******************************************
*
*
*****************************************************************************/
#define STATE1 '1'
#define STATE2 '0'

int
make_gametes (nsam, ptree, tt, newsites, ns, list)
     int nsam, newsites, ns;
     struct node *ptree;
     char **list;
     double tt;
{
  int tip, j, node;
  for (j = ns; j < ns + newsites; j++)
    {
      node = pickb (nsam, ptree, tt);
      for (tip = 0; tip < nsam; tip++)
	{
	  if (tdesn (ptree, tip, node))
	    list[tip][j] = STATE1;
	  else
	    list[tip][j] = STATE2;
	}
    }
  return 0;
}


/* 1 Nov 01 Eli changes */
/* Eli new functions */

int
makeQHgametes (int n, int start, int end, struct node *ptree, int ns,
	       struct QHparameters *QH, char **list)
{
  int ni, site, mut, segsites, i, baseint (char), segbool;
  char intbase (int);
  double muttime, tt, newttime (struct node *ptree, int n, int nodenum),
    ttthisbranch;
  double urandev, ran1 ();
  double rate, calcRate (char, struct QHparameters *QH), rategammafactor =
    0, gammadev (double alpha);
  int nlineages, baselineages[4], randomlineage;
  double baseRate[4], lineageratesum, lineageprob;
  void mutateBelow (int n, struct node *tree, int node, char *seqs);
  int biggerlist (int n, char **list);

/*   length = end-start+1;  nodenum = 2*n - 1; */
  segsites = 0;
  tt = newttime (ptree, n, 2 * n - 1);
  for (i = 0; i < 4; i++)
    baseRate[i] = calcRate (intbase (i), QH) * QH->siterate;
  site = start;
  while (site <= end)
    {
      while ((urandev = ran1 ()) == 0);
      muttime = -log (urandev);
      QH->treeseq[2 * n - 2] = QH->MRCAseq[site];
      rate = baseRate[baseint (QH->treeseq[2 * n - 2])];
      if (QH->alpha)
	rate *= (rategammafactor = gammadev (QH->alpha));
      while (muttime > (rate * tt))
	{
	  for (ni = 0; ni < n; ni++)
	    QH->treeseq[ni] = QH->treeseq[2 * n - 2];
	  QH->obsfreq[baseint (QH->treeseq[0])] += n;
	  /* update site, rate */
	  muttime -= (rate * tt);
	  if (++site > end)
	    break;
	  QH->treeseq[2 * n - 2] = QH->MRCAseq[site];
	  rate = baseRate[baseint (QH->treeseq[2 * n - 2])];
	  if (QH->alpha)
	    rate *= (rategammafactor = gammadev (QH->alpha));
	}
      if (site <= end)
	{
	  muttime /= rate;
/*       fprintf(stderr,"\nsite  %3d  ",site); */
/*       fprintf(stderr,"muttime %7.4f tt %6.4f",muttime, tt); */
	  if ((muttime <= 0.) || (muttime >= tt))
	    perror ("mkQHgamR: mutate site initial muttime calc wrong");
	  /* mutate this site */
	  nlineages = 2;
	  for (i = 0; i < 4; i++)
	    baselineages[i] = (i == baseint (QH->treeseq[2 * n - 2])) ? 2 : 0;
	  mutateBelow (n, ptree, (2 * n - 2), QH->treeseq);
/*  printf("\ntt: %10.8f", tt); */
	  for (ni = (2 * n - 3), mut = 0; ni >= 0; ni--)
	    {
/* 	fprintf(stderr,"\n          muttime: %10.8f", muttime); */
/* 	fprintf(stderr,"\nabv node time: %10.8f, node time: %10.8f, %1d lineages",(ptree+ni+1)->time,(ptree+ni)->time,nlineages); */
	      ttthisbranch =
		(((ptree + ni + 1)->time - (ptree + ni)->time) * nlineages);
/* 	fprintf(stderr,"\nsite %1d, %1d lineages, ttthisbranch: %10.8f",site, nlineages,ttthisbranch); */
	      while ((ttthisbranch -= muttime) > 0)
		{
		  /* pick a lineage at random and mutate it */
		  for (i = 0, lineageratesum = 0.; i < 4; i++)
		    lineageratesum += baseRate[i] * baselineages[i];
		  while ((urandev = ran1 ()) == 1.);
		  for (randomlineage = ni, lineageprob = 0.;
		       randomlineage >= 0; randomlineage--)
		    {
		      if ((ptree + ((ptree + randomlineage)->abv))->time >
			  (ptree + ni)->time)
			if (urandev <
			    (lineageprob +=
			     (baseRate[baseint (QH->treeseq[randomlineage])] /
			      lineageratesum)))
			  {
			    mut = 1;
			    break;
			  }
		    }
		  if (randomlineage < 0)
		    perror
		      ("\nmakeQHgametesR: mutate site pick randomlineage error");
		  baselineages[baseint (QH->treeseq[randomlineage])]--;
		  QH->treeseq[randomlineage] =
		    mutateRate (QH->treeseq[randomlineage], QH);
		  mutateBelow (n, ptree, randomlineage, QH->treeseq);
		  baselineages[baseint (QH->treeseq[randomlineage])]++;
		  /* recalculate muttime(nlineages,baselineages,baseRates) */
		  for (i = 0, rate = 0.; i < 4; i++)
		    rate += baseRate[i] * baselineages[i];
		  if (QH->alpha)
		    rate *= rategammafactor;
		  rate /= nlineages;
		  while ((urandev = ran1 ()) == 0.);
		  muttime = -log (urandev) / rate;
		}
	      if (fabs ((ptree + ni)->time) < 1e-6)
		break;
	      if ((ptree + ni)->ndes > 1)
		{
		  nlineages++;
		  baselineages[baseint (QH->treeseq[ni])]++;
		}
	      else
		{
		  nlineages--;
		  baselineages[baseint (QH->treeseq[ni])]--;
		}
	      for (i = 0, rate = 0.; i < 4; i++)
		rate += baseRate[i] * baselineages[i];
	      if (QH->alpha)
		rate *= rategammafactor;
	      rate /= nlineages;
	      if (mut)
		{
		  while ((urandev = ran1 ()) == 0.);
		  muttime = -log (urandev) / rate;
		}
	      else
		muttime = -(ttthisbranch);
	    }
	  for (i = 0; i < n; i++)
	    {
	      QH->obsfreq[baseint (QH->treeseq[i])]++;
	      if (QH->output)
		list[i][site] = QH->treeseq[i];
	    }
	  for (i = 1, segbool = 0; i < n; i++)
	    if (QH->treeseq[i] != QH->treeseq[0])
	      {
		segbool = 1;
		break;
	      }
/*   fprintf(stderr,"\nHERE "); */
/*   fprintf(stderr,"HERE2\n"); */
	  if (segbool)
	    {
	      if (((ns + segsites + 1) >= maxsites))
		{
/* 	  printf("\nREALLOC "); */
		  maxsites = ns + segsites + SITESINC;
		  posit =
		    (double *) realloc (posit, maxsites * sizeof (double));
		  if (!(QH->output))
		    biggerlist (n, list);
/* 	  printf("done\n"); */
		}
/* 	printf("- updating-"); */
	      posit[ns + segsites] = site;
	      if (!QH->output)
		for (i = 0; i < n; i++)
		  list[i][ns + segsites] = QH->treeseq[i];
	      segsites++;
/* 	printf("done "); */
	    }
	  site++;
	}
    }
  return (segsites);
}


void
assign_ndes (int n, struct node *tree)
{
  int i, ii, done;
  for (i = 0; i < (2 * n - 1); i++)
    (tree + i)->ndes = (i < n) ? 1 : 0;
  for (i = n; i < (2 * n - 1); i++)
    for (ii = (i - 1), done = 0; (ii >= 0) && (done < 2); ii--)
      if ((tree + ii)->abv == i)
	{
	  (tree + i)->ndes += (tree + ii)->ndes;
	  done++;
	}
}


/***  newttime.c : Returns the total time, tree can have extra nodes and tips with non-zero time. **/
double
newttime (ptree, nsam, nodenum)
     struct node *ptree;
     int nsam, nodenum;
{
  int i, nlineages = 2;
  double t = 0.;
  for (i = (nodenum - 2); i >= 0; i--)
    {
      t += ((ptree + (ptree + i)->abv)->time - (ptree + i)->time) * nlineages;
      if ((ptree + i)->ndes > 1)
	nlineages++;
      else if ((ptree + i)->ndes == 1)
	nlineages--;
      if (fabs ((ptree + i)->time) < 1e-6)
	break;
    }
  return (t);
}

/* end Eli new functions */
/* end Eli */


/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/
double
ttime (ptree, nsam)
     struct node *ptree;
     int nsam;
{
  double t;
  int i;
  t = (ptree + 2 * nsam - 2)->time;
  for (i = nsam; i < 2 * nsam - 1; i++)
    t += (ptree + i)->time;
  return (t);
}

/***  pickb : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   ****/
int
pickb (nsam, ptree, tt)
     int nsam;
     struct node *ptree;
     double tt;
{
  double x, y, ran1 ();
  int i;
  x = ran1 () * tt;
  for (i = 0, y = 0; i < 2 * nsam - 2; i++)
    {
      y += (ptree + (ptree + i)->abv)->time - (ptree + i)->time;
      if (y >= x)
	return (i);
    }
  return (i);
}

/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/
int
tdesn (ptree, tip, node)
     struct node *ptree;
     int tip, node;
{
  int k;
  for (k = tip; k < node; k = (ptree + k)->abv);
  if (k == node)
    return (1);
  else
    return (0);
}


/* pick2()  */
int
pick2 (n, i, j)
     int n, *i, *j;
{
  double ran1 ();
  *i = n * ran1 ();
  while ((*j = n * ran1 ()) == *i)
    ;
  return 0;
}

/**** ordran.c  ***/
int
ordran (n, pbuf)
     int n;
     double pbuf[];
{
  ranvec (n, pbuf);
  order (n, pbuf);
  return 0;
}


int
mnmial (n, nclass, p, rv)
     int n, nclass, rv[];
     double p[];
{
  double ran1 ();
  double x, s;
  int i, j = 0;

  for (i = 0; i < nclass; i++)
    rv[i] = 0;
  for (i = 0; i < n; i++)
    {
      x = ran1 ();
      j = 0;
      s = p[0];
      while ((x > s) && (j < (nclass - 1)))
	s += p[++j];
      rv[j]++;
    }
  return (j);
}

int
order (n, pbuf)
     int n;
     double pbuf[];
{
  int gap, i, j;
  double temp;
  for (gap = n / 2; gap > 0; gap /= 2)
    for (i = gap; i < n; i++)
      for (j = i - gap; j >= 0 && pbuf[j] > pbuf[j + gap]; j -= gap)
	{
	  temp = pbuf[j];
	  pbuf[j] = pbuf[j + gap];
	  pbuf[j + gap] = temp;
	}

  return 0;
}


int
ranvec (n, pbuf)
     int n;
     double pbuf[];
{
  int i;
  double ran1 ();
  for (i = 0; i < n; i++)
    pbuf[i] = ran1 ();
  return 0;
}



int
poisso (u)
     double u;
{
  double cump, ru, ran1 (), p, gasdev ();
  int i = 1;
  if (u > 30.)
    return ((int) (0.5 + gasdev (u, u)));
  ru = ran1 ();
  p = exp (-u);
  if (ru < p)
    return (0);
  cump = p;
  while (ru > (cump += (p *= u / i)))
    i++;
  return (i);
}


/* a slight modification of crecipes version */
double
gasdev (m, v)
     double m, v;
{
  static int iset = 0;
  static float gset;
  float fac, r, v1, v2;
  double ran1 ();
  if (iset == 0)
    {
      do
	{
	  v1 = 2.0 * ran1 () - 1.0;
	  v2 = 2.0 * ran1 () - 1.0;
	  r = v1 * v1 + v2 * v2;
	}
      while (r >= 1.0);
      fac = sqrt (-2.0 * log (r) / r);
      gset = v1 * fac;
      iset = 1;
      return (m + sqrt (v) * v2 * fac);
    }
  else
    {
      iset = 0;
      return (m + sqrt (v) * gset);
    }
}
