/**********  segtre_mig.c **********************************
*
*	This subroutine uses a Monte Carlo algorithm described in
*	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
*	a history of a random sample of gametes under a neutral
*	Wright-Fisher model with recombination and geographic structure. 
*	Input parameters
*	are the sample size (nsam), the number of sites between
*	which recombination can occur (nsites), and the recombination
*	rate between the ends of the gametes (r). The function returns
*	nsegs, the number of segments the gametes were broken into
*	in tracing back the history of the gametes.  The histories of
*	these segments are passed back to the calling function in the
*	array of structures seglst[]. An element of this array,  seglst[i],
* 	consists of three parts: (1) beg, the starting point of
*	of segment i, (2) ptree, which points to the first node of the
*	tree representing the history of the segment, (3) next, which
*	is the index number of the next segment.
*	     A tree is a contiguous set of 2*nsam nodes. The first nsam
*	nodes are the tips of the tree, the sampled gametes.  The other
*	nodes are the nodes ancestral to the sampled gametes. Each node
*	consists of an "abv" which is the number of the node ancestral to
*	that node, an "ndes", which is not used or assigned to in this routine,
*	and a "time", which is the time (in units of 4N generations) of the
*	node. For the tips, time equals zero.
*	Returns a pointer to an array of segments, seglst.
*
**************************************************************************/

#include <math.h>
#include "msQHarbpop.h"

#define NL putchar('\n')
#define size_t unsigned

#define MIN(x, y) ( (x)<(y) ? (x) : (y) )

#define ERROR(message) fprintf(stderr,message),NL,exit(1)

#define SEGINC 80

extern int flag;

int nchrom, begs, nsegs;
long nlinks;
static int *nnodes = NULL;
double t, cleft, pc, lnpc;

static unsigned seglimit = SEGINC;
static unsigned maxchr;

struct seg
{
  int beg;
  int end;
  int desc;
};

struct chromo
{
  int nseg;
  int pop;
  struct seg *pseg;
};

static struct chromo *chrom = NULL;

struct node *ptree1, *ptree2;

static struct segl *seglst = NULL;


int pick2_chrom (int pop, int config[], int *pc1, int *pc2);
int ca (int nsam, int nsites, int c1, int c2);
int isseg (int start, int c, int *psg);
int links (int c);
int xover (int nsam, int ic, int is);


struct segl *
segtre_mig (struct parameters P, int *pnsegs)
/* 	int npop, inconfig[], nsites, nsam, *pnsegs, nintn ; */
/* 	double r, mig_rate,  f, track_len , *nrec, *npast, *tpast ; */
{
  int i, j, dec, pop, c1, c2, ind, rchrom;
  int migrant, source_pop, *config;
  double ran1 (), ttemp, rft, clefta;
  double prec, cin, prect, ran, prob, rdum, arg;
  int re (), cinr (), cleftr ();
/* extra variables ? */
/* 	char c; */
/* 	int i, k, seg, flagint, intn ; */
/*	double alphag, trm, tcoal, mig, coal_prob, nnm1, nnm0 ; */
/* keep old Main loop extra variables */
/* 	int i, k, seg, flagint, intn ; */
  double mig_rate, trm, tcoal, mig, oldalphag, coal_prob;
/*	double alphag, trm, tcoal, mig, coal_prob, nnm1, nnm0 ; */
/* ms-arbpop CHANGE */
  int pop2, nsam, npop, nsites, Dint, lastDint, Mpattern;
  int event, mintime (int npop, double *eventtimes, double *time), migchrom;
  double tpast, prevtpast, r, f, *alphag, **M, *eventtimes;
  void Dint_updateconfig (struct parameters, int, int **);
/* END ms-arbpop CHANGE */

/* Initialization */
  nsam = P.nsam;
  nsites = P.nsites;
  r = P.r;
  f = P.f;
  Dint = 0;
  lastDint = P.Dintn - 1;
  npop = (P.Dint)->npops;
  Mpattern = (P.Dint)->Mpattern;
  M = (P.Dint)->M;
  tpast = (P.Dint)->Dtpast;
  t = prevtpast = 0.;

  if (chrom == NULL)
    {
      maxchr = nsam + 20;
      chrom =
	(struct chromo *)
	malloc ((unsigned) (maxchr * sizeof (struct chromo)));
      if (chrom == NULL)
	perror ("malloc error. segtre");
    }
  if (nnodes == NULL)
    {
      nnodes = (int *) malloc ((unsigned) (seglimit * sizeof (int)));
      if (nnodes == NULL)
	perror ("malloc error. segtre_mig");
    }
  if (seglst == NULL)
    {
      seglst =
	(struct segl *) malloc ((unsigned) (seglimit * sizeof (struct segl)));
      if (seglst == NULL)
	perror ("malloc error. segtre_mig.c 2");
    }
  config = (int *) malloc ((unsigned) (npop + 1) * sizeof (int));
  if (config == NULL)
    perror ("malloc error. segtre. config");
  for (pop = 0; pop < npop; pop++)
    config[pop] = P.config[pop];
  for (pop = ind = 0; pop < npop; pop++)
    for (j = 0; j < config[pop]; j++, ind++)
      {
	chrom[ind].nseg = 1;
	if (!
	    (chrom[ind].pseg =
	     (struct seg *) malloc ((unsigned) sizeof (struct seg))))
	  ERROR ("calloc error. se1");
	(chrom[ind].pseg)->beg = 0;
	(chrom[ind].pseg)->end = nsites - 1;
	(chrom[ind].pseg)->desc = ind;
	chrom[ind].pop = pop;
      }
  seglst[0].beg = 0;
  if (!
      (seglst[0].ptree =
       (struct node *) calloc ((unsigned) (2 * nsam), sizeof (struct node))))
    perror ("calloc error. se2");

  nnodes[0] = nsam - 1;
  nchrom = nsam;

  alphag = (double *) malloc ((unsigned) npop * sizeof (double));
  if (alphag == NULL)
    perror ("malloc error. segtre. alphag");
  eventtimes =
    (double *) malloc ((unsigned) (npop * npop + 1) * sizeof (double));
  if (eventtimes == NULL)
    perror ("malloc error. segtre. eventtimes");
  if (tpast > DBL_EPSILON)
    for (pop = 0; pop < npop; pop++)
      if (fabs ((P.Dint)->Npast[pop] - (P.Dint)->Nrecent[pop]) < DBL_EPSILON)
	alphag[pop] = 0.;
      else
	alphag[pop] =
	  -log ((P.Dint)->Npast[pop] / (P.Dint)->Nrecent[pop]) / (tpast -
								  prevtpast);
  else
    {
      if (P.D)
	Dint_updateconfig (P, ++Dint, &config);
      tpast = (P.Dint + Dint)->Dtpast;
      npop = (P.Dint + Dint)->npops;
      alphag =
	(double *) realloc (alphag, (unsigned) (npop * sizeof (double)));
      if (alphag == NULL)
	perror ("realloc error. segtre init alphag[pop]\n");
      eventtimes =
	(double *) realloc (eventtimes,
			    (unsigned) (npop * npop + 1) * sizeof (double));
      if (eventtimes == NULL)
	perror ("realloc error. segtre init eventtimes\n");
      for (pop = 0; pop < npop; pop++)
	if ((tpast - prevtpast < DBL_EPSILON)
	    || (fabs ((P.Dint + Dint)->Npast[0] - (P.Dint + Dint)->Nrecent[0])
		< DBL_EPSILON))
	  alphag[pop] = 0.;
	else
	  alphag[pop] =
	    -log ((P.Dint + Dint)->Npast[pop] / (P.Dint +
						 Dint)->Nrecent[pop]) /
	    (tpast - prevtpast);
      Mpattern = (P.Dint + Dint)->Mpattern;
      M = (P.Dint + Dint)->M;
    }

  nlinks = ((long) (nsam)) * (nsites - 1);
  nsegs = 1;
  r /= (nsites - 1);
  if (f > DBL_EPSILON)
    pc = (P.track_len - 1.0) / P.track_len;
  else
    pc = 1.0;
  lnpc = log (pc);
  cleft = nsam * (1.0 - pow (pc, (double) (nsites - 1)));
  rft = f * P.track_len;
  if (r > DBL_EPSILON)
    rft *= r;

  // Main arbpop conditional
  if (P.D)
    {
      /* Main loop */
      while (nchrom > 1)
	{
	  /* ASSIGN recombination time TO eventtimes[0] */
	  prec = nlinks * r;
	  cin = nlinks * r * f;
	  clefta = cleft * rft;
	  prect = prec + cin + clefta;
	  event = 0;
	  if (prect > 0.0)
	    {
	      while ((rdum = ran1 ()) == 0.0);
	      eventtimes[event++] = -log (rdum) / prect;
	    }
	  else
	    eventtimes[event++] = -1.0;	/* no recombination within interval */

	  /* ASSIGN COALESCENCE AND MIGRATION TIMES TO eventimes[] */
	  for (pop = 0; pop < npop; pop++)
	    if (config[pop] > 1)
	      {
		while ((rdum = ran1 ()) == 0.0);
		if (fabs (alphag[pop]) < DBL_EPSILON)
		  {
		    eventtimes[event++] =
		      -log (rdum) * (P.Dint +
				     Dint)->Nrecent[pop] / (config[pop] *
							    (config[pop] -
							     1));
		  }
		else
		  {
		    arg =
		      1. - log (rdum) * (P.Dint +
					 Dint)->Nrecent[pop] *
		      exp (-alphag[pop] * (t - prevtpast)) * alphag[pop] /
		      (config[pop] * (config[pop] - 1));
		    if (arg > 0.)
		      eventtimes[event++] = log (arg) / alphag[pop];
		    else
		      eventtimes[event++] = -1.;	/* no coalescent within interval */
		  }
	      }
	    else
	      eventtimes[event++] = -1.;
	  if (event != 1 + npop)
	    {
	      fprintf (stderr,
		       "error- segtre- coal eventtimes[] assignment, event=%1d, (1+npop)=%1d\n",
		       event, (1 + npop));
	      exit (1);
	    }
	  /* migration times */
	  if (npop > 1)
	    for (pop = 0; pop < npop; pop++)
	      if (config[pop] > 0)
		{
		  for (pop2 = 0; pop2 < npop; pop2++)
		    if (pop2 != pop)
		      if (M[pop][pop2] > DBL_EPSILON)
			{
			  while ((rdum = ran1 ()) == 0.);
			  if (Mpattern == 0)	/* migration not dependent on pop or source pop */
			    eventtimes[event++] =
			      -log (rdum) / (config[pop] * M[pop][pop2]);
			  else if (Mpattern == 1)	/* migration dependent on pop but not source pop */
			    if (fabs (alphag[pop]) < DBL_EPSILON)
			      eventtimes[event++] =
				-log (rdum) * ((P.Dint + Dint)->Nrecent[pop] *
					       exp (-alphag[pop] *
						    (t -
						     prevtpast))) /
				(config[pop] * M[pop][pop2]);
			    else
			      {
				arg =
				  1. -
				  log (rdum) *
				  ((P.Dint +
				    Dint)->Nrecent[pop] * exp (-alphag[pop] *
							       (t -
								prevtpast))) *
				  alphag[pop] / (config[pop] * M[pop][pop2]);
				if (arg > 0.)
				  eventtimes[event++] =
				    log (arg) / alphag[pop];
				else
				  eventtimes[event++] = -2.;	/* no migration[pop][pop2] within interval */
			      }
			  else	/* if ( Mpattern == 2 ) */
			  /* migration dependent on pop and source pop */
			  if (fabs (alphag[pop] - alphag[pop2]) < DBL_EPSILON)
			    eventtimes[event++] =
			      -log (rdum) * ((P.Dint + Dint)->Nrecent[pop] *
					     exp (-alphag[pop] *
						  (t -
						   prevtpast))) /
			      (config[pop] * M[pop][pop2] *
			       ((P.Dint +
				 Dint)->Nrecent[pop2] * exp (-alphag[pop2] *
							     (t -
							      prevtpast))));
			  else
			    {
			      arg =
				1. -
				log (rdum) * ((P.Dint + Dint)->Nrecent[pop] *
					      exp (-alphag[pop] *
						   (t -
						    prevtpast))) *
				(alphag[pop] -
				 alphag[pop2]) / (config[pop] * M[pop][pop2] *
						  ((P.Dint +
						    Dint)->Nrecent[pop2] *
						   exp (-alphag[pop2] *
							(t - prevtpast))));
			      if (arg > 0.)
				eventtimes[event++] =
				  log (arg) / (alphag[pop] - alphag[pop2]);
			      else
				eventtimes[event++] = -2.;	/* no migration[pop][pop2] within interval */
			    }
			}
		      else
			eventtimes[event++] = -2.;	/* else- no migration within interval */
		}
	      else
		for (pop2 = 0; pop2 < npop; pop2++)
		  if (pop2 != pop)
		    eventtimes[event++] = -2.;	/* else- no migration within interval */
	  if (event != 1 + npop * npop)
	    {
	      fprintf (stderr,
		       "error- segtre- mig eventtimes[] assignment, event=%1d, (1+npop*npop)=%1d\n",
		       event, (1 + npop * npop));
	      exit (1);
	    }

	  /* GET MINIMUM TIME AND INCREMENT  t, event=1,2,3(coal,mig,rec) */
	  // WAS event=mintime(npop,eventtimes,&t);
	  for (i = 0; i < (1 + npop * npop); i++)
	    if (eventtimes[i] > DBL_EPSILON)
	      {
		event = i;
		break;
	      }
	  for (i++; i < (1 + npop * npop); i++)
	    if ((eventtimes[i] > DBL_EPSILON)
		&& (eventtimes[i] < eventtimes[event]))
	      event = i;
	  if (event == (1 + npop * npop))
	    t = tpast + 1.;	// no event in interval
	  else
	    t += eventtimes[event];

	  /* EITHER GO TO NEXT INTERVAL, OR DO SOMETHING */
	  if ((Dint < lastDint) && (t > tpast))
	    {
	      Dint_updateconfig (P, ++Dint, &config);
	      npop = (P.Dint + Dint)->npops;
	      t = prevtpast = tpast;
	      tpast = (P.Dint + Dint)->Dtpast;
	      alphag =
		(double *) realloc (alphag,
				    (unsigned) (npop * sizeof (double)));
	      if (alphag == NULL)
		perror ("realloc error. segtre alphag[pop]\n");
	      eventtimes =
		(double *) realloc (eventtimes,
				    (unsigned) (npop * npop +
						1) * sizeof (double));
	      if (eventtimes == NULL)
		perror ("realloc error. segtre init eventtimes\n");
	      for (pop = 0; pop < npop; pop++)
		if ((tpast - prevtpast < DBL_EPSILON)
		    ||
		    (fabs
		     ((P.Dint + Dint)->Npast[0] -
		      (P.Dint + Dint)->Nrecent[0]) < DBL_EPSILON))
		  alphag[pop] = 0.;
		else
		  alphag[pop] =
		    -log ((P.Dint + Dint)->Npast[pop] / (P.Dint +
							 Dint)->
			  Nrecent[pop]) / (tpast - prevtpast);
	      Mpattern = (P.Dint + Dint)->Mpattern;
	      M = (P.Dint + Dint)->M;
	    }

	  else
	    {			/* OR DO SOMETHING */
	      //  Main loop
	      //                if (Dint==0)
/* 	            printf("\ndint: %1d tcoal %9.7f alphag %9.7f\n", Dint, nchrom, eventtimes[1], alphag[0]); */
/* 	            printf("\nnchrom=%1d nchrom[].pop=",nchrom); */
/* 		    for (i=0;i<nchrom;i++) printf("%1d,",chrom[i].pop); */
/* 		    printf("\nevent: %1d eventtimes[]=",event); */
/* 		    for (i=0;i<1+npop*npop;i++) printf("%6.4f,",eventtimes[i]); */
/* 		    printf("\n"); */


	      if (event-- == 0)
		{		/* recombination event ? */
		  if ((ran = ran1 ()) < (prec / prect))
		    {		/* Xover */
		      rchrom = re (nsam);
		      config[chrom[rchrom].pop] += 1;
		    }
		  else if (ran < (prec + clefta) / prect)
		    {		/* cleft event */
		      rchrom = cleftr (nsam);
		      config[chrom[rchrom].pop] += 1;
		    }
		  else
		    {		/* cin event */
		      rchrom = cinr (nsam, nsites);
		      if (rchrom >= 0)
			config[chrom[rchrom].pop] += 1;
		    }
		}
	      else if (event < npop)
		{		/* coalescent event */
		  pick2_chrom (event, config, &c1, &c2);	/* c1 and c2 are chrom's to coalesce */
		  dec = ca (nsam, nsites, c1, c2);
		  config[event] -= dec;
		}
	      else
		{		/* migration event */
/* 	      fprintf(stderr," \n 2got here fuck %4.2f\n", mig_rate); */
/* 	    printf(" \n got here1 config[]="); */
/* 	    for (i=0;i<npop;i++) printf("%1d,", config[i]); */
		  event -= npop;
		  pop = 0;
/* 	    fprintf(stderr,"\n outside while event=%1d>npop=%1d pop=%1d\n",event,npop,pop); */
		  while (event >= (npop - 1))
		    {
/* 	      fprintf(stderr,"\n got inside while event=%1d>npop=%1d pop=%1d\n",event,npop,pop); */
		      pop++;
		      event -= (npop - 1);
		    }
		  source_pop = (event < pop) ? event : event + 1;
/* 	    printf(" pop=%1d sourcepop=%1d ",pop,source_pop); */

		  migrant = (int) ((double) config[pop] * ran1 ());
		  if (migrant == config[pop])
		    migrant--;

/* 	    printf("migrant=%1d \n", migrant); */

		  config[pop] -= 1;
		  config[source_pop] += 1;
		  for (migchrom = 0; migchrom < nchrom; migchrom++)
		    if (chrom[migchrom].pop == pop)
		      {
			if (migrant == 0)
			  break;
			migrant--;
		      }
/* 	      fprintf(stderr," \n got here2 migrant=%1d nchrom=%1d migchrom=%1d", migrant, nchrom, migchrom); */
		  if (migchrom >= nchrom)
		    perror ("error- segtre- migchrom assignment\n");
		  chrom[migchrom].pop = source_pop;

/* 	      fprintf(stderr," \n 3got here \n"); */
		}
/* 	    fprintf(stderr, "\nend of DO SOMETHING, nchrom=%1d\n", nchrom); */
	    }

/*    fprintf(stderr, "end inside while(nchrom>1), nchrom=%1d\n", nchrom); */
	}
/*    fprintf(stderr, "end outside while(nchrom>1), nchrom=%1d\n", nchrom); */


    }
  else
    {				// !P.Dint

/* Old Main loop */
      mig_rate = (npop > 1) ? M[0][1] : 0.;
      oldalphag = alphag[0];
      while (nchrom > 1)
	{
	  prec = nlinks * r;
	  cin = nlinks * f;
	  if (r > DBL_EPSILON)
	    cin *= r;
	  clefta = cleft * rft;
	  prect = prec + cin + clefta;
	  mig = nchrom * mig_rate;
	  tcoal = -1.0;
	  trm = -1.0;
	  if ((prect + mig) > 0.0)
	    {
	      while ((rdum = ran1 ()) == 1.0);
	      trm = -log (1. - rdum) / (prect + mig);
	    }

	  coal_prob = 0.;
	  for (pop = 0; pop < npop; pop++)
	    coal_prob += ((double) config[pop]) * (config[pop] - 1.);
	  if (coal_prob > 0.0)
	    {
	      while ((rdum = ran1 ()) == 1.0);
	      if (oldalphag == 0)
		tcoal =
		  -log (1. - rdum) * (P.Dint + Dint)->Nrecent[0] / coal_prob;
	      else
		{
		  arg =
		    1. - oldalphag * (P.Dint +
				      Dint)->Nrecent[0] * exp (-oldalphag *
							       (t -
								prevtpast)) *
		    log (1. - rdum) / coal_prob;
		  if (arg <= 0.0)
		    tcoal = tpast + 1.0;	/* no coalescent within interval */
		  else
		    tcoal = log (arg) / oldalphag;
		}
	    }

	  if (tcoal >= 0.0)
	    {
	      if (trm >= 0.0)
		ttemp = MIN (trm, tcoal);
	      else
		{
		  ttemp = tcoal;
		  trm = tcoal + 1.0;
		}
	    }
	  else
	    {
	      ttemp = trm;
	      tcoal = trm + 1.0;
	    }
	  //  Old Main loop
	  //        if (Dint==0)
	  //         printf("\ndint: %1d nchrom: %1d tcoal %9.7f alphag %9.7f\n", Dint, nchrom, tcoal, oldalphag);

	  if ((Dint < lastDint) && (t + ttemp > tpast))
	    {
	      /* UPDATE */
	      t = prevtpast = tpast;
	      ++Dint;
	      tpast = (P.Dint + Dint)->Dtpast;
	      if (tpast - prevtpast < DBL_EPSILON)
		oldalphag = 0.0;
	      else
		oldalphag =
		  -log ((P.Dint + Dint)->Npast[0] / (P.Dint +
						     Dint)->Nrecent[0]) /
		  (tpast - prevtpast);
	      npop = (P.Dint + Dint)->npops;
	      mig_rate = (npop > 1) ? (P.Dint + Dint)->M[0][1] : 0.;
	    }

	  else
	    {

	      t += ttemp;
	      if (trm < tcoal)
		{
		  if ((ran = ran1 ()) < (prec / (prect + mig)))
		    {
/* 	      printf("\nevent-rec dint: %1d prevtpast %9.7f t %9.7f tpast %9.7f\n", Dint, prevtpast, t, tpast); */
		      /*recombination */
		      rchrom = re (nsam);
		      config[chrom[rchrom].pop] += 1;
		    }
		  else if (ran < (prec + clefta) / (prect + mig))
		    {		/*  cleft event */
		      rchrom = cleftr (nsam);
		      config[chrom[rchrom].pop] += 1;
		    }
		  else if (ran < prect / (prect + mig))
		    {		/* cin event */
		      rchrom = cinr (nsam, nsites);
		      if (rchrom >= 0)
			config[chrom[rchrom].pop] += 1;
		    }
		  else
		    {		/* migration event */
/* 	      printf("\nevent-mig dint: %1d prevtpast %9.7f t %9.7f tpast %9.7f", Dint, prevtpast, t, tpast); */

/* 		    fprintf(stderr," \n got here1 config[]="); */
/* 		    for (i=0;i<npop;i++) fprintf(stderr,"%1d,", config[i]); */

		      migrant = nchrom * ran1 ();
		      while ((source_pop =
			      npop * ran1 ()) == chrom[migrant].pop);

/* 		    fprintf(stderr," pop=%1d sourcepop=%1d",pop,source_pop); */
/* 		    fprintf(stderr,"migrant=%1d ", migrant); */

		      config[chrom[migrant].pop] -= 1;
		      config[source_pop] += 1;
		      chrom[migrant].pop = source_pop;


/* 	      fprintf(stderr," \n 3got here \n"); */
		    }
		}

	      else
		{		/* coalescent event */
/* 	      printf("\nevent-coal dint: %1d prevtpast %9.7f t %9.7f tpast %9.7f", Dint, prevtpast, t, tpast); */

/* 		    fprintf(stderr," \n got here1 config[]="); */
/* 		    for (i=0;i<npop;i++) fprintf(stderr,"%1d,", config[i]); */

		  /* pick the two, c1, c2  */
		  ran = ran1 ();
		  prob = 0.0;
		  for (pop = 0; pop < npop; pop++)
		    {
		      prob +=
			((double) config[pop]) * (config[pop] -
						  1.) / coal_prob;
		      if (ran < prob)
			break;
		    }
		  if (pop == npop)
		    pop = npop - 1;
		  pick2_chrom (pop, config, &c1, &c2);	/* c1 and c2 are chrom's to coalesce */
		  dec = ca (nsam, nsites, c1, c2);
		  config[pop] -= dec;
		}
	    }
	}

    }				// D or not D conditional
  *pnsegs = nsegs;
  free (config);
  return (seglst);
}


/* ms-arbpop CHANGE */
void
Dint_updateconfig (struct parameters P, int Dint, int **pconfig)
{
  int i, ind, index, pop, npop, *oldconfig, mnmial_index (int, double *);
  oldconfig =
    (int *) malloc ((unsigned) (P.Dint + Dint - 1)->npops * sizeof (int));
  for (pop = 0; pop < (P.Dint + Dint - 1)->npops; pop++)
    oldconfig[pop] = (*pconfig)[pop];
  npop = (P.Dint + Dint)->npops;
  *pconfig = (int *) realloc (*pconfig, (unsigned) npop * sizeof (int));
  for (pop = 0; pop < npop; pop++)
    (*pconfig)[pop] = 0;
  for (pop = 0; pop < (P.Dint + Dint - 1)->npops; pop++)
    if (npop > 1)
      for (i = ind = 0; i < oldconfig[pop]; i++)
	{
	  index = mnmial_index (npop, (P.Dint + Dint)->a[pop]);
	  ++(*pconfig)[index];
	  for (;; ind++)
	    if (chrom[ind].pop == pop)
	      break;
	  chrom[ind].pop = index;
	}
    else
      {
	(*pconfig)[0] += oldconfig[pop];
	for (ind = 0; ind < nchrom; ind++)
	  chrom[ind].pop = 0;
      }
  free (oldconfig);
}

// from Dick's mnmial (in ms.c)
// except returns index (class) to be multinomially incremented
int
mnmial_index (int nclass, double *p)
{
  double ran1 ();
  double x, s;
  int j;
  /* int i; */
/* 	for(i=0; i<n ; i++) { */
  x = ran1 ();
  j = 0;
  s = p[0];
  while ((x > s) && (j < (nclass - 1)))
    s += p[++j];
/* 	   } */
  return (j);
}

/* END ms-arbpop CHANGE */


/******  recombination subroutine ***************************
   Picks a chromosome and splits it in two parts. If the x-over point
   is in a new spot, a new segment is added to seglst and a tree set up
   for it.   ****/
int
re (nsam)
     int nsam;
{
  struct seg *pseg = NULL;
  int el, lsg, lsgm1, ic, is, spot;
  double ran1 ();
/* First generate a random x-over spot, then locate it as to chrom and seg. */
  spot = nlinks * ran1 () + 1.;
  /* get chromosome # (ic)  */
  for (ic = 0; ic < nchrom; ic++)
    {
      lsg = chrom[ic].nseg;
      lsgm1 = lsg - 1;
      pseg = chrom[ic].pseg;
      el = ((pseg + lsgm1)->end) - (pseg->beg);
      if (spot <= el)
	break;
      spot -= el;
    }
  is = pseg->beg + spot - 1;
  xover (nsam, ic, is);
  return (ic);
}

int
cleftr (int nsam)
{
  struct seg *pseg;
  int ic, is;
  double ran1 (), x, sum, len;
  while ((x = cleft * ran1 ()) == 0.0);
  sum = 0.0;
  ic = -1;
  while (sum < x)
    {
      sum += 1.0 - pow (pc, links (++ic));
    }
  pseg = chrom[ic].pseg;
  len = links (ic);
  is =
    pseg->beg + floor (1.0 +
		       log (1.0 - (1.0 - pow (pc, len)) * ran1 ()) / lnpc) -
    1;
  xover (nsam, ic, is);
  return (ic);
}

int
cinr (int nsam, int nsites)
{
  struct seg *pseg = NULL;
  int len, el, lsg, lsgm1 = -1, ic, is, spot, endic;
  double ran1 ();

/* First generate a random x-over spot, then locate it as to chrom and seg. */
  spot = nlinks * ran1 () + 1.;
  /* get chromosome # (ic)  */
  for (ic = 0; ic < nchrom; ic++)
    {
      lsg = chrom[ic].nseg;
      lsgm1 = lsg - 1;
      pseg = chrom[ic].pseg;
      el = ((pseg + lsgm1)->end) - (pseg->beg);
      if (spot <= el)
	break;
      spot -= el;
    }
  is = pseg->beg + spot - 1;
  endic = (pseg + lsgm1)->end;
  xover (nsam, ic, is);

  len = floor (1.0 + log (ran1 ()) / lnpc);
  if (is + len >= endic)
    return (ic);
  if (is + len < (chrom[nchrom - 1].pseg)->beg)
    {
      ca (nsam, nsites, ic, nchrom - 1);
      return (-1);
    }
  xover (nsam, nchrom - 1, is + len);
  ca (nsam, nsites, ic, nchrom - 1);
  return (ic);
}

int
xover (int nsam, int ic, int is)
{
  struct seg *pseg, *pseg2;
  int i, lsg, lsgm1, newsg, jseg, k, in;
  double ran1 (), len;

  pseg = chrom[ic].pseg;
  lsg = chrom[ic].nseg;
  len = (pseg + lsg - 1)->end - pseg->beg;
  cleft -= 1 - pow (pc, len);
  /* get seg # (jseg)  */
  for (jseg = 0; is >= (pseg + jseg)->end; jseg++);
  if (is >= (pseg + jseg)->beg)
    in = 1;
  else
    in = 0;
  newsg = lsg - jseg;

  /* copy last part of chrom to nchrom  */
  nchrom++;
  if (nchrom >= maxchr)
    {
      maxchr += 20;
      chrom =
	(struct chromo *) realloc (chrom,
				   (unsigned) (maxchr *
					       sizeof (struct chromo)));
      if (chrom == NULL)
	perror ("malloc error. segtre2");
    }
  if (!
      (pseg2 = chrom[nchrom - 1].pseg =
       (struct seg *) calloc ((unsigned) newsg, sizeof (struct seg))))
    ERROR (" alloc error. re1");
  chrom[nchrom - 1].nseg = newsg;
  chrom[nchrom - 1].pop = chrom[ic].pop;
  pseg2->end = (pseg + jseg)->end;
  if (in)
    {
      pseg2->beg = is + 1;
      (pseg + jseg)->end = is;
    }
  else
    pseg2->beg = (pseg + jseg)->beg;
  pseg2->desc = (pseg + jseg)->desc;
  for (k = 1; k < newsg; k++)
    {
      (pseg2 + k)->beg = (pseg + jseg + k)->beg;
      (pseg2 + k)->end = (pseg + jseg + k)->end;
      (pseg2 + k)->desc = (pseg + jseg + k)->desc;
    }

  lsg = chrom[ic].nseg = lsg - newsg + in;
  lsgm1 = lsg - 1;
  nlinks -= pseg2->beg - (pseg + lsgm1)->end;
  len = (pseg + lsgm1)->end - (pseg->beg);
  cleft += 1.0 - pow (pc, len);
  len = (pseg2 + newsg - 1)->end - pseg2->beg;
  cleft += 1.0 - pow (pc, len);
  if (!(chrom[ic].pseg =
	(struct seg *) realloc (chrom[ic].pseg,
				(unsigned) (lsg * sizeof (struct seg)))))
    perror (" realloc error. re2");
  if (in)
    {
      begs = pseg2->beg;
      for (i = 0, k = 0;
	   (k < nsegs - 1) && (begs > seglst[seglst[i].next].beg - 1);
	   i = seglst[i].next, k++);
      if (begs != seglst[i].beg)
	{
	  /* new tree  */
	  if (nsegs >= seglimit)
	    {
	      seglimit += SEGINC;
	      nnodes =
		(int *) realloc (nnodes,
				 (unsigned) (sizeof (int) * seglimit));
	      if (nnodes == NULL)
		perror ("realloc error. 1. segtre_mig.c");
	      seglst =
		(struct segl *) realloc (seglst,
					 (unsigned) (sizeof (struct segl) *
						     seglimit));
	      if (seglst == NULL)
		perror ("realloc error. 2. segtre_mig.c");
	      /*  printf("seglimit: %d\n",seglimit);  */
	    }
	  seglst[nsegs].next = seglst[i].next;
	  seglst[i].next = nsegs;
	  seglst[nsegs].beg = begs;
	  if (!
	      (seglst[nsegs].ptree =
	       (struct node *) calloc ((unsigned) (2 * nsam),
				       sizeof (struct node))))
	    perror ("calloc error. re3.");
	  nnodes[nsegs] = nnodes[i];
	  ptree1 = seglst[i].ptree;
	  ptree2 = seglst[nsegs].ptree;
	  nsegs++;
	  for (k = 0; k <= nnodes[i]; k++)
	    {
	      (ptree2 + k)->abv = (ptree1 + k)->abv;
	      (ptree2 + k)->time = (ptree1 + k)->time;
	    }
	}
    }
  return (ic);
}

/***** common ancestor subroutine **********************
   Pick two chromosomes and merge them. Update trees if necessary. **/
int
ca (nsam, nsites, c1, c2)
     int nsam, c1, c2;
     int nsites;
{
  int yes1, yes2, seg1, seg2, seg;
  int tseg, start, end, desc, k;
  struct seg *pseg;
  struct node *ptree;

  seg1 = 0;
  seg2 = 0;
  if (!(pseg = (struct seg *) calloc ((unsigned) nsegs, sizeof (struct seg))))
    perror ("alloc error.ca1");

  tseg = -1;
  for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++)
    {
      start = seglst[seg].beg;
      yes1 = isseg (start, c1, &seg1);
      yes2 = isseg (start, c2, &seg2);
      if (yes1 || yes2)
	{
	  tseg++;
	  (pseg + tseg)->beg = seglst[seg].beg;
	  end =
	    (k < nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : nsites - 1);
	  (pseg + tseg)->end = end;

	  if (yes1 && yes2)
	    {
	      nnodes[seg]++;
	      if (nnodes[seg] >= (2 * nsam - 2))
		tseg--;
	      else
		(pseg + tseg)->desc = nnodes[seg];
	      ptree = seglst[seg].ptree;
	      desc = (chrom[c1].pseg + seg1)->desc;
	      (ptree + desc)->abv = nnodes[seg];
	      desc = (chrom[c2].pseg + seg2)->desc;
	      (ptree + desc)->abv = nnodes[seg];
	      (ptree + nnodes[seg])->time = t;

	    }
	  else
	    {
	      (pseg + tseg)->desc = (yes1 ?
				     (chrom[c1].pseg + seg1)->desc :
				     (chrom[c2].pseg + seg2)->desc);
	    }
	}
    }
  nlinks -= links (c1);
  cleft -= 1.0 - pow (pc, (double) links (c1));
  free (chrom[c1].pseg);
  if (tseg < 0)
    {
      free (pseg);
      chrom[c1].pseg = chrom[nchrom - 1].pseg;
      chrom[c1].nseg = chrom[nchrom - 1].nseg;
      chrom[c1].pop = chrom[nchrom - 1].pop;
      if (c2 == nchrom - 1)
	c2 = c1;
      nchrom--;
    }
  else
    {
      if (!
	  (pseg =
	   (struct seg *) realloc (pseg,
				   (unsigned) ((tseg +
						1) * sizeof (struct seg)))))
	perror (" realloc error. ca1");
      chrom[c1].pseg = pseg;
      chrom[c1].nseg = tseg + 1;
      nlinks += links (c1);
      cleft += 1.0 - pow (pc, (double) links (c1));
    }
  nlinks -= links (c2);
  cleft -= 1.0 - pow (pc, (double) links (c2));
  free (chrom[c2].pseg);
  chrom[c2].pseg = chrom[nchrom - 1].pseg;
  chrom[c2].nseg = chrom[nchrom - 1].nseg;
  chrom[c2].pop = chrom[nchrom - 1].pop;
  nchrom--;
  if (tseg < 0)
    return (2);			/* decrease of nchrom is two */
  else
    return (1);
}

/*** Isseg: Does chromosome c contain the segment on seglst which starts at
	    start? *psg is the segment of chrom[c] at which one is to begin 
	    looking.  **/
int
isseg (start, c, psg)
     int start, c, *psg;
{
  int ns;
  struct seg *pseg;
  ns = chrom[c].nseg;
  pseg = chrom[c].pseg;
  for (; ((pseg + (*psg))->beg <= start) && ((*psg) < ns); ++(*psg))
    if ((pseg + (*psg))->end >= start)
      return (1);
  return (0);
}



int
pick2_chrom (pop, config, pc1, pc2)
     int pop, *pc1, *pc2, config[];
{
  int c1, c2, cs, cb, i, count;
  pick2 (config[pop], &c1, &c2);
  cs = (c1 > c2) ? c2 : c1;
  cb = (c1 > c2) ? c1 : c2;
  i = count = 0;
  for (;;)
    {
      while (chrom[i].pop != pop)
	i++;
      if (count == cs)
	break;
      count++;
      i++;
    }
  *pc1 = i;
  i++;
  count++;
  for (;;)
    {
      while (chrom[i].pop != pop)
	i++;
      if (count == cb)
	break;
      count++;
      i++;
    }
  *pc2 = i;

  return (0);
}


/****  links(c): returns the number of links between beginning and end of chrom **/
int
links (c)
     int c;
{
  int ns;

  ns = chrom[c].nseg - 1;

  return ((chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);
}
