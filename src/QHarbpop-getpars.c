/*
 * QHarbpop-getpars.c
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

#include "msQHarbpop.h"

int getpars (int argc, char *argv[], struct parameters *);
void getQpars (int argc, char *argv[], int *parg, struct QHparameters *QH);
void Qparamscheck (int argc, struct QHparameters *QH);
double fsum2Dv (int sumindex, int dim1, int dim2, double **v);
int argnoerr (const char *option, struct parameters P);
int printpars (char *arg0, struct parameters P);
int usage (void);


int
getpars (int argc, char *argv[], struct parameters *pP)
{
  int arg = 2, i, ii, iii =
    -1, sum, usage (), argnoerr (const char *, struct parameters), printflag =
    0;
  /* int Mpattflag=0 ; */
  int cased = 0, casem = 0, configsum;
  double mig_rate=0.0, fsum2Dv (int, int, int, double **);
  int *pDintn;
  struct Dintervalparams **pDint;
  struct QHparameters **pQH;

  if (argc < 4)
    usage ();
  if ((argc < 2) || (argv[arg][0] == '-'))
    argnoerr ("nsam", *pP);
  else
    pP->nsam = atoi (argv[arg++]);
  if ((argc < 3) || (argv[arg][0] == '-'))
    argnoerr ("howmany", *pP);
  else
    pP->howmany = atol (argv[arg++]);
  pDintn = &(pP->Dintn);
  *pDintn = 0;
  pDint = &(pP->Dint);
  pP->theta = pP->r = pP->f = pP->track_len = 0.;
  pP->segsites = pP->D = pP->seQmut = 0;
  pP->nsites = 1000;
  pQH = &(pP->QH);
  *pQH = (struct QHparameters *) malloc (sizeof (struct QHparameters));
  if (*pQH == NULL)
    perror ("getpars QH malloc error\n");
  (*pQH)->alpha = 0.;
  while (arg < argc)
    {
      switch (argv[arg][1])
	{
	case 'r':
	  arg++;
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("r option r", *pP);
	  else
	    pP->r = atof (argv[arg++]);
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("r option nsites", *pP);
	  else
	    pP->nsites = atoi (argv[arg++]);
	  if (pP->nsites < 2)
	    {
	      fprintf (stderr,
		       "with -r option must specify both rec_rate and nsites>1\n");
	      printpars ("ERROR", *pP);
	      usage ();
	    }
	  break;
	case 'c':
	  arg++;
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("c option f", *pP);
	  else
	    pP->f = atof (argv[arg++]);
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("c option track_len", *pP);
	  else
	    pP->track_len = atof (argv[arg++]);
	  if (pP->track_len < 1.)
	    {
	      fprintf (stderr,
		       "with -c option must specify both f and track_len>0\n");
	      printpars ("ERROR", *pP);
	      usage ();
	    }
	  break;
	case 't':
	  arg++;
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("t option", *pP);
	  else
	    pP->theta = atof (argv[arg++]);
	  break;
	case 's':
	  arg++;
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("s option", *pP);
	  else
	    pP->segsites = atoi (argv[arg++]);
	  break;

	case 'Q':
	  pP->seQmut = 1;
	  (*pQH)->output = (strchr (argv[arg], '1')) ? 1 : 0;
	  getQpars (argc, argv, &arg, *pQH);
	  Qparamscheck (argc, *pQH);
	  break;
	case 'H':
	  (*pQH)->alpha = ((++arg == argc)
			   || (argv[arg][0] ==
			       '-')) ? 1. : atof (argv[arg++]);
	  if (((*pQH)->alpha <= 0) || ((*pQH)->alpha >= 1000))
	    {
	      fprintf (stderr, "gammaHet bad alpha %4.2f\n", (*pQH)->alpha);
	      usage ();
	    }
	  break;

	case 'm':
	  casem = 1;
	  arg++;
	  if (cased == 0)
	    {
	      *pDintn = 1;
	      *pDint =
		(struct Dintervalparams *)
		malloc (sizeof (struct Dintervalparams));
	    }
	  (*pDint)->Mpattern = 0;
	  if (cased == 0)
	    (*pDint)->Dtpast = 1.;
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("m option mig_rate", *pP);
	  else
	    mig_rate = atof (argv[arg++]);
	  if (!mig_rate > 0.)
	    {
	      fprintf (stderr, " -m error mig_rate\n");
	      printpars ("ERROR", *pP);
	      usage ();
	    }
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("m option npops", *pP);
	  else
	    (*pDint)->npops = atoi (argv[arg++]);
	  if ((*pDint)->npops <= 1)
	    {
	      fprintf (stderr, " -m error npops\n");
	      printpars ("ERROR", *pP);
	      usage ();
	    }
	  if (cased == 0)
	    {
	      (*pDint)->Nrecent =
		(double *) malloc ((unsigned) (*pDint)->npops *
				   sizeof (double));
	      (*pDint)->Npast =
		(double *) malloc ((unsigned) (*pDint)->npops *
				   sizeof (double));
	      for (ii = 0; ii < (*pDint)->npops; ii++)
		{
		  (*pDint)->Nrecent[ii] = 1.;
		  (*pDint)->Npast[ii] = 1.;
		}
	    }
	  (*pDint)->M =
	    (double **) malloc ((unsigned) (*pDint)->npops *
				sizeof (double *));
	  for (ii = 0; ii < (*pDint)->npops; ii++)
	    {
	      (*pDint)->M[ii] =
		(double *) malloc ((unsigned) (*pDint)->npops *
				   sizeof (double));
	      for (iii = 0; iii < (*pDint)->npops; iii++)
		(*pDint)->M[ii][iii] = (ii == iii) ? 0 : mig_rate;
	    }
	  for (i = 1; i < *pDintn; i++)
	    {			// *pDintn>1 iff cased
	      (*pDint + i)->npops = (*pDint)->npops;
	      (*pDint + i)->Mpattern = 0;
	      (*pDint + i)->Nrecent =
		(double *) realloc ((*pDint + i)->Nrecent,
				    (unsigned) (*pDint)->npops *
				    sizeof (double));
	      (*pDint + i)->Npast =
		(double *) realloc ((*pDint + i)->Npast,
				    (unsigned) (*pDint)->npops *
				    sizeof (double));
	      (*pDint + i)->M =
		(double **) malloc ((unsigned) (*pDint)->npops *
				    sizeof (double *));
	      for (ii = 1; ii < (*pDint)->npops; ii++)
		{
		  (*pDint + i)->Nrecent[ii] = (*pDint + i)->Nrecent[0];
		  (*pDint + i)->Npast[ii] = (*pDint + i)->Npast[0];
		}
	      (*pDint + i)->M = (*pDint)->M;
	    }
	  if (cased == 1)
	    free (pP->config);
	  pP->config =
	    (int *) malloc ((unsigned) ((*pDint)->npops) * sizeof (int));
	  for (ii = configsum = 0; ii < (*pDint)->npops; ii++)
	    {
	      if ((arg == argc) || (argv[arg][0] == '-'))
		argnoerr ("m option config[]", *pP);
	      else
		pP->config[ii] = atoi (argv[arg++]);
	      configsum += pP->config[ii];
	      if (configsum == pP->nsam)
		break;
	    }
	  for (++ii; ii < (*pDint)->npops; ii++)
	    pP->config[ii] = 0;
	  break;
	case 'd':
	  cased = 1;
	  arg++;
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("d option nintn", *pP);
	  else
	    *pDintn = atoi (argv[arg++]);
	  if (casem == 0)
	    {
	      *pDint =
		(struct Dintervalparams *) malloc ((unsigned) (*pDintn) *
						   sizeof (struct
							   Dintervalparams));
	      (*pDint)->npops = 1;
	      (*pDint)->Nrecent = (double *) malloc (sizeof (double));
	      (*pDint)->Npast = (double *) malloc (sizeof (double));
	    }
	  else			// (casem==1)
	    *pDint =
	      (struct Dintervalparams *) realloc (*pDint,
						  (unsigned) (*pDintn) *
						  sizeof (struct
							  Dintervalparams));
	  // assign Nrecent[], Npast[], Dtpast for interval zero
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("d option nrec", *pP);
	  else
	    *((*pDint)->Nrecent) = atof (argv[arg++]);
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("d option npast", *pP);
	  else
	    *((*pDint)->Npast) = atof (argv[arg++]);
	  for (ii = 1; ii < (*pDint)->npops; ii++)
	    {			// npops>1 iff casem
	      (*pDint)->Nrecent[ii] = (*pDint)->Nrecent[0];
	      (*pDint)->Npast[ii] = (*pDint)->Npast[0];
	    }
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("d option tpast", *pP);
	  else
	    (*pDint)->Dtpast = atof (argv[arg++]);
	  if ((*pDintn == 1) && ((*pDint)->Dtpast < DBL_EPSILON))
	    {
	      fprintf (stderr,
		       " -d option error -- only interval has zero length.\n");
	      printpars ("ERROR", *pP);
	      usage ();
	    }
	  if (casem == 0)
	    {
	      pP->config = (int *) malloc (sizeof (int));
	      *(pP->config) = pP->nsam;
	    }
	  // assign Nrecent[], Npast[], Dtpast for the intervals >zero
	  for (i = 1; i < *pDintn; i++)
	    {
	      (*pDint + i)->npops = (*pDint)->npops;
	      (*pDint + i)->Nrecent =
		(double *) malloc ((unsigned) (*pDint)->npops *
				   sizeof (double));
	      (*pDint + i)->Npast =
		(double *) malloc ((unsigned) (*pDint)->npops *
				   sizeof (double));
	      if ((arg == argc) || (argv[arg][0] == '-'))
		argnoerr ("d option nrec", *pP);
	      else
		*((*pDint + i)->Nrecent) = atof (argv[arg++]);
	      if ((arg == argc) || (argv[arg][0] == '-'))
		argnoerr ("d option npast", *pP);
	      else
		*((*pDint + i)->Npast) = atof (argv[arg++]);
	      for (ii = 1; ii < (*pDint)->npops; ii++)
		{		// npops>1 iff casem
		  (*pDint)->Nrecent[ii] = (*pDint)->Nrecent[0];
		  (*pDint)->Npast[ii] = (*pDint)->Npast[0];
		}
	      if ((arg == argc) || (argv[arg][0] == '-'))
		argnoerr ("d option tpast", *pP);
	      else
		(*pDint + i)->Dtpast = atof (argv[arg++]);
	      // move make-Dtpast cummulative somewhere else? (into ms or stree?)
	      (*pDint + i)->Dtpast += (*pDint + i - 1)->Dtpast;
	      // relax this ??
	      if ((i < *pDintn - 1)
		  && ((*pDint + i)->Dtpast - (*pDint + i - 1)->Dtpast <
		      DBL_EPSILON))
		{		// depends on cummulative
		  fprintf (stderr,
			   " -d option error -- interval other than first with zero length.\n");
		  printpars ("ERROR", *pP);
		  usage ();
		}
	      if (casem == 1)
		(*pDint + i)->M = (*pDint)->M;
	    }
	  if ((*((*pDint + (--i))->Npast) - *((*pDint + i)->Nrecent)) >
	      DBL_EPSILON)
	    {
	      fprintf (stderr,
		       " -d option error -- growth rate negative in oldest interval.\n");
	      printpars ("ERROR", *pP);
	      usage ();
	    }
	  break;

	case 'D':
	  arg++;
	  pP->D = 1;
	  if ((arg == argc) || (argv[arg][0] == '-'))
	    argnoerr ("D option Dintn", *pP);
	  else
	    *pDintn = atoi (argv[arg++]);
	  *pDint =
	    (struct Dintervalparams *) malloc ((unsigned) (*pDintn) *
					       sizeof (struct
						       Dintervalparams));
	  // Dintn LOOP
	  for (i = 0; i < *pDintn; i++)
	    {
	      // npops
	      if ((arg == argc) || (argv[arg][0] == '-'))
		argnoerr ("D option npops", *pP);
	      else
		(*pDint + i)->npops = atoi (argv[arg++]);
	      if ((*pDint + i)->npops == 0)
		{
		  fprintf (stderr, " -D error npops=0\n");
		  printpars ("ERROR", *pP);
		  usage ();
		}
	      // config OR a
	      if (i == 0)
		{
		  pP->config =
		    (int *) malloc ((unsigned) ((*pDint + i)->npops) *
				    sizeof (int));
		  if ((*pDint + i)->npops == 1)
		    *(pP->config) = pP->nsam;
		  else
		    {
		      for (ii = configsum = 0; ii < (*pDint + i)->npops; ii++)
			{
			  if ((arg == argc) || (argv[arg][0] == '-'))
			    argnoerr ("D option config[]", *pP);
			  else
			    pP->config[ii] = atoi (argv[arg++]);
			  configsum += pP->config[ii];
			  /* when config doesn't add up to nsam, I think it should stop here, Naoki, Nov 12 2008*/
			  /*
			  if (configsum == pP->nsam)
			    break;
			  */
			}
		      if (configsum != pP->nsam) {
			fprintf(stderr, "\nERROR: numbers of samples in subpops (config[]) don't add up to nsam in option -D \n\n");
			argnoerr("D option config[]", *pP);			
		      }
		      /* The following doesn't make sense, so commented out, Naoki Nov 12 2008
		      for (++ii; ii < (*pDint)->npops; ii++)
			pP->config[ii] = 0;
		      */
		    }
		}
	      else if ((*pDint + i)->npops > 1)
		{
		  (*pDint + i)->a =
		    (double **) malloc ((unsigned) (*pDint + i - 1)->npops *
					sizeof (double *));
		  for (ii = 0; ii < (*pDint + i - 1)->npops; ii++)
		    (*pDint + i)->a[ii] =
		      (double *) malloc ((unsigned) ((*pDint + i)->npops) *
					 sizeof (double));
		  for (ii = 0; ii < (*pDint + i - 1)->npops; ii++)
		    {
		      for (iii = 0; iii < (*pDint + i)->npops; iii++)
			if ((arg == argc) || (argv[arg][0] == '-'))
			  argnoerr ("D option a[]", *pP);
			else
			  (*pDint + i)->a[ii][iii] = atof (argv[arg++]);
		    }
		  for (ii = 0; ii < (*pDint + i - 1)->npops; ii++)
		    if (fabs
			(1. -
			 fsum2Dv (2, ii, (*pDint + i)->npops,
				  (*pDint + i)->a)) > DBL_EPSILON)
		      {
			fprintf (stderr,
				 " -D error a[prevDintpop][.] should sum to 1\n");
			printpars ("ERROR", *pP);
			usage ();
		      }
		}
	      // Mpattern AND M
	      if ((*pDint + i)->npops > 1)
		{
		  if ((arg == argc) || (argv[arg][0] == '-'))
		    argnoerr ("D option Mpattern", *pP);
		  else
		    (*pDint + i)->Mpattern = atoi (argv[arg++]);
		  // Mpatt just in 1st interval?
/*        if (!Mpattflag) */
/*          if ((arg==argc)||(argv[arg][0]=='-')) argnoerr("D option Mpattern",*pP); */
/*          else { */
/*            (*pDint+i)->Mpattern = (*pDint)->Mpattern = atoi( argv[arg++] ); */
/*            if (((*pDint+i)->Mpattern<0)||((*pDint+i)->Mpattern>2)) { */
/*              fprintf(stderr," -D option error -- Mpattern not equal 0, 1 or 2.\n"); */
/*              printpars("ERROR",*pP); usage(); } */
/*            Mpattflag=1; */
/*          } */
/*        else (*pDint+i)->Mpattern = (*pDint)->Mpattern ; */
		  // M
		  (*pDint + i)->M =
		    (double **) malloc ((unsigned) (*pDint + i)->npops *
					sizeof (double *));
		  for (ii = 0; ii < (*pDint + i)->npops; ii++)
		    {
		      (*pDint + i)->M[ii] =
			(double *) malloc ((unsigned) (*pDint + i)->npops *
					   sizeof (double));
		      for (iii = 0; iii < (*pDint + i)->npops; iii++)
			(*pDint + i)->M[ii][iii] = 0.;
		    }
		  if (argv[arg][0] == 'I')
		    {
		      arg++;
		      //read M*; M[0][1]=M*/npops-1; assign all off-diag =M[0][1]
		      if ((arg == argc) || (argv[arg][0] == '-'))
			argnoerr ("D option Island M*", *pP);
		      else
			(*pDint + i)->M[0][1] =
			  atof (argv[arg++]) / ((*pDint + i)->npops - 1);
		      for (ii = 0; ii < (*pDint + i)->npops; ii++)
			for (iii = 0; iii < (*pDint + i)->npops; iii++)
			  if (((ii == 0) && (iii > 1)) || (ii > 0))
			    if (ii != iii)
			      (*pDint + i)->M[ii][iii] =
				(*pDint + i)->M[0][1];
		    }
		  else if (argv[arg][0] == 'S')
		    {
		      arg++;
		      if ((argv[arg - 1][1] == '1')
			  || (argv[arg - 1][1] == '2'))
			{
			  if ((arg == argc) || (argv[arg][0] == '-'))
			    argnoerr ("D option StepStone M1", *pP);
			  else
			    (*pDint + i)->M[0][1] = atof (argv[arg++]) / 2.;
			  // loop over all 1-off-diag and assign =M[0][1]
			  for (ii = 0; ii < (*pDint + i)->npops; ii++)
			    for (iii = 0; iii < (*pDint + i)->npops; iii++)
			      if ((ii > 0) && (abs (ii - iii) == 1))
				{
				  (*pDint + i)->M[ii][ii - 1] =
				    (*pDint + i)->M[0][1];
				  if (ii < (*pDint + i)->npops - 1)
				    (*pDint + i)->M[ii][ii + 1] =
				      (*pDint + i)->M[0][1];
				}
			  if (argv[arg - 2][1] == '2')
			    {
			      // read Minf; M[0][2]=Minf/npops; all off-diag +=M[0][2]
			      if ((arg == argc) || (argv[arg][0] == '-'))
				argnoerr ("D option StepStone Minf", *pP);
			      else
				(*pDint + i)->M[0][2] =
				  atof (argv[arg++]) / (*pDint + i)->npops;
			      for (ii = 0; ii < (*pDint + i)->npops; ii++)
				for (iii = 0; iii < (*pDint + i)->npops;
				     iii++)
				  if ((ii != iii)
				      && (!((ii == 0) && (iii == 2))))
				    (*pDint + i)->M[ii][iii] +=
				      (*pDint + i)->M[0][2];
			    }
			}
		      else
			{
			  // read in all 1-off-diag M[][]
			  for (ii = 0; ii < (*pDint + i)->npops; ii++)
			    for (iii = 0; iii < (*pDint + i)->npops; iii++)
			      if ((arg == argc) || (argv[arg][0] == '-'))
				argnoerr ("D option StepStone M[][]", *pP);
			      else if (abs (ii - iii) == 1)
				(*pDint + i)->M[ii][iii] = atof (argv[arg++]);
			}
		    }
		  else
		    // Completely arbitrary migration matrix
		    for (ii = 0; ii < (*pDint + i)->npops; ii++)
		      for (iii = 0; iii < (*pDint + i)->npops; iii++)
			if ((arg == argc) || (argv[arg][0] == '-'))
			  argnoerr ("D option M[]", *pP);
			else if (ii != iii)
			  (*pDint + i)->M[ii][iii] = atof (argv[arg++]);
		}
	      // N
	      (*pDint + i)->Nrecent =
		(double *) malloc ((unsigned) (*pDint + i)->npops *
				   sizeof (double));
	      (*pDint + i)->Npast =
		(double *) malloc ((unsigned) (*pDint + i)->npops *
				   sizeof (double));
	      if (argv[arg][0] == 'N')
		{
		  arg++;
		  if (argv[arg - 1][1] == '1')
		    for (ii = 0; ii < (*pDint + i)->npops; ii++)
		      (*pDint + i)->Nrecent[ii] = (*pDint + i)->Npast[ii] =
			1.;
		  else if (argv[arg - 1][1] == 'c')
		    for (ii = 0; ii < (*pDint + i)->npops; ii++)
		      if ((arg == argc) || (argv[arg][0] == '-'))
			argnoerr ("D option Nc N", *pP);
		      else
			(*pDint + i)->Nrecent[ii] = (*pDint + i)->Npast[ii] =
			  atof (argv[arg++]);
		}
	      else
		for (ii = 0; ii < (*pDint + i)->npops; ii++)
		  {
		    if ((arg == argc) || (argv[arg][0] == '-'))
		      argnoerr ("D option Nrecent", *pP);
		    else
		      (*pDint + i)->Nrecent[ii] = atof (argv[arg++]);
		    if ((arg == argc) || (argv[arg][0] == '-'))
		      argnoerr ("D option Npast", *pP);
		    else
		      (*pDint + i)->Npast[ii] = atof (argv[arg++]);
		  }
	      // Dtpast
	      if ((arg == argc) || (argv[arg][0] == '-'))
		argnoerr ("D option tpast", *pP);
	      else
		(*pDint + i)->Dtpast = atof (argv[arg++]);
	      if (i == 0)
		{
		  if ((*pDintn == 1) && ((*pDint)->Dtpast < DBL_EPSILON))
		    {
		      fprintf (stderr,
			       " -D option error -- only interval has zero length.\n");
		      printpars ("ERROR", *pP);
		      usage ();
		    }
		}
	      else
		{
		  (*pDint + i)->Dtpast += (*pDint + i - 1)->Dtpast;
		  if ((i < *pDintn - 1)
		      && ((*pDint + i)->Dtpast - (*pDint + i - 1)->Dtpast <
			  DBL_EPSILON))
		    {
		      fprintf (stderr,
			       " -D option error -- interval other than first with zero length.\n");
		      printpars ("ERROR", *pP);
		      usage ();
		    }
		}
	    }
	  // END OF Dintn LOOP
	  for (i--, ii = 0; ii < (*pDint + i)->npops; ii++)
	    {
	      if (((*pDint + i)->Npast[ii] - (*pDint + i)->Nrecent[ii]) >
		  DBL_EPSILON)
		{
		  fprintf (stderr,
			   " -D option error-- pop%1d growth rate negative in oldest interval.\n",
			   ii);
		  printpars ("ERROR", *pP);
		  usage ();
		}
	      // CHECK TO SEE THAT ((M[i][j]!=0)||((a[.][i]==0)||(a[.][j]==0)))
	      if ((*pDint + i)->npops == 2)
/*        for (iii=0;iii<(*pDint+i)->npops;iii++) */
/*          if ((ii!=iii)&&((*pDint+i)->M[ii][iii] < DBL_EPSILON)) */
		if (fsum2Dv (2, ii, (*pDint + i)->npops, (*pDint + i)->M) <
		    DBL_EPSILON)
		  {
		    if (i == 0)
		      {
			if ((pP->config[ii] > 0) && (pP->config[iii] > 0))
			  {
			    fprintf (stderr,
				     " -D error-- isolated pops %1d,%1d with non-zero config[i,j] in oldest interval.\n",
				     ii, iii);
			    printpars ("ERROR", *pP);
			    usage ();
			  }
		      }
		    else
		      {
			fprintf (stderr, "reached. %9.7f \n",
				 fsum2Dv (1, (*pDint + i - 1)->npops, ii,
					  (*pDint + i)->a));
			if ((fsum2Dv
			     (1, (*pDint + i - 1)->npops, ii,
			      (*pDint + i)->a) > DBL_EPSILON)
			    &&
			    (fabs
			     ((double) ((*pDint + i - 1)->npops) -
			      fsum2Dv (1, (*pDint + i - 1)->npops, ii,
				       (*pDint + i)->a)) > DBL_EPSILON))
			  {
			    fprintf (stderr,
				     " -D error-- isolated pop ii=%1d with non-zero a[.][ii] in oldest interval.\n",
				     ii);
			    printpars ("ERROR", *pP);
			    usage ();
			  }
		      }
		  }
	    }
	  break;

	case 'P':
	  arg++;
	  printflag = 1;
	  break;
	default:
	  fprintf (stderr, " option error-- extra arguments: ");
	  for (; arg < argc;)
	    fprintf (stderr, "%s ", argv[arg++]);
	  fprintf (stderr, "\n");
	  printpars ("ERROR", *pP);
	  usage ();
	}
    }

  if (*pDintn == 0)
    {
      *pDintn = 1;
      *pDint =
	(struct Dintervalparams *) malloc (sizeof (struct Dintervalparams));
      (*pDint)->Dtpast = 1.;
      (*pDint)->npops = 1;
      (*pDint)->Nrecent = (double *) malloc (sizeof (double));
      (*pDint)->Npast = (double *) malloc (sizeof (double));
      *((*pDint)->Nrecent) = 1.;
      *((*pDint)->Npast) = 1.;
      pP->config = (int *) malloc (sizeof (int));
      *(pP->config) = pP->nsam;
    }

  if (((fabs (pP->theta) < DBL_EPSILON) && (pP->segsites == 0))
      || ((fabs (pP->theta) > DBL_EPSILON) && (pP->segsites > 0)))
    {
      fprintf (stderr, " either -s or -t option must be used. \n");
      printpars ("ERROR", *pP);
      usage ();
    }
  for (i = sum = 0; i < (*pDint)->npops; i++)
    sum += pP->config[i];
  if (sum != pP->nsam)
    {
      fprintf (stderr, " sum sample sizes != nsam\n");
      printpars ("ERROR", *pP);
      usage ();
    }
  if ((pP->theta < DBL_EPSILON) && (pP->seQmut))
    {
      fprintf (stderr, " -Q and -H options must be used with -t option. \n");
      usage ();
      exit (1);
    }
  if (!(pP->seQmut) && ((*pQH)->alpha))
    {
      fprintf (stderr,
	       " -H option must be used with -Q option. alpha %8.6f\n",
	       (*pQH)->alpha);
      usage ();
      exit (1);
    }

  return printflag;
}


void
getQpars (int argc, char *argv[], int *parg, struct QHparameters *QH)
{
  int nQarg = 0;
  while (((*parg + nQarg + 1) < argc) && (argv[*parg + nQarg + 1][0] != '-'))
    nQarg++;
  (*parg)++;
  if ((nQarg == 1) || (nQarg == 2) || (nQarg == 5) || (nQarg == 6))
    {
      QH->Rtstv = atof (argv[(*parg)++]);
      QH->Ytstv = ((nQarg == 2)
		   || (nQarg == 6)) ? atof (argv[(*parg)++]) : QH->Rtstv;
    }
  else
    QH->Rtstv = QH->Ytstv = 1.;
  if ((nQarg == 4) || (nQarg == 5) || (nQarg == 6))
    {
      QH->freqA = atof (argv[(*parg)++]);
      QH->freqC = atof (argv[(*parg)++]);
      QH->freqG = atof (argv[(*parg)++]);
      QH->freqT = atof (argv[(*parg)++]);
    }
  else
    QH->freqA = QH->freqC = QH->freqG = QH->freqT = 0.25;
  if ((nQarg != 0) && (nQarg != 1) && (nQarg != 2) && (nQarg != 4)
      && (nQarg != 5) && (nQarg != 6))
    {
      fprintf (stderr, "wrong number %1d of args with -Q\n", nQarg);
      usage ();
    }
}

void
Qparamscheck (int argc, struct QHparameters *QH)
{
  if ((QH->Rtstv == 0) || (QH->Ytstv == 0))
    {
      fprintf (stderr, "\n bad %d arguments\n tstv equals zero\n", argc);
      usage ();
    }
  else if (fabs ((QH->freqA + QH->freqC + QH->freqG + QH->freqT) - 1) > 1e-6)
    {
      fprintf (stderr, "\n bad %d arguments\n freqs don't add to one\n",
	       argc);
      usage ();
    }
}


double
fsum2Dv (int sumindex, int dim1, int dim2, double **v)
{
  int i;
  double fsum = 0.;
  if (sumindex == 1)
    for (i = 0; i < dim1; i++)
      fsum += v[i][dim2];
  else
    for (i = 0; i < dim2; i++)
      fsum += v[dim1][i];
  return (fsum);
}


int
argnoerr (const char *option, struct parameters P)
{
  fprintf (stderr, "-%s error: wrong no. args\n", option);
  printpars ("ERROR", P);
  usage ();
  return 0;
}


int
printpars (char *arg0, struct parameters P)
{
  int i, ii, iii;
  printf ("%s nsam %1d howmany %1ld\n", arg0, P.nsam, P.howmany);
  printf ("  theta %1.2f segsites %1d\n", P.theta, P.segsites);
/* Aug 11, 2004 Eli CHANGES */
  if (P.seQmut)
    {
      printf ("seQmut %1d output %1d\n", P.seQmut, P.QH->output);
      printf ("tstvAG CT %4.2f %4.2f, freqACGT %4.2f %4.2f %4.2f %4.2f\n",
	      P.QH->Rtstv, P.QH->Ytstv, P.QH->freqA, P.QH->freqC, P.QH->freqG,
	      P.QH->freqT);
      printf ("gammaHet alpha %4.2f\n", P.QH->alpha);
/*     printf("QH->siterate %10.8f\n",P.QH->siterate); */
    }
/* end Eli */
  printf ("  r %1.2f f %1.2f tr_len %1.2f nsites %1d\n", P.r, P.f,
	  P.track_len, P.nsites);
  printf ("  Dintn %1d ", P.Dintn);
  for (i = 0; i < P.Dintn; i++)
    {
      printf ("\n    ");
      printf ("Dint%1d npops %1d", i, (P.Dint + i)->npops);
      if (i == 0)
	{
	  printf ("\n      config[] ");
	  for (ii = 0; ii < (P.Dint + i)->npops; ii++)
	    printf ("%1d ", P.config[ii]);
	}
      else if ((P.Dint + i)->npops > 1)
	{
	  printf ("\n      a[][] ");
	  for (ii = 0; ii < (P.Dint + i - 1)->npops; ii++)
	    for (iii = 0; iii < (P.Dint + i)->npops; iii++)
	      printf ("%1.2f ", (P.Dint + i)->a[ii][iii]);
	}
      if ((P.Dint + i)->npops > 1)
	{
	  printf ("\n      Mpattern %1d", (P.Dint + i)->Mpattern);
	  printf ("\n      M[][] ");
	  for (ii = 0; ii < (P.Dint + i)->npops; ii++)
	    for (iii = 0; iii < (P.Dint + i)->npops; iii++)
	      printf ("%1.2f ", (P.Dint + i)->M[ii][iii]);
	}
      printf ("\n      (Nrec_Npast)[] ");
      for (ii = 0; ii < (P.Dint + i)->npops; ii++)
	printf ("%1.2f %1.2f ", (P.Dint + i)->Nrecent[ii],
		(P.Dint + i)->Npast[ii]);
      printf ("\n       tpast %1.2f", (P.Dint + i)->Dtpast);
    }
  printf ("\n");
  return (0);
}


int
usage ()
{
  fprintf (stderr, "usage: msQHarbpop seed nsam howmany \n");
  fprintf (stderr, "  Options: \n");
  fprintf (stderr, "\t -t theta   (this option or the next must be used.)\n");
  fprintf (stderr, "\t -s segsites   ( fixed number of segregating sites)\n");
  fprintf (stderr,
	   "\t -Q[1] [Rtstvratio [Ytstvratio] [Afreq Cfreq Gfreq Tfreq]] (0,1,2,4,5 or 6 args)\n");
  fprintf (stderr,
	   "\t    (seQuence mutation follows JC(0), K2P(1), HKY(5) or TN(6) model (no. args) )\n");
  fprintf (stderr,
	   "\t    (-Q outputs variable sites, -Q1 outputs whole sequences)\n");
  fprintf (stderr, "\t    (use only with -t)\n");
  fprintf (stderr,
	   "\t -H gamma_alpha (Heterogeneity of mut across sites ~ gamma(alpha,mean 1) times theta)\n");
  fprintf (stderr, "\t    (use only with -Q and -t)\n");
  fprintf (stderr, "\t -r rho nsites     (rho here is 4Nc)\n");
  fprintf (stderr,
	   "\t -c f track_len   (ratio of conversion rate to r, or conv rate if r==0) \n");
  fprintf (stderr, "\t -m mig_rate npop n1 n2 ...   \n");
  fprintf (stderr, "\t -d nintn nr1  np1 t1 [nr2 np2 t2 ...]   \n");
  fprintf (stderr,
	   "\t -D Dintn npops1 [ config[] Mpattern ( 'I' M | 'S1' M1 | 'S2' M1 Minf | 'S' M1[1-offdiag] | M1[offdiag] )\n");
  fprintf (stderr,
	   "\t                   ( 'N1' | 'Nc' N[] | Nr11 Np11 [ Nr12 Np12 ... ] ) t1\n");
  fprintf (stderr,
	   "\t        [ npops2 [ a2[][] Mpattern ( M ) ( N ) ] t2 ]\n");
  fprintf (stderr, "\t        [ ... ]\n");
  fprintf (stderr,
	   "\t -P (print the parameters in human-readable format and exit(0))\n");
  fprintf (stderr,
	   " See ms.c comments at top for explanation of these parameters.\n");
  exit (1);
}
