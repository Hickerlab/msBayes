/* rejectinC.c : 

ARGUMENTS: OBS_STATS_FILE SIM_STATS_FILE TOLERANCE COLUMNS..

rejectinC, or rejeCt, does the acceptance/rejection step of Hickerson
et al's ABC algorithm in C (rather than R, which uses A LOT of
memory. This program expects a file of observed summary statistics and
a file of simulated summary statistics, in the SAME FORMAT (same
column contains same summary statistic, e.g. Tajima's D, in both
files); input filenames are the first two commandline arguments.  The
third argument is the tolerance, or the proportion of simulation
replicates that will be accepted.  The fourth and following arguments
are the columns in the summary statistic files to be considered.

rejeCt opens the observed stats file and reads the appropriate columns
into memory. Then it opens the simulated stats file and reads through
it calculating a mean and sd for each of the appropriate columns.  It
then normalizes the observed stats (in memory). It then reopens the
sim stats file and reads in the appropriate stats, and normalizes them
by mean and sd and calculates the Euclidian distance between that
simulation replicate and the observed data. This program keeps the
proportion tolerance (the number equal to tolerance times the number
of simulation replicates, i.e. lines in the sim stats file) of the
simulation replicates with the lowest Euclidean distance from the
observed. For these simulation replicates, the entire lines from the
sim stats file are written to standard out.

If simulated statistics are equal to 'nan' or 'inf', a warning message
is written to standard error (should show up on the terminal). These
statistics are not considered when calculating means and sd's, and the
component of Euclidean distance for that statistic equal to 3 (i.e. 3
sd's from the mean, or "kinda far") is arbitrarily substituted.

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

int
main (int argc, char *argv[])
{
  int i, c, ncolumns, *columns, enoughMem;
  FILE *pfin, *fopen (), *freopen ();
  char *inputstring, *definstr, *cptr = NULL, dum[30], **bestline, *strcpy ();
  double tolerance, oldTol, *obs, *vec, *x, *s, *normobs, normalizedvec, euclid,
    *besteuclid, **bestvec;
  long linecount, *count, li, lii, tolerated;
  double sqrt ();
  void usage (int);

  if ((argc > 1) && ((argv[1][0] == '-') && (argv[1][2] == 'h')))
    usage (0);
  tolerance = atof (argv[3]);
  ncolumns = argc - 4; // only the first three parameters are not sumStats to be used(??)
  columns = (int *) malloc ((unsigned) ncolumns * sizeof (int));
  if (columns == NULL)
    {
      fprintf (stderr, "No memory in msReject\n");
      exit (EXIT_FAILURE);
    }

  for (i = 4; i < argc; i++)
    {
      columns[i - 4] = atoi (argv[i]) - 1;
    }

  inputstring = (char *) malloc (100000 * sizeof (char));
  definstr = inputstring;
  obs = (double *) malloc ((unsigned) ncolumns * sizeof (double));
  vec = (double *) malloc ((unsigned) ncolumns * sizeof (double));
  count = (long *) malloc ((unsigned) ncolumns * sizeof (long));
  x = (double *) malloc ((unsigned) ncolumns * sizeof (double));
  s = (double *) malloc ((unsigned) ncolumns * sizeof (double));
  normobs = (double *) malloc ((unsigned) ncolumns * sizeof (double));
  if (inputstring == NULL || obs == NULL || vec == NULL || count == NULL ||
      x == NULL || s == NULL)
    {
      fprintf (stderr, "No memory in msReject\n");
      exit (EXIT_FAILURE);
    }
  for (i = 0; i < ncolumns; i++)
    {
      x[i] = s[i] = 0.;
      count[i] = 0;
    }
  linecount = 0;

  if ((pfin = fopen (argv[1], "r")) == NULL)
    {
      fprintf (stderr, "Cannot open the observed data stats file %s.\n",
	       argv[1]);
      exit (1);
    }
  fgets (inputstring, 99999, pfin);
  for (c = i = 0; c < ncolumns; c++)
    {
      while (i <= columns[c])
	{
	  sscanf (inputstring, "%s", dum);
	  cptr = inputstring;
	  while ((*cptr != '\0') && (*cptr != ' ') && (*cptr != '\t'))
	    cptr++;
	  while ((*cptr != '\0') && ((*cptr == ' ') || (*cptr == '\t')))
	    cptr++;
	  inputstring = cptr;
	  ++i;
	}
      obs[c] = atof (dum);
      if ((*cptr == '\0') && (c < ncolumns - 1))
	{
	  fprintf (stderr, "Not enough columns in observed stats file %s.\n",
		   argv[1]);
	  exit (1);
	}
    }
  inputstring = definstr;

  if (freopen (argv[2], "r", pfin) == NULL)
    {
      fprintf (stderr, "Cannot open the simulated data stats file %s.\n",
	       argv[2]);
      exit (1);
    }
  while (fgets (inputstring, 99999, pfin) != NULL)
    {
      for (c = i = 0; c < ncolumns; c++)
	{
	  while (i <= columns[c])
	    {
	      sscanf (inputstring, "%s", dum);
	      cptr = inputstring;
	      while ((*cptr != '\0') && (*cptr != ' ') && (*cptr != '\t'))
		cptr++;
	      while ((*cptr != '\0') && ((*cptr == ' ') || (*cptr == '\t')))
		cptr++;
	      inputstring = cptr;
	      ++i;
	    }
	  vec[c] = atof (dum);
	  if (finite (vec[c]))
	    {
	      x[c] += vec[c];
	      s[c] += vec[c] * vec[c];
	      ++count[c];
	    }
	  if ((*cptr == '\0') && (c < ncolumns - 1))
	    {
	      fprintf (stderr,
		       "Not enough columns in simulated stats file %s, line %ld.\n",
		       argv[2], linecount);
	      exit (1);
	    }
	}
      inputstring = definstr;
      ++linecount;
    }
  for (c = 0; c < ncolumns; c++)
    {
      x[c] /= (double) count[c];
      s[c] /= (double) count[c];
      s[c] -= x[c] * x[c];
      s[c] = sqrt (s[c] * (double) count[c] / ((double) count[c] - 1.0));
      normobs[c] = obs[c] - x[c];
      normobs[c] /= s[c];
    }


  /* Not enough memory, so finding a new tolerance value which will work */
  oldTol = tolerance;
  do {
    enoughMem = 1;

    tolerated = (long) (tolerance * (double) linecount + DBL_EPSILON);
    besteuclid = (double *) malloc ((unsigned) tolerated * sizeof (double));
    bestvec = (double **) malloc ((unsigned) tolerated * sizeof (double *));
    bestline = (char **) malloc ((unsigned) tolerated * sizeof (char *));

    if (besteuclid == NULL || bestvec == NULL || bestline == NULL) {
      free(besteuclid); free(bestvec); free(bestline);
      enoughMem = 0;
      tolerance /= 2;
      continue;
    }
	
    for (li = 0; li < tolerated; ++li)
      {
	bestvec[li] = (double *) malloc ((unsigned) ncolumns * sizeof (double));
	bestline[li] = (char *) malloc ((unsigned) 100000 * sizeof (char));
	if (bestvec[li] == NULL || bestline[li] == NULL) {
	  /* clean up the allocated memory */
	  long freeLine;
	  for(freeLine = 0; freeLine <= li; freeLine++) {
	    free(bestvec[freeLine]);
	    free(bestline[freeLine]);
	  }
	  free(besteuclid); free(bestvec); free(bestline);

	  enoughMem = 0;
	  tolerance /= 2;
	  break;
	}
      }
  } while (enoughMem == 0);

  if (oldTol != tolerance) {
    
    double betterTol = 500.0 / (double) linecount;
    fprintf(stderr, "ERROR: No memory in msReject.\n");
    fprintf(stderr, 
	    "ERROR: PLEASE USE A LOWER TOLERANCE. Instead of %f you "
	    "specified,\n", oldTol);
    fprintf (stderr, 
	     "ERROR: tolerance values between %g to %f should work.\n", 
	     betterTol, tolerance);
    fprintf (stderr,     
	     "ERROR: The lowest tolerance (%g) results in acceptance of "
	     "500 simulations.\n",
	     betterTol);
    fprintf(stderr,
	    "ERROR: A low tolerance value, accepting 500-1000 simulations, is\n"
	    "ERROR: recommended to improve the estimation.\n");
    exit (EXIT_FAILURE);
    /* it could continue with the new tolerance value, but I'm
       quiting to be safer */
  }


  linecount = 0;

  if (freopen (argv[2], "r", pfin) == NULL)
    {
      fprintf (stderr, "Cannot open the simulated data stats file %s (2).\n",
	       argv[2]);
      exit (1);
    }
  while (fgets (inputstring, 99999, pfin) != NULL)
    {
      euclid = 0.;
      for (c = i = 0; c < ncolumns; c++)
	{
	  while (i <= columns[c])
	    {
	      sscanf (inputstring, "%s", dum);
	      cptr = inputstring;
	      while ((*cptr != '\0') && (*cptr != ' ') && (*cptr != '\t'))
		cptr++;
	      while ((*cptr != '\0') && ((*cptr == ' ') || (*cptr == '\t')))
		cptr++;
	      inputstring = cptr;
	      ++i;
	    }
	  vec[c] = atof (dum);
	  if (finite (vec[c]))
	    {
	      normalizedvec = (vec[c] - x[c]) / s[c];
	    }
	  else
	    {
	      fprintf (stderr,
		       "WARNING: Nan or inf in simulated stats file %s, line %ld, column %d. Arbitrarily, the Euclidean distance of this simulation replicate to the observed has been incremented by 3 for this statistic.\n",
		       argv[2], linecount, columns[c]);
	      normalizedvec = 3.;
	    }
	  euclid +=
	    (normalizedvec - normobs[c]) * (normalizedvec - normobs[c]);
	}
      inputstring = definstr;
      euclid = sqrt (euclid);

      if (linecount == 0)
	{
	  besteuclid[0] = euclid;
	  for (c = 0; c < ncolumns; c++)
	    bestvec[0][c] = vec[c];
	  strcpy (bestline[0], inputstring);
	}
      else
	{
	  for (li = 0; li < (linecount < tolerated ? linecount : tolerated);
	       li++)
	    {
	      if (euclid < besteuclid[li])
		break;
	    }
	  
	  for (lii = (linecount < tolerated ? linecount : (tolerated - 1));
	       lii > li; lii--)
	    {
	      besteuclid[lii] = besteuclid[lii - 1];
	      for (c = 0; c < ncolumns; c++)
		bestvec[lii][c] = bestvec[lii - 1][c];
	      strcpy (bestline[lii], bestline[lii - 1]);
	    }

	  if (li < tolerated)
	    {
	      besteuclid[li] = euclid;
	      for (c = 0; c < ncolumns; c++)
		{
		  bestvec[li][c] = vec[c];
		}
	      strcpy (bestline[li], inputstring);
	    }
	}
      ++linecount;
    }

  for (li = 0; li < tolerated; li++)
    {
      printf ("%s", bestline[li]);
    }

  exit (0);
}



void
usage (int exitstatus)
{
  printf
    ("\nUsage: rejeCt [--help] Obsered_Stats_file Simulated_Stats_file Tolerance Columns[]\n");
  printf
    ("\nrejectinC, or rejeCt, does the acceptance/rejection step of Hickerson et al's ABC algorithm in C\n");
  printf
    ("  (rather than R, which uses A LOT of memory. This program expects a file of observed summary statistics\n");
  printf
    ("  and a file of simulated summary statistics, in the SAME FORMAT (same column contains same summary\n");
  printf
    ("  statistic, e.g. Tajima's D, in both files); input filenames are the first two commandline arguments. \n");
  printf
    ("  The third argument is the tolerance, or the proportion of simulation replicates that will be accepted. \n");
  printf
    ("  The fourth and following arguments are the columns in the summary statistic files to be considered. \n");
  printf
    ("\nrejeCt opens the observed stats file and reads the appropriate columns into memory. Then it opens the \n");
  printf
    ("  simulated stats file and reads through it calculating a mean and sd for each of the appropriate columns.\n");
  printf
    ("  It then normalizes the observed stats (in memory). It then reopens the sim stats file and reads in the\n");
  printf
    ("  appropriate stats, and normalizes them by mean and sd and calculates the Euclidian distance between that\n");
  printf
    ("  simulation replicate and the observed data. This program keeps the proportion tolerance (the number equal\n");
  printf
    ("  to tolerance times the number of simulation replicates, i.e. lines in the sim stats file) of the \n");
  printf
    ("  simulation replicates with the lowest Euclidean distance from the observed. For these simulation\n");
  printf
    ("  replicates, the entire lines from the sim stats file are written to standard out. \n");
  printf
    ("\nIf simulated statistics are equal to 'nan' or 'inf', a warning message is written to standard error (should\n");
  printf
    ("  show up on the terminal). These statistics are not considered when calculating means and sd's, and the\n");
  printf
    ("  component of Euclidean distance for that statistic equal to 2 (i.e. 2 sd's from the mean, or \"kinda far\")\n");
  printf ("  is arbitrarily substituted.\n\n");
  exit (exitstatus);
}
