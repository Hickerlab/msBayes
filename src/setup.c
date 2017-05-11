/*
 * setup.c
 *
 * Copyright (C) 2006   Naoki Takebayashi and Michael Hickerson
 *
 * This file is a part of msprior, distributed with msBayes.
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

/* for strcasestr */
#define _GNU_SOURCE
#include <string.h>

#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "setup.h"
#include "msprior.h"
#include "whiteSpaces.h"
#include "initvars.h"
#include "stringUtils.h"	/* for cmatrix, uniqueStrings */

#define LNSZ 256		/* max length of a line */

int setup_done = 0;
int debug_level = 0;

static void ParseCommandLine (int argc, char *agv[]);
static int InteractiveSetupParams (runParameters * paramPtr);
static void SetupParams (FILE * fp, runParameters * paramPtr);
static int SetConfigFile (char *fName);
void PrintParam (FILE *fp);
static int InitMutPara (FILE * fp, mutParameterArray * mutParaArray);
static int ProcessTaxonLociInMutParaArray (mutParameterArray * mpaPtr);
static int CheckMutParaArray (mutParameterArray * mpaPtr, int index);
static int ReadMutLine (mutParameter * mpp, char *line, int ncol);
static int InitConPara (FILE * fp, constrainedParameterArray * conParaArray);
static int CheckConParaArray (constrainedParameterArray * cpaPtr, int index);
static int ReadConLine (constrainedParameter * cpa, char *line, int ncol);

/*
 * The functions in this file is to set global parameters gParam and gMutParam.
 * To test the functionality of the functions in this file,
 *
 * gcc -g -DTEST_SETUP setup.c initvars.c whiteSpaces.c
 *
 * There are several ways to set the parameter values: command line options,
 * config files, and interactive way.
 *
 * If config file is not specified with -c command line option, it
 * will use interactive mode.  If some parameters are set by
 * command-line options, these values will be used as the default
 * values for the interactive mode.  Otherwise, default values are set
 * in msprior.h.  In interactive mode, pushing return will accept the
 * default values.
 *
 * If config file is specified, and parameters are set by command line
 * options, the command line options overwrite the settings of config
 * file.
 */


/*
 * If you want to make other parameters configuarable, you need to do this:
 * 1. add the parameter struct runParameters (in msprior.h)
 * 2. define the constant for the default value (in msprior.h)
 * 3. Add a statement to assign the default value to parameter in 
 *    SetDefaultParams().  All modifications for step 3-6 are in setup.c.
 * 4. Add a for loop in InteractiveSetupParams(), pay attention to var. type
 * 5. Add the entry to init_globals call inside of SetupParams()
 * 6. Add an fprintf() in PrintParam()
 */

/*
 * Loads all the appropriate configuration data to gParam and
 * gMutParam exit on failure.
 */
void
LoadConfiguration (int argc, char *argv[])
{
  FILE *fp;

  int rc;

  /* initialize the parameters with 0 */
  memset (&gParam, 0, sizeof (runParameters));
  memset (&gMutParam, 0, sizeof (gMutParam));
  memset (&gConParam, 0, sizeof (gConParam));

  /* Process the options */
  ParseCommandLine (argc, argv);

  /*
   * Load the parameteters for this simulation.  
   * We could load it from a config file.
   * If a file isn't specified (with -c), then we'll try to
   *   get the parameters interactively.
   */
  if (gParam.prngSeed == 0)
    gParam.prngSeed = -1;
  if (gParam.configFile[0])
    {
      fp = fopen (gParam.configFile, "r");
      if (!fp)
	{
	  fprintf (stderr, "Unable to open %s for reading\n",
		   gParam.configFile);
	  exit (EXIT_FAILURE);
	}

      SetupParams (fp, &gParam);
      fclose (fp);

    }
  else
    {
      InteractiveSetupParams (&gParam);
    }

  /* rand num seed */
  if (gParam.prngSeed < 0)
    (void) time (&gParam.prngSeed);

  /* read in mutation model and populate gMutParam */
  fp = fopen (gParam.configFile, "r");
  if (!fp)
    {
      fprintf (stderr, "Unable to open %s for reading\n", gParam.configFile);
      exit (EXIT_FAILURE);
    }

  rc = InitMutPara (fp, &gMutParam);
  if (rc != 0)
    {
      fprintf (stderr, "Unable to read in muataion parameters from %s\n",
	       gParam.configFile);
      exit (EXIT_FAILURE);
    }

  fclose (fp);


  if (gParam.constrain > 0)
    {
      //reopen file so it will be read from beginning
      fp = fopen (gParam.configFile, "r");

      rc = InitConPara (fp, &gConParam);
      fclose (fp);

      if (rc != 0)
	{
	  fprintf (stderr, "Unable to read in constrain paramters from %s\n",
		   gParam.configFile);
	}
    }
}

/*
 * Set the configuration file if it hasn't already been set.
 *
 * Return -1 on error,
 *         0 on success.
 */
static int
SetConfigFile (char *fName)
{
  int len;

  if (!fName || !fName[0])
    return -1;

  len = strlen (fName);
  if ((len + 1) >= MAX_FILENAME_LEN)
    return -1;

  strncpy (gParam.configFile, fName, MAX_FILENAME_LEN);

  setup_done = 1;
  return 0;
}


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
	   "\nUsage: %s [--help] [--reps N] [--debug N] [--seed N] "
	   "[--config <filename>] [--info]\n\n"
	   "       help: Print this usage function (-h)\n"
	   "       reps: Run N replications (-r)\n"
	   "      debug: Set the debug level to N (larger N => more info) (-d)\n"
	   "       seed: Set the pseudo-random number generator seed to N (-s)\n"
	   "     config: Specify a configuration file (-c)\n"
	   "       info: Print config and exit (-i)\n\n"
	   "If config file is not given, the parameter values are set "
	   "interactively.\n\n", p);
  exit (EXIT_FAILURE);
}

/*
 * Parse the command line.
 *
 */
static struct option sim_opts[] = {
  {"help", 0, NULL, 'h'},	/* list available options */
  {"debug", 1, NULL, 'd'},	/* specify debug level */
  {"seed", 1, NULL, 's'},	/* specify the prng seed */
  {"reps", 1, NULL, 'r'},	/* how many reps */
  {"config", 1, NULL, 'c'},	/* which config file to use */
  {"info", 0, NULL, 'i'},       /* print config info and exit */
  {NULL, 0, NULL, 0}
};

static void
ParseCommandLine (int argc, char *argv[])
{
  while (1)
    {
      int rc;
      long rcl;
      int opt;
      int opt_index;

      opt = getopt_long (argc, argv, "hd:s:r:c:p:i", sim_opts, &opt_index);
      if (opt < 0)
	break;

      switch (opt)
	{
	case 'h':		/* Print usage and exit */
	  PrintUsage (argv[0]);	/* This function will exit */
	  break;
	case 'd':		/* set the debug level */
	  if (!optarg)
	    {
	      fprintf (stderr, "Must select debug level\n");
	      PrintUsage (argv[0]);
	    }
	  rc = (int) strtol (optarg, NULL, 10);
	  if (errno || (rc < 0))
	    {
	      fprintf (stderr, "Invalid debug level: %s\n", optarg);
	      PrintUsage (argv[0]);
	    }
	  debug_level = rc;
	  break;
	case 's':		/* set the PRNG seed */
	  if (!optarg)
	    {
	      fprintf (stderr,
		       "Must select a pseudo-ransom number seed with -s\n");
	      PrintUsage (argv[0]);
	    }
	  rcl = strtol (optarg, NULL, 10);
	  if (errno)
	    {
	      fprintf (stderr, "Invalid pseudo-ransom number seed: %s\n",
		       optarg);
	      PrintUsage (argv[0]);
	    }
	  gParam.prngSeed = rcl;
	  break;
	case 'r':		/* Specify the number of repetitions */
	  if (!optarg)
	    {
	      fprintf (stderr, "Must select number of repetitions\n");
	      PrintUsage (argv[0]);
	    }
	  gParam.reps = strtoull (optarg, NULL, 10);
	  if (errno || (gParam.reps < 0))
	    {
	      fprintf (stderr, "Invalid number of repetitions: %s\n", optarg);
	      PrintUsage (argv[0]);
	    }
	  break;
	case 'c':		/* Specify the configuration file */
	  rc = SetConfigFile (optarg);
	  if (rc < 0)
	    {
	      fprintf (stderr, "Could not set config filename as '%s'\n",
		       optarg);
	      PrintUsage (argv[0]);
	    }
	  break;
	case 'i':              /* Print config info and exit */
	  gParam.printConf = 1;
	  break;
	default:
	  PrintUsage (argv[0]);	/* This function will exit */
	  break;
	}
    }
}

/* 
 * read a line from stdin, the string will be assigned to line[]
 * Return the length of the string
 */
int
GetLine (char *line, int max)
{
  if (fgets (line, max, stdin) == NULL)
    return 0;
  else
    return strlen (line);
}

/* 
 * If the config file is not specified with the command line option
 * (-c), this function is called from InteractiveSetupParams().  If
 * some parameter values are not set by command line options, this
 * function set gParam to the default values.  The default constants
 * are in msPriors.h
 */
static void
SetDefaultParams (runParameters * paramPtr)
{
  if (!paramPtr->upperTheta)
    {
      paramPtr->upperTheta = DEFAULT_UPPER_THETA;
    }
  if (!paramPtr->lowerTheta)
    {
      paramPtr->lowerTheta = DEFAULT_LOWER_THETA;
    }
  if (!paramPtr->upperTau)
    {
      paramPtr->upperTau = DEFAULT_UPPER_TAU;
    }
  if (!paramPtr->numTauClasses)
    {
      paramPtr->numTauClasses = 0;
    }
  if (!paramPtr->bufferTauClasses)
	{
	  paramPtr->bufferTauClasses = 0;
	}
  if (!paramPtr->concentrationShape)
	{
	  paramPtr->concentrationShape = 0;
	}
  if (!paramPtr->concentrationScale)
	{
	  paramPtr->concentrationScale= 0;
	}
  if (!paramPtr->upperMig)
    {
      paramPtr->upperMig = DEFAULT_UPPER_MIG;
    }
  if (!paramPtr->upperRec)
    {
      paramPtr->upperRec = DEFAULT_UPPER_REC;
    }
  if (!paramPtr->upperAncPopSize)
    {
      paramPtr->upperAncPopSize = DEFAULT_UPPER_ANC_POPSIZE;
    }
  if (!paramPtr->reps)
    {
      paramPtr->reps = DEFAULT_REPS;
    }
  if (!paramPtr->configFile || !paramPtr->configFile[0])
    {
      strncpy (paramPtr->configFile, DEFAULT_MUT_FILE, MAX_FILENAME_LEN);
    }
  if (!paramPtr->subParamConstrain || !paramPtr->subParamConstrain[0])
    {
      strncpy (paramPtr->subParamConstrain, DEFAULT_SUBPARAMCONSTRAIN,
	       NUMBER_OF_CONPARAM);
    }
  return;
}


/* 
 * If the config file is not specified with the command line option
 * (-c), this function set gParam interactively after setting them to
 * default values (with SetDefaultParams()).  The default constants
 * are in msPriors.h
 */
#define MAX_INPUT_LINE_LENGTH  1024
static int
InteractiveSetupParams (runParameters * paramPtr)
{
  unsigned long long tempValULL;
  unsigned int tempValUI;
  double tempValDouble;
  char line[MAX_INPUT_LINE_LENGTH], fn[MAX_FILENAME_LEN];
  int lineLen, badInput = 1;

  if (!paramPtr)
    return -1;

  fprintf (stderr, "\nPress [return] to accept the default value in [ ]\n\n");

  SetDefaultParams (paramPtr);

  /* theta */
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Lower limit of uniform prior distribution for theta "
	       "[%lf]: \n", paramPtr->lowerTheta);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%lf", &tempValDouble) == 1))
	{
	  badInput = 0;
	  paramPtr->lowerTheta = tempValDouble;
	}
    }
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Upper limit of uniform prior distribution for theta "
	       "[%lf]: \n", paramPtr->upperTheta);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%lf", &tempValDouble) == 1))
	{
	  badInput = 0;
	  paramPtr->upperTheta = tempValDouble;
	}
    }

  /* tau */
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Upper limit of uniform prior distribution for tau, time of divergence (tau-max) "
	       "[%lf]: \n", paramPtr->upperTau);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%lf", &tempValDouble) == 1))
	{
	  badInput = 0;
	  paramPtr->upperTau = tempValDouble;
	}
    }

  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Number of divergenece times across all Y of the taxon-pairs (Psi).  The hyper-parameter value should be between 1 and #taxon pairs (Y).  For example, 2 means that the model is constrained to have two divergence times (Psi=2), and each taxon pair formed at either of the two divergence times (each of which are drawn from  a uniform distribution).  Specify 0 (default) if you do not want to constrain Psi (number of divergence times), and want to draw it from the discrete uniform distribution of [1, #taxon pairs]. "
	       " [%u]: \n", paramPtr->numTauClasses);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%u", &tempValUI) == 1))
	{
	  badInput = 0;
	  paramPtr->numTauClasses = tempValUI;
	}
    }

  /* mig */
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Upper limit of uniform prior distribution for migration rate [%lf]: \n",
	       paramPtr->upperMig);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%lf", &tempValDouble) == 1))
	{
	  badInput = 0;
	  paramPtr->upperMig = tempValDouble;
	}
    }

  /* rec */
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Upper limit of uniform prior distribution for recombination rate: "
	       "[%lf]: \n", paramPtr->upperRec);

      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%lf", &tempValDouble) == 1))
	{
	  badInput = 0;
	  paramPtr->upperRec = tempValDouble;
	}
    }

  /* ancPop */
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Coefficient for the upper limit of uniform prior distribution for ancestral theta "
	       ": [%lf]\n", paramPtr->upperAncPopSize);
      fprintf (stderr,
	       "  The upper limit for ancestral theta is determined by "
	       "this coefficient\n  multiplied by the upper limit for (current) theta (%lf) : \n",
	       paramPtr->upperTheta);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%lf", &tempValDouble) == 1))
	{
	  badInput = 0;
	  paramPtr->upperAncPopSize = tempValDouble;
	}
    }

  /* rep */
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Number of draws from the Hyperprior (#simulations) [%llu]: \n",
	       paramPtr->reps);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%llu", &tempValULL) == 1))
	{
	  badInput = 0;
	  paramPtr->reps = tempValULL;
	}
    }

  /* mutation file */
  for (badInput = 1; badInput;)
    {
      fprintf (stderr,
	       "Filename of master infile (sample sizes, #bp and mutation parameters) [%s]: \n",
	       paramPtr->configFile);
      lineLen = GetLine (line, MAX_INPUT_LINE_LENGTH);
      if (lineLen == 1)
	badInput = 0;		/* use default value */
      else if ((lineLen > 1) && (sscanf (line, "%s", fn) == 1))
	{
	  if (SetConfigFile (fn) < 0)
	    {
	      fprintf (stderr, "Bad filename\n");
	      strncpy (paramPtr->configFile, DEFAULT_MUT_FILE,
		       MAX_FILENAME_LEN);
	      continue;
	    }
	  badInput = 0;
	}
    }

  return (0);
}

/*
 * Take the file pointer to the config file, and set up the
 * parameter (usually gParam).
 */
static void
SetupParams (FILE * fp, runParameters * paramPtr)
{
  int retVal;

  unsigned long r;

  r = paramPtr->reps;		/* save this in case this is set by command line */

  SetDefaultParams (paramPtr);

  retVal = init_globals (fp,
			 "lowerTheta upperTheta upperTau bufferTauClasses concentrationShape concentrationScale upperMig upperRec upperAncPopSize numTauClasses reps constrain subParamConstrain",
			 "ddddddddduVus",
			 &paramPtr->lowerTheta, &paramPtr->upperTheta,
			 &paramPtr->upperTau, &paramPtr->bufferTauClasses, 
             &paramPtr->concentrationShape, &paramPtr->concentrationScale, &paramPtr->upperMig, 
			 &paramPtr->upperRec, &paramPtr->upperAncPopSize, 
			 &paramPtr->numTauClasses, &paramPtr->reps,
			 &paramPtr->constrain, &paramPtr->subParamConstrain);

  if (retVal != 0)
    {
      if (debug_level)
	{
	  PrintParam (stderr);
	}
      fprintf (stderr, "Error reading in the parameter config file\n");
      exit (EXIT_FAILURE);
    }

  /* paramPtr->reps = 0;  Wen put this, but this will screw up */

  if (r > 0)
    {				/* over-ride with the command line option */
      paramPtr->reps = r;
    }

  return;
}

/*
 * Read in the config file.
 * A non-comment line containing "BEGIN SMAPLE_TBL" (case insensitive) indicates
 * the beginning of the table, and "END SAMPLE_TBL" indicates the end.
 *
 * Old method (kept for backward compat., works if constraints are not used):
 * It ignores empty lines, comments line (start with #), and a line 
 * with '=' sign.
 * The first line which doesn't contain '=' is considered to be the 1st
 * line with mutation parameter data.  Then initialize the mutParaArray
 * with these values.
 */

static int
InitMutPara (FILE * fp, mutParameterArray * mpaPtr)
{
  char ln[LNSZ];
  char *p;
  int tmpCol, numColumns = 0;
  int index, rc;
  mutParameter *mpp;

  while (fgets (ln, LNSZ, fp))
    {				/* read init file */
      RmLeadingSpaces (ln);
      RmExtraWhiteSpaces (ln);

      if (ln[0] == 0 || ln[0] == '#')	/* skip if blank line and */
	continue;		/* comments starting with # */

      if (strcasestr (ln, "BEGIN SAMPLE_TBL"))
	{			/* found the sample table */
	  fgets (ln, LNSZ, fp);
	  break;
	}
      /* keeping this for backward compatibility */
      p = strchr (ln, '=');	/* find equal sign */
      if (p == NULL)		/* beginning of mut data */
	break;
    }

  /* process the mutation data table */
  index = 0;
  do
    {
      /* clean up spaces */
      RmLeadingSpaces (ln);
      if (ln[0] == 0 || ln[0] == '#')	/* skip if blank line and */
	continue;		/* comments starting with # */
      RmTrailingSpaces (ln);

      tmpCol = RmExtraWhiteSpaces (ln) + 1;

      /* make sure the number of columns are correct */
#ifdef W_GAMMA
      if (tmpCol != 13)		/* should be 13 columns */
#else
      if (tmpCol != 12)		/* should be 12 columns */
#endif
	{
	  fprintf (stderr,
		   "WARN: row with incorrect number of columns encountered "
		   "in sample size & mutation model table.  Ignoring the row.\n");
	  continue;
	}
      if (!numColumns)
	numColumns = tmpCol;
      else if (numColumns != tmpCol)
	{			/* config file screwed */
	  fprintf (stderr, "The number of columns for mutation data should be"
		   "%d, but following line has %d columns\n%s",
		   numColumns, tmpCol, ln);
	  exit (EXIT_FAILURE);
	}

      /* check if mem allocated */
      rc = CheckMutParaArray (mpaPtr, index);
      if (rc < 0)
	{
	  fprintf (stderr, "Error found in CheckMutParaArray\n");
	}			/* error */

      mpp = &(mpaPtr->data[index]);

      rc = ReadMutLine (mpp, ln, numColumns);

      if (rc < 0)
	{
	  fprintf (stderr, "WARN: The following is weird, ignoring\n%s\n",
		   ln);
	}
      else
	{
	  index++;
	}

    }
  while (fgets (ln, LNSZ, fp) && (!strcasestr (ln, "END SAMPLE_TBL")));

  rc = ProcessTaxonLociInMutParaArray (mpaPtr);
  if (rc < 0)			/* error */
    {
      fprintf (stderr,
	       "In InitMutPara of msprior, Trying to Process the SAMPLE_TBL, and error encountered: Err code = %d\n",
	       rc);
      return -1;
    }

  gParam.numTaxonLocusPairs = mpaPtr->numElements;

  if (gParam.numTaxonLocusPairs < 1)
    {
      /* didn't find any sample size, mut para entries */
      return -1;
    }

  return 0;
}

/* primary goal is to create locus table */
static int
ProcessTaxonLociInMutParaArray (mutParameterArray * mpaPtr)
{
  int row, col, index, numUniqLoc, numUniqTaxon;
  mutParameter *mpp;

  char **locusNames, **taxonNames, **uniqLocNames, **uniqTaxonNames;

  /* cmatrix in sumStatsVector.c */
  locusNames = cmatrix (mpaPtr->numElements, MAX_NAME_CHAR_LEN);
  uniqLocNames = cmatrix (mpaPtr->numElements, MAX_NAME_CHAR_LEN);
  taxonNames = cmatrix (mpaPtr->numElements, MAX_NAME_CHAR_LEN);
  uniqTaxonNames = cmatrix (mpaPtr->numElements, MAX_NAME_CHAR_LEN);

  for (index = 0; index < mpaPtr->numElements; index++)
    {
      mpp = &mpaPtr->data[index];
      /* copy all locusNames to an array */
      strcpy (locusNames[index], mpp->locusName);
      /* copy all species Names to an array */
      strcpy (taxonNames[index], mpp->taxonName);
    }

  /* Extract Unique locus, spNames */
  numUniqTaxon =
    UniqueStrings (taxonNames, uniqTaxonNames, mpaPtr->numElements);
  numUniqLoc = UniqueStrings (locusNames, uniqLocNames, mpaPtr->numElements);

  /* fill out the hashTbl in mpPtr */
  mpaPtr->taxonIDTbl = ht_init (numUniqTaxon, NULL);
  mpaPtr->locusIDTbl = ht_init (numUniqLoc, NULL);
  if (mpaPtr->taxonIDTbl == NULL || mpaPtr->locusIDTbl == NULL)
    {
      fprintf (stderr,
	       "In ProcessTaxonLociInMutParaArray, no mem for HashTbls\n");
      exit (EXIT_FAILURE);
    }

  /* fill up a hash table, key is taxon name, value is index */
  for (index = 0; index < numUniqTaxon; index++)
    {
      int *val;
      val = (int *) malloc (sizeof (int));
      if (val == NULL)
	{
	  fprintf (stderr,
		   "ERROR: Not enough memory in ProcessTRaxonLociInMutParaArray\n");
	  exit (EXIT_FAILURE);
	}

      *val = index;
      if (ht_insert (mpaPtr->taxonIDTbl, uniqTaxonNames[index],
		     strlen (uniqTaxonNames[index]),
		     val, sizeof (int)) == NULL)
	{
	  fprintf (stderr,
		   "ERROR: Not enough memory in ProcessTRaxonLociInMutParaArray, inserting to hash table\n");
	  exit (EXIT_FAILURE);
	}
    }

  /* fill up a hash table, key is locus name, value is index */
  for (index = 0; index < numUniqLoc; index++)
    {
      int *val;
      val = (int *) malloc (sizeof (int));
      if (val == NULL)
	{
	  fprintf (stderr,
		   "ERROR: Not enough memory in ProcessTaxonLociInMutParaArray\n");
	  exit (EXIT_FAILURE);
	}

      *val = index;
      if (ht_insert (mpaPtr->locusIDTbl, uniqLocNames[index],
		     strlen (uniqLocNames[index]), val, sizeof (int)) == NULL)
	{
	  fprintf (stderr,
		   "ERROR: Not enough memory in ProcessTaxonLociInMutParaArray, inserting to hash table\n");
	  exit (EXIT_FAILURE);
	}
    }

  /* go through each line and assign the locus/taxonID, value is index */
  /* initialize locTbl, fill with -1 */
  mpaPtr->locTbl = malloc (sizeof (lociTbl));
  if (mpaPtr->locTbl == NULL)
    {
      fprintf (stderr,
	       "ERROR: not enough memory for lociTbl in ProcessTaxonLociInMutParaArray\n");
      exit (EXIT_FAILURE);
    }
  mpaPtr->locTbl->numTaxon = numUniqTaxon;
  mpaPtr->locTbl->numLoci = numUniqLoc;
  /* These two values need to be copied to numTaxonPairs, numTaxonPairs of gParam */
  
  gParam.numTaxonPairs = numUniqTaxon;
  gParam.numLoci = numUniqLoc;

  mpaPtr->locTbl->tbl = (int **) calloc (numUniqTaxon, sizeof (int *));
  if (mpaPtr->locTbl->tbl == NULL)
    {
      fprintf (stderr,
	       "ERROR: not enough memory (2) for lociTbl in ProcessTaxonLociInMutParaArray\n");
      exit (EXIT_FAILURE);
    }
  for (row = 0; row < numUniqTaxon; row++)
    {
      mpaPtr->locTbl->tbl[row] = (int *) calloc (numUniqLoc, sizeof (int));
      if (mpaPtr->locTbl->tbl[row] == NULL)
	{
	  fprintf (stderr,
		   "ERROR: not enough memory (3) for lociTbl in ProcessTaxonLociInMutParaArray\n");
	  exit (EXIT_FAILURE);
	}
      for (col = 0; col < numUniqLoc; col++)
	{
	  mpaPtr->locTbl->tbl[row][col] = -1;
	}
    }

  /* fill up the locTbl, with the correct 0-offset index */
  for (index = 0; index < mpaPtr->numElements; index++)
    {
      int *locIndex, *taxonIndex;
      mpp = &mpaPtr->data[index];

      taxonIndex = (int *) ht_search (mpaPtr->taxonIDTbl, mpp->taxonName,
				      strlen (mpp->taxonName));
      locIndex = (int *) ht_search (mpaPtr->locusIDTbl, mpp->locusName,
				    strlen (mpp->locusName));
      /* update the numerical IDs in the mutParameter */
      mpp->taxonID = (unsigned int) (*taxonIndex);
      mpp->locusID = (unsigned int) (*locIndex);

      if (mpaPtr->locTbl->tbl[*taxonIndex][*locIndex] >= 0)
	{
	  fprintf (stderr,
		   "WARN: In the configuration file, there is multiple lines for the following taxon:locus = %s : %s.  The first info is used.\n",
		   mpp->taxonName, mpp->locusName);
	  mpp->taxonID = mpp->locusID = -1;	/* note it is unsigned int */
	}
      else
	{
	  mpaPtr->locTbl->tbl[*taxonIndex][*locIndex] = index;
	}
    }

  /* clean the momory of charMats */
  freeCMatrix (mpaPtr->numElements, locusNames);
  freeCMatrix (mpaPtr->numElements, uniqLocNames);
  freeCMatrix (mpaPtr->numElements, taxonNames);
  freeCMatrix (mpaPtr->numElements, uniqTaxonNames);

  return 0;
}

static int
InitConPara (FILE * fp, constrainedParameterArray * cpaPtr)
{
  char ln[LNSZ];		// char array with length = 256(LNSZ)

  int index, rc, foundConstraints = 0;
  int tmpCol, numColumns = 0;
  constrainedParameter *cpp;


  while (fgets (ln, LNSZ, fp))
    {
      RmLeadingSpaces (ln);
      RmExtraWhiteSpaces (ln);

      if (ln[0] == 0 || ln[0] == '#')	/* skip if blank line and */
	continue;		/* comments starting with # */

      /* found the beginning of constrain table */
      if (strcasestr (ln, "BEGIN CONSTRAIN"))
	{
	  foundConstraints = 1;
	  fgets (ln, LNSZ, fp);
	  break;
	}
    }

  if (!foundConstraints)
    {				/* no constraints */
      return -1;
    }

  index = 0;
  do
    {
      RmLeadingSpaces (ln);
      RmTrailingSpaces (ln);

      if (ln[0] == 0 || ln[0] == '#')
	continue;

      tmpCol = RmExtraWhiteSpaces (ln) + 1;

      // printout tmpCol to see the value

      if (!numColumns)
	numColumns = tmpCol;
      else if (numColumns != tmpCol)
	{
	  fprintf (stderr,
		   "The number of columns for mutation data should be %d, but following line has %d columns\n%s",
		   numColumns, tmpCol, ln);
	  exit (EXIT_FAILURE);
	}

      // check might not be used
      rc = CheckConParaArray (cpaPtr, index);
      if (rc < 0)
	{
	  fprintf (stderr, "Hmm, error found in CheckConParaArray\n");
	}

      cpp = &cpaPtr->conData[index];

      rc = ReadConLine (cpp, ln, numColumns);
      if (rc < 0)
	{
	  fprintf (stderr, "Hey, error in reading constrain parameters\n%s\n",
		   ln);
	}
      else
	index++;

    }
  while (fgets (ln, LNSZ, fp) && (!strcasestr (ln, "END CONSTRAIN")));

  return 0;
}


/* 
 * Check if mpp->data[index] is allocated.
 * If the array mpp->data isn't big enough, it will increase
 * the size of the array by increment of growBy until the array
 * is big enough to use index-th element.
 * It will update the mpaPtr->numElements and mpaPtr->numAllocated if needed.
 * Returns:
 *  0 for success
 *  1 for error
 */
static int
CheckMutParaArray (mutParameterArray * mpaPtr, int index)
{
  const int growBy = 10;

  if (index < 0)
    {
      fprintf (stderr,
	       "In CheckMutParaArray(), 2nd arg should be non-negative.\n");
      return -1;
    }

  if (!mpaPtr)
    {
      fprintf (stderr, "In CheckMutParaArray(), NULL pointer is given\n");
      return -1;
    }

  if (index + 1 <= mpaPtr->numElements)
    {
      return 0;			/* the element is allocated, in use */
    }
  else if (index + 1 <= mpaPtr->numAllocated)
    {
      mpaPtr->numElements = index + 1;
      return 0;
    }
  else
    {				/* all alocated ones are used up, so grow the memory */
      size_t newNumEle = (1 + (index + 1) / growBy) * growBy;
      mutParameter *rc;
      rc =
	(mutParameter *) realloc (mpaPtr->data,
				  sizeof (mutParameter) * newNumEle);
      if (!rc)
	{
	  fprintf (stderr, "Not enough memory in CheckMutParaArray()\n");
	  exit (EXIT_FAILURE);
	}

      mpaPtr->data = rc;
      mpaPtr->numAllocated = newNumEle;
      mpaPtr->numElements = index + 1;

      return 0;
    }
  return -1;			/* It should never come here */
}


static int
CheckConParaArray (constrainedParameterArray * cpaPtr, int index)
{
  int growBy = 10;


  if (index < 0)
    {
      fprintf (stderr,
	       "In CheckConArray(), 2nd arg should be non-negative. \n");
      return -1;
    }

  if (!cpaPtr)
    {
      fprintf (stderr, "In CheckConArray(), NULL pointer is givin\n");
    }

  if (index + 1 <= cpaPtr->conNumElements)
    {
      return 0;
    }
  else if (index + 1 <= cpaPtr->conNumAllocated)
    {
      cpaPtr->conNumElements = index + 1;
      return 0;
    }
  else
    {
      size_t newNumEle = (1 + (index + 1) / growBy) * growBy;
      constrainedParameter *cc;

      cc =
	(constrainedParameter *) realloc (cpaPtr->conData,
					  sizeof (constrainedParameter) *
					  newNumEle);
      if (!cc)
	{
	  fprintf (stderr, "Not enough memory in CheckConParaArray()\n");
	}

      cpaPtr->conData = cc;
      cpaPtr->conNumAllocated = newNumEle;
      cpaPtr->conNumElements = index + 1;

      return 0;
    }
  return -1;
}



/*
 * Take a line which specify the mutation model and sample sizes of a
 * pair of taxon, and assign the values to mpp.
 */
static int
ReadMutLine (mutParameter * mpp, char *line, int ncol)
{
  /* probably we should check integer, double is correct in the file */
  int rc;
  /* char dummyTaxonPairName[1024]; */


#ifdef W_GAMMA

  if (ncol == 13)
    {				/* 13 column */
      rc = sscanf (line, "%s %s %lf %lf %u %u %lf %lf %u %lf %lf %lf %s",
		   mpp->taxonName, mpp->locusName, &mpp->NScaler, 
		   &mpp->mutScaler,
		   &mpp->sample[0], &mpp->sample[1],
		   &mpp->tstv[0], &mpp->gamma, &mpp->seqLen,
		   &mpp->freqA, &mpp->freqC, &mpp->freqG, mpp->filename);
      mpp->numPerTaxa = mpp->sample[0] + mpp->sample[1];
    }
  else
    {
      return (-1);
    }

#else /* basically sscanf and ncol == are different */

  if (ncol == 12)
    {
      rc = sscanf (line, "%s %s %lf %lf %u %u %lf %u %lf %lf %lf %s",
		   mpp->taxonName, mpp->locusName, &mpp->NScaler,
		   &mpp->mutScaler,
		   &mpp->sample[0], &mpp->sample[1],
		   &mpp->tstv[0], &mpp->seqLen,
		   &mpp->freqA, &mpp->freqC, &mpp->freqG, mpp->filename);
      mpp->gamma = 999;
      mpp->numPerTaxa = mpp->sample[0] + mpp->sample[1];
    }
  else
    {
      return (-1);
    }

#endif

  mpp->freqT = 1 - mpp->freqA - mpp->freqC - mpp->freqG;


  return (rc);
}


static int
ReadConLine (constrainedParameter * cpp, char *line, int ncol)
{

  int cc;


  if (ncol == 9)		// so far there are nine values
    {
      cc = sscanf (line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &cpp->conTau, &cpp->conBottPop1, &cpp->conBottPop2,
		   &cpp->conBottleTime, &cpp->conMig, &cpp->conTheta,
		   &cpp->conN1, &cpp->conNanc, &cpp->conRec);
    }

  return (cc);
}


/*
 * print out the value of parameters, useful for debugging
 */
void
PrintParam (FILE *fp)
{
  int i;

  fprintf (fp, "## gParam ##\n");
  fprintf (fp, "lowerTheta =\t%.17lf\n", gParam.lowerTheta);
  fprintf (fp, "upperTheta =\t%.17lf\n", gParam.upperTheta);
  fprintf (fp, "upperTau =\t%.17lf\n", gParam.upperTau);
  fprintf (fp, "upperMig =\t%.17lf\n", gParam.upperMig);
  fprintf (fp, "upperRec =\t%.17lf\n", gParam.upperRec);
  fprintf (fp, "upperAncPopSize =\t%.17lf\n", gParam.upperAncPopSize);
  fprintf (fp, "reps =\t%llu\n", gParam.reps);
  fprintf (fp, "numTaxonLocusPairs =\t%u\n", gParam.numTaxonLocusPairs);
  fprintf (fp, "numTaxonPairs =\t%u\n", gParam.numTaxonPairs);
  fprintf (fp, "numLoci =\t%u\n", gParam.numLoci);
  fprintf (fp, "numTauClasses =\t%u\n", gParam.numTauClasses);
  fprintf (fp, "bufferTauClasses =\t%.17lf\n", gParam.bufferTauClasses);
  fprintf (fp, "concentrationShape =\t%.17lf\n", gParam.concentrationShape);
  fprintf (fp, "concentrationScale =\t%.17lf\n", gParam.concentrationScale);
  fprintf (fp, "prngSeed =\t%ld\n", gParam.prngSeed);
  fprintf (fp, "configFile =\t%s\n", gParam.configFile);
  fprintf (fp, "scratchFile =\t%s\n", gParam.scratchFile);
  fprintf (fp, "constrain =\t%u\n", gParam.constrain);
  fprintf (fp, "subParamConstrain =\t%s\n", gParam.subParamConstrain);
  fprintf (fp, "printConf =\t%u\n", gParam.printConf);
  
  fprintf (fp, "## gMutParam ##\n");
  for (i = 0; i < gMutParam.numElements; i++)
    {
      fprintf (fp, "### taxon:locus pair ID %d taxonID %u (%s) locusID %u (%s) NScaler %lf mutScaler %lf ###\n",
	       i + 1, gMutParam.data[i].taxonID+1, gMutParam.data[i].taxonName,
	       gMutParam.data[i].locusID+1, gMutParam.data[i].locusName,
	       gMutParam.data[i].NScaler, gMutParam.data[i].mutScaler);
      fprintf (fp, "numPerTaxa =\t%u\n", gMutParam.data[i].numPerTaxa);
      fprintf (fp, "sample =\t%u %u\n",
	       gMutParam.data[i].sample[0], gMutParam.data[i].sample[1]);
      fprintf (fp, "tstv =\t%lf  %lf\n",
	       gMutParam.data[i].tstv[0], gMutParam.data[i].tstv[1]);
      fprintf (fp, "gamma =\t%lf\n", gMutParam.data[i].gamma);
      fprintf (fp, "seqLen =\t%u\n", gMutParam.data[i].seqLen);
      fprintf (fp, "freq:A, C, G, T = %f, %f %f %f\n",
	       gMutParam.data[i].freqA, gMutParam.data[i].freqC,
	       gMutParam.data[i].freqG, gMutParam.data[i].freqT);
      fprintf (fp, "fileName =\t%s\n", gMutParam.data[i].filename);
    }

  fprintf (fp, "## gConParam, for contraint only ##\n");
  for (i = 0; i < gConParam.conNumElements; i++)
    {
      fprintf (fp, "### taxon pair %d ###\n", i + 1);
      fprintf (fp, "tau = \t%lf\n", gConParam.conData[i].conTau);
      fprintf (fp, "bottle neck populations = \t%lf %lf\n",
	       gConParam.conData[i].conBottPop1,
	       gConParam.conData[i].conBottPop2);
      fprintf (fp, "bottle neck time = \t%lf\n",
	       gConParam.conData[i].conBottleTime);
      fprintf (fp, "migration rate= \t%lf\n",
	       gConParam.conData[i].conMig);
      fprintf (fp, "theta = \t%lf\n", gConParam.conData[i].conTheta);
      fprintf (fp, "current population sizes = \t%lf\n",
	       gConParam.conData[i].conN1);
      fprintf (fp, "ancestral population size = \t%lf\n",
	       gConParam.conData[i].conNanc);
      fprintf (fp, "recombination rate = \t%lf\n",
	       gConParam.conData[i].conRec);
    }
}

#ifdef TEST_SETUP

runParameters gParam;
mutParameterArray gMutParam;

int
main (int argc, char *argv[])
{
  LoadConfiguration (argc, argv);
  PrintParam (stderr);

}

#endif
