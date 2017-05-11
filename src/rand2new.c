/*
 * rand2new.c
 *
 * Copyright (C) 2006  Michael Hickerson
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
/* modified by Naoki Takebayashi */


#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>	/* for base rand num gen's */

const gsl_rng *gBaseRand;	/* global rand number generator */

double
ran1 ()
{
  return (gsl_rng_uniform (gBaseRand));
}


int
seedit (unsigned long seed)
{
  gBaseRand = gsl_rng_alloc (gsl_rng_mt19937);	/* set the base PRNG to
						   Mersenne Twister */
  gsl_rng_set (gBaseRand, seed);	/* seed the PRNG */
  return 0;
}

int
cleanPRNG (void)
{
  gsl_rng_free ((gsl_rng *) gBaseRand);
  return 0;
}

/* unsigned int */
/* seedit (const char *flag) */
/* { */
/*   FILE *fopen (), *pfseed; */
/*   unsigned int seed2; */
/*   if (flag[0] == 't') */
/*     { */
/*       pfseed = fopen ("seedms", "r"); */
/*       if (pfseed == NULL) */
/* 	{ */
/* 	  seed2 = 59243; */
/* 	} */
/*       else */
/* 	{ */
/* 	  fscanf (pfseed, " %d", &seed2); */
/* 	  fclose (pfseed); */
/* 	} */
/*       srandom (seed2); */
/*       return (seed2); */
/*     } */
/*   else */
/*     { */
/*       pfseed = fopen ("seedms", "w"); */
/*       fprintf (pfseed, "%ld \n", random ()); */
/*       return (0); */
/*     } */
/* } */
