/*
 * QHsubs.c
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

#include "msQHarbpop.h"
#include "QHsubs.h"

void
assign_QHsiterate (double theta, int nsites, struct QHparameters *QH)
{
  double ucoeff;
  ucoeff =
    2 * (QH->freqA + QH->freqG) * (QH->freqC + QH->freqT) +
    2 * (QH->Rtstv) * (QH->freqA * QH->freqG) +
    2 * (QH->Ytstv) * (QH->freqC * QH->freqT);
  QH->siterate = theta / nsites / ucoeff;
  return;
}


void
createrandomseQ (char *seq, int length, struct QHparameters *QH)
{
  int base;
  double ran1 (), x, prob;
  for (base = 0; base < length; base++)
    {
      while ((x = ran1 ()) == 1);
      if (x < (prob = QH->freqA))
	seq[base] = 'A';
      else if (x < (prob += QH->freqC))
	seq[base] = 'C';
      else if (x < (prob += QH->freqG))
	seq[base] = 'G';
      else
	seq[base] = 'T';
    }
}


char
mutateProb (char ancbase, double ut, struct QHparameters *QH)
{
  char newbase = '?';
  double PR, PY;
  double eut, eutRA, eutYA;
  double ran1 (), x, prob;
  PR = (QH->freqA) + (QH->freqG);
  PY = (QH->freqC) + (QH->freqT);
  eut = exp (-1. * ut);
  eutRA = exp (-1. * ut * (PR * (QH->Rtstv) + PY));
  eutYA = exp (-1. * ut * (PY * (QH->Ytstv) + PR));
  while ((x = ran1 ()) == 1.);
  switch (ancbase)
    {
    case 'A':
      if (x <
	  (prob =
	   (QH->freqA) + (QH->freqA) * PY * eut / PR +
	   (QH->freqG) * eutRA / PR))
	newbase = 'A';
      else if (x <
	       (prob +=
		(QH->freqG) + (QH->freqG) * PY * eut / PR -
		(QH->freqG) * eutRA / PR))
	newbase = 'G';
      else if (x < (prob += (QH->freqC) * (1. - eut)))
	newbase = 'C';
      else
	newbase = 'T';
      break;
    case 'C':
      if (x <
	  (prob =
	   (QH->freqC) + (QH->freqC) * PR * eut / PY +
	   (QH->freqT) * eutYA / PY))
	newbase = 'C';
      else if (x <
	       (prob +=
		(QH->freqT) + (QH->freqT) * PR * eut / PY -
		(QH->freqT) * eutYA / PY))
	newbase = 'T';
      else if (x < (prob += (QH->freqA) * (1. - eut)))
	newbase = 'A';
      else
	newbase = 'G';
      break;
    case 'G':
      if (x <
	  (prob =
	   (QH->freqG) + (QH->freqG) * PY * eut / PR +
	   (QH->freqA) * eutRA / PR))
	newbase = 'G';
      else if (x <
	       (prob +=
		(QH->freqA) + (QH->freqA) * PY * eut / PR -
		(QH->freqA) * eutRA / PR))
	newbase = 'A';
      else if (x < (prob += (QH->freqC) * (1. - eut)))
	newbase = 'C';
      else
	newbase = 'T';
      break;
    case 'T':
      if (x <
	  (prob =
	   (QH->freqT) + (QH->freqT) * PR * eut / PY +
	   (QH->freqC) * eutYA / PY))
	newbase = 'T';
      else if (x <
	       (prob +=
		(QH->freqC) + (QH->freqC) * PR * eut / PY -
		(QH->freqC) * eutYA / PY))
	newbase = 'C';
      else if (x < (prob += (QH->freqA) * (1. - eut)))
	newbase = 'A';
      else
	newbase = 'G';
      break;
    default:
      perror ("mutateProb switch (ancbase): no case satisfied\n");
    }
  return newbase;
}


void
calcPij (double Pij[4][4], double ut, struct QHparameters *QH)
{
  double PA, PC, PG, PT, PR, PY;
  double eut, eutRA, eutYA;

  PA = QH->freqA;
  PC = QH->freqC;
  PG = QH->freqG;
  PT = QH->freqT;
  PR = PA + PG;
  PY = PC + PT;
  eut = exp (-1. * ut);
  eutRA = exp (-1. * ut * (PR * (QH->Rtstv) + PY));
  eutYA = exp (-1. * ut * (PY * (QH->Ytstv) + PR));

  Pij[0][0] = PA + PA * PY * eut / PR + PG * eutRA / PR;
  Pij[0][1] = PC * (1. - eut);
  Pij[0][2] = PG + PG * PY * eut / PR - PG * eutRA / PR;
  Pij[0][3] = PT * (1. - eut);
  if (fabs (Pij[0][0] + Pij[0][1] + Pij[0][2] + Pij[0][3] - 1) > 1e-6)
    perror ("calcPij: row0 doesn't sum to 1");
  Pij[1][0] = PA * (1. - eut);
  Pij[1][1] = PC + PC * PR * eut / PY + PT * eutYA / PY;
  Pij[1][2] = PG * (1. - eut);
  Pij[1][3] = PT + PT * PR * eut / PY - PT * eutYA / PY;
  if (fabs (Pij[1][0] + Pij[1][1] + Pij[1][2] + Pij[1][3] - 1) > 1e-6)
    perror ("calcPij: row1 doesn't sum to 1");
  Pij[2][0] = PA + PA * PY * eut / PR - PA * eutRA / PR;
  Pij[2][1] = PC * (1. - eut);
  Pij[2][2] = PG + PG * PY * eut / PR + PA * eutRA / PR;
  Pij[2][3] = PT * (1. - eut);
  if (fabs (Pij[2][0] + Pij[2][1] + Pij[2][2] + Pij[2][3] - 1) > 1e-6)
    perror ("calcPij: row2 doesn't sum to 1");
  Pij[3][0] = PA * (1. - eut);
  Pij[3][1] = PC + PC * PR * eut / PY - PC * eutYA / PY;
  Pij[3][2] = PG * (1. - eut);
  Pij[3][3] = PT + PT * PR * eut / PY + PC * eutYA / PY;
  if (fabs (Pij[3][0] + Pij[3][1] + Pij[3][2] + Pij[3][3] - 1) > 1e-6)
    perror ("calcPij: row3 doesn't sum to 1");
}

char
mutatePij (char ancbase, double Pij[4][4])
{
  char newbase = '?';
  double ran1 (), x, prob;

  while ((x = ran1 ()) == 1.);
  switch (ancbase)
    {
    case 'A':
      if (x < (prob = Pij[0][0]))
	newbase = 'A';
      else if (x < (prob += Pij[0][2]))
	newbase = 'G';
      else if (x < (prob += Pij[0][1]))
	newbase = 'C';
      else
	newbase = 'T';
      break;
    case 'C':
      if (x < (prob = Pij[1][1]))
	newbase = 'C';
      else if (x < (prob += Pij[1][3]))
	newbase = 'T';
      else if (x < (prob += Pij[1][0]))
	newbase = 'A';
      else
	newbase = 'G';
      break;
    case 'G':
      if (x < (prob = Pij[2][2]))
	newbase = 'G';
      else if (x < (prob += Pij[2][0]))
	newbase = 'A';
      else if (x < (prob += Pij[2][1]))
	newbase = 'C';
      else
	newbase = 'T';
      break;
    case 'T':
      if (x < (prob = Pij[3][3]))
	newbase = 'T';
      else if (x < (prob += Pij[3][1]))
	newbase = 'C';
      else if (x < (prob += Pij[3][0]))
	newbase = 'A';
      else
	newbase = 'G';
      break;
    default:
      perror ("mutatePij switch (ancbase): no case satisfied\n");
    }
  return newbase;
}


double
calcRate (char base, struct QHparameters *QH)
{
  double r = 0.0;
  switch (base)
    {
    case 'A':
      r = (QH->freqC + (QH->Rtstv) * (QH->freqG) + QH->freqT);
      break;
    case 'C':
      r = (QH->freqA + QH->freqG + (QH->Ytstv) * (QH->freqT));
      break;
    case 'G':
      r = ((QH->Rtstv) * (QH->freqA) + QH->freqC + QH->freqT);
      break;
    case 'T':
      r = (QH->freqA + (QH->Ytstv) * (QH->freqC) + QH->freqG);
      break;
    default:
      perror ("calcRate: base failed all cases\n");
    }
  return r;
}


char
mutateRate (char ancbase, struct QHparameters *QH)
{
  char newbase = '?';
  double ran1 (), x, prob;

  while ((x = ran1 ()) == 1.);
  switch (ancbase)
    {
    case 'A':
      x *= QH->freqC + (QH->Rtstv) * (QH->freqG) + QH->freqT;
      if (x < (prob = (QH->Rtstv) * (QH->freqG)))
	newbase = 'G';
      else if (x < (prob += QH->freqC))
	newbase = 'C';
      else
	newbase = 'T';
      break;
    case 'C':
      x *= QH->freqA + QH->freqG + (QH->Ytstv) * (QH->freqT);
      if (x < (prob = (QH->Ytstv) * (QH->freqT)))
	newbase = 'T';
      else if (x < (prob += QH->freqA))
	newbase = 'A';
      else
	newbase = 'G';
      break;
    case 'G':
      x *= (QH->Rtstv) * (QH->freqA) + QH->freqC + QH->freqT;
      if (x < (prob = (QH->Rtstv) * (QH->freqA)))
	newbase = 'A';
      else if (x < (prob += QH->freqC))
	newbase = 'C';
      else
	newbase = 'T';
      break;
    case 'T':
      x *= QH->freqA + (QH->Ytstv) * (QH->freqC) + QH->freqG;
      if (x < (prob = (QH->Ytstv) * (QH->freqC)))
	newbase = 'C';
      else if (x < (prob += QH->freqA))
	newbase = 'A';
      else
	newbase = 'G';
      break;
    default:
      perror ("\nmutateRate switch (ancbase): no case satisfied");
    }
  return newbase;
}


void
mutateBelow (int n, struct node *tree, int nd, char *seqs)
{
  int ni, done;
  for (ni = (nd - 1), done = 0; (ni >= 0) && (done < 2); ni--)
    {
      if (((tree + ni)->abv == nd))
	{
	  seqs[ni] = seqs[nd];
	  if ((tree + ni)->ndes > 1)
	    mutateBelow (n, tree, ni, seqs);
	  done++;
	}
    }
  return;
}


int
baseint (char base)
{
  int value = -1;
  switch (base)
    {
    case 'A':
      value = 0;
      break;
    case 'C':
      value = 1;
      break;
    case 'G':
      value = 2;
      break;
    case 'T':
      value = 3;
      break;
    }
  return (value);
}


char
intbase (int base)
{
  char value = '?';
  switch (base)
    {
    case 0:
      value = 'A';
      break;
    case 1:
      value = 'C';
      break;
    case 2:
      value = 'G';
      break;
    case 3:
      value = 'T';
      break;
    }
  return (value);
}

double
rndgamma1 (double s)
{
/* random standard gamma for s<1
   switching method
*/
  double ran1 (), r, x = 0., small = 1e-37, w;
  static double a, p, uf, ss = 10., d;
  if (s != ss)
    {
      a = 1 - s;
      p = a / (a + s * exp (-a));
      uf = p * pow (small / a, s);
      d = a * log (a);
      ss = s;
    }
  for (;;)
    {
      while (((r = ran1 ()) == 0) || (r == 1));
      if (r > p)
	x = a - log ((1. - r) / (1. - p)), w = a * log (x) - d;
      else if (r > uf)
	x = a * pow (r / p, 1. / s), w = x;
      else
	return (0);
      while (((r = ran1 ()) == 0) || (r == 1));
      if (1. - r <= w && r > 0)
	if (r * (w + 1.) >= 1 || -log (r) <= w)
	  continue;
      break;
    }
  return (x);
}

double
rndgamma2 (double s)
{
/* random standard gamma for s>1
   Best's (1978) t distribution method */
  double ran1 (), r, d, f, g, x;
  static double b, h, ss = 0;
  if (s != ss)
    {
      b = s - 1;
      h = sqrt (3 * s - 0.75);
      ss = s;
    }
  for (;;)
    {
      while (((r = ran1 ()) == 0) || (r == 1));
      g = r - r * r;
      f = (r - 0.5) * h / sqrt (g);
      x = b + f;
      if (x <= 0)
	continue;
      while (((r = ran1 ()) == 0) || (r == 1));
      d = 64 * r * r * g * g * g;
      if (d * x < x - 2 * f * f || log (d) < 2 * (b * log (x / b) - f))
	break;
    }
  return (x);
}

double
gammadev (double s)
{
/* random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
      r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
           rndgamma2() if s>1 or
           exponential if s=1 */
  double ran1 (), r = 0, x;
  if (s <= 0)
    return (-1);
  else if (s < 1)
    r = rndgamma1 (s);
  else if (s > 1)
    r = rndgamma2 (s);
  else
    {
      while (((x = ran1 ()) == 0) || (x == 1));
      r = -log (ran1 ());
    }
/* 1Nov01 Eli Stahl- divide by s before returning. */
  r /= s;
  return (r);
}
