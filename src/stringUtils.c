/*
 * uniqueStrings.c
 *
 * Copyright (C) 2007  Wen Huang, Naoki Takebayashi and Michael Hickerson
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "stringUtils.h"

#ifdef TEST_UNIQUE_STRINGS
int
main ()
{
  char **inputArr;
  char **returnArr;
  int i, j, NumOfString;

  if (!(inputArr = (char **) malloc ((unsigned) 6 * sizeof (char *))))
    {
      fprintf (stderr, "malloc error 1\n");
      exit (1);
    }

  for (i = 0; i < 6; i++)
    {
      if (!(inputArr[i] = (char *) malloc ((unsigned) 10 * sizeof (char))))
	{
	  //cout<<"malloc error 2"<<endl;
	  fprintf (stderr, "malloc error 2\n");
	  exit (1);
	}
    }

  inputArr[0] = "Tom", inputArr[1] = "Sam", inputArr[2] = "Amy", inputArr[3] =
    "Tom", inputArr[4] = "Sam", inputArr[5] = "Tracey";

  NumOfString = UniqueStrings (inputArr, returnArr, 6);

  printf ("inputArr\n");
  for (i = 0; i < 6; i++)
    printf (inputArr[i]);
  printf ("\n");

  printf ("Number of unique strings is:%d\n", NumOfString);

  return 0;

}
#endif /* TEST_UNIQUE_STRINGS */

/*
 * Count the number of unique strings in array of string
 * with duplicated elements, the program will terminates 
 * if failed to allocate memory for objects
 *
 * Arguments:
 *  inputCharArr: array of string to be analysed
 *  returnCharArr: array contains only unique string of inputCharArr
 *                 memory should be pre allocated.
 *  size: number of strings in inputCharArr
 *
 * Returns: the number of strings in returnCharArr
 */

int
UniqueStrings (char **inputCharArr, char **returnCharArr, int size)
{
  int count, indexGrow;
  int i, j;
  int copied;

  if (inputCharArr == NULL)
    {
      fprintf (stderr, "ERROR: inputCharArr empty in UniqueStrings\n");
      return -1;
    }

  /* Get count: number of unique strings in the array */
  count = 0;
  for (i = 0; i < size; i++)
    {
      int uniqQ = 1;
      for (j = i + 1; j < size; j++)
	{
	  if (strcmp (inputCharArr[i], inputCharArr[j]) == 0)
	    {			/* match */
	      uniqQ = 0;
	      break;
	    }
	}
      count += uniqQ;
    }

  /* extracting and copying from inputCharArr */
  /* Assume that memory is already allocated to returnCharArr by cmatrix() */
  indexGrow = 0;
  for (i = 0; i < size; i++)
    {
      copied = 0;
      for (j = 0; j < indexGrow; j++)
	{
	  if (!strcmp (returnCharArr[j], inputCharArr[i]))
	    {
	      copied = 1;
	      break;
	    }
	}

      if (!copied)
	{
	  strcpy (returnCharArr[indexGrow], inputCharArr[i]);
	  indexGrow++;
	}
    }

  return count;

}

char **
cmatrix (int nsam, int len)
{
  int i;
  char **m;
  if (!(m = (char **) malloc ((unsigned) (nsam * sizeof (char *)))))
    {
      perror ("alloc error in cmatrix");
      return NULL;
    }
  for (i = 0; i < nsam; i++)
    {
      if (!(m[i] = (char *) malloc ((unsigned) (len * sizeof (char)))))
	{
	  perror ("alloc error in cmatric. 2");
	  return NULL;
	}
    }
  return (m);
}

void
freeCMatrix (int nsam, char **list)
{
  int i;

  for (i = 0; i < nsam; i++)
    free (list[i]);

  free (list);
  return;
}
