/* main.c - Hashtable test
 * Copyright (C) 2007 Christopher Wellons 
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hashtab.h"

int
main (int argc, char **argv)
{
  hashtab_t *test_ht = ht_init (2, NULL);

  /* create some data */
  char *ak = "Chris";
  char *av = "Wellons";

  char *bk = "Kelsey";
  char *bv = "Grove";

  char *ck = "Matt";
  char *cv = "Takach";

  char *dk = "Gabe";
  char *dv = "Murray";

  char *ek = "Abbie";
  char *ev = "Spinella";

  /* stick the data in the table */
  ht_insert (test_ht, ak, strlen (ak), av, strlen (av));
  ht_insert (test_ht, bk, strlen (bk), bv, strlen (bv));
  ht_insert (test_ht, ck, strlen (ak), cv, strlen (cv));
  ht_insert (test_ht, dk, strlen (ak), dv, strlen (dv));
  ht_insert (test_ht, ek, strlen (ak), ev, strlen (ev));

  /* display table data */
  hashtab_iter_t ii;
  ht_iter_init (test_ht, &ii);
  for (; ii.key != NULL; ht_iter_inc (&ii))
    {
      printf ("%s => %s\n", (char *) ii.key, (char *) ii.value);
    }

  /* grow the table */
  printf ("---\nGROW!\n---\n");
  test_ht = ht_grow (test_ht, 10);

  /* print the table contents again */
  ht_iter_init (test_ht, &ii);
  for (; ii.key != NULL; ht_iter_inc (&ii))
    {
      printf ("%s => %s\n", (char *) ii.key, (char *) ii.value);
    }

  /* free the hashtable */
  ht_destroy (test_ht);

  return EXIT_SUCCESS;
}
