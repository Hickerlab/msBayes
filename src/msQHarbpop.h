#ifndef MS_QH_ARBPOP_H
#define MS_QH_ARBPOP_H
#include <stdlib.h>
#include <stdio.h>
#include <float.h>



/* msQH.c:main(),*, QHsubs.c:*   */
struct QHparameters
{
  int output;			/* dumb internal flag for outputing sequence versus polymorphic sites */
  double Rtstv, Ytstv, freqA, freqC, freqG, freqT, obsfreq[4];
  double alpha, siterate;	/* if alpha (gamma parameter) ==0, then no rate heterogeneity */
  char *MRCAseq, *treeseq;
  //  int *segbool; /* for QHPms.c */
};
struct Dintervalparams
{
  double Dtpast;
  int npops, Mpattern;
  double *Nrecent, *Npast;
  double **M, **a;
};
struct parameters
{
  int nsam, segsites, nsites, *config;
  int D, Dintn, seQmut;
  long howmany;
  double theta, r, f, track_len;
  struct Dintervalparams *Dint;
  struct QHparameters *QH;
};


/* in msQH.c:gensam(),makegametes(), streec.c:streec(),*   */
struct node
{
  int abv;
  int ndes;
  float time;
};

struct segl
{
  int beg;
  struct node *ptree;
  int next;
};


int pick2 (int n, int *i, int *j);

#endif /* MS_QH_ARBPOP_H */
