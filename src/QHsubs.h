#ifndef QH_SUBS_H
#define QH_SUBS_H


/* function prototypes */
extern void assign_QHsiterate (double theta, int nsites,
			       struct QHparameters *QH);
extern void createrandomseQ (char *seq, int length, struct QHparameters *QH);
extern char mutateProb (char ancbase, double ut, struct QHparameters *QH);
extern void calcPij (double Pij[4][4], double ut, struct QHparameters *QH);
extern char mutatePij (char ancbase, double Pij[4][4]);
extern double calcRate (char base, struct QHparameters *QH);
extern char mutateRate (char ancbase, struct QHparameters *QH);
extern void mutateBelow (int n, struct node *tree, int nd, char *seqs);
extern int baseint (char base);
extern char intbase (int base);
/* *gamma* functions come from tools.c in Yang's PAML package */
extern double gammadev (double s);
extern double rndgamma1 (double s);
extern double rndgamma2 (double s);

#endif /* QH_SUBS_H */
