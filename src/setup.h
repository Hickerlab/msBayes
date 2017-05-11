#ifndef SIM_SETUP_H
#define SIM_SETUP_H

#include <stdio.h>

#include "msprior.h"

extern void LoadConfiguration (int argc, char *argv[]);
extern void PrintParam (FILE *fp);

extern int setup_done;
extern int debug_level;

#endif /* SIM_SETUP_H */
