/* +++Date last modified: 05-Jul-1997 */

/**** init_globals(fp, names, types, ...)
**
** public domain by Raymond Gardner     Sept. 1991
** Ability to handle double, long, and comment lines "#" are added 
** by Naoki Takebayashi Jan. 2002
** unsigned, and long long were added by Naoki Takebayashi Jan 2006
**
** fp is a FILE * to the (already fopen'ed) file containing
**      initialization data
**  names is a space-separated list of names of globals (as they
**      are to appear in the data file)
**  types is a list of datatype characters, corresponding to names
**      i  for a pointer to integer
**      l  for a pointer to long
**      L  for a pointer to long long (int64_t)
**      u  for a pointer to unsigned integer
**      v  for a pointer to unsigned long integer
**      V  for a pointer to unsigned long long integer
**      d  for a pointer to double
**      s  for a pointer to string (already allocated char array)
**      p  for a pointer to pointer to char (space will be malloc'd)
**    (NOTE: no whitespace allowed in types !!)
**  followed by var arg list of pointers to variables to init
**
**  Make sure the very last line in the config file has the "\n" at the end.
**  Some unused variables in the config file will be ignored if variable
**  names are not listed in "names".
**  A line starting with "#" are comments and will be ignored.  But unlike
**  shell script, comments can't be put after an assignment, i.e., comment
**  lines has be in separate lines from assignment lines.
*/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "initvars.h"

int init_globals(FILE *fp, char *names, char *types, ...)
{
    char ln[LNSZ];
    char *p;
    va_list arglist;
    char *namep, *typep, name[40], *e;
    void *argp;
    int k;
    int numleft;

    numleft = strlen(types);
    while ( numleft && fgets(ln, LNSZ, fp) ) {  /* read init file */
        while ( isspace((int)ln[0]) )           /* drop leading whitespace */
            memmove(ln, ln+1, strlen(ln));
        if ( ln[0] == 0 || ln[0] == '#')        /* skip if blank line and */
	  continue;                             /* comments starting with # */
        p = strchr(ln, '=');                    /* find equal sign */
        if ( p == NULL )                        /* continue if none */
            continue;
        while ( p > ln && isspace((int)(p[-1])) ) { /* remove whitespace */
            memmove(p-1, p, strlen(p-1));           /*   before = sign */
            --p;
        }
        *p++ = 0;                               /* plug EOS over = sign */
        while ( isspace((int)(p[0])) )          /* remove leading space on */
            memmove(p, p+1, strlen(p));         /*    init string */
        k = strlen(p) - 1;                      /* init string length */
        if ( k < 1 )
            return -1;

        if ( p[k] != '\n' )             /* if '\n' is missing, input */
            return -1;                  /*   exceeded buffer; error return */
        p[k] = 0;                       /* plug EOS over newline */

        va_start(arglist, types);       /* setup for arglist search */

        namep = names;                  /* init ptr to var names */
        typep = types;                  /* init ptr to var types */
        while ( *namep == ' ' )         /* skip blanks before namelist */
            ++namep;
        while ( *typep ) {          /* while any typelist items left...*/

            argp = (void *)va_arg(arglist, void *); /* get var arg */

            k = strcspn(namep, " ");        /* length of namelist entry */
            memmove(name, namep, k);        /* put into name hold area */
            name[k] = 0;                    /* terminate it */
            if ( strcmp(name, ln) != 0 ) { /* if it doesn't match... */
                namep += k;                 /* get next name */
                while ( *namep == ' ' )
                    ++namep;
                ++typep;                    /* get next type */
            } else {                        /* else name is found... */
                if ( *typep == 'i' ) {      /* if it's an int, init it */
                    *(int *)argp = atoi(p);
		} else if ( *typep == 'd' ) { /* double */
		  *(double *)argp = strtod(p, NULL);
		  /* overflow or underflow check may be nice here */
		} else if ( *typep == 'l' ) { /* long */
		  /*                    *(long *)argp = strtol(p, NULL, 10); */
		  *(long *)argp = atol(p); 
		} else if ( *typep == 'L' ) { /* long long */
		  /* *(long long *)argp = atoll(p); no atoll in Apple*/
		  *(long long *)argp = strtoll(p, NULL, 10);
		} else if ( *typep == 'u' ) { /* unsigned int */
		  *(unsigned int *)argp = strtoul(p, NULL, 10);
		  /* I wonder if there is another function for this, Naoki */
		} else if ( *typep == 'v' ) { /* unsigned long */
		  *(unsigned long *)argp = strtoul(p, NULL, 10);
		} else if ( *typep == 'V' ) { /* unsigned long long */
		  *(unsigned long long *)argp = strtoull(p, NULL, 10);
                } else if ( *typep == 's' || *typep == 'p' ) {
                    if ( *p == '"' ) {      /* is string in quotes? */
                        ++p;                /* skip leading quote, and */
                        e = strchr(p, '"'); /* look for trailing quote */
                        if ( e )            /* terminate string if found */
                            *e = 0;
                    }
                    if ( *typep == 'p' ) {      /* if it's a char *ptr */
                        e = malloc(strlen(p) + 1); /* get space */
                        if ( e == 0 ) {         /* error if no space */
                            return -1; /* call va_end(arglist); first? */
                        }
                        *(char **)argp = e;
                        strcpy(*(char **)argp, p);  /* copy in string */
                    } else                          /* must be char array */
                        strcpy(argp, p);            /* copy in string */
                } else {
                    return -1;                      /* bad type */
                }
                numleft--;
                break;              /* break search; get next line */
            }
        }
        va_end(arglist);
    }
    return 0;
}

#ifdef TEST

/* use this test code with some config file like:
foo=123
bar = "a test"
dVar =  2.414
#  INT_MAX       2147483647 
 # LONG_MAX    9223372036854775807L (in alpha, but same as INT_MAX for others)
 
   baz             =   123456789012  
quux=       what is this
*/
#include <limits.h>
int foo;
char bar[80];
long baz;
char *quux;
double dVar;

int main(int argc, char **argv)
{
    FILE *fp;
    int k;
    unsigned long longVar;

    if ( argc < 2 ) {
        fprintf(stderr, "missing arg\n");
        exit(1);
    }

    if((fp = fopen(argv[1], "r") ) == NULL) {
      fprintf(stderr, "Can't open the file\n");
      exit(1);
    }

    longVar=INT_MAX;
    longVar=LONG_MAX;

    k = init_globals(fp, "foo bar baz quux dVar", "islpd",
        &foo, bar, &baz, &quux, &dVar);
    printf("k: %d\n", k);
    printf("foo: %d\nbar: <%s>\nbaz: %ld\nquux: <%s>\ndVar: <%f>\n",
        foo, bar, baz, quux, dVar);
    fclose(fp);

    return 0;
}

#endif /* TEST */
