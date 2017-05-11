/*
 * The source codes are derived from public domain snippets,
 * and slightly modified by Naoki Takebayashi
 * License: Public Domain
 */
/* collection of functions which remove white spaces */
#include <ctype.h>		/* for isspace() */
#include <string.h>		/* for strlen() */

/* Macro to move string from s to d */
#define strMove(d,s) memmove(d,s,strlen(s)+1)

/*
 * Convert multiple white space characters to a single space.
 * Returns the number of spaces in the new converted string.
 */
int
RmExtraWhiteSpaces (char *str)
{
  char *ibuf, *obuf;
  int i, cnt;
  int wsCnt = 0;

  if (str)
    {
      ibuf = obuf = str;
      i = cnt = 0;

      while (*ibuf)
	{
	  if (isspace (*ibuf) && cnt)
	    ibuf++;
	  else
	    {
	      if (!isspace (*ibuf))
		cnt = 0;
	      else
		{
		  *ibuf = ' ';
		  cnt = 1;
		  wsCnt++;
		}
	      obuf[i++] = *ibuf++;
	    }
	}
      obuf[i] = '\0';
    }

  return (wsCnt);
}

/* 
 *  Remove leading white spaces from a string 
 */
char *
RmLeadingSpaces (char *str)
{
  char *obuf;

  if (str)
    {
      for (obuf = str; *obuf && isspace (*obuf); ++obuf)
	;
      if (str != obuf)
	strMove (str, obuf);
    }
  return str;
}

/* 
 *  Remove trailing white spaces from a string 
 */
char *
RmTrailingSpaces (char *str)
{
  if (str)
    {
      int i;

      for (i = strlen (str) - 1; i >= 0 && isspace (str[i]); --i)
	;
      str[++i] = '\0';
    }
  return str;
}

/* Find the first space character and return the pointer to it 
 * Returns NULL if no space is found
 */
char *
FindFirstSpace (char *str)
{
  char *cPtr;
  if (str)
    {
      for (cPtr = str; *cPtr && !(isspace (*cPtr)); cPtr++)
	{
	  ;
	}
      if (cPtr != str + strlen (str))
	{
	  return cPtr;
	}
    }

  return NULL;
}

/* 
 *  Check if blank char string
 *    Return value:
 *      1 blank or empty string.
 *      0 Otherwise
 */
int
BlankCharStringQ (char *str)
{
  char *cPtr;
  if (str)
    {
      if (str[0] == '\0')
	return 0;

      for (cPtr = str; *cPtr != '\0' && isspace(*cPtr); cPtr++)
	{
	  ;
	}

      if (cPtr != str + strlen (str))
	  return 0;
    }

  return 1;
}
