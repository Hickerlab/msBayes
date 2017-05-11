#ifndef WHITE_SPACEH_H
#define WHITE_SPACEH_H
/*
 * Convert multiple white space characters to a single space.
 * Returns the number of spaces in the new converted string.
 */
int RmExtraWhiteSpaces (char *str);

/* 
 *  Remove leading white spaces from a string 
 */
char *RmLeadingSpaces (char *str);

/* 
 *  Remove trailing white spaces from a string 
 */
char *RmTrailingSpaces (char *str);

/* 
 * Find the first space character and return the pointer to it 
 * Returns NULL if no space is found
 */
char *FindFirstSpace (char *str);


/* 
 *  Check if it is blank string, Retuns 1 for empty or blank string.
 */
int BlankCharStringQ (char *str);
#endif /* WHITE_SPACEH_H */
