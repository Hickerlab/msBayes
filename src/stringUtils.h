#ifndef STRING_UTILS_H
#define STRING_UTILS_H

/* Character matrix */
char **cmatrix (int nsam, int len);
void freeCMatrix (int nsam, char **list);

/*
 * Count the number of unique strings in array of string
 * with duplicated elements, the program will terminates 
 * if failed to allocate memory for objects
 *
 * Arguments:
 *  inputCharArr: array of string to be analysed
 *  returnCharArr: array contains only unique string of inputCharArr
 *               memory equivalent to inputCharArr should be pre-allocated
 *  size: number of strings in inputCharArr
 *
 * Returns: the number of strings in returnCharArr
 *
 * Note: returnCharArr NOT available from outside of the function
 */
int UniqueStrings (char **inputCharArr, char ** returnCharArr, int size);

#endif /* STRING_UTILS_H  */
