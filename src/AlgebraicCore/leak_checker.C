//   Copyright (c)  2001-2009  John Abbott and Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include <stdio.h>
#include <stdlib.h>
// using exit

// This is pretty clearly a C program... maybe I'll rewrite it in C++ one day.

enum word_enum { NEITHER, UNMATCHED_ALLOC, UNMATCHED_FREE, MATCHED_ALLOC, MATCHED_FREE };

enum action_enum { COUNT, READ, UPDATE };

/* It is VITAL that neither word contains the first letter of the other */
/* First come the letters of the word, then zero.                       */

int ALLOC[] = { 'A', 'L', 'L', 'O', 'C', 0};
int FREED[] = { 'F', 'R', 'E', 'E', 'D', 0 };

char *type;
unsigned int *ptr;
int error_count;

enum word_enum scan_word(FILE *file, int *word)
{
  int i, ch;

  for (i=1; word[i]; ++i)
  {
    ch = getc(file);
    if ((ch == EOF) || ch != word[i])  /* match failed */
    {
      ungetc(ch, file);
      return NEITHER;
    }
  }
  if (word == ALLOC) return UNMATCHED_ALLOC;
  if (word == FREED) return UNMATCHED_FREE;
  i = i/word[i];  /* Should never get here; NB word[i] == 0 */
  return NEITHER; /* Keep the compiler quiet */
}

int scan_file(FILE *file, enum action_enum action)
{
  int n, ch;
  enum word_enum word_value;
  unsigned int ptr_val;
  fpos_t posn;

  n = -1;
  while ((ch = getc(file)) != EOF)
  {
    if (ch != ALLOC[0] && ch != FREED[0]) continue;
    if (ch == ALLOC[0]) word_value = scan_word(file, ALLOC);
    else word_value = scan_word(file, FREED);
    if (word_value == NEITHER) continue;
    ++n;
    if (action == COUNT) continue;
    if (action == READ)
    {
      getc(file); /* skip one character, it will be a space */
      const int NumMatches = fscanf(file, "%x", &ptr_val);
      if (NumMatches != 1) { fprintf(stderr, "\n***ERROR: Bad input file***\n"); exit(1); }
      type[n] = word_value;
      ptr[n] = ptr_val;
      continue;
    }
    /* action == UPDATE */
    if (type[n] == MATCHED_FREE || type[n] == MATCHED_ALLOC) continue;
    /* the next 4 lines merely put a ! in the file at the current posn */
    fgetpos(file, &posn);
    fsetpos(file, &posn);
    putc('!', file);
    ++error_count;
    fflush(file);
  }
  return n;
}

void match(int imax)
{
  for (int i = 0; i < imax; ++i)
  {
    if (type[i] != UNMATCHED_FREE) continue;
    int j;
    for (j=i; j >= 0; --j)
      if (ptr[j] == ptr[i] && type[j] == UNMATCHED_ALLOC) break;
    if (j < 0) continue;
    type[i] = MATCHED_ALLOC;
    type[j] = MATCHED_FREE;
  }
}

int main(int argc, char **argv)
{
  FILE *file;
  int n;
  if (argc != 2)
  {
    fprintf(stderr, "%s: requires 1 arg, a filename for a file containing alloc/free reports.\n", argv[0]);
    exit(1);
  }
  file = fopen(argv[1], "r+");
  n = 1 + scan_file(file, COUNT);
  printf("Found %d reports of calls to alloc or free.\n", n);
  printf("Searching for unmatched calls.\n");
  rewind(file);
  type = (char*)malloc(n*sizeof(char));
  ptr = (unsigned int*)malloc(n*sizeof(unsigned int));
  scan_file(file, READ);
  match(n);
  rewind(file);
  error_count = 0;
  scan_file(file, UPDATE);
  printf("Number of unmatched calls to alloc or free = %d\n", error_count);
  fclose(file);
  free(type);
  free(ptr);
  return 0;
}
