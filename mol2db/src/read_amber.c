/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

int read_amber(void) {
   char delim[] = " ";
   int ct;
   ct = 0;
   while ( fgets (buffer, LINE_LENGTH, ambf) != NULL) {
      if ( strncmp("!", buffer,1)) {
         strcpy (amb[++ct].atom, substr(buffer,1,4));
         strcpy (amb[ct].atom, strtok(amb[ct].atom, delim));
         strcpy (amb[ct].res, substr(buffer,8,3));
         amb[ct].crg = (float)atof(substr(buffer,18,6));
         amb[ct].type = atoi(substr(buffer,25,2));
      }
   }
   return ct;
}
