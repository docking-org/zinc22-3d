/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

void read_solv(void) {
   char delim[] = " ";
   int pnum, ct;
   ct = 0;
   fgetpos(solvtable, &solvptr);
   while ( fgets (buffer, LINE_LENGTH, solvtable) != NULL) {
      pnum = 0;
/* April 2014. change from 4+8 to 4+16 , thus 4,9  to 4,13 */
/* if buffer[13] is ' ' then we are reading old format else new format */
      strcpy (solv[++ct].mfcd, substr(strtok(buffer, delim),4,13));
      solv[ct].cnt = atoi(strtok(NULL, delim));
      solv[ct].crg[pnum] = (float)atof(strtok(NULL, delim));
      solv[ct].pol[pnum] = (float)atof(strtok(NULL, delim));
      solv[ct].surf[pnum] = (float)atof(strtok(NULL, delim));
      solv[ct].apol[pnum] = (float)atof(strtok(NULL, delim));
      solv[ct].tot[pnum] = (float)atof(strtok(NULL, delim));
      for (pnum=1; pnum <= solv[ct].cnt; pnum++) {
         fgets (buffer, LINE_LENGTH, solvtable);
         solv[ct].crg[pnum] = (int)roundD2F(atof(strtok(buffer, delim)));
         solv[ct].pol[pnum] = (float)atof(strtok(NULL, delim));
         solv[ct].surf[pnum] = (float)atof(strtok(NULL, delim));
         solv[ct].apol[pnum] = (float)atof(strtok(NULL, delim));
         solv[ct].tot[pnum] = (float)atof(strtok(NULL, delim));
      }
   }
   nsolv = ct;
   return;
}
