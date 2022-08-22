/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */

#include "mol2hier.h"

/* April 2014. 4+8 -> 4+16, thus 13 -> 21 */
void skip_mol(char mfcd[17]) {
   char delim[] = " ";
   
   printf("%s has too many atoms -- discarded\n",mol[1].mfcd);
   while ( fgets (buffer, LINE_LENGTH, mol2file) != NULL) {
      if (!strncmp (buffer, "@<TRIPOS>MOLECULE", 17)) {
         fgets (buffer, LINE_LENGTH, mol2file);
         strcpy(mol[1].mfcd, substr(buffer,4,13));
         fgets (buffer, LINE_LENGTH, mol2file);
         numats = atoi( strtok(buffer, delim) );
         numbnd = atoi( strtok(NULL, delim) );
         if (numats > MAXATOMS) {
            printf("%s has too many atoms -- discarded\n",mol[1].mfcd);
         }
         else {
            return;
         }
      }
   } 
}
