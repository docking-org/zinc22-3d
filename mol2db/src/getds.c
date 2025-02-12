/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

/* this gets the next downstream (ds) atoms for each atom if input array.
It puts the connected atoms in output */

void getds(int input[200], int output[200]) {
   int i, j, k, cnum;
   i = j = k = cnum = 0;
  
   reset(output);
   while (input[++i] != 0) {  /*for each input atom */
      if (input[i] == 999) continue; /*so some branches wont be followed */
      for (j=1; j<=numbnd; j++) {
         if (bond[j].atom1 == input[i] ) {
            src = bond[j].atom1;
            dest = bond[j].atom2;
         }
         else if (bond[j].atom2 == input[i] ) {
            src = bond[j].atom2;
            dest = bond[j].atom1;
         }
         else continue;
         k = 0;
         skip = FALSE;
         while (set_a[++k] != 0) if (dest == set_a[k]) skip = TRUE;
         k = 0;
         while (input[++k] != 0) if (dest == input[k]) skip = TRUE;
         k = 0;
         while (output[++k] != 0) if (dest == output[k]) skip = TRUE;
         if (!skip) {
            output[k] = dest;
            cnum++;
         }
      }
   }
}


