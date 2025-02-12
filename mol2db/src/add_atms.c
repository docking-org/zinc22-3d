/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"
/* this adds all confs to the hierarchy */

void add_atms(int keynum, int conf, int molnum) {
   int i, j, h[MAXATOMS];
   int num;
   j = num = 0;
   reset(set_c);
   for (i=1; i<=numats; i++) /* get the atoms in the level */
      if (atlvl[i] == key[keynum][0]) set_c[++j] = numats * (conf-1) + i;  
   j = 0;
   for (i=1; i <= key[keynum][2]; i++) {  /* grab heavy atoms from key*/
      num = set_c[i];
      if (atom[num].vdw == 6 || atom[num].vdw == 7 || atom[num].vdw == 25) {
         h[++j] = num;
         continue;
      }
      hindex++;
      hier[hindex][0] = key[keynum][0];
      hier[hindex][1] = num;
      hier[hindex][2] = num - (numats * (conf-1));
      hier[hindex][3] = conf;
   }
   for (i=1; i <= j; i++) {  /* grab hydrogens from key*/
      hindex++;
      hier[hindex][0] = key[keynum][0];
      hier[hindex][1] = h[i];
      hier[hindex][2] = h[i] - (numats * (conf-1));
      hier[hindex][3] = conf;
   }
}
