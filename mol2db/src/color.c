/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

bool con_right(int src, int bndcnt, char* dest);

void colorit(void) {
   int cnum, atnum;
   int i;
   for (atnum=1; atnum<=numats; atnum++)
      for (i=0; i < nconf; i++) atom[numats *i + atnum].color = def_col;

   for (cnum=1; cnum <= ncolor; cnum++) {
/*      printf ("checking color %d %s %s %d %s\n", color[cnum].num, 
         color[cnum].name, color[cnum].a1, color[cnum].c1, color[cnum].a2);*/
      for (atnum=1; atnum<=numats; atnum++) {
         if (atom[atnum].vdw == 999) continue;
         if (strcmp(atom[atnum].name,color[cnum].a1)) continue;
         skip = FALSE;
         if (color[cnum].c1 != 0) {
            skip = !con_right(atnum, color[cnum].c1, color[cnum].a2);
         }
         if (!skip) {
            atom[atnum].color = color[cnum].num;
/*            printf ("found: %d %d\n",atom[atnum].num,atnum); */
            for (i=0; i < nconf; i++) 
               atom[numats *i + atnum].color = color[cnum].num;
         }
      }
   }
}

bool con_right(int src, int bndcnt, char* dest) {
   int i;
   set_a[1] = set_b[1] = set_c[2] = 0;
   set_c[1] = src;
   for (i=1; i<= abs(bndcnt); i++) {
      moveset(set_b, set_a);
      moveset(set_c, set_b);
      getds(set_b, set_c);
   }
   i=0;
   while (set_c[++i] != 0) {
      if (bndcnt > 0 && !strcmp(atom[set_c[i]].name, dest)) return TRUE;
      if (bndcnt < 0 && !strcmp(atom[set_c[i]].name, dest)) return FALSE;
      }
   if (bndcnt > 0) return FALSE;
   if (bndcnt < 0) return TRUE;
   return FALSE;
}
