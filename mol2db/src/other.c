/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"
#include <math.h>

int unique_positions(int atm) {
/* Loop over atm from each conf.  While atm is not in "seen" list
(i.e., temp[500]) add it to the "seen" list and add to npos */  
   int i, j, k, npos;
   int temp[MAXCONF];
   i = j = k = 0;
   npos = 1;
   for (i = 0; i < MAXCONF; i++) temp[i] = 0;
   temp[1] = atm;  /* first element of temp has conf #1 values */
   for (i=0; i<nconf; i++) {
      k = (i * numats) + atm;
      j = 0;
      skip = FALSE;
      while (++j <= npos) 
         if ( proximal(k,temp[j]) ) {
            skip = TRUE;
            break; //can quit now, already found match
         }
      if (!skip) {
         temp[j] = k;
         npos++;
      }
   }
   return npos;
}

//distol is actually necessary. use 2 or more. differences in thousands column
bool proximal(int a, int b) {
  if (atom[a].x==atom[b].x && atom[a].y==atom[b].y && atom[a].z==atom[b].z)
    return TRUE; //easiest case
  int diffX = atom[a].x-atom[b].x;
  int diffY = atom[a].y-atom[b].y;
  int diffZ = atom[a].z-atom[b].z;
  if (diffX*diffX+diffY*diffY+diffZ*diffZ < distol) {
    //printf("%d\n", diffX*diffX+diffY*diffY+diffZ*diffZ);    //debugging
    return TRUE;
  }
  return FALSE; 
}

void reset(int temp[MAXATOMS]) {  /* clears the set */
   int i = 0;
   for (i = 1; i < MAXATOMS; i++) temp[i] = 0;
}

/* replaces set2 with set1.  Should be done by swapping pointers */
void moveset(int set1[MAXATOMS], int set2[MAXATOMS]) {
   int i = 0;
   while (set1[++i] !=0 ) 
     set2[i] = set1[i];
   reset(set1);
}

bool is_member(int at, int temp[MAXATOMS], int natm) {
   int i;
   for (i=1; i<=natm; i++) 
     if (temp[i] == at) 
       return TRUE;
   return FALSE;
}

int gen_ctab(void) {
   int cnum, i, ncolors;
   ncolors = 0;
   for (cnum=1; cnum <= ncolor; cnum++) {
      skip = FALSE;
      for (i=1; i< cnum; i++) {
         if (!strcmp(color[cnum].name,color[i].name)) {
            skip = TRUE;
            color[cnum].num = color[i].num;
         }
      }
      if (!skip) {
         ncolors++;
         color[cnum].num = ncolors;
      }
/*      if (color[cnum].c1 == 0) printf("%d %s %s\n", color[cnum].num, 
         color[cnum].name, color[cnum].a1);
      else
      printf("%d %s %s %d %s\n", color[cnum].num, color[cnum].name,
         color[cnum].a1, color[cnum].c1, color[cnum].a2); */
   }
   def_col = ncolors+1;
   printf("%d %s (DEFAULT COLOR)\n",def_col,def_color);
   return def_col; 
}

void translate(void) {
   int i;
   for (i=1; i <= tot_atm; i++) {
      minx < 0 ? (atom[i].x -= minx) : (atom[i].x += minx); 
      miny < 0 ? (atom[i].y -= miny) : (atom[i].y += miny);
      minz < 0 ? (atom[i].z -= minz) : (atom[i].z += minz);
   }
}

//new method that puts all the atoms in a level into save
//assume save is initialized, but put a 0 at the end anyway
int lvl_cnt_save(int lvl, int save[MAXATOMS]) {
  int i, cnt;
  cnt = 0;
  for (i=1; i<=numats; i++) 
    if (atlvl[i] == lvl) {
      save[cnt] = i;
      cnt++;
    }
  save[cnt] = 0; //flag end of list
  return cnt;
}

int lvl_cnt(int lvl, int *lead) {
   int i, cnt;
   cnt = 0;
   for (i=1; i<=numats; i++) if (atlvl[i] == lvl) {
     if (cnt == 0) *lead = i;
     cnt++;
   }
   return cnt;
}

void countat(int br, int *brhvy, int *brhyd, 
             int *torh, int *cth, int brnum[MAXGRP]) {
   int i,j;
   int hv, hy, th, ct;
   hv = hy = th = 0;
   ct = 1;
   for (i=0; i<=knum; i++) if (brnum[i] == br+1) 
      for (j=1; j<=numats; j++) if (atlvl[j] == key[i][0]) {
         if (atom[j].vdw == 6 || atom[j].vdw == 7 ||
             atom[j].vdw == 25) {
            hy++;
            if (torque[j].hyd) {
               ct = ct * 360/torque[j].dihe[0];
               th = th + upos[j] - unique_positions(j);
            }
         }
         else {
            hv++;
         }
      }
   fflush(outf);
   *brhyd = hy;
   *brhvy = hv;
   *torh = th;
   *cth = ct;
   return;
}
