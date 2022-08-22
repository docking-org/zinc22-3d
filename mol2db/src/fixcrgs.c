/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

/* this routine also counts hydrogens*/
bool is_term_o(int atomnum, int *dest);

void fixcrgs(void) {
   int i, j;
   int atnum, eqnum, dest, sum, charge;
   char sybtype[2];
   int whichat[5][2]; /* atom number, charge */
   bool cpres;

  	for (i=0; i < MAXATOMS; done[++i] = FALSE);
 	for (i=0; i < 5; whichat[++i][1] = 0); /* which O's get averaged */
   i = j = dest = sum = numhyd = charge = 0;
   for (atnum=1; atnum<=numats; atnum++) {
      if (!strncmp("H",atom[atnum].name,1)) numhyd++;
      /* equalize: carboxylate, sulfate, phosphate, and nitro groups */
      if (!strcmp("O.3",atom[atnum].name)) {
         skip = TRUE;
         if (is_term_o(atnum, &dest)) continue;
         if (done[dest]) continue;
         if (!strcmp("C.2",atom[dest].name)) skip = FALSE;
         if (!strcmp("P.3",atom[dest].name)) skip = FALSE;
         if (!strcmp("S.o2",atom[dest].name)) skip = FALSE;
         if (!strcmp("S.O2",atom[dest].name)) skip = FALSE;
         if (!strcmp("N.2",atom[dest].name)) skip = FALSE;
         if (skip) continue;
         done[atom[dest].num] = TRUE;
         j = eqnum = 0;
         set_a[1] = whichat[++eqnum][1] = atom[atnum].num;
         set_a[2] = 0;
         set_b[1] = dest;
         set_b[2] = 0;
         getds(set_b, set_c);
         j=0;
         sum = whichat[eqnum][2] = atom[atnum].charge;
         while (set_c[++j] != 0) {
            if (!strcmp("O.2",atom[set_c[j]].name) ||
                !strcmp("O.3",atom[set_c[j]].name)) {
               if (is_term_o(set_c[j], &dest)) continue;
               whichat[++eqnum][1] = atom[set_c[j]].num;
               whichat[eqnum][2] = atom[set_c[j]].charge;
               sum = sum + atom[set_c[j]].charge;
            }
         }
         if (eqnum == 1) continue;
         charge = (int)sum/eqnum;
         printf("equalizing charges for %d conformation(s) of %s.\n",
            nconf,mol[1].mfcd);
         for (j=1; j <= eqnum; j++) {
            printf("  atom #%3d %-5s (%5.3f) changed to %5.3f\n", 
               whichat[j][1], atom[whichat[j][1]].name,  
               (float)whichat[j][2]/1000, (float) charge/1000);
            for (i=0; i < nconf; i++) {
               atom[(numats * (i)) + whichat[j][1] ].charge = charge;
               strcpy(atom[(numats * (i)) + whichat[j][1] ].name,"O.X");
            }
         }
      }
/* Charge guanidinium groups, (add +(1/3) to each nitrogen */
      if (!strcmp("C.2",atom[atnum].name) && !done[atnum]) {
         cpres = FALSE; /* was addchrg */
         set_a[1] = atom[atnum].num;
         set_a[2] = 0;
         getds(set_a, set_b);
         j = eqnum = 0;
         while (set_b[++j] > 0) {
            if (!strncmp("C",atom[set_b[j]].name,1)) {
               cpres = TRUE;
               set_b[j] = 999;
            }
            else if (!strncmp("N", atom[set_b[j]].name, 1)) {
               whichat[++eqnum][1] = atom[set_b[j]].num;
               whichat[eqnum][2] = atom[set_b[j]].charge;
            }
         }
         getds(set_b, set_c);
         skip = FALSE;
         j = 0;
         while (set_c[++j] != 0) {
            strcpy(sybtype, substr(atom[set_c[j]].name,3,1));
            if (!strcmp("2", sybtype) || !strcmp("a", sybtype)
               || !strcmp("p", sybtype)) skip = TRUE;
         }
         j=0;
         if (!skip && eqnum == 3) {
            done[atom[atnum].num] = TRUE;
            printf("charging guanidinium N in %d conformation(s) of %s.\n",
               nconf, mol[1].mfcd);
            for (j=1; j <= eqnum; j++) {
               printf("  atom #%3d %-5s (%5.3f) + 0.333 = %5.3f\n", 
                  whichat[j][1], atom[whichat[j][1]].name,  
                  (float) whichat[j][2]/1000, 
                  (float) ((whichat[j][2])+333)/1000);
               for (i=0; i < nconf; i++)
                  atom[(numats * (i)) + whichat[j][1] ].charge = 
                  whichat[j][2] + 333; 
            }
            continue;
         }
/* Charge amidinium groups, (add +(1/2) to each nitrogen */
         if (!skip && eqnum == 2 && cpres) {
            done[atom[atnum].num] = TRUE;
            printf("charging amidinium N in %d conformation(s) of %s.\n",
               nconf, mol[1].mfcd);
            for (j=1; j <= eqnum; j++) {
               printf("  atom #%3d %-5s (%5.3f) + 0.500 = %5.3f\n", 
                  whichat[j][1], atom[whichat[j][1]].name,  
                  (float) whichat[j][2]/1000, 
                  (float) ((whichat[j][2])+500)/1000);
               for (i=0; i < nconf; i++) 
                  atom[(numats * (i)) + whichat[j][1] ].charge = 
                  whichat[j][2] + 500; 
            }
         } 
      }
   }
}


bool is_term_o(int atomnum, int *dest) {
   int j, nlinks;
   nlinks = 0;
   for (j=1; j<=numbnd; j++) {
      if (bond[j].atom1  == atom[atomnum].num) {
         *dest = bond[j].atom2;
         nlinks++;
      }
      if (bond[j].atom2  == atom[atomnum].num) {
         *dest = bond[j].atom1;
         nlinks++;
      }
   }
   return nlinks > 1 ? TRUE : FALSE;
} 
