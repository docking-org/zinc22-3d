/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"
/* this writes out the hierarchy array and rotates hydrogens */

void rotateH(float* A,float* B,float* C,float* D, float ang);

void write_hierarchy(int br) {
   int i, k, j, at1, at2, conf, pr;
   float di_a[3], di_b[3], di_c[3], di_d[3], di_st[3], angleTo;
/* pr is a protein offset.  for proteins the hierarchy should be numbered
   starting with 29 to make room for CB (the new 19)
*/
   if (protein) {
      pr = 10; 
   } else {
      pr = 0;
   }
   if (br==0) brid[0] = 1;
      for (i=brid[br]; i <= brid[br+1]-1; i++) {
         if (protein && hier[i][0] == 9 &&
             strncmp(atom[hier[i][1]].atname,"CB",2)) continue; 
         if (torque[hier[i][2]].hyd == 1) {
            conf = hier[i][3];
            /* get the four atoms, get conformation, get coords */
            /* at2 = atom in first mol (has dihe vals associated) */
            /* torque[at2].dihe[1]  atoms for dihe */
            /* at1  moves to correct instance of atom */
/*
            hier[i][1]  <-- the coords I want
            hier[i][2]  <-- the atom id (for getting dihe)
            atom[tba].x <-- original data
*/
            at2 = hier[i][2];
            if (conf == 1) {
               at1 = 0;
            } else {
               at1 = numats * (conf - 1);
            }
            di_a[0] = (float)atom[at1+torque[at2].dihe[1]].x/1000;
            di_a[1] = (float)atom[at1+torque[at2].dihe[1]].y/1000;
            di_a[2] = (float)atom[at1+torque[at2].dihe[1]].z/1000;
            di_b[0] = (float)atom[at1+torque[at2].dihe[2]].x/1000;
            di_b[1] = (float)atom[at1+torque[at2].dihe[2]].y/1000;
            di_b[2] = (float)atom[at1+torque[at2].dihe[2]].z/1000;
            di_c[0] = (float)atom[at1+torque[at2].dihe[3]].x/1000;
            di_c[1] = (float)atom[at1+torque[at2].dihe[3]].y/1000;
            di_c[2] = (float)atom[at1+torque[at2].dihe[3]].z/1000;
            di_d[0] = (float)atom[at1+torque[at2].dihe[4]].x/1000;
            di_d[1] = (float)atom[at1+torque[at2].dihe[4]].y/1000;
            di_d[2] = (float)atom[at1+torque[at2].dihe[4]].z/1000; 
            for (j=0; j < 360; j=j+5) {
               angleTo = (float)j;

               /* store original H position */
               if (j>0) for (k=0; k <= 2; k++) di_st[k] = di_d[k];

               if (j>0) rotateH(di_a, di_b, di_c, di_d, angleTo);
               if (hier_sp) {
                  for(k=1; k<=tens(hier[i][0]); k++) fprintf(outf," ");
                  fprintf (outf,"%3d%6d%6d%6d\n", hier[i][0]+pr,
                  (int)(di_d[0]*1000), (int)(di_d[1]*1000),
                  (int)(di_d[2]*1000));
               } else {
                  fprintf (outf,"%3d%6d%6d%6d\n", hier[i][0]+pr,
                  (int)(di_d[0]*1000), (int)(di_d[1]*1000),
                  (int)(di_d[2]*1000));
               }

               /* restore original H position */
               if (j>0) for (k=0; k <= 2; k++) di_d[k] = di_st[k];
            }
         } else {
            if (hier_sp) {
               for(k=1; k<=tens(hier[i][0]); k++) fprintf(outf," ");
               fprintf (outf,"%3d%6d%6d%6d\n", hier[i][0]+pr,
               atom[hier[i][1]].x, atom[hier[i][1]].y,
               atom[hier[i][1]].z);
            } else {
               fprintf (outf,"%3d%6d%6d%6d\n", hier[i][0]+pr,
               atom[hier[i][1]].x, atom[hier[i][1]].y,
               atom[hier[i][1]].z);
            }
         }
      }
} 
