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

void write_charges(int keynum, int conf, int molnum) {
   int i, j, k, h[MAXATOMS];
   int num, hydlvl, crg;
   j = num = 0;
   reset(set_c);
   for (i=1; i<=numats; i++) /* get the atoms in the level */
      if (atlvl[i] == key[keynum][0]) set_c[++j] = numats * (conf-1) + i;  
   j = 0;
   for (i=1; i <= key[keynum][2]; i++) {  /* grab heavy atoms from key*/
      num = set_c[i];
      if (atom[num].vdw == 6 || atom[num].vdw == 7 || 
         atom[num].vdw == 25) {
         h[++j] = num;
         continue;
      }
      if (molnum == 0) {
         crg = atom[num].charge;
      } else {
         crg = solv[molnum].crg[num];
      }
      if (!protein) 
         fprintf (outf,"%3d%2d%5d%d%2d%9.3f%9.3f\n",key[keynum][0],
         atom[num].vdw, crg, 0, atom[num].color,
         solv[molnum].pol[num], solv[molnum].apol[num]);
      else {
         if (key[keynum][0] == 9 && !strncmp(atom[num].atname,"CB",2)) {
            fprintf(outf, "%3d%2d%5d%3d %-4s %3s\n",19,atom[num].vdw, crg,
            atom[num].subst, atom[num].atname, atom[num].resname);
         }
         if (key[keynum][0] > 9) {
            fprintf(outf, "%3d%2d%5d%3d %-4s %3s\n",
            key[keynum][0]+10, atom[num].vdw, crg, 
            atom[num].subst, atom[num].atname, atom[num].resname);
         }
      }
   }
   for (i=1; i <= j; i++) {  /* grab hydrogens from key*/
      hydlvl = key[keynum][0];
      if (molnum == 0) {
         crg = atom[h[i]].charge;
      } else {
         crg = solv[molnum].crg[h[i]];
      }

      if (torque_H) if (torque[h[i]].hyd == 1) {
         /* find atom in hierarchy, get the level */
         k=1;
         while (h[i] != hier[k][2]) k++;
         hydlvl = hier[k][0];
      }
      if (!protein) {
         fprintf (outf,"%3d%2d%5d%d%2d%9.3f%9.3f\n",hydlvl,
         atom[h[i]].vdw, crg, 0, atom[h[i]].color,
         solv[molnum].pol[h[i]],solv[molnum].apol[h[i]]);
      } else {
         if (key[keynum][0] > 9) {
            fprintf(outf, "%3d%2d%5d%3d %-4s %3s\n",
            hydlvl+10, atom[h[i]].vdw, crg,
            atom[h[i]].subst, atom[h[i]].atname, atom[h[i]].resname);
         }
      }
   }
}
