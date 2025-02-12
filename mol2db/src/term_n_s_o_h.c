/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

bool term_n_s_o_h(int dihe[4]) {
   int n1[6], n1c, i, bo, nb;
   /*  get atom connected to H called nos */
   get_neighbors(dihe[4], n1, &n1c);
   dihe[3] = n1[1];
   dihe[0] = dihe[1] = dihe[2] = 0;

   /* check to see if atom is N, O, or S */
   if (!strncmp("S", atom[dihe[3]].atname, 1) ||
       !strncmp("N", atom[dihe[3]].atname, 1) ||
       !strncmp("O", atom[dihe[3]].atname, 1)) {
      /* H is connected to nos, does nos have only two neighbors? */
      get_neighbors(dihe[3], n1, &n1c);
      if (n1c == 2) {
         /* check bond order to dihe[3] & get dihe[2] */
         i = 1;
         bo = 1;
         nb = 0;
         while (i <= numbnd && nb < 2 ) {
            if (bond[i].atom1 == dihe[3] ||
                bond[i].atom2 == dihe[3]) {
               nb++;
               if (strncmp("H", atom[bond[i].atom1].atname, 1) &&
                   strncmp("H", atom[bond[i].atom2].atname, 1)) {
                  if (bond[i].atom1 == dihe[3]) {
                     dihe[2] = bond[i].atom2;
                  } else {
                     dihe[2] = bond[i].atom1;
                  }
               }
               if (bond[i].order == 2) {
                  bo = 2; nb = 2;
               }
            }
            i++;
         }
         get_neighbors(dihe[2], n1, &n1c);
         if (n1[1] == dihe[3]) {
            dihe[1] = n1[2];
         } else {
            dihe[1] = n1[1];
         } 

         /* if it's N it must have a double bond */
         if (!strncmp("N", atom[dihe[3]].atname, 1)) {
            if (bo == 1) return FALSE;
            nb = 0;
            bo = 0;
            while (i <= numbnd && nb < 2 ) {
               if (bond[i].atom1 == dihe[2] ||
                   bond[i].atom2 == dihe[2]) {
                  nb++;
                  if (bond[i].order == 2) bo++;
               }
               i++;
            }
            /* if N has 2 double bonds, linear is hard to calc rot */
            if (bo == 2) {
               dihe[0] = 0;
               return FALSE;
            } else {
               dihe[0] = 180;
               return TRUE;
            }
         } 
         if (!strncmp("O", atom[dihe[3]].atname, 1) ||
             !strncmp("S", atom[dihe[3]].atname, 1)) {
            if (bo == 2) return FALSE;
            if (!strncmp("C.ar", atom[dihe[2]].name, 4)) {
               dihe[0] = 180;
            } else {
               if (!strncmp("C.1", atom[dihe[2]].name, 3)) {
                  dihe[0] = 0;
                  return FALSE;
               } else {
                  dihe[0] = 120;
               }
            }
            return TRUE;
         }
      }
   }
   return FALSE;
}
