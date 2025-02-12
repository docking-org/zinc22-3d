/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

int lookup_amb(char *res, char *at) {
   int i;
   for (i=1; i<=amb_cnt; i++) {
      if (!strncmp(amb[i].res, res, 3)) {
         if (!strncmp(amb[i].atom, at, 3)) return amb[i].type;
      }
   }
   printf ("Unable to find residue %3s, atom %s in table\n",res, at);
   return 0;
}
