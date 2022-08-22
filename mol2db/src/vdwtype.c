/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

bool carbon_linked(int hnum);

 /* Adapted from fixid.f and fixhc.f from DOCK3.5 mol2db */
int vdwtype(int num, char name[] ) {  
   if ( !strcmp("H", name)) {
      if (carbon_linked(num)) return 7;
      return 6;
   }
   if ( !strncmp("C.", name, 2)) {
      if ( !strncmp("C.3", name, 3)) return 5;
      return 1;
   }
   if (!strncmp("N.", name, 2)) {
      if ( !strncmp("N.4", name, 3)) return 9;
      if ( !strncmp("N.3", name, 3)) return 10;
      return 8;
   }
   if ( !strncmp("O.", name, 2)) {
      if ( !strncmp("O.3", name, 3)) return 12;
      return 11;
   }
   if ( !strcmp("F", name)) return 15;
   if ( !strcmp("I", name)) return 18;
   if ( !strncmp("S.", name, 2)) return 14;
   if ( !strncmp("P.", name, 2)) return 13;
   if ( !strncmp("P2", name, 2)) return 13; /* tmp fix, remove */
   if ( !strncmp("Cl", name, 2)) return 16;
   if ( !strncmp("Br", name, 2)) return 17;
   if ( !strncmp("Si", name, 2)) return 24;
   if ( !strncmp("Ca", name, 2)) return 21;
   if ( !strncmp("Na", name, 2) || !strcmp("K", name)) return 19;
   if ( !strncmp("Li", name, 2) || !strncmp("Al", name, 2)) return 20;
   if ( !strcmp("B", name)) return 20;
   if ( !strncmp("Du", name, 2) || !strncmp("LP", name, 2)) return 25;
   printf("id number not found for atom %s\n",name);
   printf("Program stops in vdwtype.\n");
   exit(1);
   return 0;
}

bool carbon_linked(int hnum) {
   int j;
   for (j=1; j<=numbnd; j++) {
      if (bond[j].atom1 == hnum) 
         if (!strncmp("C", atom[bond[j].atom2].name,1)) return TRUE;
      if (bond[j].atom2 == hnum) 
         if (!strncmp("C", atom[bond[j].atom1].name,1)) return TRUE;
   }
   return FALSE;
}
