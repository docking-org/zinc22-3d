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

void write_atoms(int lvl, int conf);
void countat(int br, int *brhvy, int *brhyd, 
             int *torh, int *cth, int brnum[MAXGRP]);
char *substr( char *line, int start, int length);
int create_index(int bend[50], int bbeg[50]);

int tens(int x);
int lvl2k(int lvl); 


int conf_hier(int tpos) {
   int bend[50], bbeg[50]; //index into the keys for each branch 
                           //       (9=rigid, 19=flex1, 18=flex2, 17etc=etc)
                           //bbeg goes to first (beginning)
                           //bend goes to last (end)
   int conf, atm, k, i, nlist, curr, next;
   int j, crr, nxt, lead, max, tconf, molnum, confrm;
   int next_list, mark_lvl, list_atm;
   int br, brnum[50], brhvy, brhyd, torh, cth;
   hindex = conf = atm = k = i = nlist = curr = next = 0;
   for (i = 1; i < 200 ; i++) keyorder[i] = 0;
   //rootlist[index][0] is either 0 or the root of a level's atom index
   //  this is sometimes incorrect, but is fixed when converting to keys later
   //rootlist[index][1] is either 0 or the output name (9, 19, 18, etc) of level
   if (rootlist[1][0] == 0) { //i don't think this loop is ever run
      rootlist[1][0] = 1;     //unless these are not set up somehow
      rootlist[1][1] = 9;
   } 
   //next 2 lines debug the rootlist if necessary   
   //i=0;
   //while (rootlist[++i][0] != 0) printf("rootlist[%d][0] = %d [1] = %d\n", i, rootlist[i][0], rootlist[i][1]);
   i=0;
   while (rootlist[++i][0] != 0) keyorder[i] = i;
   /* convert back to keys so I don't have to rework everything downstream.
     Whoever designed this using the keys in numerical order was an idiot.  This
     should be rewritten so as to not use keys, then we don't have to do the
     conversion k=keyorder[keycntr] all over.  DL  */
   //key conversion is necessary because the rootlist[x][0] is sometimes wrong
   knum = i - 1; //knum is the number of groups?
   for (i=1; i<=knum; i++) {
      key[i][0] = rootlist[i][1]; //still the output group id
      key[i][2] = lvl_cnt(key[i][0], &lead);//lvl_cnt is the number in the level
      key[i][1] = lead;//lead is set in the previous method call, an atom id
      //next line debugs what is going on here
      //printf("%d %d %d\n", key[i][0], key[i][1], key[i][2]);
      int savelvllength = lvl_cnt_save(key[i][0], keyWholeAtom[i]);
      //keyWholeAtom is indexed the same as key but has all the atom info
      //each one ends in a 0 when there are no more atoms in that level
   }
   /* add first molecule */
   for (keycntr = 1; keycntr <= knum; keycntr++) {
      k = keyorder[keycntr];
      //printf("%d %d\n", k, keycntr); //confirm keys are useless 1=1 2=2 etc
      nlist++;
      list[nlist][0] = key[k][0]; //the output group number (9, 19, 18, 29, etc)
      list[nlist][1] = 1; //[1] is the number of conformations here?
      list[nlist][2] = nlist + 1; //links to elsewhere in the tree
      //next line debugs what is going on here
      //printf("%d %d %d %d\n", nlist, list[nlist][0], list[nlist][1], list[nlist][2]);
   }
   keyorder[knum+1] = 0;
   list[nlist][2] = 0; //links nowhere since it is the end
   /* integrate rest of confs into hierarchy */
   for (conf = 2; conf <= nconf; conf++) {
      mark_lvl = 9; //idea is to set this to the level above the one being added
      i = 0;
      keycntr = next = 1; //start at 1 since 0 is rigid component
      while (++keycntr <= knum) {
         k = keyorder[keycntr]; //k = keycntr
         while (next != 0) {
            curr = next;
            next = list[curr][2];
            next_list = list[next][0];
            //printf("%d\n", atm); //confirming not every atm is compared
            //printf("curr list = %d %d\n",list[curr][0], atom[atm].x);  
            if (list[curr][0] == key[k][0]) { //found the place to add it??
               bool allMatch = TRUE;
               //keyWholeAtom[k] is the complete key list
               int currentKWAindex = 0;
               while (allMatch && keyWholeAtom[k][currentKWAindex] != 0) {
                 //atm = (conf-1) * numats + key[k][1]; //new atom id to add
                 //list_atm = (list[curr][1]-1) * numats + key[k][1];
                 int currentAtom = keyWholeAtom[k][currentKWAindex]; 
                 atm = (conf-1) * numats + currentAtom;
                 list_atm = (list[curr][1]-1) * numats + currentAtom;
                 if (!proximal(atm,list_atm)) //if atoms overlap
                   allMatch = FALSE;
                 currentKWAindex++; //go to next atom
               }
               //instead of just checking atm and list_atm, want to check
               //each keyWholeAtom[k] and make sure they are ALL the same
               //putting this off until i really understand what happens here
               if (allMatch) { //if atoms overlap
                  k=keyorder[++keycntr];     //advance the key (down the tree)
                  if (keycntr > knum)        //if we're done break out
                     break; 
                  mark_lvl = key[k][0];      //otherwise continue
                  continue;
               }
            }
            if (key[k][0] < mark_lvl) {
               mark_lvl = 9; //found an earlier element, move back up the list
            } 
            //printf("check %d %d %d\n",next_list,mark_lvl,key[k][0]); 
            if ( next_list < mark_lvl ) { 
               //printf("\nbreak1 %d %d %d\n",next_list,mark_lvl,key[k][0]);
               break;
            }  
            if ( tens(next_list) == tens(key[k][0]) && 
                 next_list < key[k][0] ) {
                 //printf("\nbreak2 %d %d\n",next_list,key[k][0]);
               break;
            }  
            mark_lvl = key[k][0];

         }
         if (keycntr > knum)  
            break;
         //printf ("%d Adding %d m=%d, %d, %d\n",conf,key[k][0],mark_lvl,list[curr][2], nlist);
         nlist++;
         list[nlist][0] = key[k][0];
         list[nlist][1] = conf;
         list[nlist][2] = list[curr][2];
         list[curr][2] = nlist;
         next = curr;
      }

   }
   //prints out the final list[][] structure
   //for (i = 0; i <= nlist; i++) 
   //   printf("%d %d %d %d\n", i, list[i][0], list[i][1], list[i][2]);
   
   /* header info */
   nhyd = 0;
   for (k=1; k <= numats; k++)
      if (atom[k].vdw == 6 || atom[k].vdw == 7 || atom[k].vdw == 25) nhyd++;
   k = 0;

   /*  lookup solvation data */
   i = molnum = 0;
   while (++i <= nsolv ) {
      if (!strcmp( mol[1].mfcd, solv[i].mfcd)) {
         molnum = i;
         i = nsolv + 1;
      }
   }
   (molnum == 0) ? (confrm = -3) : (confrm = 1);
   if (protein) confrm = 1;
   /* write hierarchy to hier array */
   nxt = 1;
   br = 0;
   k = 20;
   for (j=1; j <= nlist; j++) {
      crr = nxt;
      nxt = list[crr][2];
      if (tens(list[crr][0]) == 1) {
         if (list[crr][0] < k) {
            k = list[crr][0];
            br++;
            brid[br] = hindex + 1; //i think this is the hierarchy index
         }
      }
      //printf("crr=%d, level=%d, conf=%d\n",crr,list[crr][0],list[crr][1]); 
      add_atms(lvl2k(list[crr][0]), list[crr][1], molnum); 
   }
   brid[++br] = hindex+1;
   /* this routine assigns a branch number to each group */
   br = 1;
   bbeg[1] = 1;
   for (i=1; i<=knum; i++) {
      if (tens(list[i][0]) == 1) {
         br++;
         bend[br-1] = i - 1;
         bbeg[br] = i;
      }
      brnum[i] = br;
   }
   br++;
   bend[br-1] = i - 1;
   bbeg[br] = i;
   //for (i=1; i<=br; i++)   //these lines debug the bbeg and bend data
   //  printf("%d %d %d\n", i, bbeg[i], bend[i]);
   /* this establishes pointer order for extracting */
   max = 0;
   for (i=1; i<=knum; i++) 
     if (key[i][0] > max) 
       max = key[i][0];
   max = tens(max);
   nxt = 1;
   tconf = create_index(bend, bbeg); 
   if (tconf == -1) {
      printf("Macrocycle error on %s\n", mol[1].mfcd);
      exit(11);
   }
   if (tconf == -2) 
      return -2;
   //printf("tconf:%d nconf:%d\n", tconf, nconf); //debugging
   if (tconf < nconf) 
      return -3;

   /* generate output for rigid fragment*/
   if (!protein) { 
       countat(0, &brhvy, &brhyd, &torh, &cth, brnum);
       fprintf (outf,"Family  %5d %3d %3d %3d\n",1,1,nbr,nbr);
/* April 2014, 4+8 -> 4+16 , thus 9s -> 17s */
       fprintf (outf,"%-47.47s%13s ",mol[1].name,mol[1].mfcd);
       for (i=2; i<=nbr+1; i++) fprintf (outf,"%2d",i);
       fprintf (outf,"\n");
       fprintf (outf,"%7d%7d%6d%6d%12.4f%6d%9.2f %4d %10d\n", brid[1]-1,
                brid[1]-1, brhvy, brhyd, solv[molnum].pol[0], confrm, 
                solv[molnum].apol[0],ind, 1);
       nwritten = brid[1]-1;
       write_charges(1, 1, molnum);
       write_hierarchy(0);
   } else {
      if (!strncmp(atom[1].resname,"PRO",3)) {
         fprintf (outf,"%7d%7d%6d%6d%12.4f%6d%9.2f %4d %10d\n", 3,
         3, 3, 0, solv[molnum].pol[0], confrm, 
         solv[molnum].apol[0],ind, 1);
         for (k=5; k<=7; k++) fprintf(outf, 
            "%3d%2d%5d%3d %-4s %3s\n",19,atom[k].vdw,atom[k].charge,
            atom[k].subst, atom[k].atname, atom[k].resname);
         for (k=5; k<=7; k++) fprintf (outf,"%3d%6d%6d%6d\n",
            19, atom[hier[k][1]].x, atom[hier[k][1]].y,
            atom[hier[k][1]].z);
      }
      else if (!strncmp(atom[1].resname,"ALA",3)) {
         fprintf (outf,"%7d%7d%6d%6d%12.4f%6d%9.2f %4d %10d\n", 1,
         1, 1, 0, solv[molnum].pol[0], confrm,
         solv[molnum].apol[0],ind, 1);
         fprintf(outf, "%3d%2d%5d%3d %-4s %3s\n",19,atom[5].vdw, 
         atom[5].charge, atom[5].subst, atom[5].atname, atom[5].resname);
         fprintf (outf,"%3d%6d%6d%6d\n",19, atom[hier[5][1]].x,
         atom[hier[5][1]].y, atom[hier[5][1]].z);
      }
      else if (!strncmp(atom[1].resname,"GLY",3)) {
         fprintf (outf,"%7d%7d%6d%6d%12.4f%6d%9.2f %4d %10d\n", 0,
         0, 0, 0, solv[molnum].pol[0], confrm,
         solv[molnum].apol[0],ind, 1);
      }
      else {}
   }

   /* generate output for flexible branches*/
   oconf = 1;
   for (br=1; br<=nbr; br++) {       
      countat(br, &brhvy, &brhyd, &torh, &cth, brnum); 
      if (bconf[br+1] == 0) bconf[br+1] = 1;
      if (protein) {
         fprintf (outf,"%7d%7d%6d%6d%12.4f%6d%9.2f %4d %10d\n",
                  brid[br+1] - brid[br] + torh+1, brhvy + brhyd+1, brhvy+1,
                  brhyd, solv[molnum].pol[0], confrm,
                  solv[molnum].apol[0],ind, bconf[br+1]*cth);
         write_charges(1, 1, molnum);
      } else {
         fprintf (outf,"%7d%7d%6d%6d%12.4f%6d%9.2f %4d %10d\n",  
                  brid[br+1] - brid[br] + torh, brhvy + brhyd, brhvy,
                  brhyd, solv[molnum].pol[0], confrm, solv[molnum].apol[0],ind, 
                  bconf[br+1]*cth); 
      }
      nwritten = nwritten + brid[br+1] - brid[br] + torh;
      oconf = oconf * (bconf[br+1]*cth);
      for (i=1; i<=max+1; i++) for (j=9; j>=0; j--) 
           for (k=1; k<=knum; k++) {
               if ((i*10)+j == key[k][0]) { 
                  if (brnum[k]==br+1) write_charges(k, 1, molnum);
               }
            } 
      if (br == 1 && protein) write_hierarchy(0);
      write_hierarchy(br);
   } 

/* April 2014. 4+8 -> 4+16, thus 8s -> 16s */
   printf("%12s %4d %5d %7d     %5d %7d   %9d %14.0f\n", mol[1].mfcd,
     key[1][2], numats - key[1][2], key[1][2]+(numats - key[1][2]) * nconf,
     nconf, nwritten, tconf, oconf);
   fflush(NULL); 
}

int tens(int x) {
   return x = (int)(x/10);
}

int lvl2k(int lvl) {
   int i = 0;
   for (i=1; i<=knum; i++) if (key[i][0] == lvl) return i;
   return 0;
} 

char *substr( char *line, int start, int length) {
   char *retstr;
   int i;

   i = strlen(line) + 1;
   if ( i < start + length) 
     return "";
   retstr = (char*) malloc( length * sizeof( char ) );

   for( i=0; i<length; i++ ) 
      retstr[i]=line[start+i-1];
   retstr[i] = '\0';   
   return retstr; 
}

