/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

int w[200];  /* i believe this is a pointer into the atom list, updated on
                the fly to reflect updating for each moving branch */
int stream_order[100][2];
int p_list[100][3]; //no idea what this is, rename later

void update_mol(int ptr);
int base(int ptr, int br, int bbeg[50], int bend[50]);

int create_index(int bend[50], int bbeg[50]) {
   int i, j, k, ndone, atm;
   int br, tconf; //br is the branch id (0 for rigid, 1,2,etc for flexible)
   ndone = j = 0;
   atm = 1;
   for (i=0; i < MAXCONF; done[i++] = FALSE);
   for (i=0; i < 100; i++) stream_order[i][0] = stream_order[i][1] = 0;
   for (i=0; i < 10; i++) bconf[i] = 0;
  
   /* set initial pointers, no idea what they do */
   p_list[knum-1][2] = 0;  //i believe p_list[x][2] is a list of indices for 
                      //each branch start point
   key[0][2] = 0;
   for (k=1; k <= knum; k++) 
     p_list[k][0] = p_list[k][1] = -1; //initialize now
     p_list[k][2] = p_list[k-1][2] + key[k-1][2]; //initialize to last plus next group
   if (key[2][0] == 0) { //only ever executed when just rigid
       key[1][1] = 1;
       key[1][2] = numats;
   } 
   //this sets p_list[x][0-1], still not sure how used later
   //also sets the stream order which seems to always be in order
   //but is the reason for the extremely strange loop logic
   while (ndone < knum) {
      //printf ("%d %d %d\n",atm, hier[atm][0], hier[atm-1][0]);
      if (hier[atm][0] != hier[atm-1][0] || atm == 1) { 
         for (i=knum; i >= 1; i--) 
            if (hier[atm][0] == key[i][0] && !done[i]) {
              done[i] = TRUE;
              //printf ("key %d of %d done level %d\n",i,knum,hier[atm][0]);
              ndone++;
              p_list[i][0] = atm;
              p_list[i][1] = atm - 1 + key[i][2];
              j++;
              stream_order[j][0] = key[i][0];
              stream_order[j][1] = i;
              //printf("i:%d j:%d key[i][0]:%d\n", i, j, key[i][0]);
              break;
            }
      }
      atm++;
   }
   p_list[knum+1][0] = hindex + 1;
   oconf = 1;
   if (nconf == 1) {
      ind = 1;
      return 1; //only 1 conformation--rigid case?
   } 
   if (key[1][2] == numats) 
      return 1; //only 1 conformation--rigid case?
   /* find position of keys in input data */
   i = 0;
   j = 1;
   for (k=1; k <= knum; k++) { 
      j = j + key[k-1][2];
      for (atm = 1; atm <= numats; atm++) { 
         if (proximal(atm, hier[w[j]][1])) {
            set_a[++i] = atm;
            break;
         }
      }
   }
   /* extract molecules (1/2 ark style -- 1 by 1)  */
   tconf = 1;
   for (br = 2; br <= nbr+1; br++) { //iterate through each flexible branch
      i = bend[br]; 
      bconf[br] = 1; //initialize branch conformations to 1
      while (i > bend[br-1]) {
         //printf("i bend[br-1]: %d %d\n", i, bend[br-1]);
         //printf("i hier hier: %d %d %d\n", i,hier[p_list[i][1]+1][0], hier[p_list[i][0]][0]);
         if (hier[p_list[i][1]+1][0] == hier[p_list[i][0]][0]) {
            j = i;
            //printf("j p_list p_list: %d %d %d\n", j, p_list[j][1] ,p_list[bbeg[br+1]][0]);
            if (p_list[j][1] >= p_list[bbeg[br+1]][0])
               break;
            while (j <= bend[br]) {
               if (stream_order[j][1] != stream_order[i][1]) {
                  if (stream_order[j][1] < stream_order[i][1])  {
                     p_list[j][1] = base(j,br,bbeg,bend) - 1; //i don't think this is ever executed
                  }
                  else {
                     p_list[j][1] = base(j,br,bbeg,bend);
                  }
                  if ((p_list[j][1] == -999999999) || 
                                 (p_list[j][1] + 1 == -999999999))
                     return -1; //this should really be an exception
                  while (hier[p_list[j][1]+1][0] != hier[p_list[j][0]][0] &&
                                   p_list[j][1] < p_list[bbeg[br+1]][0]) 
                     p_list[j][1]++;
                  p_list[j][0] = p_list[j][1] + 1;
               }
               if (p_list[j][1] >= p_list[bbeg[br+1]][0]) 
                  break;
               update_mol(j); 
               j++;
            }
            bconf[br]++;
            //printf("%d %d \n", br, bconf[br]);
            i = bend[br];
         }
         else if (hier[p_list[i][1]+1][0] < hier[p_list[i][0]][0]) {
            i--;
         }
         else {
            while (hier[p_list[i][1]+1][0] > hier[p_list[i][0]][0]) {
               p_list[i][1]++;
               if(p_list[i][1] > hindex) 
                  return 0;
            }
         }
      }
      tconf = tconf * bconf[br];
      //printf ("tconf = %d %d %d\n",tconf, br, bconf[br]);
   }
   return tconf;
}

//i think this moves to the next conformation, updating w from p
//w & p are the worst names ever
void update_mol(int ptr) {
   int i;
   //printf("update_mol %d ", ptr);
   for (i=1; i <= key[ptr][2]; i++) {
      p_list[ptr][1]++;
      w[p_list[ptr][2]+i] = p_list[ptr][1];
      //printf("%d,%d ", p_list[ptr][2]+i, w[p_list[ptr][2]+i]);
   }
   //printf("\n");
}

//i think this function finds the first atom in the previous branch 
//as a way to walk back up the tree, i think the error is when it has
//walked all the way back up. there was a bug because it used 9999
//which sometimes the return could reasonably be (since it is a pointer into
//an array of all the atom coordinates. now it uses a still-bad but not
//possible to reach under normal circumstances value of -99999999
int base(int ptr, int br, int bbeg[50], int bend[50]) {
   int i;
   i = bbeg[br];
   while (stream_order[i][1] != ptr && i<=bend[br]) 
     i++;
   if (tens(stream_order[i][0]) == 0) 
     return -999999999; //new 'flag'. really should be an exception
   while (tens(stream_order[i][0]) >= tens(hier[p_list[ptr][0]][0])) 
     i--;
   return p_list[stream_order[i][1]][1];
}

