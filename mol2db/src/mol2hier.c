/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

void ui(void);
void write_col_tab(void);
void read_solv(void);
int read_amber(void);
/* April 2014. 4+8 -> 4+12, thus 13 -> 21 */
void skip_mol(char mfcd[17]);

FILE* inhier;
FILE* outf;
FILE* mol2file;
FILE* filelist;
FILE* solvtable;
FILE* ambf;
char *flist, inputmol2[80], *outputdb, *solvtab, *ambtab;

Atm atom[MAXCONF * MAXATOMS];
Bnd bond[MAXATOMS];
Molecule mol[MAXCONF];
Solv solv[MAXMOLS];
hyd_tor torque[MAXATOMS];
Colortab color[MAXCOLOR];
Ambertable amb[MAXAMB]; 

int nconf, numats, numbnd, numhyd, cnum;
int ncolor, tot_atm, ens, def_col, nsolv;
int set_a[MAXATOMS], set_b[MAXATOMS];
int set_c[MAXATOMS];
int atlvl[MAXATOMS], rootlist[MAXATOMS][2];
int src, dest, knum, nwritten, nhyd, hindex;
int minx, miny, minz, ind;
int hier[MAXFATOMS][4]; /* 1: level, 2: atom number */
int upos[MAXATOMS];
int nbr;
int list[MAXFGRP][3], brid[MAXGRP];
int bconf[10], amb_cnt;
int distol;

int key[200][3];  /* 0: level, 1: lead atom, 2: atoms in level */
int keyWholeAtom[200][MAXATOMS];
int keyorder[200], keycntr;
double oconf;

bool skip, protein, torque_H, amber;
bool done[MAXCONF];

fpos_t solvptr, solvbeg;
char c;
char buffer[LINE_LENGTH], tmp[LINE_LENGTH];
char fam[5];
char def_color[31]; /*name of default color */


int main( int argc, char* argv[] )
{

   int i = 0;
   int j = 0;
   int atm_start, id, slen; 
   int dihe[4];
   char delim[] = " ";
   char gz[4];
   char gzip[80], origfile[80];
   char *tok;
   char tempchar[7];
   char inputtemp[80];
   char lenchk[80];
   bool skip_flag;        
   int  rconf0 = 0;       // index of rigid fragment of previous conformer
   int  rconf = 0;       // index of rigid fragment of current conformer
   char seps[] = " _\n";   
   char buf[LINE_LENGTH]; 
   distol = 7; //actually needs to be 7 or higher. do NOT change.
               //this is the difference allowed in the THOUSANDTHS column
               //which can be slightly different from omega 

   i = nconf = ens = 0;
   atm_start = cnum = 1;

/* Check calling parms and open files  */
   if ( argc != 2 && argc !=3 ) {
      printf("USAGE: mol2hier <INHIER>\n");
      exit(1);
   }
   if ( argc == 3) {
      if (!strcmp(argv[1], "-i")) {
         strcpy(inputtemp, argv[2]);
         if ((inhier = fopen ("inhier", "r") ) == NULL) {
            printf("Error opening parameter (%s).\n",argv[1]);
            exit(1);
         }
      }
      else {
         exit(1);
      }
   }
   else if ((inhier = fopen (argv[1], "r") ) == NULL) {
      printf("Error opening parameter (%s).\n",argv[1]); 
      exit(1);
   }
   ui();
   if (col_tab) write_col_tab();
   printf ("D  MFCD   rigid  flex   I_ats  I_confs");
   printf ("    O_ats    O_confs     O_hconfs\n");
   if (flist != NULL) { 
      if ( fgets (buffer, LINE_LENGTH, filelist) == NULL) {
         printf("Empty file list (%s).  Please play again.\n",flist);
         exit(1);
      }
      strcpy(inputmol2, buffer);
   }
   if (solv_t) {
      fgetpos(solvtable, &solvbeg);
      read_solv();
   }

/* read amber table */
   if (ambtab != NULL) {
      amb_cnt = read_amber();
   }
fileloop:
   inputmol2[strlen(inputmol2)-1] =  '\0';
   if (argc == 3) strcpy(inputmol2, inputtemp);
   strcpy(origfile, inputmol2);
   gz[0] = inputmol2[strlen(inputmol2)-3];
   gz[1] = inputmol2[strlen(inputmol2)-2];
   gz[2] = inputmol2[strlen(inputmol2)-1];
   gz[3] = inputmol2[strlen(inputmol2)];
   if (!strncmp(".gz", gz, 3)) {
      strcpy(gzip, "/usr/bin/gzip -cd ");
      strcat(gzip, inputmol2);
      strcat(gzip, " > mol2hierfile.tmp");
      printf("Decompressing: %s\n",inputmol2);
      system (gzip);
      strcpy(inputmol2, "mol2hierfile.tmp");
   }
   if ((mol2file = fopen( inputmol2, "r" ) ) == NULL) {
      printf("Error opening mol2 file %s from list.\n",inputmol2);
      exit(1);
   }
/* Read in list all atoms of all confs of a given molecule  */
   while ( fgets (buffer, LINE_LENGTH, mol2file) != NULL) {
      if (!strncmp (buffer, "@<TRIPOS>MOLECULE", 17)) {
         reset(set_a);
         reset(set_b);
         reset(set_c);
         nconf++;
         if (nconf > MAXCONF) {
            printf ("You want I should get a hernia?!  Please use");
            printf (" fewer than\n %d conformations at once\n",MAXCONF);
            exit(1);
         }
         fgets (buffer, LINE_LENGTH, mol2file); 
         strcpy(fam, substr(buffer,14,5));
/* if buffer[13] is ' ' then old format else new format */
         strcpy(mol[nconf].mfcd, substr(buffer,4,13));
// BW, to identify multiple conformations of the rigid fragment
         strcpy(buf, buffer);
         tok = strtok(buf, seps);
//         printf("ID1:%s\n", substr(tok,4,9));
         while (tok != NULL) {
               if (!strcmp(tok, "rconf")) {
                     rconf = atoi(strtok(NULL, seps));
//                     printf("%s\n", buffer);
               }
               tok = strtok(NULL, seps);
         }
         if (nconf>1) {
            if (strcmp(mol[nconf].mfcd, mol[nconf-1].mfcd) ){  // new mfcd 
               if (!skip_flag) {
                  nconf--;
                  tot_atm = atm_start - 1;
                  if (equal_q) fixcrgs();
                  if (coloring) colorit();
                  if (trans_coord) translate();
                  hiergen(protein);
                  ens++;               // ensemble counter ++
                  nconf = atm_start = 1;          // beg of new molecule
                  minx = miny = minz = 9999;
                  ind = 0;
                  strcpy(mol[nconf].mfcd, substr(buffer,4,13));
                  rconf0 = rconf;    
               }
            }
            else if ( rconf != rconf0 && rconf != 1) {  // new rconf
//         printf("previous rconf:%d ; current rconf:%d\n",rconf0, rconf);
//               printf("nconf1:%d\n ",nconf);
               if (!skip_flag) {
                  nconf--;
                  tot_atm = atm_start - 1;
                  if (equal_q) fixcrgs();
                  if (coloring) colorit();
                  if (trans_coord) translate();
                  hiergen(protein);
                  ens++;               // ensemble counter ++
                  nconf = atm_start = 1;          // beg of new molecule
                  minx = miny = minz = 9999;
                  ind = 0;
                  strcpy(mol[nconf].mfcd, substr(buffer,4,13));
               }
               rconf0 = rconf;
            }
            else {}
         } 
 
         fgets (buffer, LINE_LENGTH, mol2file);

         if (nconf == 1) {   
            numats = atoi( strtok(buffer, delim) );
            numbnd = atoi( strtok(NULL, delim) );
            if (numats > MAXATOMS) skip_mol(mol[nconf].mfcd); 
            id = atoi(substr(mol[1].mfcd,2,12));
         }
         for (i=0; i<4; i++) fgets (buffer, LINE_LENGTH, mol2file);
         i = strlen(buffer) - 1;
         strncpy (mol[nconf].name, buffer,i);
         mol[nconf].name[strlen(buffer)-1] = '\0';
         //the following loop has been changed 13/Nov/2009 by rgc
         //instead of processing mol2 files that have a @<TRIPOS>MOLECULE
         //section of only 6 lines exactly, it now works with any number 5 or 
         //higher by reading until it finds the @<TRIPOS>ATOM section. 
         //this removes the necessity to add a 6th line (since other .mol2 
         //file tools expect exactly 5 lines in the MOLECULE section
         while ( strncmp(buffer, "@<TRIPOS>ATOM", 13) ) {            
            fgets (buffer, LINE_LENGTH, mol2file);
         }
         skip_flag = FALSE;
         for (i = atm_start; i < atm_start + numats; i++) {
            fgets (buffer, LINE_LENGTH, mol2file);
            atom[i].num = atoi(strtok(buffer, delim) );
            strcpy(atom[i].atname, (strtok(NULL, delim)) );

/*          get the xyz coordinates if they're too big, set skip */
            strcpy(lenchk, strtok(NULL, delim));
            if (strlen(lenchk) < 10) {
              atom[i].x = (int)roundD2F(atof(lenchk));
            }
            else {
              skip_flag = TRUE;
            }
            strcpy(lenchk, strtok(NULL, delim));
            if (strlen(lenchk) < 10) {
              atom[i].y = (int)roundD2F(atof(lenchk));
            }
            else {
              skip_flag = TRUE;
            }
            strcpy(lenchk, strtok(NULL, delim));
            if (strlen(lenchk) < 10) {
              atom[i].z = (int)roundD2F(atof(lenchk));
            }
            else {
              skip_flag = TRUE;
            }

            if (skip_flag) break;
            if (atom[i].x < minx) minx = atom[i].x;
            if (atom[i].y < miny) miny = atom[i].y;
            if (atom[i].z < minz) minz = atom[i].z;
            strcpy(atom[i].name, (strtok(NULL, delim)) );
            atom[i].subst = atoi(strtok(NULL, delim) );
            if (protein) {
               strcpy(tempchar, strtok(NULL, delim));
               strncpy(atom[i].resname, tempchar, 3);
               slen = strlen(tempchar) - 3;
               atom[i].subst = atoi(substr(tempchar, 4,slen));
            }
            else {
               strtok(NULL, delim);
            }
            atom[i].charge = (int)roundD2F(atof(strtok(NULL, delim)));
            if (amber) {
/*               lookup_amb(atom[i].num, atom[i].atname);*/
            }
         }
         if (skip_flag) {
            nconf = 0;
            printf("coordinates too big %d, skipping\n",id);
         }
         else {
            tot_atm = i - 1;
            atm_start = i;
            if (nconf == 1) {
               fgets (buffer, LINE_LENGTH, mol2file);
               if ( strncmp(buffer, "@<TRIPOS>BOND", 13) ) {
                  printf ("Read error in bonds (mol %d). Exiting\n",nconf);
                  exit(1);
               }
               for (i=1; i<=numbnd; i++) {
                  fgets (buffer, LINE_LENGTH, mol2file);
                  bond[i].num = atoi(strtok(buffer, delim) );
                  bond[i].atom1 = atoi(strtok(NULL, delim) );
                  bond[i].atom2 = atoi(strtok(NULL, delim) );
                  tok = strtok(NULL, delim);
                  if (atoi(tok) > 0) bond[i].order = atoi(tok);
                  else if( !strncmp( tok, "am", 2)) bond[i].order = 5;
                  else bond[i].order = 4;
               }
               if (torque_H) { 
                  /* determine which H's to rotate */
                  for (i=1; i<=numats; i++) {
                     if (!strncmp("H", atom[i].atname, 1)) {
                        dihe[4] = i;
                        torque[i].hyd = term_n_s_o_h(dihe);
                        for (j=0; j<=4; j++) torque[i].dihe[j] = dihe[j];
                     } else {
                        torque[i].hyd = FALSE;
                        for (j=0; j<=4; j++) torque[i].dihe[j] = 0;
                     }
                  }
               }
            }
         }
      }
   }
   if (equal_q) fixcrgs();
   if (coloring) colorit();
   if (trans_coord) translate();
   hindex = 0;
   hiergen(protein);
   fclose(mol2file); 
   ens++;
   nconf = tot_atm = 0;
   atm_start = 1;
   printf("Finished processing %d ensembles from %s\n",ens,origfile);

   if (flist != NULL) {
      strcpy( buffer, "xtemp");
      if ( fgets (buffer, LINE_LENGTH, filelist) == NULL) exit(1);
      if (!strcmp(buffer, "xtemp") ) exit(1);
      strcpy(inputmol2, buffer);
      ens = 0;
      goto fileloop;
   }
   return 0;
}


float roundD2F(double x) {
   return x > 0 ? (float)(1000 * x + 0.5) : (float)(1000 * x - 0.5);
}

void write_col_tab(void) {
   int i, j;
   fprintf(outf,"DOCK 5.2 ligand_atoms\n");
   for (i=1; i <= def_col; i++) {
      for (j=1; j <= ncolor; j++)
         if (color[j].num == i) {
            fprintf(outf,"%-30s (%d)\n",color[j].name,i);
            break;
         }
   }
   fprintf(outf,"%-30s (%d)\n",def_color, def_col);
}

