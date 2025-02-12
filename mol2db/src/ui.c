/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

bool col_tab, coloring, solv_t, equal_q, hier_sp, trans_coord;

int gen_ctab(void);

bool tf(char *yn);

void ui(void) {
   /* int def_col;   integer value of default color */
   int b = 0;
   bool openfile;
   char delim[] = " ";
   char *field, *value;
   openfile = col_tab = coloring = solv_t = equal_q = hier_sp = FALSE;
   amber = protein = trans_coord = torque_H = FALSE;
   while ( fgets (buffer, LINE_LENGTH, inhier) != NULL) {
      field = strtok(buffer, delim);
      if (!strncmp("#", field, 1)) continue;
      if (field[1] == '\0') continue;
      value = strtok(NULL, delim);
      if (!strcmp("protein",field ))  {
         if (tf(value)) protein = TRUE;
      }
      else if (!strcmp("equalize_charges",field ))  {
         if (tf(value)) equal_q = TRUE;
      }
      else if (!strcmp("solvation_correction",field )) {
         if (tf(value)) solv_t = TRUE;
      }
      else if (!strcmp("color_atoms",field )) {
         if (tf(value)) coloring = TRUE;
      }
      else if (!strncmp("output_color_table",field,18 )) {
         if (tf(value)) col_tab = TRUE;
      }
      else if (!strcmp("hierarchy_spacing",field )) {
         if (tf(value)) hier_sp = TRUE;
      }
      else if (!strcmp("translate_coordinates",field )) {
         if (tf(value)) trans_coord = TRUE;
      }
      else if (!strcmp("torque_hydrogens",field )) {
         if (tf(value)) torque_H = TRUE;
      }
      else if (!strcmp("mol2_file",field)) {
         if (openfile == TRUE) {
            printf("You have selected both mol2_file and mol2_file_list.");
            printf("  One to a customer please.\n");
            exit(1);
         }
         openfile = TRUE;
         strcpy (inputmol2, value);
         continue;
      }
      else if (!strncmp("mol2_file_list", field, 14)) {
         if (openfile == TRUE) {
            printf("You have selected both mol2_file and mol2_file_list.");
            printf("  One to a customer please.\n");
            exit(1);
         }
         openfile = TRUE;
         flist = value;
         flist[strlen(flist)-1] = '\0';
         if ((filelist = fopen( flist, "r" ) ) == NULL) {
            printf("Error opening list of mol2 files (%s).\n",flist);
            exit(1);
         }
      }
      else if (!strcmp("db_file",field )) {
         outputdb = value;
         outputdb[strlen(outputdb)-1] = '\0';
         if ((outf = fopen( outputdb, "w" ) ) == NULL) {
            printf("Error opening output file (%s).\n",outputdb);
            exit(1);
         }
         else printf("Open output file %s . . . . .Succeeded.\n",outputdb);
         continue;
      }
      else if (solv_t && !strcmp("solvation_table",field )) {
         solvtab = value;
         solvtab[strlen(solvtab)-1] = '\0';
         if ((solvtable = fopen( solvtab, "r" ) ) == NULL) {
            printf("Error opening solvation table (%s).\n",solvtab);
            exit(1);
         }
         else printf("Open solvation file %s . . . Succeeded\n",solvtab);
         continue;
      }
      else if (!strcmp("amber_table",field )) {
         amber = TRUE;
         ambtab = value;
         ambtab[strlen(ambtab)-1] = '\0';
         if (!strcmp("default", ambtab)) {
            strcpy(ambtab, "/usr/people/dock/v3.5b/parms/prot.table.ambcrg.ambH");
         }
         if ((ambf = fopen( ambtab, "r" ) ) == NULL) {
            printf("Error opening AMBER atom parameters (%s).\n",ambtab);
            exit(1);
         }
         else printf("Open AMBER Parameters %s . . . Succeeded\n",ambtab);
         continue;
      } 
      else if (!strcmp("flag",field )) printf("Not yet implimented\n");
      else if (!strncmp("color_table",field, 11 )) {
         while ( fgets (buffer, LINE_LENGTH, inhier) != NULL) {
            field = strtok(buffer, delim);
            if (!strcmp("default_color",field )) {
               strcpy(def_color, strtok(NULL, delim));
               def_color[strlen(def_color)-1] = '\0';
               break;
            }
            ncolor++;
            strcpy(color[ncolor].name, field);
            strcpy(color[ncolor].a1, strtok(NULL, delim));
            color[ncolor].c1 = 0;
            if ( (field=strtok(NULL, delim)) != NULL) {
               b = atoi(field);
               color[ncolor].c1 = b;
               value = strtok(NULL, delim);
               value[strlen(value)-1] = '\0';
               strcpy(color[ncolor].a2, value);
            }
            else {
               color[ncolor].a1[strlen(color[ncolor].a1)-1] = '\0';
            }
         }
         def_col = gen_ctab();
      }
      else {
         printf("field name unknown:  %s\n",field);
         printf("If this should be a color, make sure coloring is on\n");
         exit(1);
      }
   }
   printf("\n   Procedures\n");
   if (protein) {
      printf("Generate protein hierarchy\n");
   }
   if (equal_q) printf("Charge Equalization\n");
   if (solv_t) {
      printf("Solvation Correction\n");
      if (solvtab == NULL) {
         printf("No solvation file specifide\n");
         exit(1);
      }
   }
   if (amber) {
      printf("AMBER van der Waals types will be used\n");
   }
   if (torque_H) {
      printf("Terminal hydrogens will be rotated\n");
   } else {
      printf("Terminal hydrogens will NOT be rotated\n");
   }
   if (coloring) printf("Atoms will be colored\n");
   if (col_tab) printf("Print Color Table  (need the answer key, huh?)\n\n");
   if (!openfile) {
      printf(" *** Hey! I need an input file (mol2_file or mol2_file_list.");
      exit(1);
   }
}
 
bool tf(char *yn) {
   if (!strncmp("Y",yn,1) || !strncmp("y",yn,1)
      || !strncmp("T",yn,1) || !strncmp("t",yn,1)) return TRUE;
   return FALSE;
}

