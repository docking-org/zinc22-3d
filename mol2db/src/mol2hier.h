/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#ifndef mol2hier_h_
#define mol2hier_h_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LINE_LENGTH 150
#define MAXCONF 5010
#define MAXATOMS 210
#define MAXCOLOR 200
#define MAXMOLS 1500
#define MAXGRP 1000
#define MAXFGRP 25000
#define MAXFATOMS 300000
#define MAXAMB 400

typedef enum {FALSE, TRUE} bool;

/* Atm struct stores all atoms of all confs for a given molecule */
typedef struct _Atm{ int num,subst, x, y, z, charge, vdw, color;
                     char name[6], atname[5], resname[4];} Atm;
extern Atm atom[MAXCONF * MAXATOMS];

/* atom[i].name actually means atom type because I'm an idiot */


/* Bnd struct stores a single copy of the connect information for a mol */
typedef struct _Bnd{ int num, atom1, atom2, order; } Bnd;
extern Bnd bond[MAXATOMS];

/* Molecule stores FCD and Name.  One copy per mol  */
/* April 2014. was 4+8 now 4+16, thus 13 ->17 */
typedef struct _Molecule{ char mfcd[17]; char name[80]; } Molecule;
extern Molecule mol[MAXCONF];

/* Solv stores the solvation table */
typedef struct Solv{ char mfcd[17]; int cnt, crg[MAXATOMS]; 
     float pol[MAXATOMS], surf[MAXATOMS], apol[MAXATOMS], tot[MAXATOMS]; } Solv;
extern Solv solv[MAXMOLS];

/* hyd_tor stores info about dihedrals of H's to be torqued */
typedef struct hyd_tor{ bool hyd; int ang, cnt, dihe[5]; } hyd_tor;
extern hyd_tor torque[MAXATOMS];

/* Color table.  One copy per mol  */
typedef struct Colortable{ int num; char name[30]; char a1[5]; 
          int c1; char a2[5]; } Colortab;
extern Colortab color[MAXCOLOR];

/* AMBER parameters */
typedef struct Ambertable{char atom[5]; char res[4]; float crg; int
type; } Ambertable;
extern Ambertable amb[MAXAMB];

extern FILE* inhier;
extern FILE* mol2file;
extern FILE* filelist;
extern FILE* outf;
extern FILE* solvtable;
extern FILE* ambf;
extern fpos_t solvptr, solvbeg;
extern char *flist, inputmol2[80], *outputdb, *solvtab, *ambtab;

extern int distol; //needs to be 2 since omega actually outputs slightly
                   //different thousands of a decimal point angstroms
extern int nconf; //nconf is number of input conformatios (1 for each mol2)
extern int numats, numbnd, numhyd, cnum;
extern int ncolor, tot_atm, ens, def_col, nsolv;
extern int set_a[MAXATOMS], set_b[MAXATOMS]; //what are these???
extern int set_c[MAXATOMS]; //set_c used to keep track of atoms in current level
extern int atlvl[MAXATOMS]; //atlvl is the map from atom -> level name
extern int rootlist[MAXATOMS][2];
//rootlist[index][0] is either 0 or the root of a level's atom index
//rootlist[index][1] is either 0 or the output name (9, 19, 18, etc) of level
extern int src, dest, knum, nwritten, nhyd, hindex;
extern int minx, miny, minz, ind;
extern int key[200][3];  /* 0: level, 1: lead atom, 2: atoms in level */
extern int keyWholeAtom[200][MAXATOMS]; //saves the whole list of atom indxs
extern int keyorder[200], keycntr;  /* order or travers branches */
extern int hier[MAXFATOMS][4]; /* 0: level, 1: atom number, 2: mol. ptr */
extern int upos[MAXATOMS];  /* unique positions for atom */
extern int nbr, amb_cnt; 
extern int list[MAXFGRP][3]; //0 is output group id, 1 is # of confs, 2 is links
extern int brid[MAXGRP];
extern int bconf[10]; /* number of confs for branch n */
extern double oconf;

extern bool skip, trans_coord, protein, torque_H, amber;
extern bool log_file, col_tab, coloring, solv_t, equal_q, hier_sp;
extern bool done[MAXCONF];

extern char c;
extern char buffer[LINE_LENGTH];
extern char fam[5];
extern char def_color[31]; /*name of default color */
float roundD2F(double);

int vdwtype(int num, char name[] );

void hiergen(bool protein);

void reset(int temp[MAXATOMS]);

void moveset(int set1[MAXATOMS], int set2[MAXATOMS]);

void getds(int input[MAXATOMS], int output[MAXATOMS]);

bool is_member(int at, int temp[MAXATOMS], int natm);

bool proximal(int a, int b);

void add_confs(int tpos);

int unique_positions(int atm);

void colorit(void);

char *substr( char *line, int start, int length);

void fixcrgs(void);

void list_branches(int root, int lead);

void add_atms(int lvl, int conf, int molnum);

void write_hierarchy(int br);

void translate(void);

int tens(int x);

int lvl_cnt(int lvl, int *lead);

int lvl_cnt_save(int lvl, int save[MAXATOMS]);

bool term_n_s_o_h(int dihe[4]);

#endif
