/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
#include "mol2hier.h"

void list_branches(int root, int lead);
int complete_level(int root, int lead, int branch[50], int hl);
void get_neighbors(int atm, int neighbors[6], int *n1c);
void assign_level(int root, int hl, int current_level[MAXATOMS]);
int reset_lead(int atm);
int conf_hier(int tpos);
int lookup_amb(char *res, char *at); 
int find_rigid_component(void);

/* identify rigid fragment and atoms connected to it */
void hiergen(bool protein) {
   int b, i, j, k, hl, hel, n1[6], n1c, atm, tmp; 
   int tpos, branch[50];
   tpos = 0;
   hel = 0;
   nbr = 0; //number of branches
   reset(set_a); //set_a is the a map from possible branch id to last atom
                 //in the rigid fragment
   for (i=0; i < MAXCONF; done[i++] = FALSE); 
   for (i=0; i < MAXATOMS; upos[i++] = 0); 
   for (i=0; i < MAXATOMS; atlvl[i++] = 0); 
   for (i=0; i < MAXATOMS; i++) rootlist[i][0] = rootlist[i][1] = 0;
   printf("distol = %d\n", distol);
   /* identify rigid atoms */
   for (j=1; j<=numats; j++) {
      //one way (amber) or another (dock 3.5 adapted vdwtype function)
      //get and set the vdw type
      if (amber && !strncmp("C", atom[j].atname,1)) {
         atom[j].vdw = lookup_amb(atom[j].resname, atom[j].atname);
      } else {
         atom[j].vdw = vdwtype(j, atom[j].name);
      }
      for (i=1; i<nconf; i++) {
         //this copies the vdw type to all conformations
         k = j + (i * numats);
         atom[k].vdw = atom[j].vdw;
      }
      upos[j] = unique_positions(j);
      //printf("%d %d\n", j, upos[j]); //debug the unique positions
      if ((torque_H) && (torque[j].hyd == 1)) {
         upos[j] = upos[j] * (int)(360/torque[j].dihe[0]);
      } 
      //printf("unique positions for atom %d are %d\n", j, upos[j]);
      tpos = tpos + upos[j]; //the total positions flag tpos
   }
   i = find_rigid_component();//new way of finding first atom in rigid component
   hl = 9;
   hel = complete_level(0, i, branch, hl); //this call makes the rigid component
   if (hel != numats) { //numats is the total number of atoms read in
      hel = numats;
      b = 0;
      while (branch[++b] != 0) { //look through the branches found
         n1c = 0;
         atm = branch[b]; //atm is potential branch start point
         for (i=1; i <= numbnd; i++) {
           if (bond[i].atom1 == atm ) 
              n1[++n1c] = bond[i].atom2;
           else if (bond[i].atom2 == atm ) 
              n1[++n1c] = bond[i].atom1;
         }
         //n1 is the list of things bound to the branch start point
         for (i=1; i<=n1c; i++) 
           if (upos[n1[i]] == 1) { //unique positions == 1, then
             set_a[b] = n1[i];    //put this branch into set_a  (root list)
           }
      }
      //branch is the list of neighboring atoms of the rigid component
      b = 0;
      while (branch[++b] != 0) {
         //printf("calling list_branches %d %d %d\n", b, set_a[b], branch[b]);
         list_branches(set_a[b], branch[b]);
         if (set_a[b] > 0) {
            nbr++; //number of branches, used later (global variable)
            int count = 0;
            for (count = 1; count < b; count++) {
               if (set_a[count] == set_a[b]) {
                  nbr--; //actually not a branch, already seen same rigid start
                  break; //can quit with the for loop now
               }   
            }
         }
      }
      //printf("nbr:%d\n", nbr);
      //debugging output that prints input atom#->db group assignments 
      //for (j=1; j<=numats; j++) printf("%d = %d\n",j,atlvl[j]);  
   }
   //all atoms have been assigned to a group/level
   tmp = conf_hier(tpos);
   for (i=1; i<= hindex; i++) 
      hier[i][0] = hier[i][1] = 0;
   if ((tmp == -1 || tmp == -3)) {
      //this used to trigger a lowering of the distance tolerance
      //and running hiergen again. now it just exits with the error
      //i'm not sure anyone knows what/when/how these errors occur
      printf("Fatal error %d from conf_hier()\n", tmp);
      exit(tmp);
   } 
}

// recursive explore down branch, root is in the rigid component and lead
// is the first atom in the branch, root and lead are bonded to each other
void list_branches(int root, int lead) {
   int branch[50];
   int i = 0;
   complete_level(root, lead, branch, 9); //finishes this branch
   i = 0;
   //now branch has data about the neighbors of this component that haven't 
   //been processed, so do them
   while (branch[++i] != 0) {
      lead = reset_lead(branch[i]);
      list_branches(lead, branch[i]);
   }
}

//root is the atom index of the atom in the previous group connected to this 1
// special case 0 for rigid component
//lead is the first atom in this level (group, component, etc)
//branch is used to keep track of the branches and is used by the calling funct
//hl is always 9??, basically the root output name 
//branch[50] is where potential branches are stored
int complete_level(int root, int lead, int branch[50], int hl) {
   //printf("root:%d lead:%d hl:%d\n", root, lead, hl);
   int neighs[6], neighs_count, bc, i, j, level_upos;
   int current_level[MAXATOMS]; //replacing set_c with this local, list of 
                                //the atoms in the current level
   reset(current_level);
   cnum = bc = i = 0;
   level_upos = upos[lead];
   current_level[++cnum] = lead; //put the first atom index in the set 
                         //representing the current level
   //quit this loop if nothing new was added to the current_level in last run
   while (current_level[++i] != 0) {
      //get neighbors of the last thing added to the set
      get_neighbors(current_level[i], neighs, &neighs_count); 
      //neighs has been updated by the call to contain up to 6 bonded atom indxs
      for (j=1; j<=neighs_count; j++) { //j is a neighbor index
         //if this neighbor has more unique positions then it is a branch
         if (upos[neighs[j]] > level_upos)  
            branch[++bc] = neighs[j]; 
         //otherwise if it is not already in a level and it hasn't been seen
         //already then put it in this level
         else if (!is_member(neighs[j],current_level,cnum)&& !done[neighs[j]]) {
            current_level[++cnum] = neighs[j];  
         }                         
      }                            
   }
   branch[++bc] = 0; //make sure the branch list ends with 0
   //set up the rootlist to be a list of atom indices and level output names
   //also uses current_level (list of atom indices in this level) to set atlvl (
   // the map from atom index to output level) and update the done flag
   if (root != 0) {
      assign_level(root, hl, current_level);
      return 0; //added but return not important 
   } else { //stuff to do for root case
      i = 0;
      int total_atoms_assigned = 0;
      while (current_level[++i] != 0) {
         atlvl[current_level[i]] = hl; //hl is always 9 so this sets the rigid output group
         done[current_level[i]] = TRUE; //sets done since not done in assign_level for rig
         total_atoms_assigned++; //counter of total number of atoms assigned
      }
      return total_atoms_assigned; //number of atoms in the rigid component
   }
}

//this function is called with the atm first (input) and a 6 length int array
// and a pass by reference int, neighbors and n1c respectively (output)
//it uses the information in bond to find the neighboring atoms (up to 6)
//and returns the atm reference id numbers in neighbors and the number of them
//in n1c (neighbors1count basically).
void get_neighbors(int atm, int neighbors[6], int *n1c) {
   int i;
   int c = 0;
   for (i=1; i <= 5; i++) neighbors[i] = 0;
   for (i=1; i <= numbnd; i++) {
      if (bond[i].atom1 == atm ) neighbors[++c] = bond[i].atom2;
      else if (bond[i].atom2 == atm ) neighbors[++c] = bond[i].atom1;
   }
   *n1c = c;
}

//rootlist[index][0] is either 0 or the root of a level's atom index
//rootlist[index][1] is either 0 or the output name (9, 19, 18, etc) of level
//assign_level sets these up after each level is found in complete_level
void assign_level(int root, int hl, int current_level[MAXATOMS]) {
   int i, lvl;
   i = lvl = 0;
   while (rootlist[++i][0] != 0)
      if (rootlist[i][0] == root) 
         lvl = rootlist[i][1];
   if (lvl == 0) {
      if (rootlist[1][0] == 0) {
         rootlist[1][0] = root;
         rootlist[1][1] = hl;
      }
      lvl = tens(atlvl[root]) * 10 + 19; 
      i = 0;
      //to find a new level name, make sure it isn't in use already
      while (rootlist[++i][0] != 0) 
         if (rootlist[i][1] == lvl) {
            lvl--;
            i=0;
         }
      //do the assignment for this level. each index of rootlist is a new level
      rootlist[i][0] = root;
      rootlist[i][1] = lvl;
   }
   //also want to set the atlvl which is used during output, use current_level
   //as the list of atoms in this level, also set the done flag
   i=0;
   while (current_level[++i] !=0) {
      atlvl[current_level[i]] = lvl;
      //printf("%d root = %d (%d)\n",current_level[i],root,lvl); 
      done[current_level[i]] = TRUE;
   }
}

//pick a new 'lead' atom for the group??
int reset_lead(int atm) {
   int i, tmp;
   tmp = 0;   
   for (i=1; i <= numbnd; i++) {
      if (bond[i].atom1 == atm ) 
         tmp = bond[i].atom2;
      else if (bond[i].atom2 == atm ) 
         tmp = bond[i].atom1;
      else {
         continue;
      }
      if (atlvl[tmp] > 0) return tmp;
   }
   //printf("failure in reset_lead %d\n", tmp);
   return tmp; //i think this is a failure case
}

//finds the largest connected group of atoms with a single unique position
int find_rigid_component(void) {
  int atom_count = 0;
  int possibly_rigid[MAXATOMS]; //0 if no way, 1 if 1 unique position
                                //2 or greater during checks (takes others down)
  for (atom_count = 1; atom_count < MAXATOMS; atom_count++) {
    if (1 == upos[atom_count])
      possibly_rigid[atom_count] = 1;
    else
      possibly_rigid[atom_count] = 0;
  }
  for (atom_count = 1; atom_count < MAXATOMS; atom_count++) {
    if (1 == possibly_rigid[atom_count]) {
      int bonded_neighbors[6];
      int bond_neigh_total = 0;
      int neighbor_stack[MAXATOMS];
      int neighbor_count = 0;
      for (neighbor_count = 0; neighbor_count < MAXATOMS; neighbor_count++)
        neighbor_stack[neighbor_count] = -1;
      neighbor_stack[0] = atom_count; //everything on this stack needs processed
      neighbor_count = 1; //keeps track of size of stack
      while (neighbor_count > 0) {
        neighbor_count--;
        int temp_atom = neighbor_stack[neighbor_count];
        get_neighbors(temp_atom, bonded_neighbors, &bond_neigh_total);
        int bonded_count = 1; //1 -indexed return from get_neighbors
        for (bonded_count=1;bonded_count<=bond_neigh_total;bonded_count++) {
          int current_neighbor = bonded_neighbors[bonded_count];
          if (1 == possibly_rigid[current_neighbor]) {
            possibly_rigid[atom_count]++;
            possibly_rigid[current_neighbor] = 0; //reset since added to prev line
            neighbor_stack[neighbor_count] = current_neighbor; //add to stack
            neighbor_count++; //increase counter
          }
        }
      }
    }
  }  
  int biggest_component = 0;
  int biggest_count = 0;
  for (atom_count = 1; atom_count < MAXATOMS; atom_count++) {
    if (biggest_count < possibly_rigid[atom_count]) {
      biggest_count = possibly_rigid[atom_count];
      biggest_component = atom_count;  
      //printf("%d %d\n", atom_count, possibly_rigid[atom_count]);
    }
  }
  return biggest_component;
}


