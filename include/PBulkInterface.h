// Author: Ingo Froehlich
// Written: 10/07/2007
// Modified: 
// PBulkInterface Class Header

#ifndef _PBULKINTERFACE_H_
#define _PBULKINTERFACE_H_

#include "PParticle.h"
#include "TTree.h"

#define PPROJECTOR_PRIORITY    99
#define FILTER_PRIORITY        60
#define DECAY_PRIORITY         50
#define FILEINPUT_PRIORITY     1

class PBulkInterface;
PBulkInterface *makeGlobalBulk();
PBulkInterface &fPBulkInterface();

class PBulkInterface: public TObject {

 private:

    static Int_t gBulkCounter; 
    
 protected:
    
    Double_t current_weight;
    Int_t    bulk_id;      //Unique bulk ID
    Int_t    fPriority;    //order when adding the bulk in PReaction
    TTree   *tree;         //Pointer to storage tree
    Int_t   *size_branches;
    Int_t   *key_branches;
    Int_t   *current_size_branches[MAX_NUM_BRANCHES];
    PParticle **particle_array_branches[MAX_NUM_BRANCHES];

 public:

    PBulkInterface();
    
    virtual bool Modify(PParticle ** array, int *decay_done, int *num, int maxnum);  //Modify particle bulk

    void SetWeight(Double_t c) {current_weight = c;};
    Double_t GetWeight(void) {return current_weight;};

    void SetPriority(Int_t p) {
	fPriority = p;
    };
    Int_t GetPriority() {return fPriority;};

    void SetTree(TTree *my_tree) {tree = my_tree;};
    TTree *GetTree(void)         {return tree;};

    void SetSizeBranches(Int_t *my_size_branches)  {size_branches = my_size_branches;};
    void SetKeysBranches(Int_t *my_key_branches)   {key_branches = my_key_branches;};
    void SetBranchNum   (Int_t i, Int_t *my_size)  {current_size_branches[i] = my_size;};
    void SetBranchArray (Int_t i, PParticle **my_particle) {particle_array_branches[i] = my_particle;};

    Int_t *GetSizeBranches()     {return size_branches;};
    Int_t *GetKeysBranches()     {return key_branches;};
    Int_t *GetBranchNum(Int_t i) {return current_size_branches[i];};
    PParticle **GetBranchArray (Int_t i) {return particle_array_branches[i];};

ClassDef(PBulkInterface,0) // Pluto bulk interface base class
};


#endif 

















