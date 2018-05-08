// Author: I. Froehlich
// Written: 14.02.2008
// 

#ifndef _PBATCH_H_
#define _PBATCH_H_

#include "TInterpreter.h"
#include "TMethodCall.h"
#include "TMethodArg.h"

#include "PFormula.h"
#include "PValues.h"
#include "PParticle.h"

#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include <iostream>
#include "PDefines.h"

using namespace std;

#define MAX_COMMAND_POINTER    1000
#define MAX_COMMAND_TMETHODS   100
#define MAX_COMMAND_OPTIONS    10
#define MAX_STACK_GOSUB        20
#define MAX_READLINE_ARGS      20


#define IS_OBJECT 1
#define IS_DOUBLE 2

#define kFOREACH    2
#define kGOTO       4
#define kFOREACHEND 8
#define kUPDATE     16
#define kEOF        32
#define kPUSH       64
#define kPUSHBRANCH 128
#define kELSE       256

#define COMMAND_MASS  '1'
#define COMMAND_MASS2 '2'
#define COMMAND_GETBEAM   '3'
#define COMMAND_GETTARGET '4'
#define COMMAND_IF    '?'
#define COMMAND_ELSE  '%'
#define COMMAND_PLUS  '+'
#define COMMAND_MINUS '-'
#define COMMAND_DIV   '/'
#define COMMAND_MULT  '*'
#define COMMAND_IS    '='
#define COMMAND_EQUAL '~'


#define COMMAND_BOOST     'b'
#define COMMAND_COS       'c'
#define COMMAND_EVAL      'e'
#define COMMAND_FABS      'f'
#define COMMAND_INTERNAL  'i'
#define COMMAND_PRINT     'p'
#define COMMAND_ROT       'r'
#define COMMAND_THETA     't'
#define COMMAND_GETRANDOM 'v'
#define COMMAND_GETRANDOMX 'x'
#define COMMAND_GETRANDOMY 'y'


#define COMMAND_ANGLE    'A'
#define COMMAND_BRANCH   'B'
#define COMMAND_P3E      'E'
#define COMMAND_GOTO     'G'
#define COMMAND_PUSH     'H'
#define COMMAND_LABEL    'L'
#define COMMAND_PFORMULA 'P'
#define COMMAND_P3M      'Q'
#define COMMAND_RETURN   'R'
#define COMMAND_GOSUB    'S'
#define COMMAND_PVALUE   'V'
#define COMMAND_EXIT     'X'
#define COMMAND_ECHO     'Y'
#define COMMAND_READLINE 'Z'
#define COMMAND_FORMORE  '>'
#define COMMAND_FOREACH  '<'



#define METHOD_RETURN_DOUBLE 1
#define METHOD_RETURN_INT    2
#define METHOD_RETURN_VOID   4
#define METHOD_RETURN_PPARTICLE   8

#define STREAM_MAX_POS       "npar"
#define STREAM_DEFAULT_POS   "cpos"

class PBatch;
PBatch& fBatch();
PBatch* makeGlobalBatch();


class PBatch : public TObject {

 private:

    Bool_t CheckAndSplit(char * arg, char delim, int *key1, int *key2);
    
    Bool_t GetArguments(const char *a, /* const char *b, */ 
			char *name, char **function, char **arg1, char **arg2);
    Int_t  CheckObjectType(Int_t key);
    Int_t  GetKey(char *name, int fl, int makeflag);
    Int_t  GetDelimPosition(char *arg, char delim, Int_t *yes=NULL);

    Int_t  GetMethodHandle(char *name, Int_t flag=0);
    Int_t  CrackMethodArgs(char *name);
    Int_t  arg1, arg2, arg3, arg4;

    Int_t command_pointer, method_pointer, last_command_pointer;
    TMethodCall *methods[MAX_COMMAND_TMETHODS];
    char *method_name[MAX_COMMAND_TMETHODS];
    Int_t methods_flags[MAX_COMMAND_TMETHODS],
	methods_arg_flags[4][MAX_COMMAND_TMETHODS];

    char lst_command[MAX_COMMAND_POINTER];
    Int_t lst_command_int[MAX_COMMAND_POINTER];
    Int_t flag_command_int[MAX_COMMAND_POINTER];  //=0: PParticle, =1: PUtilsRengine
    Int_t lst_key_a[MAX_COMMAND_POINTER];
    Int_t error_flag[MAX_COMMAND_POINTER];
    Int_t lst_key[MAX_COMMAND_OPTIONS][MAX_COMMAND_POINTER];
    Int_t  lst_options_counter[MAX_COMMAND_POINTER];
    PFormula *lst_form[MAX_COMMAND_POINTER];
    //    Int_t lst_key2[MAX_COMMAND_POINTER];
    //    Int_t lst_key3[MAX_COMMAND_POINTER];
    //    Int_t lst_key4[MAX_COMMAND_POINTER];

    char *echo_string[MAX_COMMAND_POINTER];
    char *readline_format_string[MAX_COMMAND_POINTER];
    char *readline_string[MAX_COMMAND_POINTER];
    Double_t **readline_args[MAX_COMMAND_POINTER];
    int readline_num_args[MAX_COMMAND_POINTER];
    char *varlist;   //Allowed command for new variables

    Int_t batch_particle_param, batch_value_param, pid_param, num_command_param, num_batch_param, num_bulk_param,
	stream_default_pos_param, stream_max_pos_param, batch_update_param, batch_models_param;
    Int_t batch_histogram_param;
    Int_t num_batch,num_bulk, locnum_batch, locnum_bulk, locnum_command, locnum_old_command;  //Position in the projector & reaction loop
    Int_t locnum_branch; //current branch number
    Int_t else_position;
    Int_t current_position;

    static Int_t stack_num_batch[MAX_STACK_GOSUB], stack_num_bulk[MAX_STACK_GOSUB], 
	stack_num_command[MAX_STACK_GOSUB], stack_num_pos;

    void  AddSpacePlaceholder(char *command);
    void  RemoveSpacePlaceholder(char *command);
    Int_t EvalPFormula(char *command);
    

    PValues pdummy;

    TH1       *fHisto1;
    TH2       *fHisto2;
    TH3       *fHisto3;
    TGraph    *fGraph;
    TGraph2D  *fGraph2D;
    Int_t      is_readonly;
    TH1D      **slicesx, **slicesy;

    FILE      *file, *tmp_file;
    PParticle *current_particle; //For push command

    Double_t *x, *y, *z;  //Axis values for _eval
    Int_t eval_err_dumped;

    Int_t status;

    TTree *tree;         //Pointer to storage tree
    Int_t *size_branches;
    Int_t *key_branches;
    

 public:

    //constructor
    PBatch();

    void Print(const Option_t *delme=NULL) const ;

    Bool_t AddCommand(const char *command); //adds a command line to batch
    Bool_t AddCommand(char command, int key_a, int key1, int key2, int key3=-1, int key4=-1, int key5=-1);

    void SetPosition(Int_t my_num_batch, Int_t my_num_bulk) {
	num_batch = my_num_batch; 
	num_bulk  = my_num_bulk;
    };
    Int_t GetNewBatch()   {return locnum_batch;};
    Int_t GetNewBulk()    {return locnum_bulk;};
    Int_t GetNewCommand() {return locnum_command;};
    Int_t GetOldCommand() {return locnum_old_command;};
    Int_t GetBranch()     {return locnum_branch;};
    Int_t GetElsePosition() {return else_position;};
    Int_t GetCurrentPosition() {return current_position;};

    using TObject::Execute;
    Int_t Execute(Int_t command_pos=0, Int_t retval=kFALSE);
    Int_t ExecuteLastLine(void) {
	Bool_t retval = Execute(last_command_pointer);
	last_command_pointer = command_pointer;
	return retval;
    }
    Bool_t Execute(const char *command) {
	if (AddCommand(command)) {
	    return (Bool_t)ExecuteLastLine();
	}
	last_command_pointer = command_pointer;
	Error("Execute", "Command not executed");
	return kFALSE;
    }

    void SetToolObject(TH1 *Histo1) {
	fHisto1 = Histo1;
    }
    void SetToolObject(TH2 *Histo2) {
	fHisto2 = Histo2;
    }
    void SetToolObject(TH3 *Histo3) {
	fHisto3 = Histo3;
    }
    void SetToolObject(TGraph *Graph) {
	fGraph = Graph;
    }
    void SetToolObject(TGraph2D *Graph2D) {
	fGraph2D = Graph2D;
    }
    void SetToolObject(FILE *f) {
	file = f;
    }
    void SetToolObjectTmp(FILE *f) {
	tmp_file = f;
    }
    Int_t IsReadonly() {
	return is_readonly;
    }

    PParticle *GetCurrentParticle(void) {
	return current_particle;
    }

    Int_t Status(void) {
	//workaround to check for a "[+]"
	//=1 for such a case
	//=0 in all other cases
	return status;
    }

    void SetVarList(char *my_x) {
	//Allowed commands for new variables
	//Format must be "a;b;c;" with trailing ;'s 
	//If NULL don't care!
	varlist = my_x;
    };

    void ReplaceAll(TString *op, const char *oldstring, const char *newstring);

    void SetTree(TTree *my_tree) {tree = my_tree;};
    TTree *GetTree(void)         {return tree;};

    void SetSizeBranches(Int_t *my_size_branches) {size_branches = my_size_branches;};
    void SetKeysBranches(Int_t *my_key_branches)  {key_branches = my_key_branches;};

    ClassDef(PBatch, 0)  //Batch commands
};

#endif
