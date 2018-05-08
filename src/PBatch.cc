////////////////////////////////////////////////////////
// Batch commands for particle and value operations
//
// This class allows for particle operations for
// analysing angles, masses, or making boost and
// rotations of TLorentzVectors/PParticles
//
// It is used for the PProjector to make histograms
// _inside_ the event loop without the need of an
// analysis macro
//
// The syntax is inspired by C++, but only a few
// rules of C++ have been implemented. Furthermore, the syntax
// has been simplified to make it more readable, and an interface
// to the particle input (and maybe to the database in a later step)
// has been foreseen
//  
// As usual we have assignments, operators, functions and methods
// Each operation can be seperated by a semicolon (;)
//
// A simple case is outlined in the following command
// _x = [p,1]->M();
// This means that the double _x (which is automatically constructed)
// is assigned to the mass of the particle [p,1]
// Particle in brackets are the only allowed use without instantiantion
// as they are (hopefully) filled by the input. The additional number 
// can be used for ordering information
//
// + and - operation
// This adds/subtracts PParticles or doubles:
// pp = [p,1] + [p,2];
// ...nothing else then the pp composite
// _x = pp->M();
// gives the invariant mass of the pp pair
// CAVEAT: Do _not_ use object name which are
// already assigned to particles in the data base, e.g.
// pi0=[pi0] is not allowed, use _pi0=[pi0] instead
// 
// Brackets: Objects can be nested:
// _x = ([p,1] + [p,2])->M();
// ...makes the same as the 2 lines above
//
// Build-in methods:
// M()              : Invariant mass of the object
// M2()             : Invariant mass^2 of the object
// Boost(obj)       : Boost into rest frame of "obj"
//   Example:
//   pp = p1->Boost([p,1] + [p,2]);
//   ...boost p1 in the rest frame of [p,1]+[p,2]
// Rot(obj)         : Rotate object such that obj
//                    would point to z-Direction
// Angle(obj)       : Opening angle between obj
// Theta()          : Theta of momentum
// Print()          : Dump 4momentum or double
//
// In addition, all browsable methods of the class PParticle
// can be used (functions without arguments)
//
// Example:
// obj->Px();
//
// Build-in functions:
// cos()
// fabs()
//
// In addition, the TFormula syntax can be used
//
// Example:
// val_new = (val + 1.2);
// val = (val * TMath::Pi());
//
// Conditions:
// if(arg)          : Interrupts the current chain if arg == 0
//
// A very complex example: Look for the helicity distribution
// of the e+e- pair in the eta Dalitz decay
// _eta=[eta]; _ep=[e+]; _em=[e-]; 
// _eta->Boost([p + p]); _em->Boost([p + p]); _ep->Boost([p + p]); 
// _ep->Rot(_eta); _em->Rot(_eta); _eta->Rot(_eta) ; 
// _ep->Boost(_eta); _em->Boost(_eta); dil=_ep+_em; 
// _ep->Rot(dil); _em->Rot(dil); dil->Rot(dil); _ep->Boost(dil); _em->Boost(dil) ; 
// s1= _ep->Theta(); _x = cos(s1)
//
//
//                    Author: I. Froehlich
//                    Written: 14.02.2008
//                    Revised: 
//
////////////////////////////////////////////////////////

//makeDataBase()->ListEntries(-1,1,"name,*batch_particle,*pid,*batch_value");

#include "PBatch.h"
#include "TString.h"
#include "PUtils.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "PParticle.h"
#include "PStaticData.h"
#include "PChannelModel.h"

#include "TClass.h"
#include "TMethod.h"

#include <cmath>

Int_t PBatch::stack_num_pos = 0;
Int_t PBatch::stack_num_batch[MAX_STACK_GOSUB], 
    PBatch::stack_num_bulk[MAX_STACK_GOSUB], 
    PBatch::stack_num_command[MAX_STACK_GOSUB];

PBatch &fBatch() {
    static PBatch *ans = new PBatch();
    return *ans;
}

PBatch *makeGlobalBatch() {
    return &fBatch();
}

PBatch::PBatch() {

    makeStaticData();
    
    command_pointer = last_command_pointer = 0;
    method_pointer  = 0;

    fHisto1          = NULL;
    fHisto2          = NULL;
    fHisto3          = NULL;
    fGraph           = NULL;
    fGraph2D         = NULL;
    is_readonly      = 0;
    slicesx = slicesy = NULL;
    current_particle = NULL;

    batch_particle_param = makeDataBase()->GetParamTObj("batch_particle");
    batch_models_param   = makeDataBase()->GetParamTObj("batch_models");
    batch_value_param    = makeDataBase()->GetParamDouble("batch_value");
    pid_param            = makeDataBase()->GetParamInt("batch_pid");
    Int_t batch_position_param = makeDataBase()->GetParamInt("batch_position");

    if (batch_particle_param < 0) 
	batch_particle_param = makeDataBase()->MakeParamTObj("batch_particle", 
							     "PParticle storage for batch");

    if (batch_value_param < 0) 
	batch_value_param = makeDataBase()->MakeParamDouble("batch_value", 
							    "Value storage for batch");

    if (batch_models_param < 0) 
	batch_models_param = makeDataBase()->MakeParamTObj("batch_models",
							   "Storage for distribution objects");

    if (pid_param < 0) 
	pid_param = makeDataBase()->MakeParamInt("batch_pid", "PID for batch");

    if (batch_position_param < 0) 
	makeDataBase()->MakeParamInt("batch_position", "PID position for batch");

    batch_histogram_param = makeDataBase()->GetParamTObj("batch_histogram");
    if (batch_histogram_param < 0) 
	batch_histogram_param = makeDataBase()->MakeParamTObj("batch_histogram", "Histogram storage for batch");
    

    //This is used for the "goto"
    num_command_param = makeDataBase()->GetParamInt("num_command");    
    if (num_command_param < 0) 
	num_command_param = makeDataBase()->MakeParamInt("num_command", "Number of command for a label");
    num_batch_param = makeDataBase()->GetParamInt("num_batch");    
    if (num_batch_param < 0) 
	num_batch_param = makeDataBase()->MakeParamInt("num_batch", "Number of batch object for a label");
    num_bulk_param = makeDataBase()->GetParamInt("num_bulk");    
    if (num_bulk_param < 0) 
	num_bulk_param = makeDataBase()->MakeParamInt("num_bulk", "Number of bulk for a label");
    num_batch = num_bulk = -1;

    //for the branching
    if (makeDataBase()->GetParamInt("branch_idx") < 0) {
	makeDataBase()->MakeParamInt("branch_idx", "Index of the tree branch in PReaction");
    }
    
    //Used for "formore"
    stream_default_pos_param = makeDataBase()->GetParamInt(STREAM_DEFAULT_POS);
    if (stream_default_pos_param < 0)  
	stream_default_pos_param = makeDataBase()->MakeParamInt(STREAM_DEFAULT_POS, "Default position");
    stream_max_pos_param = makeDataBase()->GetParamInt(STREAM_MAX_POS);
    if (stream_max_pos_param < 0)  
	stream_max_pos_param = makeDataBase()->MakeParamInt(STREAM_MAX_POS, "Max position in stream");

    batch_update_param = makeDataBase()->GetParamInt("batch_update");
    if (batch_update_param < 0)
	makeDataBase()->MakeParamInt("batch_update", 
				     "If set this is a variable which must trigger an update"); 

    for (int i=0; i<MAX_COMMAND_POINTER; i++) error_flag[i] = 0;

    x = makeStaticData()->GetBatchValue("_x");
    y = makeStaticData()->GetBatchValue("_y");
    z = makeStaticData()->GetBatchValue("_z");
    eval_err_dumped = status = 0;
    varlist = NULL;
    file    = tmp_file = NULL;

    tree          = NULL;
    size_branches = NULL;
    key_branches  = NULL;
    locnum_branch = 0;
    else_position = -1;
}

Int_t PBatch::Execute(Int_t command_pos, Int_t retval) {

    //retval = kFALSE; //standard, can be overwritten for a foreach
    //cout << "command:" << command_pos << ":" << command_pointer << ":" << retval << endl;

    current_position = -1;
    locnum_old_command = -1;

    for (int i=command_pos; i<command_pointer; i++) {
	TObject  *res  = NULL;
	TObject  *res2 = NULL;
	TObject  *res3 = NULL;
	Double_t *val  = NULL;
	Double_t *val2 = NULL;
	Double_t *val3 = NULL;
	Double_t *val4 = NULL;
	
	current_position = i;

	//cout << "command:" << lst_command[i] << " num "<< i <<  endl;

	if (lst_command[i] == COMMAND_PLUS) {
	    //
	    //+ of pparticles or values
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    makeDataBase()->GetParamTObj(lst_key[2][i], batch_particle_param, &res2);
	    makeDataBase()->GetParamTObj(lst_key_a[i],  batch_particle_param, &res3);
	
	    if (res && res2 && res3) {
		PParticle *a1 = (PParticle *) res;
		PParticle *a2 = (PParticle *) res2;
		PParticle *r  = (PParticle *) res3;
		*r = *a1;
		r->AddTmp(*a2);
		found = kTRUE;
	    } 
	
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key_a[i],  batch_value_param, &val3);
	
	    if (val && val2 && val3) {
		*val3 = *val2 + *val;
		found = kTRUE;
	    }
	    if (!found) return retval;
	} else if (lst_command[i] == COMMAND_MINUS) {
	    //
	    //+ of pparticles or values
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    makeDataBase()->GetParamTObj(lst_key[2][i], batch_particle_param, &res2);
	    makeDataBase()->GetParamTObj(lst_key_a[i],  batch_particle_param, &res3);
	
	    if (res && res2 && res3) {
		PParticle *a1 = (PParticle *) res;
		PParticle *a2 = (PParticle *) res2;
		PParticle *r  = (PParticle *) res3;
		*r = *a1 - *a2;
		found = kTRUE;
	    } 
	
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key_a[i],  batch_value_param, &val3);
	
	    if (val && val2 && val3) {
		*val3 = *val - *val2;
		found = kTRUE;
	    }
	    if (!found) return found;
	} else if (lst_command[i] == COMMAND_MULT) { 
	    //
	    // Internal *=
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key_a[i],  batch_value_param, &val3);	
	    if (val && val2 && val3) {
		*val3 = *val * *val2;
		found = kTRUE;
	    } 
	    if (!found) return found;
	} else if (lst_command[i] == COMMAND_DIV) { 
	    //
	    // Internal *=
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key_a[i],  batch_value_param, &val3);	
	    if (val && val2 && val3) {
		*val3 = *val / *val2;
		found = kTRUE;
	    } 
	    if (!found) return found;
	} else if (lst_command[i] == COMMAND_EQUAL) {
	    //
	    // compare values (+/- 0.5)
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key_a[i],  batch_value_param, &val3);
	
	    if (val && val2 && val3) {
		if (fabs(*val - *val2) < 0.5) 
		    *val3 = 1.;
		else 
		    *val3 = 0.;
		found = kTRUE;
	    }
	    if (!found) return retval;
	} else if (lst_command[i] == COMMAND_BOOST) { 
	    //
	    //boost pparticles
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a1  = (PParticle *) res;
	    makeDataBase()->GetParamTObj(lst_key[2][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a2  = (PParticle *) res;
	
	    a1->Boost(-a2->BoostVector()) ;
	} else if (lst_command[i] == COMMAND_GETBEAM) { 
	    //
	    //read the beam from composite
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) {
		Warning("GetBeam()", "Composite not found");
		return retval;
	    }
	    PParticle *a1  = (PParticle *) res;
	    makeDataBase()->GetParamTObj(lst_key_a[i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a2  = (PParticle *) res;
	    PParticle *a3  = a1->GetBeam();
	    if (!a3) {
		Warning("GetBeam()", "Beam not found");
		return retval;
	    }
	    *a2 = *a3;
	} else if (lst_command[i] == COMMAND_GETTARGET) { 
	    //
	    //read the target from composite
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a1  = (PParticle *) res;
	    makeDataBase()->GetParamTObj(lst_key_a[i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a2  = (PParticle *) res;
	    PParticle *a3  = a1->GetTarget();
	    if (!a3) {
		Warning("GetTarget()", "Target not found");
		return retval;
	    }
	    *a2 = *a3;
	} else if (lst_command[i] == COMMAND_ANGLE) { 
	    //
	    //angle between particle tracks
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a1  = (PParticle *) res;
	    makeDataBase()->GetParamTObj (lst_key[2][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a2  = (PParticle *) res;
	    
	    Double_t angle = a1->Vect().Angle(a2->Vect());
	    
	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
	    if (!val) {
		Warning("Angle()", "Result value not found");
		return kFALSE;
	    }
	    *val = angle;
	} else if (lst_command[i] == COMMAND_ROT) { 
	    //
	    //rotate
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a1  = (PParticle *) res;
	    makeDataBase()->GetParamTObj(lst_key[2][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a2  = (PParticle *) res;
	    double tmp_phi   = a2->Phi();
	    double tmp_theta = a2->Theta();

	    a1->RotateZ(-tmp_phi);
	    a1->RotateY(-tmp_theta);	 
	} else if (lst_command[i] == COMMAND_IS) {
	    //
	    // '='
	    //
	    Bool_t found = kFALSE;
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    makeDataBase()->GetParamTObj(lst_key_a[i], batch_particle_param, &res2);
	    if (res && res2) {
		PParticle *r   = (PParticle *) res2;
		PParticle *a1  = (PParticle *) res;
		//store PID, if old PParticle had already an pid, keep it
		Int_t pid  = r->ID();
		Double_t w = r->W();
		*r = *a1;
		if (pid) {
		    r->SetID(pid);
		    r->SetW(w);
		}
		found = kTRUE;
	    } 

	    makeDataBase()->GetParamDouble(lst_key_a[i],  batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);

	    if (val && val2) {
		*val2 = *val;
		found = kTRUE;
	    }
	    
	    if (!found) return retval;
	    Int_t *update;
	    if (makeDataBase()->GetParamInt(lst_key_a[i], batch_update_param, &update)) {
		if (*update == 1) {
		    locnum_old_command = i + 1;
		    return kUPDATE;
		}
	    }

	} else if (lst_command[i] == COMMAND_PUSH) {
	    //
	    // push()
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (res) {
		current_particle = (PParticle *)res;
		//current_particle->Print();
		if (locnum_old_command < 0) locnum_old_command = i + 1;
		if (lst_key[2][i] == -1) {
		    //old version (function or empty method)
		    return retval | kPUSH;
		}
		locnum_branch = 0;
		makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val);
		if (!val) return retval;
		locnum_branch = *val;
		return retval | kPUSHBRANCH;
	    }
	} else if (lst_command[i] == COMMAND_MASS2) {
	    //
	    // M2()
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a1  = (PParticle *) res;
	    
	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
	    if (!val) return retval;
	    *val = a1 -> M2();
	} else if (lst_command[i] == COMMAND_MASS) {
	    //
	    // M()
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) {
		return retval;
	    }
	    PParticle *a1  = (PParticle *) res;

	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
	    if (!val) {
		Warning("M()", "Result value %s not found",
			makeDataBase()->GetName(lst_key_a[i]));
		return retval;
	    }
	    *val = a1 -> M();
	} else if (lst_command[i] == COMMAND_THETA) {
	    //
	    // Theta()
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (!res) return retval;
	    PParticle *a1  = (PParticle *) res;
	    
	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
	    if (!val) return retval;
	    *val = (a1 -> Theta());
	} else if (lst_command[i] == COMMAND_INTERNAL) {
	    //
	    // Pointer to PUtils/PParticle/PChannelModel
	    //
	    TObject *a1;
	    if (flag_command_int[i] == 0) {
		makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
		if (!res) return retval;
		a1 = (PParticle *) res;
	    }
	    else if (flag_command_int[i] == 1) {
		a1 = makePUtilsREngine();
	    } else {
		makeDataBase()->GetParamTObj(lst_key[1][i], batch_models_param, &res);
		if (!res) return retval;
		a1  = (PChannelModel *) res;
	    }

 	    Double_t *argval;

	    //First, we have to set the variables
	    methods[lst_command_int[i]]->ResetParam();
	    for (int j=0; j<4; j++) {
		//cout << lst_command_int[i] << ":" << j << endl;
		if (methods_arg_flags[j][lst_command_int[i]] == METHOD_RETURN_DOUBLE) {
		    //cout << lst_key[j+2][i] << ":" << batch_value_param << endl;
		    if (makeDataBase()->GetParamDouble(lst_key[j+2][i], batch_value_param, &argval)) { //j=1 is object
			if (!argval) {
			    Warning("Execute", "Double argument for key %i is NULL", lst_key[j+2][i]);
			    return retval;
			}
			//cout << argval << endl;
			methods[lst_command_int[i]]->SetParam(*argval);
		    } else {
			Warning("Execute","Double argument for key %i not found", lst_key[j+2][i]);
			return retval;
		    }
		} else if (methods_arg_flags[j][lst_command_int[i]] == METHOD_RETURN_INT) {
		    if (makeDataBase()->GetParamDouble(lst_key[j+2][i], batch_value_param, &argval)) { //j=1 is object
			if (!argval) {
			    Warning("Execute", "Int argument for key %i is NULL", lst_key[j+2][i]);
			    return retval;
			}
		    } else {
			Warning("Execute", "Int argument for key %i not found", lst_key[j+2][i]);
			return retval;
		    }
		    methods[lst_command_int[i]]->SetParam((Long_t)*argval);
		} 
	    }

	    //Execute internal
	    if (methods_flags[lst_command_int[i]] == METHOD_RETURN_DOUBLE) {
		makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
		if (!val) return retval;
		methods[lst_command_int[i]]->Execute(a1, *val);
		//cout << "Double_t meth " << method_name[lst_command_int[i]] << " called, result " << *val << endl;
	    } else if (methods_flags[lst_command_int[i]] == METHOD_RETURN_INT) { //INT
		Long_t ret;
		makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
		if (!val) return retval;
		methods[lst_command_int[i]]->Execute(a1, ret);
		//cout << "Int_t meth " << method_name[i] << " called, result " << ret << endl;
		(*val) = (Double_t) ret;
	    } else if (methods_flags[lst_command_int[i]] == METHOD_RETURN_PPARTICLE) { 
		//PParticle* ret2;//, *ret;
		makeDataBase()->GetParamTObj(lst_key_a[i], batch_particle_param, &res);
		if (!res) return retval;
		//ret2  = (PParticle * ) res;
		//methods[lst_command_int[i]]->Execute(a1,ret); //BUGBUG
		Fatal("Execute", "Return type PParticle* not yet working");
		//if (ret) *ret2 = * ret;
	    } else { //VOID
		methods[lst_command_int[i]]->Execute(a1);
		//cout << "void meth " << method_name[lst_command_int[i]] << " called" << endl;
		makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
		if (val) (*val) = 0;
	    }

	} else if (lst_command[i] == COMMAND_PVALUE) {
	    //
	    // Reads the PValue ("PParticle.val")
	    // or the database entry
	    //
	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
	    if (!val) {
		if (!error_flag[i]) {
		    error_flag[i] = 1;
		    Error("Execute", "Result object %s not found",
			      makeDataBase()->GetName(lst_key_a[i]));
		}
		return retval;
	    }	    
	    if (lst_key[3][i] == 0) { //PValue
		makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
		if (!res) return retval;
		PParticle *a1  = (PParticle *) res;
		if (!a1->GetValue(lst_key[2][i], val)) {
		    if (!error_flag[i]) {
			error_flag[i] = 1;
			Error("Execute", "PValue %i not set", lst_key[2][i]);
		    }
		    return retval;
		} 
	    } else { //database entry
		if (lst_key[3][i] < 0) {
		    makeDataBase()->GetParamDouble(lst_key[1][i], (-lst_key[3][i]) - 1, &val2);
		    if (!val2) {
			Error("Execute", "Connection to data base for double %i failed", (-lst_key[3][i]) - 1);
			return retval;
		    }
		    *val = *val2;		    
		} else {
		    Int_t *intval;
		    makeDataBase()->GetParamInt(lst_key[1][i], lst_key[3][i] - 1, &intval);
		    *val = (Double_t)*intval;
		} 
	    }
	} else if (lst_command[i] == COMMAND_BRANCH) {
	    //
	    // Returns the branch index
	    // 
	    //
	    makeDataBase()->GetParamDouble (lst_key_a[i], batch_value_param, &val);
	    if (!val) {
		if (!error_flag[i]) {
		    error_flag[i] = 1;
		    Error("Execute", "Result object %s not found",
			  makeDataBase()->GetName(lst_key_a[i]));
		}
		return retval;
	    }	    
	    *val = lst_key[2][i];
	} else if (lst_command[i] == COMMAND_PFORMULA) {
	    //
	    // Pointer to PFormula
	    //
	    Double_t pars[MAX_COMMAND_OPTIONS];
	    for (int j=0; j<lst_options_counter[i]; j++) {
		makeDataBase()->GetParamDouble(lst_key[j+1][i], batch_value_param, &val);
		if (!val) return retval;
		pars[j] = *val;
	    }
	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
	    if (!val) return retval;
	    lst_form[i]->SetParameters(pars);
	    *val = lst_form[i]->Eval(1.);
	} else if (lst_command[i] == COMMAND_PRINT) {
	    //
	    // Dumps the value(s)
	    //
	    makeDataBase()->GetParamTObj(lst_key[1][i], batch_particle_param, &res);
	    if (res) {
		cout << "*******: " << makeDataBase()->GetName(lst_key[1][i]) << endl;
		PParticle *a1  = (PParticle *) res;
		//a1->Dump();
		a1->Print("scatter");
	    }
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    if (val) {
		cout << "*******: " << makeDataBase()->GetName(lst_key[1][i]) << endl;
		cout << *val << endl;
	    }
	} else if (lst_command[i] == COMMAND_COS) {
	    //
	    // Internal cos()
	    //
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    if (!val) {
		return retval;
	    }
 	    Double_t myres = cos(*val);
	    
 	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);

 	    if (!val) return retval;
 	    *val = myres;
	} else if (lst_command[i] == COMMAND_FABS) {
	    //
	    // Internal fabs()
	    //
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    if (!val) {
		return retval;
	    }
 	    Double_t myres = fabs(*val);
	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
 	    if (!val) return retval;
 	    *val = myres;
	} else if (lst_command[i] == COMMAND_EVAL) {
	    //
	    // Eval the attached histogram and returns the value
	    //
	    Double_t myres  = 0.;
	    TH1 *my_fHisto1 = NULL;
	    TH2 *my_fHisto2 = NULL;
	    TH3 *my_fHisto3 = NULL;
	    TGraph   *my_fGraph   = NULL;
	    TGraph2D *my_fGraph2D = NULL;
	    
	    if (lst_key[1][i] == -1) {
		//no external histo
		my_fHisto1  = fHisto1;
		my_fHisto2  = fHisto2;
		my_fHisto3  = fHisto3;
		my_fGraph   = fGraph;
		my_fGraph2D = fGraph2D;
	    } else {
		TObject *ret;
		if (makeDataBase()->GetParamTObj(lst_key[1][i], batch_histogram_param, &ret)) {
		    if (((TH1*)ret)->GetDimension() == 1)
			my_fHisto1 = (TH1*)ret;
		    if (((TH1*)ret)->GetDimension() == 2)
			my_fHisto2 = (TH2*)ret;
		    if (((TH1*)ret)->GetDimension() == 3)
			my_fHisto3 = (TH3*)ret;
		}
	    }

	    //Eval()
	    if (lst_key[2][i] == -1 && lst_key[3][i] == -1 && lst_key[4][i] == -1) {
		if (!x) {
		    return retval;
		    }
		if (my_fHisto1) {
		    int bin = my_fHisto1->FindBin(*x);
		    myres = my_fHisto1->GetBinContent(bin);
		} else if (my_fHisto2) {
		    int bin = my_fHisto2->FindBin(*x,*y);
		    myres = my_fHisto2->GetBinContent(bin);
		} else if (my_fHisto3) {
		    int bin = my_fHisto3->FindBin(*x,*y,*z);
		    myres = my_fHisto3->GetBinContent(bin);
		} else if (my_fGraph) {
		    myres = my_fGraph->Eval(*x);
		} else if (my_fGraph2D) {
		    myres = my_fGraph2D->Interpolate(*x,*y);
		} else {
		    if (!eval_err_dumped) {
			eval_err_dumped = 1;
			Error("Execute", "Eval() called, but no object present");
		    }
		    return retval;
		}
	    } else  { //Eval(x,y,z)
		if (lst_key[2][i] != -1  && lst_key[3][i] == -1) { //1dim
		    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
		    if (!val2) {
			return retval;
		    }
		    if (my_fHisto1) {
			int bin = my_fHisto1->FindBin(*val2);
			myres = my_fHisto1->GetBinContent(bin);
		    } else if (my_fGraph) {
			myres = my_fGraph->Eval(*x); 
		    } else {
			if (!eval_err_dumped) {
			    eval_err_dumped = 1;
			    Error("Execute", "Eval(x) called, but no TH1 or TGraph object present");
			}
			return retval;
		    }
		} else if (lst_key[2][i] != -1  && lst_key[3][i] != -1 && lst_key[4][i] == -1) { //2dim
		    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
		    makeDataBase()->GetParamDouble(lst_key[3][i], batch_value_param, &val3);
		    if (!val2 || !val3) {
			return retval;
		    }
		    if (my_fHisto2) {
			int bin = my_fHisto2->FindBin(*val2, *val3);
			myres = my_fHisto2->GetBinContent(bin);
		    } else if (my_fGraph2D) {
		    myres = my_fGraph2D->Interpolate(*x, *y);
		    } else {
			if (!eval_err_dumped) {
			    eval_err_dumped = 1;
			    Error("Execute", "Eval(x,y) called, but no TH2 or TGraph2D object present");
			}
			return retval;
		    }
		} else  { //3dim
		    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
		    makeDataBase()->GetParamDouble(lst_key[3][i], batch_value_param, &val3);
		    makeDataBase()->GetParamDouble(lst_key[4][i], batch_value_param, &val4);
		    
		    if (!val2 || !val3 || !val4) {
			return retval;
		    }
		    if (my_fHisto3) {
			int bin = my_fHisto3->FindBin(*val2, *val3, *val4);
			myres = my_fHisto3->GetBinContent(bin);
		    } else {
			if (!eval_err_dumped) {
			    eval_err_dumped = 1;
			    Error("Execute", "Eval(x,y,z) called, but no TH3 object present");
			}
			return retval;
		    }
		}
	    }

 	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);

 	    if (!val) return retval;
 	    *val = myres;
	} else if (lst_command[i] == COMMAND_GETRANDOM) {
	    //
	    // get a random number from the attached histogram
	    //
	    Double_t myres  = 0.;
	    TH1 *my_fHisto1 = NULL;
	    TH2 *my_fHisto2 = NULL;
	    TH3 *my_fHisto3 = NULL;
	    
	    if (lst_key[1][i] == -1) {
		//no external histo
		my_fHisto1 = fHisto1;
		my_fHisto2 = fHisto2;
		my_fHisto3 = fHisto3;
	    } else {
		TObject *ret;
		if (makeDataBase()->GetParamTObj(lst_key[1][i], batch_histogram_param, &ret)) {
		    if (((TH1*)ret)->GetDimension() == 1 )
			my_fHisto1 = (TH1*)ret;
		    if (((TH1*)ret)->GetDimension() == 2 )
			my_fHisto2 = (TH2*)ret;
		    if (((TH1*)ret)->GetDimension() == 3 )
			my_fHisto3 = (TH3*)ret;
		}
	    }

	    //GetRandom()
	    if (lst_key[2][i] == -1 && lst_key[3][i] == -1 && lst_key[4][i] == -1) {
		if (!x) {
		    return retval;
		    }
		if (my_fHisto1) {
		    myres = my_fHisto1->GetRandom();
		} else {
		    if (!eval_err_dumped) {
			eval_err_dumped = 1;
			Error("Execute", "GetRandom() called, but no 1-dimensional histogram object present");
		    }
		    return retval;
		}
	    } else  { //GetRandom(x,y,z)
		if (lst_key[2][i] != -1 && lst_key[3][i] == -1) { //1dim
		    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
		    if (!val2) {
			return retval;
		    }
		    if (my_fHisto1) {
			myres = my_fHisto1->GetRandom();
			*val2 = myres;
		    } else {
			if (!eval_err_dumped) {
			    eval_err_dumped = 1;
			    Error("Execute", "GetRandom(x) called, but no TH1 object present");
			}
			return retval;
		    }
		} else if (lst_key[2][i] != -1 && lst_key[3][i] != -1 && lst_key[4][i] == -1) { //2dim
		    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
		    makeDataBase()->GetParamDouble(lst_key[3][i], batch_value_param, &val3);
		    if (!val2 || !val3) {
			return retval;
		    }
		    if (my_fHisto2) {
			my_fHisto2->GetRandom2(*val2, *val3);
			myres = 0;
		    } else {
			if (!eval_err_dumped) {
			    eval_err_dumped = 1;
			    Error("Execute", "GetRandom(x,y) called, but no TH2 object present");
			}
			return retval;
		    }
		} else  { //3dim
		    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
		    makeDataBase()->GetParamDouble(lst_key[3][i], batch_value_param, &val3);
		    makeDataBase()->GetParamDouble(lst_key[4][i], batch_value_param, &val4);
		    
		    if (!val2 || !val3 || !val4) {
			return retval;
		    }
		    if (my_fHisto3) {
			my_fHisto3->GetRandom3(*val2, *val3, *val4);
			myres = 0;
		    } else {
			if (!eval_err_dumped) {
			    eval_err_dumped = 1;
			    Error("Execute", "GetRandom(x,y,z) called, but no TH3 object present");
			}
			return retval;
		    }
		}
	    }

 	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);

 	    if (!val) return retval;
 	    *val = myres;
	} else if (lst_command[i] == COMMAND_GETRANDOMX) {
	    //
	    // get a random number from the attached histogram
	    //
	    Double_t myres=0.;
	    TH2 * my_fHisto2 = NULL;
	    
	    if (lst_key[1][i] == -1) {
		//no external histo
		my_fHisto2  = fHisto2;
	    } else {
		TObject *ret = NULL;
		if (makeDataBase()->GetParamTObj(lst_key[1][i], batch_histogram_param, &ret)) {
		    if (((TH1*)ret)->GetDimension() == 2)
			my_fHisto2 = (TH2*)ret;
		}
	    }
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    if (!val2) {
		return retval;
	    }
	    if (my_fHisto2) {
		int n_x = my_fHisto2->GetNbinsX();
		int n_y = my_fHisto2->GetNbinsY();
		
		if (!slicesy) { //first time call
		    slicesy = new TH1D*[n_y];
		    for (int i=0; i<n_y; i++) {
			slicesy[i] = new TH1D(Form("%s_y_%i",my_fHisto2->GetName(), i), "Projected", n_x, 
					      my_fHisto2->GetXaxis()->GetXmin(), 
					      my_fHisto2->GetXaxis()->GetXmax());
			slicesy[i]->SetDirectory(0);
			for (int j=0; j<n_x; j++) {
			    slicesy[i]->SetBinContent(j+1, my_fHisto2->GetBinContent(j+1, i+1));
			}
		    }
		}
		int ypos = my_fHisto2->GetYaxis()->FindBin(*val2);
		makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
		if (!val) return retval;
		if (ypos > 0 && ypos <= n_y) {
		    myres = slicesy[ypos-1]->GetRandom();
		    *val  = myres;
		} else {
		    *val  = 0;
		}
	    } else {
		if (!eval_err_dumped) {
		    eval_err_dumped = 1;
		    Error("Execute", "GetRandomX(y) called, but no TH2 object present");
		}
		return retval;
	    }
	} else if (lst_command[i] == COMMAND_GETRANDOMY) {
	    //
	    // get a random number from the attached histogram
	    //
	    Double_t myres  = 0.;
	    TH2 *my_fHisto2 = NULL;
	    
	    if (lst_key[1][i] == -1) {
		//no external histo
		my_fHisto2 = fHisto2;
	    } else {
		TObject *ret = NULL;
		if (makeDataBase()->GetParamTObj(lst_key[1][i], batch_histogram_param, &ret)) {
		    if (((TH1*)ret)->GetDimension() == 2)
			my_fHisto2 = (TH2*)ret;
		}
	    }
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    if (!val2) {
		return retval;
	    }
	    if (my_fHisto2) {
		int n_x = my_fHisto2->GetNbinsX();
		int n_y = my_fHisto2->GetNbinsY();
		
		if (!slicesx) { //first time call
		    slicesx = new TH1D*[n_x];
		    for (int i=0; i<n_x; i++) {
			slicesx[i] = new TH1D(Form("%s_x_%i",my_fHisto2->GetName(), i), "Projected", n_y, 
					      my_fHisto2->GetYaxis()->GetXmin(), 
					      my_fHisto2->GetYaxis()->GetXmax());
			slicesx[i]->SetDirectory(0);
			for (int j=0; j<n_y; j++) {
			    slicesx[i]->SetBinContent(j+1, my_fHisto2->GetBinContent(i+1, j+1));
			}
		    }
		}
		int xpos = my_fHisto2->GetXaxis()->FindBin(*val2);
		makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);
		if (!val) return retval;
		if (xpos > 0 && xpos <= n_x) {
		    myres = slicesx[xpos-1]->GetRandom();
		    *val  = myres;
		} else {
		    *val  = 0;
		}
	    } else {
		if (!eval_err_dumped) {
		    eval_err_dumped = 1;
		    Error("Execute", "GetRandomY(x) called, but no TH2 object present");
		}
		return retval;
	    }
	} else if (lst_command[i] == COMMAND_IF) {
	    //
	    // if (....)
	    //
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    if (!val) {
		Warning("Execute", "if: argument not found");
		return retval;
	    }
 	    Double_t myres = *val;
	    
 	    makeDataBase()->GetParamDouble(lst_key_a[i], batch_value_param, &val);

 	    if (!val) return retval;
 	    *val = myres;
	    if (fabs(*val) == 0) return retval;
	} else if (lst_command[i] == COMMAND_P3M) {
	    //
	    // Constructor for PParticles (with mass)
	    //
	    makeDataBase()->GetParamTObj(lst_key_a[i],    batch_particle_param, &res);
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key[3][i], batch_value_param, &val3);
	    makeDataBase()->GetParamDouble(lst_key[4][i], batch_value_param, &val4);

	    if (!res || !val || !val2 || !val3 || !val4) {
		return kTRUE;
	    }
	    
	    ((PParticle *) res)->SetPxPyPzE((*val),(*val2),(*val3),
					    sqrt((*val)*(*val)+(*val2)*(*val2)+
						 (*val3)*(*val3)+(*val4)*(*val4)));
	    
	} else if (lst_command[i] == COMMAND_P3E) {
	    //
	    // Constructor for PParticles (with energy)
	    //
	    makeDataBase()->GetParamTObj(lst_key_a[i],    batch_particle_param, &res);
	    makeDataBase()->GetParamDouble(lst_key[1][i], batch_value_param, &val);
	    makeDataBase()->GetParamDouble(lst_key[2][i], batch_value_param, &val2);
	    makeDataBase()->GetParamDouble(lst_key[3][i], batch_value_param, &val3);
	    makeDataBase()->GetParamDouble(lst_key[4][i], batch_value_param, &val4);

	    if (!res || !val || !val2 || !val3 || !val4) {
		return kTRUE;
	    }
	    
	    ((PParticle *) res)->SetPxPyPzE((*val),(*val2),(*val3),(*val4));
	    
	} else if (lst_command[i] == COMMAND_GOTO) {
	    //
	    // Goto
	    //
	    Bool_t labelfound = kTRUE;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i], num_command_param, &locnum_command)) 
		labelfound = kFALSE;
	    makeDataBase()->GetParamInt(lst_key[1][i], num_batch_param, &locnum_batch);
	    makeDataBase()->GetParamInt(lst_key[1][i], num_bulk_param,  &locnum_bulk);
	    if ((locnum_command>=0) && (labelfound) && (locnum_batch < 0)) { //local batch
		i = locnum_command-1;
	    }
	    else if((locnum_command>=0) && (labelfound)) { //go back to PProjector
		return kGOTO;
	    }
	    if (!labelfound)
		Warning("Execute", "Label '%s' not found in goto command", makeDataBase()->GetName(lst_key[1][i]));
	} else if (lst_command[i] == COMMAND_GOSUB) {
	    //
	    // Gosub
	    //
	    //space left on stack?
	    if (stack_num_pos == MAX_STACK_GOSUB) {
		Warning("Execute", "Cannot call '%s': stack full (loop?)", makeDataBase()->GetName(lst_key[1][i]));
		return kFALSE;
	    }

	    Bool_t labelfound = kTRUE;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i], num_command_param, &locnum_command)) 
		labelfound = kFALSE;
	    makeDataBase()->GetParamInt(lst_key[1][i], num_batch_param, &locnum_batch);
	    makeDataBase()->GetParamInt(lst_key[1][i], num_bulk_param,  &locnum_bulk);
	    if((locnum_command>=0) && (labelfound)) {
		stack_num_batch[stack_num_pos] = num_batch;
		stack_num_bulk[stack_num_pos]  = num_bulk;
		stack_num_command[stack_num_pos] = i;
		stack_num_pos++;
		if (locnum_batch < 0) { //local batch
		    i = locnum_command-1;
		}
		else  { //go back to PProjector
		    return kGOTO;
		}
	    }
	    if (!labelfound)
		Warning("Execute", "Label '%s' not found in goto command", makeDataBase()->GetName(lst_key[1][i]));
	} else if (lst_command[i] == COMMAND_RETURN) {
	    //
	    // Returns from gosub
	    //
	    //something on stack?
	    if (stack_num_pos == 0) {
		Warning("Execute", "Cannot return: stack empty");
		return kFALSE;
	    }
	    stack_num_pos--;
	    locnum_command = stack_num_command[stack_num_pos]+1;
	    locnum_batch   = stack_num_batch[stack_num_pos];
	    locnum_bulk    = stack_num_bulk[stack_num_pos];
	    if (locnum_batch < 0) { //local batch
		i = locnum_command-1;
	    }
	    else  { //go back to PProjector
		return kGOTO;
	    }
	} else if (lst_command[i] == COMMAND_ELSE) {
	    //
	    // Just return...
	    //
	    return retval | kELSE;
	} else if (lst_command[i] == COMMAND_EXIT) {
	    //
	    // Exit
	    //
	    locnum_command     = 999999999;
	    locnum_old_command = 999999999;
	    locnum_batch       = 999999999;
	    locnum_bulk        = num_bulk;
	    return kGOTO;
	} else if ((lst_command[i] == COMMAND_FORMORE) || (lst_command[i] == COMMAND_FOREACH)) {
	    //
	    // Build-in loops for PProjector
	    //
	    Int_t  stream_max;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i], stream_max_pos_param, &stream_max)) {
		//cout << "particle not in stream" << endl;
		return kFALSE; //particle not in stream, abort
	    }

	    Int_t *stream_default;
	    if (!makeDataBase()->GetParamInt(lst_key[1][i], stream_default_pos_param, &stream_default)) {
		Int_t *dummy = new Int_t(1);
		stream_default = dummy;
		makeDataBase()->SetParamInt(lst_key[1][i], STREAM_DEFAULT_POS, dummy);
	    } else {
		(*stream_default)++;
	    }

	    if (lst_command[i] == COMMAND_FORMORE) {
		if ( (*stream_default) >= stream_max) { //Maximum reached, abort
		    //reset value:
		    (*stream_default) = 0;
		    return kFOREACHEND; 
		} 
	    } else { //FOREACH
		if ((*stream_default) > stream_max) { //Maximum reached, abort
 		    //reset value:
 		    (*stream_default) = 0;
 		    return kFOREACHEND; 
 		} else {
 		    retval = kFOREACH;
 		}
		if (locnum_old_command < 0) locnum_old_command = i; //Jump to the "foreach" command
	    }

	} else if (lst_command[i] == COMMAND_READLINE) {

	    if (readline_format_string[i]) {
		if (!file) {
		    Error("Execute", "readline called, but no file open");
		    return kFALSE;
		}
		char line[1000];
		if (!fgets(line, 1000, file)) return kEOF;
		Double_t my_readline_args[MAX_READLINE_ARGS];
		int r = sscanf(line, readline_format_string[i],
			       my_readline_args,   my_readline_args+1,
			       my_readline_args+2, my_readline_args+3,
			       my_readline_args+4, my_readline_args+5,
			       my_readline_args+6);

		if (r != readline_num_args[i]) {
		    Warning("Execute", "%i arguments expected, but only %i read", readline_num_args[i], r);
		    cout << "< " << readline_string[i] << endl;
		    cout << "> " << line;
		}
		for (int j=0; j<readline_num_args[i]; j++) {
		    *(*(readline_args[i]+j)) = my_readline_args[j];
		}
	    }

	} else if (lst_command[i] == COMMAND_ECHO) {
	    //
	    // echo blabla
	    //
	    char puffer[1000]; //I hope there will be never a string longer...
	    char format_puffer[1000]; //I hope there will be never a string longer...
	    unsigned int puffer_pointer = 0, format_puffer_pointer = 0;

	    if (!file)
		printf("<PBatch> "); 
	    Int_t seek_mode = 0;
	    for (unsigned int j=0; j<=strlen(echo_string[i]); j++) {
		char current = (*(echo_string[i]+j));
		if (current == '$')
		    seek_mode = 1;
		else {
		    //if (seek_mode && (((*(echo_string[i]+j)) == ' ')  ||  (*(echo_string[i]+j)) == '\0') ) {
		    if ((seek_mode==1) && (current!='%') && 
			(current!='#') && (current!='_') && (!isalnum(current))) {
			//end
			seek_mode = 0;
			puffer[puffer_pointer] = '\0';
			Double_t *x = makeStaticData()->GetBatchValue(puffer, 0);
			if (file) {
			    if (x) {
				if (format_puffer_pointer<2) 
				    fprintf(file, "%#g", *x);
				else {
				    format_puffer[format_puffer_pointer] = '\0';
				    char type = format_puffer[format_puffer_pointer-1];
				    if ((type == 'd') ||
					(type == 'i') ||
					(type == 'o') ||
					(type == 'x') ||
					(type == 'X')||
					(type == 'u'))
					fprintf(file, format_puffer, (int)*x);
				    else if (type == 'c')
					fprintf(file, format_puffer, (char)*x);
				    else
					fprintf(file, format_puffer, *x);			    
				}
			    }
			    else fprintf(file, "-");
			} else { 
			    if (x) {
				if (format_puffer_pointer<2) printf("%#g", *x);
				else {
				    format_puffer[format_puffer_pointer] = '\0';
				    char type = format_puffer[format_puffer_pointer-1];
				    if ((type == 'd') ||
					(type == 'i') ||
					(type == 'o') ||
					(type == 'x') ||
					(type == 'X')||
					(type == 'u'))
					printf(format_puffer, (int)*x);
				    else if (type == 'c')
					printf(format_puffer, (char)*x);
				    else
					printf(format_puffer, *x);
				}
			    }
			    else printf("[Error: %s not found]", puffer);
			}
			
			puffer_pointer = 0;
			format_puffer_pointer = 0;
		    }
		    if ((seek_mode==2) && (current=='%')) { //end format string
			seek_mode = 1;
		    }
		    else if ((seek_mode==1) && (current=='%')) { //found format string
			seek_mode = 2;
		    }
		    else if (seek_mode == 1) {
			puffer[puffer_pointer] = current;
			if (puffer_pointer < 1000) puffer_pointer++;
			else {
			    if (!file)
				printf("[Error: variable too large]");
			    seek_mode = 0;
			}
		    }
		    if (seek_mode == 2) {
			format_puffer[format_puffer_pointer] = current;
			
			if (format_puffer_pointer < 1000) format_puffer_pointer++;
			else {
			    if (!file)
				printf("[Error: format too large]");
			    seek_mode = 1;
			}
		    }
		    
		    if (!seek_mode) {
			if (file && (current != '\0')) fprintf(file,"%c",current);
			else printf("%c",current);
		    }
		}
	    }
	    if (file) 
		fprintf(file, "\n");
	    else {
		if (seek_mode == 2) printf("[Error: format string not closed]");
		printf("\n");
	    }
	}
    }

    return kTRUE | retval ;
}

Bool_t PBatch::AddCommand(const char *_command) {

    //cout << "AddCommand1:" << command << endl;
    char *command = PUtils::NewString(_command);
    PUtils::remove_spaces(&command);
    if (strlen(command) == 0) return kTRUE;
    Bool_t has_something  = kFALSE;
    int    curly_brackets = 0;

    for (unsigned int i=0; i<strlen(command); i++) {
	if ((command[i] != ' ') && (command[i] != ';'))
	    has_something = kTRUE;
	if (command[i]  == '{') curly_brackets++;
	if (command[i]  == '}') curly_brackets--;
	if ((command[i] == ';') && curly_brackets) command[i] = '\"'; //placeholder dummy
    }
    if (!has_something) return kTRUE;

    //cout << "AddCommand:" << command << endl;
    //First check if we have a composite command

    int is_composite = 0;

    for (UInt_t i=0; i<strlen(command); i++) {
	if (command[i]==';') is_composite = 1;
    }

    if (is_composite) {
	char *array[200];
	Int_t array_s = 200; //max products
	PUtils::Tokenize(command, ";", array, &array_s);

	for (int i=0; i<array_s; i++) {
	    AddCommand(array[i]);
	}
	return kTRUE;
    }

    //single command:

    //echo has the highest priority:
    if (!strncmp(command, "echo", 4)) {
	// 	    key_a = makeStaticData()->
	// 		MakeDirectoryEntry("batch_objects",command3);
	if (strlen(command) > 4) {
	    char *dummy = new char[strlen(command)-2];
	    strncpy(dummy, command+4, strlen(command)-4);
	    dummy[strlen(command)-4] = '\0';
	    if (tmp_file) file = tmp_file;
	    PUtils::remove_spaces(&dummy);
	    echo_string[command_pointer] = dummy;
	} else 
	    echo_string[command_pointer] = new char('\0');
	AddCommand(COMMAND_ECHO, -1, -1, -1);
	return kTRUE;
    }    
    
    //readline
    if (!strncmp(command, "readline", 8)) {
	if (strlen(command) > 8) {
	    char *dummy = new char[strlen(command)-6];
	    strncpy(dummy, command+8, strlen(command)-8);
	    dummy[strlen(command)-8] = '\0';
	    PUtils::remove_brackets(&dummy, '{', '}');
	    for (unsigned int i=0; i<strlen(dummy); i++) {
		if (dummy[i] == '\"') dummy[i] = ';';
	    }
	    if (tmp_file) file = tmp_file;
	    char puffer[1000]; //I hope there will be never a string longer...
	    unsigned int puffer_pointer = 0, dummy2_pointer = 0, args_pointer = 0;
	    char *dummy2 = new char[strlen(dummy)+20];
	    readline_args[command_pointer] = new Double_t*[MAX_READLINE_ARGS];
	    
	    Int_t seek_mode = 0;
	    for (unsigned int j=0; j<=strlen(dummy); j++) {
		char current = (*(dummy+j));
		if (current == '@') {
		    seek_mode = 1;
		    dummy2[dummy2_pointer]   = '%';
		    dummy2[dummy2_pointer+1] = 'l';
		    dummy2[dummy2_pointer+2] = 'e';
		    //dummy2[dummy2_pointer+3] = ' ';
		    dummy2_pointer += 3; 
		    if (dummy2_pointer > (strlen(dummy)+20)) {
			Fatal("AddCommand", "strlen for tmp string too short");
		    }
		} else {
		    if (seek_mode && ((!(current=='#') and !(current=='_') and !isalnum(current)))) {
			//end
			seek_mode = 0;
			puffer[puffer_pointer] = '\0';
			Info("Readline", "Write to batch variable '%s'", puffer);
			char *dummy3 = new char[strlen(puffer)+1];
			strcpy(dummy3, puffer);
			
			Double_t *x = makeStaticData()->GetBatchValue(dummy3, 1);
			if (x) {
			    if (args_pointer < MAX_READLINE_ARGS) {
				//push pointer on the list
				*((readline_args[command_pointer])+args_pointer) = x;
				args_pointer++;
			    } else {
				Error("AddCommand", "Readline: %s too many variables (>MAX_READLINE_ARGS)", dummy);
			    }
			} else 
			    Error("AddCommand", "Readline: puffer %s not initialized", puffer);
			    
			puffer_pointer = 0;
		    }
		    if (seek_mode) {
			puffer[puffer_pointer] = current;
			if (puffer_pointer < 1000) puffer_pointer++;
			else {
			    Error("AddCommand", "Readline: %s too large", dummy);
			    seek_mode = 0;
			}
		    }
		    
		    if (!seek_mode) {
			dummy2[dummy2_pointer] = dummy[j];
			dummy2_pointer++;
		    }
		}
	    }

	    readline_num_args[command_pointer] = args_pointer;
	    readline_string[command_pointer]   = dummy;
	    readline_format_string[command_pointer] = dummy2;
	} else {
	    readline_format_string[command_pointer] = new char('\0');
	}
	AddCommand(COMMAND_READLINE, -1, -1, -1);
	return kTRUE;
    }

    //first, check for a label
    for (UInt_t i=1; i<strlen(command); i++) {
	
	Bool_t ret = kFALSE;
	if (command[i] == ':') {
	    if (((i == (strlen(command))) ||  (command[i+1] != ':'))
		&&  ((i == 0) ||  (command[i-1] != ':'))) {
		//add new label
		char *label_name = new char[i+1];
		strncpy(label_name, command, i);
		label_name[i] = '\0';
		int key_a = makeStaticData()->
		    MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command);
		if (key_a >= 0) {
		    Int_t *delme = new Int_t(command_pointer);
		    AddCommand(COMMAND_LABEL, key_a, -1, -1);
		    ret = AddCommand(command+i+1);
		    int key_l = GetKey(label_name, 1, 1);
		    if (key_l > 0) {
			makeDataBase()->SetParamInt(key_l, "num_command", delme);
			delme = new Int_t(num_batch);
			makeDataBase()->SetParamInt(key_l, "num_batch", delme);
			delme = new Int_t(num_bulk);
			makeDataBase()->SetParamInt(key_l, "num_bulk", delme);
		    }
		}
		return ret;
	    }
	}
    }

    //second thing: I check for a "=" operator

    int is_operator = 0;
    Bool_t found = kFALSE;
    char *prod[2];
    Int_t prod_s = 2; //max 2 products
    int key_a=-1, key1=-1, key2=-1;

    int numbrack = 0;
    for (UInt_t i=0; i<strlen(command); i++) {
	if (command[i]=='(' || command[i]==')') numbrack++;
	if (command[i]=='=' && !numbrack) {
	    if (i>0 && i<(strlen(command)-1) 
		&& command[i-1]!='<' && command[i-1]!='>' 
		&& command[i-1]!='*' && command[i-1]!='/' 
		&& command[i-1]!='-' && command[i-1]!='+' 
		&& command[i-1]!='=' && command[i+1]!='=') is_operator++;
	}
    }

    //cout << is_operator << endl;

    //let us copy the command first, because it might be modified
    //after the pointer is given to the data base
    char *command2 = new char[strlen(command)+1];
    strcpy(command2, command);
    command = command2;
    char *command3 = new char[strlen(command)+1];
    strcpy(command3, command);
    Int_t copy_ctor = 0;
    if (!is_operator) {
	//first the functions
	if (!strncmp(command,"goto",4) || !strncmp(command,"Goto",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    char *dummy = new char[strlen(command)-3];
	    strncpy(dummy, command+4, strlen(command)-4);
	    dummy[strlen(command)-4] = '\0';
	    key1 = GetKey(dummy, 1, 1);
	    if (AddCommand(COMMAND_GOTO, key_a, key1, -1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"gosub",5) || !strncmp(command,"Gosub",5)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    char *dummy = new char[strlen(command)-3];
	    strncpy(dummy, command+5, strlen(command)-5);
	    dummy[strlen(command)-5] = '\0';
	    key1 = GetKey(dummy, 1, 1);
	    if (AddCommand(COMMAND_GOSUB, key_a, key1, -1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"return",6) || !strncmp(command,"Return",6)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    if (AddCommand(COMMAND_RETURN, key_a, key1, -1)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"else",4) || !strncmp(command,"Else",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    if (AddCommand(COMMAND_ELSE, key_a, -1, -1)) {
		if (else_position > -1) {
		    Error("AddCommand", "Multiple use of 'else'");
		    return kFALSE;
		}
		else_position = command_pointer;
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		AddCommand(command+4); //just to make sure is the ';' is forgotten...
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"exit",4) || !strncmp(command,"Exit",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    if (AddCommand(COMMAND_EXIT, key_a, key1, -1)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"branch",6) || !strncmp(command,"Branch",6)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    char *dummy = new char[strlen(command)-5];
	    strncpy(dummy, command+6, strlen(command)-6);
	    dummy[strlen(command)-6] = '\0';
	    key1 = GetKey(dummy, 1, 1);
	    if (AddCommand(COMMAND_BRANCH, key_a, key1, -1)) {
		//check if we use this key already
		int used_id = -1;
		if (size_branches && key_branches) {
		    for (int br=0; br<*size_branches; br++) {
			if (key_branches[br] == key1) used_id = br;
		    }
		} else {
		    Error("AddCommand", "[%s] No access to branches (batch not attached to PReaction)", command);
		}
		if (used_id == -1) {
		    if (*size_branches != MAX_NUM_BRANCHES) {
			key_branches[*size_branches] = key1;
			used_id = (*size_branches);
			(*size_branches)++;
		    } else {
			Error("AddCommand", "[%s] MAX_NUM_BRANCHES reached", command);
		    }
		}
		lst_key[2][command_pointer-1] = used_id+1;  //new branches start with 1
		//    lst_key[3][command_pointer]=key3;
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"formore",7) || !strncmp(command,"Formore",7)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    char *dummy = new char[strlen(command)-6];
	    strncpy(dummy, command+7, strlen(command)-7);
	    dummy[strlen(command)-7] = '\0';
	    if (strcmp(dummy, "(*)"))
		key1 = GetKey(dummy, 1, 1);
	    else
		key1 = GetKey((char *)"dummy", 1, 1);
	    if (AddCommand(COMMAND_FORMORE, key_a, key1, -1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"foreach",7) || !strncmp(command,"Foreach",7)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    char * dummy = new char[strlen(command)-6];
	    strncpy(dummy,command+7,strlen(command)-7);
	    dummy[strlen(command)-7]='\0';
	    if (strcmp(dummy,"(*)"))
		key1 = GetKey(dummy,1,1);
	    else
		key1 = GetKey((char *)"dummy",1,1);
	    if (AddCommand(COMMAND_FOREACH,key_a,key1,-1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"push",4) || !strncmp(command,"Push",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    char *dummy = new char[strlen(command)-3];
	    strncpy(dummy, command+4, strlen(command)-4);
	    dummy[strlen(command)-4] = '\0';
	    key1 = GetKey((char *)dummy, 1, 1);
	    TObject *delme2 = NULL;
	    makeDataBase()->GetParamTObj(key1, batch_particle_param, &delme2);
	    if (!delme2) {
		Warning("AddCommand", "[%s]: the object you want to push is not a particle", command);
		return kFALSE;
	    }
	    if (AddCommand(COMMAND_PUSH, key_a, key1, -1)) 
		return kTRUE;
	    return kFALSE;
	}
	if (!strncmp(command,"eval()",6) || !strncmp(command,"Eval()",6)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    if (AddCommand(COMMAND_EVAL, key_a, -1, -1)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		is_readonly = 1;
		return kTRUE;
	    }
	    return kFALSE;
	} else if (!strncmp(command,"eval(",5) || !strncmp(command,"Eval(",5)) {
	    //Exists also as method!!!
	    char *prodx[10];
	    Int_t prodx_s = 10; //max 10 products (but problems catched below)

	    PUtils::Tokenize(command+5, ",", prodx, &prodx_s);
	    if (prodx_s > 3) {
		Warning("AddCommand", "Maximum 3 arguments for Eval()");
	    } 
	    *(prodx[prodx_s-1]+strlen(prodx[prodx_s-1])-1) = '\0'; //remove trailing )
	    
	    int key3=-1, key4=-1;
	    if (prodx_s > 0) 
		key2 = GetKey(prodx[0], 1, -1);
	    if (prodx_s > 1) 
		key3 = GetKey(prodx[1], 1, -1);
	    if (prodx_s > 2) 
		key4 = GetKey(prodx[2], 1, -1);
	    
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    
	    if (AddCommand(COMMAND_EVAL, key_a, -1, key2, key3, key4)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		is_readonly = 1;
		return kTRUE;
	    }
	    return kFALSE;
	} 

	if (!strncmp(command,"GetRandom()",11) || !strncmp(command,"getrandom()",11)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    if (AddCommand(COMMAND_GETRANDOM, key_a, -1, -1)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		is_readonly = 1;
		return kTRUE;
	    }
	    return kFALSE;
	} else if (!strncmp(command,"GetRandom(",10) || !strncmp(command,"getrandom(",10)) {
	    //Exists also as method!!!
	    char *prodx[10];
	    Int_t prodx_s = 10; //max 10 products (but problems catched below)
	    
	    PUtils::Tokenize(command+10, ",", prodx, &prodx_s);
	    if (prodx_s > 3) {
		Warning("AddCommand", "Maximum 3 arguments for GetRandom()");
	    } 
	    *(prodx[prodx_s-1]+strlen(prodx[prodx_s-1])-1) = '\0'; //remove trailing )
	    
	    int key3=-1, key4=-1;
	    if (prodx_s > 0) 
		key2 = GetKey(prodx[0], 1, -1);
	    if (prodx_s > 1) 
		key3 = GetKey(prodx[1], 1, -1);
	    if (prodx_s > 2) 
		key4 = GetKey(prodx[2], 1, -1);
	    	    
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    
	    if (AddCommand(COMMAND_GETRANDOM, key_a, -1, key2, key3, key4)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		is_readonly = 1;
		return kTRUE;
	    }
	    return kFALSE;
	} 

	if (!strncmp(command,"GetRandomX",10) || !strncmp(command,"getrandomx",10)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    key2 = GetKey(command+10, 1, -1);
	    if (AddCommand(COMMAND_GETRANDOMX, key_a, -1, key2)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		is_readonly = 1;
		return kTRUE;
	    }
	    return kFALSE;
	} 

	if (!strncmp(command,"GetRandomY",10) || !strncmp(command,"getrandomy",10)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    key2 = GetKey(command+10, 1, -1);
	    if (AddCommand(COMMAND_GETRANDOMY, key_a, -1, key2)) {
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		is_readonly = 1;
		return kTRUE;
	    }
	    return kFALSE;
	} 

	if (!strncmp(command,"cos(",4)) {
	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    key1 = GetKey(command+3, 1, -1);
	    if (AddCommand(COMMAND_COS, key_a, key1, -1)) {
		
		//final thing is to create the result
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	if (!strncmp(command,"fabs(",5)) {

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);

	    key1 = GetKey(command+4, 1, -1);
	    if (AddCommand(COMMAND_FABS, key_a, key1, -1)) {
		//final thing is to create the result
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}
	

	if (!strncmp(command,"if(",3) || !strncmp(command,"if ",3) || !strncmp(command,"ifx",3)) {

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	   
	    //check if the next char after () is a ';':
	    int bracket_pos=-1, bracket_pos2=-1, bracket_num=0, found_trailer=0;
	    for (unsigned int i=2; i<strlen(command); i++) {
		//cout << command[i] << ":" << bracket_num  << ":" << bracket_pos<< endl;
		if (command[i] != '(' && command[i] != ' ' && !bracket_num) i=strlen(command); //exit
		if (command[i] == '(' && bracket_pos == -1) {
		    bracket_pos = i;
		    bracket_num++;
		} else if (command[i] == '(') bracket_num++;

		if (command[i] == ')' && bracket_pos != -1 && bracket_pos2 == -1 && bracket_num==1) {
		    bracket_pos2 = i;
		} else if (bracket_pos2 != -1 && command[i] != ' ') found_trailer=1;
		else if (command[i] == ')') bracket_num--;
	    }

	    //cout << command << endl;
	    if (bracket_pos != -1 && bracket_pos2 != -1 && found_trailer) {
		command[bracket_pos2] = ';';
		command[bracket_pos] = ' ';
		command[2] = 'x';
		return AddCommand(command);
	    }
	    if (!strncmp(command,"if(",3))
		key1 = GetKey(command+2, 1, -1);
	    else
		key1 = GetKey(command+3, 1, -1);
	    if (AddCommand(COMMAND_IF, key_a, key1, -1)) {
		//final thing is to create the result
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		return kTRUE;
	    }
	    return kFALSE;
	}

	if (!strncmp(command,"P3M(",4)) {
	    char *prodx[4];
	    Int_t prodx_s = 4; //max 4 products

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);

	    PUtils::Tokenize(command+4, ",", prodx, &prodx_s);
	    if (prodx_s < 4) {
		Warning("AddCommand", "P3M needs 4 arguments");
	    }

	    *(prodx[3]+strlen(prodx[3])-1) = '\0'; //remove trailing )

	    key1 = GetKey(prodx[0], 1, -1);
	    key2 = GetKey(prodx[1], 1, -1);
	    int key3 = GetKey(prodx[2], 1, -1);
	    int key4 = GetKey(prodx[3], 1, -1);

	    if (AddCommand(COMMAND_P3M, key_a, key1, key2, key3, key4)) {
		//adding objects, if needed
		TObject *delme = NULL;
		makeDataBase()->GetParamTObj(key_a, batch_particle_param, &delme);
		if (!delme) {
		    delme = (TObject *) (new PParticle(0,0,0,0));
		    makeDataBase()->SetParamTObj (key_a, "batch_particle", delme);
		}
		Double_t *val = NULL;
		makeDataBase()->GetParamDouble(key1, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key1, "batch_value", val);
		}
		val = NULL;
		makeDataBase()->GetParamDouble(key2, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key2, "batch_value", val);
		}
		val = NULL;
		makeDataBase()->GetParamDouble(key3, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key3, "batch_value", val);
		}
		val = NULL;
		makeDataBase()->GetParamDouble(key4, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key4, "batch_value", val);
		}
		
		return kTRUE;
	    }
	    return kFALSE;
	}


	if (!strncmp(command,"P3E(",4)) {
	    char *prodx[4];
	    Int_t prodx_s = 4; //max 4 products

	    key_a = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);

	    PUtils::Tokenize(command+4, ",", prodx, &prodx_s);
	    if (prodx_s < 4) {
		Warning("AddCommand", "P3E needs 4 arguments");
	    }

	    *(prodx[3]+strlen(prodx[3])-1) = '\0'; //remove trailing )

	    key1 = GetKey(prodx[0], 1, -1);
	    key2 = GetKey(prodx[1], 1, -1);
	    int key3 = GetKey(prodx[2], 1, -1);
	    int key4 = GetKey(prodx[3], 1, -1);

	    if (AddCommand(COMMAND_P3E, key_a, key1, key2, key3, key4)) {
		//adding objects, if needed
		TObject *delme = NULL;
		makeDataBase()->GetParamTObj(key_a, batch_particle_param, &delme);
		if (!delme) {
		    delme = (TObject *) (new PParticle(0,0,0,0));
		    makeDataBase()->SetParamTObj(key_a, "batch_particle", delme);
		}
		Double_t *val = NULL;
		makeDataBase()->GetParamDouble(key1, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key1, "batch_value", val);
		}
		val = NULL;
		makeDataBase()->GetParamDouble(key2, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key2, "batch_value", val);
		}
		val = NULL;
		makeDataBase()->GetParamDouble(key3, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key3, "batch_value", val);
		}
		val = NULL;
		makeDataBase()->GetParamDouble(key4, batch_value_param, &val);
		if (!val) {
		    val = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key4, "batch_value", val);
		}

		return kTRUE;
	    }
	    return kFALSE;
	}

	//Try operators like *= */ *- *+
	if (strstr(command, "*=")) {
	    PUtils::Tokenize(command, "*=", prod, &prod_s);
	    key1 = GetKey(prod[0], 1, -1);
	    key2 = GetKey(prod[1], 1, -1);
	    if (AddCommand(COMMAND_MULT, key1, key1, key2)) 
		return kTRUE;
	    return kFALSE;
	}
	if (strstr(command, "/=")) {
	    PUtils::Tokenize(command, "/=", prod, &prod_s);
	    key1 = GetKey(prod[0], 1, -1);
	    key2 = GetKey(prod[1], 1, -1);
	    if (AddCommand(COMMAND_DIV, key1, key1, key2)) 
		return kTRUE;
	    return kFALSE;
	}
	if (strstr(command, "+=")) {
	    PUtils::Tokenize(command, "+=", prod, &prod_s);
	    key1 = GetKey(prod[0], 1, -1);
	    key2 = GetKey(prod[1], 1, -1);
	    if (AddCommand(COMMAND_PLUS, key1, key1, key2)) 
		return kTRUE;
	    return kFALSE;
	}
	if (strstr(command, "-=")) {
	    PUtils::Tokenize(command, "-=", prod, &prod_s);
	    key1 = GetKey(prod[0], 1, -1);
	    key2 = GetKey(prod[1], 1, -1);
	    if (AddCommand(COMMAND_MINUS, key1, key1, key2)) 
		return kTRUE;
	    return kFALSE;
	}

	//Methods below here...
	//Only methods which are
	//not embedded in brackets are taken into account
	Int_t dot_version = 0, prev_dot_version = 0;
	Int_t method_position = -1;
	Int_t brack_counter = 0, total_brack_counter = 0;
	for (unsigned int i=0; i<(strlen(command)-1); i++) {
	    if (command[i] == '(') {
		brack_counter++;
		if (prev_dot_version)
		    total_brack_counter++;
	    }
	    if (command[i] == ')') brack_counter--;
	    if (brack_counter == 0) {
		if (command[i]=='.' && isalpha(command[i+1])) {
		    method_position = i;
		    dot_version     = 1;
		}
		if (command[i]=='-' && command[i+1]=='>') {
		    method_position = i;
		    dot_version     = 2;
		}
		if (total_brack_counter && prev_dot_version && ((command[i]!=' ') || (command[i]!=';')))
		    dot_version = 0; //must be something different
		prev_dot_version = dot_version;
	    }
	}

	brack_counter = 0;
	if (dot_version > 0) {
	    for (int i = method_position-1; i>=0; i--) {
		if (command[i]=='(' || command[i]=='[' || command[i]=='{') 
		    brack_counter--;
		if (command[i]==')' || command[i]==']' || command[i]=='}') 
		    brack_counter++;
		if (brack_counter==0 && ((command[i]==' ') || (command[i]=='+') || (command[i]=='-')))
		    dot_version = 0; //cases like "a + b.x"
	    }
	}

	if ((dot_version>0) && (!found)) {
	    //catch the cases where comparators are used
	    for (unsigned int i=0; i< (unsigned int)method_position; i++) {
		if ((command[i]=='~') || (command[i]=='=') || (command[i]=='<')  || (command[i]=='>'))
		    dot_version = 0;
	    }
	}

	if ((dot_version>0) && (!found)) {
	    char *prodx[2];
	    prodx[0] = new char[method_position+1];
	    strncpy(prodx[0], command, method_position);
	    prodx[0][method_position] = '\0';
	    
 	    key_a = makeStaticData()->
 		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command);

	    if (dot_version == 1) {
		prodx[1] = new char[strlen(command)-method_position];
		strncpy(prodx[1], command+method_position+1, strlen(command)-method_position-1);
		prodx[1][strlen(command)-method_position-1] = '\0';
	    } else {
		prodx[1] = new char[strlen(command)-method_position-1];
		strncpy(prodx[1], command+method_position+2, strlen(command)-method_position-2);
		prodx[1][strlen(command)-method_position-2] = '\0';
		
	    } 

	    if (!strcmp(prodx[1], "M2()")) {
		key1 = GetKey(prodx[0], 1, -1);
		if (key1 < 0) return kFALSE;
		if (AddCommand(COMMAND_MASS2, key_a, key1, -1)) 
		    found = kTRUE;
	    }
	    if (!strcmp(prodx[1], "M()")) {
		key1 = GetKey(prodx[0], 1, -1);
		if (AddCommand(COMMAND_MASS, key_a, key1, -1)) 
		    found = kTRUE;
	    }
	    if (!strcmp(prodx[1], "Theta()")) {
		key1 = GetKey(prodx[0], 1, -1);
		if (AddCommand(COMMAND_THETA, key_a, key1, -1)) 
		    found = kTRUE;
	    }
	    if (!strcmp(prodx[1], "Print()")) {
		key1 = GetKey(prodx[0], 1, -1);
		if (AddCommand(COMMAND_PRINT, key_a, key1, -1)) 
		    found = kTRUE;
	    }
// 	    if (!strcmp(prodx[1],"Parent()")) {
// 		key1 = GetKey(prodx[0],1,0);
// 		if (AddCommand(COMMAND_PARENT,key_a,key1,-1)) 
// 		    found=kTRUE;
// 	    }
// 	    if (!strcmp(prodx[1],"ID()")) {
// 		key1 = GetKey(prodx[0],1,0);
// 		if (AddCommand(COMMAND_ID,key_a,key1,-1)) 
// 		    found=kTRUE;
// 	    }

	    //with arguments:
	    if (!strncmp(prodx[1], "Boost(", 6)) {
		key1 = GetKey(prodx[0],   1, -1);
		key2 = GetKey(prodx[1]+5, 1, -1);
		if (AddCommand(COMMAND_BOOST, key_a, key1, key2)) 
		    found = kTRUE;
	    }
	    if (!strncmp(prodx[1], "Angle(", 6)) {
		key1 = GetKey(prodx[0],   1, -1);
		key2 = GetKey(prodx[1]+5, 1, -1);
		if (AddCommand(COMMAND_ANGLE, key_a, key1, key2)) 
		    found = kTRUE;
	    }
	    if (!strncmp(prodx[1], "Rot(", 4)) {
		key1 = GetKey(prodx[0],   1, -1);
		key2 = GetKey(prodx[1]+3, 1, -1);
		if (AddCommand(COMMAND_ROT, key_a, key1, key2)) 
		    found = kTRUE;
	    }
	    if (!strncmp(prodx[1], "GetBeam(", 8)) {
		key1 = GetKey(prodx[0], 1, -1);
		if (AddCommand(COMMAND_GETBEAM, key_a, key1, -1)) 
		    found = kTRUE;
		TObject *delme = NULL;
		makeDataBase()->GetParamTObj(key_a, batch_particle_param, &delme);
		if (!delme) {
		    delme = (TObject *) (new PParticle(0,0,0,0));
		    makeDataBase()->SetParamTObj(key_a, "batch_particle", delme);
		}
	    }
	    if (!strncmp(prodx[1], "GetTarget(", 10)) {
		key1 = GetKey(prodx[0], 1, -1);
		if (AddCommand(COMMAND_GETTARGET, key_a, key1, -1)) 
		    found = kTRUE;
		TObject *delme = NULL;
		makeDataBase()->GetParamTObj(key_a, batch_particle_param, &delme);
		if (!delme) {
		    delme = (TObject *) (new PParticle(0,0,0,0));
		    makeDataBase()->SetParamTObj(key_a, "batch_particle", delme);
		}
	    }
	    if (!strncmp(prodx[1],"Push(",5) || !strncmp(prodx[1],"push(",5)) {
		key1 = GetKey(prodx[0],   1, -1);
		key2 = GetKey(prodx[1]+4, 1, -1);
		if (AddCommand(COMMAND_PUSH, key_a, key1, key2)) 
		    found = kTRUE;
	    }
	    if (!strncmp(prodx[1],"Eval(",5) || !strncmp(prodx[1],"eval(",5)) {
		key1 = GetKey(prodx[0], 1, -1);
		char *prodxx[4];
		Int_t prodxx_s = 4; //max 4 products (see test below)

		PUtils::Tokenize(prodx[1]+5, ",", prodxx, &prodxx_s);
		if (prodxx_s > 3) {
		    Warning("AddCommand", "[%s] Too many arguments for ->Eval()", command);
		} 

		*(prodxx[prodxx_s-1]+strlen(prodxx[prodxx_s-1])-1) = '\0'; //remove trailing )
		
		int key3=-1, key4=-1;
		if (prodxx_s > 0) 
		    key2 = GetKey(prodxx[0], 1, -1);
		if (prodxx_s > 1) 
		    key3 = GetKey(prodxx[1], 1, -1);
		if (prodxx_s > 2) 
		    key4 = GetKey(prodxx[2], 1, -1);
		//check if histogram is there...
		TObject *ret;
		if (!makeDataBase()->GetParamTObj(key1, batch_histogram_param, &ret)) {
		    Error("AddCommand", "[%s] Histogram '%s' not found in ->Eval()", command, prodx[0]);
		}
		if (((TH1*)ret)->GetDimension() !=prodxx_s) {
		    Error("AddCommand", "[%s] Histogram '%s' in ->Eval() has dimension %i, but it should have %i",
			  command, prodx[0], ((TH1*)ret)->GetDimension(), prodxx_s);
		}

		if (AddCommand(COMMAND_EVAL, key_a, key1, key2, key3, key4)) {
		    Double_t *delme = new Double_t(0.);
		    makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		}
		found = kTRUE;
	    }

	    //If there are no (), it could be an internal PValue
	    Int_t found_brackets = 0;
	    for (unsigned int j=0; j<strlen(prodx[1]); j++) {
		if (((prodx[1])[j] == ')') || ((prodx[1])[j] == '('))
		    found_brackets++;
	    }

	    if (!found_brackets) {
		int val_id = pdummy.StringToValueID(prodx[1]);
		int database_id = 0;
		if (val_id < 0) {
		    //check for database entry
		    //positive: param_int (+1)
		    //negative: param_double (-1)
		    database_id = makeDataBase()->GetParamInt(prodx[1]) -
			makeDataBase()->GetParamDouble(prodx[1]);
		    if (database_id == 0)
			Error("AddCommand", "[%s] The value %s is unknown", command, prodx[1]);
		} 
		if ((val_id>=0) || database_id){		    
		    if (strcmp(prodx[0],"*") == 0) { //dummy *
			key1 = GetKey((char *)"dummy", 1, 0);
		    } else {
			key1 = GetKey(prodx[0], 1, 0);
		    }
		    if (key1 >= 0) {
			AddCommand(COMMAND_PVALUE, key_a, key1, val_id, database_id);
			found = kTRUE;
		    } 
		}
	    }

	    //Try to get build-in method
	    if (!found) {
		if (*(prodx[0]) == '{') { //for pchannelmodels:
		    PUtils::remove_brackets(&prodx[0], '{', '}');
		    key1 = makeDataBase()->GetEntry(prodx[0]);
		    if (key1 >= 0) {
			Int_t handle = GetMethodHandle(prodx[1], 2);
			if (handle>=0) {	
			    //cout << "key1:" << key1 << " arg1:" << arg1 << " handle:" << handle << endl;
			    if (AddCommand(COMMAND_INTERNAL, key_a, key1, arg1, arg2, arg3, arg4))  {
				found = kTRUE;
				lst_command_int[command_pointer-1]  = handle;
				flag_command_int[command_pointer-1] = 2;
			    }
			} else {
			    Error("AddCommand", "[%s] The method [%s] is unknown", command, prodx[1]);
			    return kFALSE;
			}
		    } else {
			Error("AddCommand", "[%s] Model [%s] not found", command, prodx[0]);
		    }
		} else { //for pparticle:
		    Int_t handle = GetMethodHandle(prodx[1]);
		    if (handle >= 0) {		
			key1 = GetKey(prodx[0], 1, -1);
			if (AddCommand(COMMAND_INTERNAL, key_a, key1, arg1, arg2, arg3, arg4))  {
			    found = kTRUE;
			    lst_command_int[command_pointer-1]  = handle;
			    flag_command_int[command_pointer-1] = 0;
			}
		    }
		}
	    }

	    //Make another try using the PFormula
	    if (!found) {
		found = EvalPFormula(command);
		if (found) return kTRUE;
	    }
	    
	    if (!found) {
		Error("AddCommand", "[%s] The method ->%s is unknown", command, prodx[1]);
		return kFALSE;
	    } else {
		//final thing is to create the result
		Double_t *delme =  new Double_t(0.);
		makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
		
		return kTRUE;
	    }
	} else {
	    //Warning("AddCommand","[%s] Unknown single command",command);
	    key_a=makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command);
	    prod[1] = command; //continue and try the best...
	    //return kFALSE; //not yet done
	}
    } else if (is_operator > 1) {
	//equal
	if (strstr(command, "==")) {
	    for (UInt_t i=0; i<(strlen(command)-1); i++) {
#if 0
		if (command[i]=='=' && command[i+1]=='=' ) {
		    Info("AddCommand","[%s] the '==' is not exact for doubles - lets use '~'",
			 command);
		    command[i] = '~';
		    command[i+1] = ' ';

		}
#endif
	    }
	} else {
	    Error("AddCommand", "[%s] Too many ='s", command);
	    return kFALSE;
	}
    } else {
	PUtils::Tokenize(command, "=", prod, &prod_s);

	key_a = GetKey(prod[0], 2, 2);
	copy_ctor = 1;
    }
    
    found = kFALSE;
    
    if (!is_operator) {
	//look for the arguments
	if (CheckAndSplit(prod[1], '~', &key1, &key2)) {
	    if (AddCommand(COMMAND_EQUAL, key_a, key1, key2))
		found = kTRUE;
	} else if (CheckAndSplit(prod[1], '+', &key1, &key2)) {
	    if (AddCommand('+', key_a, key1, key2))
		found = kTRUE;
	} else if (CheckAndSplit(prod[1], '-', &key1, &key2)) {
	    if (AddCommand('-', key_a, key1, key2))
		found = kTRUE;
	}
    }

    if (!found && copy_ctor) {
	//copy ctor
	key1 = GetKey(prod[1],1,-1);
	key2 = -1;
	if (is_operator == 1) {
	    if (AddCommand(COMMAND_IS, key_a, key1, key2))
		found = kTRUE;
	} else
	    if (AddCommand(COMMAND_EQUAL, key_a, key1, key2))
		found = kTRUE;
//	GetKey(prod[0],0,1);

	Bool_t makenew = kFALSE;

	//is lvalue one of the "switching" functions?
	for (UInt_t i=0; i<strlen(prod[0]); i++) 
	    if (*(prod[0]+i) == ')' || *(prod[0]+i) == '(') {
		makenew = kTRUE;
	    }
	if (makenew) AddCommand(prod[0]);
    } else if (!found) {
	//if nothing helps, try to use the PFormula....
	found = EvalPFormula(command);
	if (found) return kTRUE;
	//....or the wrapper to PUtils
	Int_t handle = GetMethodHandle(command, 1);
	if (handle >= 0) {		
	    //cout << "found PUtils with " << command3<< endl;
	    key1 = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command3);
	    Double_t *delme = new Double_t(0.);
	    makeDataBase()->SetParamDouble(key1, "batch_value", delme);
	    //key1 = GetKey(command3,1,1);
	    if (AddCommand(COMMAND_INTERNAL, key_a, key1, arg1, arg2, arg3, arg4))  {
		found = kTRUE;
		lst_command_int[command_pointer-1]  = handle;
		flag_command_int[command_pointer-1] = 1;
	    }
	    return kTRUE;
	}
    }

    if (found) {
	
	//Now look to the arguments, in order to know
	//if the result object is a PParticle or a double

	Int_t is_obj_result_type   = 0;
	Int_t is_value_result_type = 0;

	Int_t num_args = 1;
	if (key2 > -1) num_args = 2;

	for (int pos=0; pos<num_args; pos++) {
	    //loop over 2 args
	    Int_t key = key1;
	    if (pos == 1) key = key2;
	    
	    if (key > -1) {
		Double_t *delme;
		if (makeDataBase()->GetParamDouble(key, "batch_value", &delme))
		    is_value_result_type++;
		TObject *delme2;
		if (makeDataBase()->GetParamTObj(key, "batch_particle", &delme2))
		    is_obj_result_type++;
		if (!is_value_result_type && !is_obj_result_type) {
		    //have to check if objects are there and what kind they are
		    //4momentum?
		    
		    Int_t type = CheckObjectType(key);
		    if (type < 0) {
			Error("AddCommand", "[%s] Object not identified: %s", command,
			      makeDataBase()->GetName(key));
			return kFALSE;
		    }
		    
		    if (type == IS_OBJECT) is_obj_result_type++;
		    else is_value_result_type++;
		}
	    }
	}

	//cout << key_a << ":" << is_value_result_type << ":"<< is_obj_result_type << endl;

	if (is_obj_result_type) {
	    TObject *delme;
	    if (!makeDataBase()->GetParamTObj(key_a, "batch_particle", &delme)) {
		PParticle *newparticle = new PParticle(0,0,0,0);
		//pid existing? If yes, use it
		Int_t oldpid = makeStaticData()->GetParticleIDByKey(key_a);
		if (oldpid >= 0) newparticle->SetID(oldpid);
		TObject *delme = (TObject *) newparticle;
		makeDataBase()->SetParamTObj(key_a, "batch_particle", delme);
	    }
	}
	if (is_value_result_type) {
	    Double_t *delme;
	    if (!makeDataBase()->GetParamDouble(key_a, "batch_value", &delme)) {
		delme = new Double_t (0.);
		makeDataBase()->SetParamDouble (key_a, "batch_value", delme);
		//some warning for the access to PFormula
		if (is_operator && (strcmp(prod[0],"t") == 0 || strcmp(prod[0],"x") == 0 
				    || strcmp(prod[0],"y") == 0 || strcmp(prod[0],"z") == 0)) {
		    Warning("AddCommand",
			    "[%s] The '%s' is a reserved keyword in TFormula. You might run into problems.",
			    command, prod[0]);
		}
	    }
	}
	return kTRUE;
    } //end found operator

    Error("AddCommand", "[%s] Syntax error", command);

    return kFALSE;
}

void PBatch::AddSpacePlaceholder(char *command) {
    int numbra=0, numcbra=0;
    for (unsigned int i=0; i<strlen(command); i++) {
	if (command[i] == '[') numbra++;
	if (command[i] == '{') numcbra++;
	if (command[i] == ']') numbra--;
	if (command[i] == '}') numcbra--;
	if (numbra && (command[i] == ' '))  command[i] ='#';
	if (numcbra && (command[i] == '+')) command[i] ='\\';
	if (numcbra && (command[i] == '-')) command[i] ='`';
	if (numcbra && (command[i] == '/')) command[i] ='@';
    }
}

void PBatch::RemoveSpacePlaceholder(char *command) {
    int numbra=0, numcbra=0;
    for (unsigned int i=0; i<strlen(command); i++) {
	if (command[i] == '[') numbra++;
	if (command[i] == '{') numcbra++;
	if (command[i] == ']') numbra--;
	if (command[i] == '}') numcbra--;
	if (numbra && (command[i] == '#'))   command[i] =' ';
	if (numcbra && (command[i] == '\\')) command[i] ='+';
	if (numcbra && (command[i] == '`'))  command[i] ='-';
	if (numcbra && (command[i] == '@'))  command[i] ='/';
    }
}

Int_t PBatch::EvalPFormula(char *command) {
    //This helper function evaluates a (possible)
    //PFormula command, and put the arguments on the stack
    //if this did not worked out, the stack is resetted
    AddSpacePlaceholder(command);

    //cout << "PFormula called:" << command << endl;
    PFormula *tmp = new PFormula(command,command);    
    Int_t worked_out = 1, num_params = 0;

    const char *mod_command = command;
    Int_t lst_key_tmp[MAX_COMMAND_OPTIONS]; //can be overwritten by "daughter" commands
    Int_t lst_options_counter_tmp = 0;

    for (int i=0; i<MAX_COMMAND_OPTIONS; i++) lst_key_tmp[i] = -1;

    while (tmp->error_code && num_params<(MAX_COMMAND_OPTIONS-1)) {
	
	//try to get the ugly guy
	//copy the error string first
	char *internal_command = new char[strlen(tmp->error_string.Data())+1];
	strcpy(internal_command, tmp->error_string.Data());

	//cout << "Internal command error: "<< internal_command << endl;
	//BUGBUG
	//if the error string is *exactly* the input command, something is wrong....
	if (strcmp(internal_command,command) == 0) {
	    delete(tmp);
	    RemoveSpacePlaceholder(command);
	    return 0;
	}

	//	mod_command = ...aus PFormula, after replacing	
	TString *op = new TString(tmp->chaine);
	//cout << "Chaine from TFormula: "<< op->Data() << endl;
	char opt[5];
	sprintf(opt, "[%i]", num_params);
	//op->ReplaceAll(tmp->error_string,opt);
	ReplaceAll(op, tmp->error_string, opt);
	mod_command = op->Data();
	//cout << "Modified command: "<< mod_command  << endl;

	delete(tmp);
	RemoveSpacePlaceholder(internal_command);
	Int_t key = GetKey(internal_command, 1, -1);
	if (key >= 0) lst_key_tmp[num_params+1] = key;
	lst_options_counter_tmp++;
	tmp = new PFormula(mod_command, mod_command);   
	num_params++;
    }

    if (tmp->error_code) worked_out = 0;

    if (worked_out) {
	lst_form[command_pointer] = tmp;
	lst_options_counter[command_pointer] = lst_options_counter_tmp;
	for (int i=0; i<MAX_COMMAND_OPTIONS; i++)
	    lst_key[i][command_pointer] = lst_key_tmp[i];
	RemoveSpacePlaceholder(command);	
	Int_t key_a = makeStaticData()->
	    MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, command);
	Double_t *delme =  new Double_t(0.);
	makeDataBase()->SetParamDouble(key_a, "batch_value", delme);
	lst_key_a[command_pointer] = key_a;
	//	lst_key[1][command_pointer]=key1;
	lst_command[command_pointer] = COMMAND_PFORMULA;
	command_pointer++;
    } else RemoveSpacePlaceholder(command);
    
    return worked_out;
}

void PBatch::ReplaceAll(TString *op, const char *oldstring, const char *newstring) {
    //Similar to TString::ReplaceAll, but checks also the following
    //char

    bool doloop = kTRUE;

    while (doloop) {
	const char *mystring = op->Data();
	doloop = kFALSE;

	for (int i=0; i<((int)strlen(mystring)-(int)strlen(oldstring)+1); i++) {
	    //loop over possible matching positions

	    Bool_t varname = kFALSE;

	    if (strncmp(oldstring, mystring+i, strlen(oldstring))==0) {

		if (i) { //not the first position
		    if (PUtils::ValidVariableName(mystring+i-1, strlen(oldstring)+1))
			varname = kTRUE;
		}
		
		if (mystring[i+strlen(oldstring)]) { //not the last
		    if (PUtils::ValidVariableName(mystring+i, strlen(oldstring)+1))
			varname = kTRUE;
		}
		
		// 	    if (strncmp(oldstring,mystring+i,strlen(oldstring))==0) {
		// 		//match
		// 		//cout << "match: " << mystring << ":" << oldstring <<  ":" << newstring << endl;
		// 		char nextchar = mystring[i+strlen(oldstring)];
		
		// 		if ((nextchar == '\0') || 
		// 		    (!isalpha(nextchar) && 
		// 		     (nextchar != '.') &&
		// 		     (nextchar != '-') &&
		// 		     (nextchar != '_')) || 
		// 		    ((nextchar == '-') && (mystring[i+strlen(oldstring)+1] != '>'))
		// 		    ) {
		if (!varname) {
		    op->Replace(i, strlen(oldstring), newstring, strlen(newstring));
		    i = strlen(mystring)+1;
		    doloop = kTRUE;  //yet another try
		}
	    }
	}
    }
}

Int_t PBatch::GetMethodHandle(char *name, Int_t flag) {
    //seeks for the internal method
    //returns a -1 if not found or invalid

    //BUGBUG: have to include overloading

    const TList *list = NULL;
    if (flag == 0)
	list = gROOT->GetClass("PParticle")->GetListOfAllPublicMethods();
    else if (flag == 1)
	list = gROOT->GetClass("PUtilsREngine")->GetListOfAllPublicMethods();
    else if (flag == 2)
	list = gROOT->GetClass("PChannelModel")->GetListOfAllPublicMethods();
    else Fatal("GetMethodHandle", "Unsupported flag");

    TIterator *iter = list->MakeIterator();
    TMethod *meth = NULL;
    int error = 0;

    while (error < 2) {
	if (error == 1) error = 2; //already 2nd turn
	while ((meth = (TMethod *) iter->Next())) {

	    UInt_t realname_len = 0;
	    for (;realname_len<strlen(name); realname_len++) 
		if (name[realname_len] == '(') break;
	    
	    if ((strncmp(meth->GetName(),name,realname_len) == 0) && 
		(strlen(meth->GetName()) == realname_len)) {
		//found name
		//cout << meth->GetReturnTypeName() << endl;
		if ((strcmp(meth->GetReturnTypeName(), "Double_t")   ==0 ) || 
		    (strcmp(meth->GetReturnTypeName(), "Int_t")      ==0 ) || 
		    (strcmp(meth->GetReturnTypeName(), "PParticle*") ==0 ) || 
		    (strcmp(meth->GetReturnTypeName(), "void")       ==0 )) {
		    
		    //Set handle and prepare the pointers		 
		    //is MethodCall existing?

		    for (int i=0; i<method_pointer; i++) {
			if (strcmp(method_name[i],meth->GetName()) == 0) {
			    //found old entry			 
			    CrackMethodArgs(name+realname_len+1);
			    return i;
			}
		    }
		    
		    if (method_pointer == MAX_COMMAND_TMETHODS) {
			Error("GetMethodHandle", "MAX_COMMAND_TMETHODS reached");
			return -1;
		    }
		    
		    arg1 = arg2 = arg3 = arg4 = -1;
		    Int_t numargs = 0;
		    TString argstring("");

		    methods_arg_flags[0][method_pointer]
			= methods_arg_flags[1][method_pointer]
			= methods_arg_flags[2][method_pointer]
			= methods_arg_flags[3][method_pointer] = -1;

		    if (meth->GetNargs() != 0) {
			TList      *arg_list = meth->GetListOfMethodArgs();
			TIterator  *arg_iter = arg_list->MakeIterator();
			TMethodArg *arg_meth;
			
			int j = 0;
			while ((arg_meth=(TMethodArg *) arg_iter->Next())) {
			    //set up the argument list
			    //cout << arg_meth->GetTypeName() << endl;
			    if (j == 4 && error == 2) {
				Error("GetMethodHandle", "More then 4 arguments, not yet supported");
				return -1;
			    }
			    if (strlen(argstring.Data()) > 0) argstring += TString(",");
			    if (strcmp(arg_meth->GetTypeName(),"Double_t") ==0 ) {
				methods_arg_flags[j][method_pointer] = METHOD_RETURN_DOUBLE;
				argstring += TString("Double_t");
			    }
			    else if (strcmp(arg_meth->GetTypeName(),"Int_t") ==0 ) {
				methods_arg_flags[j][method_pointer] = METHOD_RETURN_INT;
				argstring += TString("Int_t");
			    }
			    else if (strcmp(arg_meth->GetTypeName(),"void") ==0 )
				methods_arg_flags[j][method_pointer] = METHOD_RETURN_VOID;
			    else {
				Error("GetMethodHandle", "Method %s has a %s as an argument, not yet supported",
				      name,arg_meth->GetTypeName());
				return -1;
			    }
			    j++;
			}
			
			numargs = CrackMethodArgs(name+realname_len+1);						
		    }

		    methods[method_pointer] = new TMethodCall();
		    
		    //cout << numargs << ":" << meth->GetNargs() << endl;

		    if (numargs == meth->GetNargs()) {
			
			if (strcmp(meth->GetReturnTypeName(), "Double_t") ==0 )
			    methods_flags[method_pointer] = METHOD_RETURN_DOUBLE;
			else if (strcmp(meth->GetReturnTypeName(), "Int_t") ==0 )
			    methods_flags[method_pointer] = METHOD_RETURN_INT;
			else if (strcmp(meth->GetReturnTypeName(), "PParticle*") ==0 )
			    methods_flags[method_pointer] = METHOD_RETURN_PPARTICLE;
			else 
			    methods_flags[method_pointer] = METHOD_RETURN_VOID;
			
			char *new_name = new char[strlen(argstring.Data())+1];
			strcpy(new_name, argstring.Data());
			
			if (flag == 0)
			    methods[method_pointer]->InitWithPrototype(gROOT->GetClass("PParticle"),
								       meth->GetName(), new_name);
			else if (flag == 1)
			    methods[method_pointer]->InitWithPrototype(gROOT->GetClass("PUtilsREngine"),
								       meth->GetName(), new_name);
			else if (flag == 2)
			    methods[method_pointer]->InitWithPrototype(gROOT->GetClass("PChannelModel"),
								       meth->GetName(), new_name);
			new_name = new char[strlen(meth->GetName())+1];
			strcpy(new_name, meth->GetName());
			method_name[method_pointer] = new_name;
			method_pointer++;
			return method_pointer-1;
		    } else { 
			error = 3;
		    }
		  
		}  else {
		    if (!error) {
			error = 1; 
		    } else if (error == 2){
			//unsupported return type
			Error("GetMethodHandle", "Method %s has a return type %s, not yet supported",
			      name, meth->GetReturnTypeName());
		    }
		} // return is matching
	    } // found meth name

	} //meth iter
	if (!error) error = 2;
    } //while error

    if (error == 3) Error("GetMethodHandle", "Method %s does not match any prototype",
			  name);
    return -1;
}

Int_t PBatch::CrackMethodArgs(char *name) {
    //Input: the arg string including the trailing ")"
    arg1 = arg2 = arg3 = arg4 = -1;

    //cout << "CrackMethodArgs:" << name << endl;
    //nested objects should not be cracked!
    //workaround: replace , by "

    Int_t numbrack = 0;
    for (UInt_t i=0; i<strlen(name); i++) {
	if ((name[i]=='(') || (name[i]=='[')) numbrack++;
	if ((name[i]==')') || (name[i]==']')) numbrack--;
	if ((name[i]==',') && (numbrack>0)) name[i]='"';
    }

    char *prodx[4];
    Int_t prodx_s = 4; //max 4 args
    
    PUtils::Tokenize(name, ",", prodx, &prodx_s);
    *(prodx[prodx_s-1]+strlen(prodx[prodx_s-1])-1) = '\0'; //remove trailing )

    if (!prodx_s) return 0;    

    numbrack = 0;
    for (UInt_t i=0; i<strlen(name); i++) {
	if ((name[i]=='(') || (name[i]=='[')) numbrack++;
	if ((name[i]==')') || (name[i]==']')) numbrack--;
	if ((name[i]=='"') && (numbrack>0)) name[i]=',';
    }
    
    arg1 = GetKey(prodx[0], 1, -1);
    if (arg1 < 0) return 0;
    if (prodx_s > 1) arg2 = GetKey(prodx[1], 1, -1);
    if (prodx_s > 2) arg3 = GetKey(prodx[2], 1, -1);
    if (prodx_s > 3) arg4 = GetKey(prodx[3], 1, -1);

    return prodx_s;
}

Int_t PBatch::CheckObjectType(Int_t key) {

    Int_t *ii;
    if (makeDataBase()->GetParamInt (key, pid_param, &ii))
	return IS_OBJECT;

    return -1;
}

Bool_t PBatch::AddCommand(char command, int key_a, int key1, int key2, int key3, int key4, int key5) {

    if (command == '=' && key_a==key1)
	return kTRUE; //filter nonsense

    if (command_pointer == MAX_COMMAND_POINTER) {
	Error ("AddCommand","MAX_COMMAND_POINTER reached");
	return kFALSE;
    }

    lst_command[command_pointer] = command;
    lst_key_a[command_pointer]   = key_a;
    lst_key[1][command_pointer]  = key1;
    lst_key[2][command_pointer]  = key2;
    lst_key[3][command_pointer]  = key3;
    lst_key[4][command_pointer]  = key4;
    lst_key[5][command_pointer]  = key5;

    command_pointer++;

    return kTRUE;
}

Bool_t PBatch::GetArguments(const char *a, char *name, char **function, char **arg1, char **arg2) {
    //looks for syntax like f(a,b);
    
    char *prod[2];
    Int_t prod_s = 2; //max 2 products

    PUtils::Tokenize(name, a, prod, &prod_s);

    if (prod[1] == NULL) {
	*function = NULL;
	prod[1] = prod[0];	
    }

    *(prod[1]+strlen(prod[1])-1) = '\0';
    //kill ")"

    int komma = 0;

    for (UInt_t i=0; i<strlen(prod[1]); i++) {
	if (*(prod[1]+i) == ',') komma++;
    }
    
    if (komma > 1) {
	Error("AddCommand", "[%s] Too many kommas", prod[1]);
	return kFALSE;
    }

    if (komma) {
	char *prodx[2];
	Int_t prodx_s = 2; //max 2 products
	
	PUtils::Tokenize(prod[1], ",", prodx, &prodx_s);
	*arg1 = prodx[0];
	*arg2 = prodx[1];
    } else {
	*arg1 = prod[1];
	*arg2 = NULL;
    }

    return kTRUE;
}

Int_t PBatch::GetDelimPosition(char *arg, char delim, Int_t *yes) {
    
    Int_t brackets = 0, 
	split_pos = -1,
	found_something = 0;
    for (UInt_t i=0; i<strlen(arg); i++) {
	if ((arg[i]=='(') || (arg[i]=='[')) brackets++;
	if ((arg[i]==')') || (arg[i]==']')) brackets--;
	if ((arg[i]==delim) && (brackets==0)) {
	    split_pos = i;
	    if (arg[i+1]=='>' || !found_something) split_pos=-1;
	    //...skip -> and - in the beginning
	    else if (yes) (*yes)++;
	}
	if (arg[i]!=' ') found_something++;
    }

    return split_pos;
}

Bool_t PBatch::CheckAndSplit(char *arg, char delim, int *key1, int *key2) {

    int yes = 0, 
	split_pos = -1;
    
    char *arg2 = new char[strlen(arg)+1];
    strcpy(arg2, arg);
    arg = arg2;

    if (strstr(arg, "->")) return kFALSE;
    if (strstr(arg, "{"))  return kFALSE;

    split_pos = GetDelimPosition(arg, delim, &yes);
    if (!split_pos) return kFALSE;  //nothing to split

    //look if the "-" is a fake minus, e.g. "a < -0.5"
    for (int i=split_pos-1; i>=0; i--) {
	if (arg[i] != ' ') {
	    if (arg[i] == '<' || arg[i] == '>' || arg[i] == '=' 
		|| arg[i] == ':'|| arg[i] == ',') {
		return kFALSE; 
	    }
	    i = -1;	    
	}
    }

    if (!yes) {
	return kFALSE; 
    }

    char *prod[2];
    
    prod[0] = new char[split_pos+1];
    prod[1] = new char[strlen(arg)-split_pos+1];
    strncpy(prod[0], arg, split_pos);
    *(prod[0]+split_pos) = '\0';
    strncpy(prod[1], arg+split_pos+1, strlen(arg)-split_pos-1);
    *(prod[1]+strlen(arg)-split_pos-1) = '\0';

    *key1 = GetKey(prod[0], 1, -1);
    *key2 = GetKey(prod[1], 1, -1);

    return kTRUE;
}

Int_t PBatch::GetKey(char *name, int fl, int makeflag) {

    //makeflag=0  : Pure GetKey, just take the key from DB
    //              no AddCommand()
    //makeflag=1  : put our object into the DB in any case
    //              this is the case for lvalues
    //              e.g. obj = ....
    //              (but see also SetVarList)
    //makeflag=2  : Same as above, but it must be a valid variable name
    //makeflag=-1 : put into DB only if we have a composite object
    //              This is the case if we find () or + or - 
    //              e.g. obj->...
    //              make AddCommand
    //fl=2        : Take SetVarList into account

    //cout << "getKey" << name << ":" << makeflag <<  endl;
  
    //do we contain ()?
    //int br = 0;
    //int found_br = 0;

    //remove brackets if fl==1
    PUtils::remove_spaces(&name);
    if (fl) {
     	//found_br =  PUtils::remove_brackets(&name,'(',')' );
	PUtils::remove_brackets(&name, '(', ')');
    }
    PUtils::remove_spaces(&name);

    Int_t key = -1; 

    if (strlen(name) == 0) return -1;

    //first check for file-input
    if (name[0] == '[' && name[strlen(name)-1]==']' && strlen(name)>2) {
	//found [] at the very begin & end
		
	//further look to the content
	//if we find any additional brackets, it is a composite object!
	Int_t num_brackets = 0;
	for (UInt_t i=1; i<(strlen(name)-1); i++) {
	    if (name[i] == '[' || name[i]==']') num_brackets++;
	}

	if (!num_brackets) {

	    key = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, name);
	    
	    char *function, *arg1, *arg2;
	    GetArguments("[", /* "]", */ name, &function, &arg1, &arg2);
	   	    
	    if (!arg1) { 
		Error("AddCommand", "[%s] Argument not found", name);
		return kFALSE;
	    }
	   
	    Int_t pid = 0; //DUMMY

	    if (strcmp(arg1, "*")) {
		//no Joker
		pid = makeStaticData()->GetParticleID(arg1);
		if (pid==0)  { 
		    Error("AddCommand", "[%s] Unknown particle %s", name, arg1);
		    return kFALSE;
		}
	    }
	    
	    Int_t *ii = new int(pid);  //never destructed, but called only once!
	   	
	    if (!makeDataBase()->SetParamInt (key, "batch_pid", ii)) {
		delete ii;
		return kFALSE;
	    }
	    
	    Int_t number = -999;
	    
	    if (arg2) { //try to get number
		Int_t *delme = new Int_t(number);
		Int_t numkey = -999;
		//First, let's check for a variable
		if (arg2[0] == '$') {
		    if (PUtils::ValidVariableName(arg2+1)) {
			numkey = makeDataBase()->GetEntry(arg2+1);
			if (numkey > 0) {
			    Double_t *val = NULL;
			    if (!makeDataBase()->GetParamDouble(numkey, batch_value_param, &val)) {
				numkey = -1;
			    } else {
				Int_t *result = new Int_t(1);
				makeDataBase()->SetParamInt(numkey, "batch_update", result);   
			    }
			}
		    } 
		    if (numkey < 0) {
			Error("AddCommand", "[%s] Unknown variable '%s'", name, arg2);
		    } 
		    *delme = -1000 - numkey; //-999 means not found!
		} else if (PUtils::IsInt(arg2)) {
		    sscanf(arg2, "%i", &number);
		    *delme = number;
		    if (strcmp(arg2,"+") == 0) { 
			*delme = -111;
			status = 1;
		    }
		} else 
		    Error("AddCommand", "[%s] Unknown value '%s'", name, arg2);
		makeDataBase()->SetParamInt(key, "batch_position", delme);
	    }
	    makeflag = 0;
	}
    } //END file input

    if (makeflag==-1 && PUtils::ValidVariableName(name)) {
	//cancel AddCommand
	makeflag = 0;
    }

    if (makeflag==2 && !PUtils::ValidVariableName(name)) {
	//cancel AddCommand
	makeflag = 0;
	Error("GetKey", "[%s] is not a valid variable name", name);
    }
    if (makeflag == 2) 
	makeflag = 1;
    if (makeflag == -1) 
	AddCommand(name);

    if (makeflag) {
	if ((varlist == NULL) || (fl!=2)) {
	    key = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, name);
	} else {
	    for (unsigned int j=0; j<((strlen(varlist) - strlen(name))); j++) {
		if ((strncmp(varlist+j,name,strlen(name)) == 0) && 
		    varlist[j+strlen(name)] == ';') {
		    key = makeStaticData()->
			MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, name);
		}
	    }
	    if (key < 0) 
		Error("AddCommand",
		      "[%s] is not an allowed object name, the list in this context is [%s] ",
		      name, varlist);
	}
    } else if (key < 0) {
	key = makeDataBase()->GetEntry(name);
    }
    if (key < 0) {
	Error("AddCommand", "[%s] Unknown object", name);
    }
    
    return key;
}

void PBatch::Print(const Option_t *) const {

    cout << "Command list: res <=== arg{key,name}" << endl;

    for (int i=0;i<command_pointer;i++) {

// 	if (lst_command[i] == COMMAND_PFORMULA) {

// 	    cout << makeDataBase()->GetName(lst_key_a[i]) 
// 		 << " <" << lst_command[i] << "> ";
// 	    for (int j=0;j<lst_options_counter[i];j++) {
// 		cout << makeDataBase()->GetName(lst_key[j+1][i])  << " ";
// 	    }
// 	    cout << endl;
// 	} else 
	{
	    
	    if (lst_command[i] == COMMAND_PFORMULA) 
		cout << "[PFormula]  :";
	    if (lst_command[i] == COMMAND_MASS) 
		cout << "->M()       :";
	    if (lst_command[i] == COMMAND_MASS2) 
		cout << "->M2()      :";
	    if (lst_command[i] == COMMAND_IF) 
		cout << "if          :";
	    if (lst_command[i] == COMMAND_PLUS) 
		cout << "+           :";
	    if (lst_command[i] == COMMAND_MINUS) 
		cout << "-           :";
	    if (lst_command[i] == COMMAND_EQUAL) 
		cout << "~           :";
	    if (lst_command[i] == COMMAND_IS) 
		cout << "=           :";
	    if (lst_command[i] == COMMAND_COS) 
		cout << "cos()       :";
	    if (lst_command[i] == COMMAND_BOOST) 
		cout << "->Boost()   :";
	    if (lst_command[i] == COMMAND_GETBEAM) 
		cout << "->GetBeam() :";
	    if (lst_command[i] == COMMAND_GETTARGET) 
		cout << "->GetTarget():";
	    if (lst_command[i] == COMMAND_FABS) 
		cout << "fabs()      :";
	    if (lst_command[i] == COMMAND_PRINT) 
		cout << "->Print()   :";
	    if (lst_command[i] == COMMAND_ROT) 
		cout << "->Rot()     :";
	    if (lst_command[i] == COMMAND_THETA) 
		cout << "->Theta()   :";
	    if (lst_command[i] == COMMAND_INTERNAL) 
		cout << "[internal]  :";
	    if (lst_command[i] == COMMAND_EVAL) 
		cout << "Eval()      :";
	    if (lst_command[i] == COMMAND_GETRANDOM) 
		cout << "GetRandom() :";
	    if (lst_command[i] == COMMAND_ANGLE) 
		cout << "->Angle()   :";
	    if (lst_command[i] == COMMAND_P3M) 
		cout << "P3M()       :";
	    if (lst_command[i] == COMMAND_P3E) 
		cout << "P3E()       :";
	    if (lst_command[i] == COMMAND_PVALUE) 
		cout << ".val        :";
	    if (lst_command[i] == COMMAND_LABEL) 
		cout << "label:      :";
	    if (lst_command[i] == COMMAND_GOTO) 
		cout << "goto        :";
	    if (lst_command[i] == COMMAND_GOSUB) 
		cout << "gosub       :";
	    if (lst_command[i] == COMMAND_RETURN) 
		cout << "return      :";
	    if (lst_command[i] == COMMAND_EXIT) 
		cout << "exit        :";
	    if (lst_command[i] == COMMAND_FORMORE) 
		cout << "formore     :";
	    if (lst_command[i] == COMMAND_FOREACH) 
		cout << "foreach     :";
	    if (lst_command[i] == COMMAND_ECHO) 
		cout << "echo        :";
	    if (lst_command[i] == COMMAND_READLINE) 
		cout << "readline    :";
	    if (lst_command[i] == COMMAND_PUSH) 
		cout << "push        :";
	    if (lst_command[i] == COMMAND_BRANCH) 
		cout << "branch      :";
	    if (lst_command[i] == COMMAND_ELSE) 
		cout << "else        :";

	    if (lst_key_a[i] > -1) {
		cout << makeDataBase()->GetName(lst_key_a[i]);
	    }

	    if (lst_key[1][i] > -1) {
 		cout << " <=== {" << lst_key[1][i] 
		     << "=" << makeDataBase()->GetName(lst_key[1][i]) 
		     << "}";
	    }

	    if (lst_key[2][i] > -1) {
 		cout << ",{" << lst_key[2][i] 
		     << "=" << makeDataBase()->GetName(lst_key[2][i]) 
		     << "}";
	    }

	    if (lst_key[3][i] > -1) {
 		cout << ",{" << lst_key[3][i] 
		     << "=" << makeDataBase()->GetName(lst_key[3][i]) 
		     << "}";
	    }

	    if (lst_key[4][i] > -1) {
 		cout << ",{" << lst_key[4][i] 
		     << "=" << makeDataBase()->GetName(lst_key[4][i]) 
		     << "}";
	    }

	    cout << endl;
	}
    }
}


ClassImp(PBatch)

