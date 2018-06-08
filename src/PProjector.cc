////////////////////////////////////////////////////////
//  Pluto projector interface
//
//  The main idea of this class is to provide a very
//  simple, fast and easy-to-use analysis tool, which
//  make the writing of an analysis macro unnecessary.
//  It uses the Bulk interface and can be plugged into
//  the reaction before or after the decay
//
//  The syntax of the commands are based on the PBatch syntax,
//  for details look into the documentation of this class.
//
//  Commands for the batch can be added in 2 ways: Either
//  with AddCommand(char * command) or
//  in one step together with a histogram which will 
//  be filled after the command has been executed.
//
//  Input: The input for the commands must be particles
//  from the data stream in the fillowing format:
//  [pid,num] where pid is the pid-string (e.g. "pi0")
//  and num is the number of the particle (can be omitted if only
//  one particle of this type present)
//
//  Output: _x,_y as doubles. The latter one only for 2dim histograms.
//
//                    Author:  Ingo Froehlich
//                    Written: 14/02/2008
//                    Revised: 
//
////////////////////////////////////////////////////////



#include "PProjector.h"
#include "PChannel.h"


PProjector::PProjector() {
    batch_pos = 0;
    pid_param = makeDataBase()->GetParamInt("batch_pid");
    if (pid_param < 0) 
	pid_param = makeDataBase()->MakeParamInt("batch_pid", "PID for batch");
    
    link_param = makeDataBase()->GetParamInt("batch_position");
    if (link_param < 0) 
	link_param = makeDataBase()->MakeParamInt("batch_position", "PID position for batch");
    
    batch_particle_param = makeDataBase()->GetParamTObj("batch_particle");
    w = makeStaticData()->GetBatchValue("_w");
    
    if (batch_particle_param < 0) 
	batch_particle_param = makeDataBase()->MakeParamTObj("batch_particle", "PParticle storage for batch");

    batch_value_param = makeDataBase()->GetParamDouble("batch_value");

    stream_default_pos_param = makeDataBase()->GetParamInt(STREAM_DEFAULT_POS);
    if (stream_default_pos_param < 0)  
	stream_default_pos_param = makeDataBase()->MakeParamInt(STREAM_DEFAULT_POS, "Default position");
    stream_max_pos_param = makeDataBase()->GetParamInt(STREAM_MAX_POS);
    if (stream_max_pos_param < 0)  
	stream_max_pos_param = makeDataBase()->MakeParamInt(STREAM_MAX_POS, "Max position in stream");


    for (int i=0; i<PROJECTOR_MAX_BATCH; i++) { 
	hist3[i]       = NULL;
	hist2[i]       = NULL;
	hist1[i]       = NULL;
	fp_out[i]      = NULL;
	fp_in[i]       = NULL;
	key_pos_out[i] = 0;
	key_pos_in[i]  = 0;
    }

    batch_pos    = 0;
    force_weight = 1;
    current_ascii_file = NULL;

    fPriority = PPROJECTOR_PRIORITY;
    proj_nr   = makeStaticData()->GetBatchValue("_system_embedded_particle_projector");

    for (int i=0; i<MAX_NUM_BRANCHES; i++) {
	current_size_branches[i]   = NULL;
	particle_array_branches[i] = NULL;
    }
}

PProjector::~PProjector() {
    for (int i=0; i<batch_pos; i++) 
	delete batch[i];
}

Bool_t PProjector::AddCommand(const char *command) {

    //adds a command line to batch
    if (batch_pos == PROJECTOR_MAX_BATCH) {
	Error("AddCommand", "PROJECTOR_MAX_BATCH reached");
	return kFALSE;
    }
    
    batch[batch_pos] = new PBatch();
    batch[batch_pos]->SetSizeBranches(size_branches);
    batch[batch_pos]->SetKeysBranches(key_branches);

    if (current_ascii_file) {
	batch[batch_pos]->SetToolObjectTmp(current_ascii_file);
    }

    batch[batch_pos]->SetPosition(batch_pos, bulk_id);  //Set absolute adress
    if (!batch[batch_pos]->AddCommand(PUtils::NewString(command)))  {
	delete batch[batch_pos];
	return kFALSE;
    }

    key = makeDataBase()->GetEntry("batch_objects");
    if (batch[batch_pos]->Status()) {
	//[+]
	if (*(proj_nr)) {
	    Error("AddCommand", "Embedded particle ([+]) used 2 times");
	}
	*(proj_nr) = bulk_id;
    }

    batch_pos++;
    return kTRUE;
}

Bool_t PProjector::AddHistogram(TH3 *histo, const char *command, Int_t fillflag) {
    
    if (!AddCommand(command)) return kFALSE;
    hist3[batch_pos-1] = histo;

    //get the result
    key   = makeDataBase()->GetEntry("batch_objects");
    key_x = makeDataBase()->GetEntry("_x");
    key_y = makeDataBase()->GetEntry("_y");
    key_z = makeDataBase()->GetEntry("_z");

    if (key_x < 0) {
	Error ("AddHistogram", "result _x not found");
	return kFALSE;
    }
    if (key_y < 0) {
	Error ("AddHistogram", "result _y not found");
	return kFALSE;
    }
    if (key_z < 0) {
	Error ("AddHistogram", "result _z not found");
	return kFALSE;
    }
    if (!makeDataBase()->GetParamDouble(key_x, batch_value_param, &x)) {
	Error ("AddHistogram", "Double _x not found");
	return kFALSE;
    }
    if (!makeDataBase()->GetParamDouble(key_y, batch_value_param, &y)) {
	Error ("AddHistogram", "Double _y not found");
	return kFALSE;
    }
    if (!makeDataBase()->GetParamDouble(key_z, batch_value_param, &z)) {
	Error ("AddHistogram", "Double _z not found");
	return kFALSE;
    }

    if (batch[batch_pos-1]->IsReadonly()) 
	fillflag = 0;
    batch[batch_pos-1]->SetToolObject(histo);
    fFillFlag[batch_pos-1] = fillflag;

    return kTRUE;
}

Bool_t PProjector::AddHistogram(TH2 *histo, const char *command, Int_t fillflag) {
    
    if (!AddCommand(command)) return kFALSE;
    hist2[batch_pos-1] = histo;

    //get the result
    key   = makeDataBase()->GetEntry("batch_objects");
    key_x = makeDataBase()->GetEntry("_x");
    key_y = makeDataBase()->GetEntry("_y");

    if (key_x < 0) {
	Error ("AddHistogram", "result _x not found");
	return kFALSE;
    }
    if (key_y < 0) {
	Error ("AddHistogram", "result _y not found");
	return kFALSE;
    }
    if (!makeDataBase()->GetParamDouble(key_x, batch_value_param, &x)) {
	Error ("AddHistogram", "Double _x not found");
	return kFALSE;
    }
    if (!makeDataBase()->GetParamDouble(key_y, batch_value_param, &y)) {
	Error ("AddHistogram", "Double _y not found");
	return kFALSE;
    }

    if (batch[batch_pos-1]->IsReadonly()) 
	fillflag = 0;
    batch[batch_pos-1]->SetToolObject(histo);
    fFillFlag[batch_pos-1] = fillflag;

    return kTRUE;
}

Bool_t PProjector::AddHistogram(TH1 *histo, const char *command, Int_t fillflag) {
    
    if (!AddCommand(command)) return kFALSE;
    hist1[batch_pos-1] = histo;

    //get the result
    key   = makeDataBase()->GetEntry("batch_objects");
    key_x = makeDataBase()->GetEntry("_x");

    if (key_x < 0) {
	Error ("AddHistogram", "result _x not found");
	return kFALSE;
    }

    if (!makeDataBase()->GetParamDouble(key_x, batch_value_param, &x)) {
	Error ("AddHistogram", "Double _x not found");
	return kFALSE;
    }
    
    if (batch[batch_pos-1]->IsReadonly()) 
	fillflag = 0;
    batch[batch_pos-1]->SetToolObject(histo);
    fFillFlag[batch_pos-1] = fillflag;

    return kTRUE;
}

Bool_t PProjector::AddTGraph(TGraph *graph, const char *command) {
    
    if (!AddCommand(command)) return kFALSE;

    //get the result
    key   = makeDataBase()->GetEntry("batch_objects");
    key_x = makeDataBase()->GetEntry("_x");

    if (key_x < 0) {
	Error ("AddTGraph", "result _x not found");
	return kFALSE;
    }

    if (!makeDataBase()->GetParamDouble(key_x, batch_value_param, &x)) {
	Error ("AddTGraph", "Double _x not found");
	return kFALSE;
    }

    batch[batch_pos-1]->SetToolObject(graph);
    return kTRUE;
}

Bool_t PProjector::AddTGraph2D(TGraph2D *graph, const char *command) {
    
    if (!AddCommand(command)) return kFALSE;

    //get the result
    key   = makeDataBase()->GetEntry("batch_objects");
    key_x = makeDataBase()->GetEntry("_x");
    key_y = makeDataBase()->GetEntry("_y");

    if (key_x < 0) {
	Error ("AddTGraph2D", "result _x not found");
	return kFALSE;
    }
    if (key_y < 0) {
	Error ("AddTGraph2D", "result _y not found");
	return kFALSE;
    }
    if (!makeDataBase()->GetParamDouble(key_x, batch_value_param, &x)) {
	Error ("AddTGraph2D", "Double _x not found");
	return kFALSE;
    }
    if (!makeDataBase()->GetParamDouble(key_y, batch_value_param, &y)) {
	Error ("AddTGraph2D", "Double _y not found");
	return kFALSE;
    }

    batch[batch_pos-1]->SetToolObject(graph);
    return kTRUE;
}

Bool_t PProjector::AddOutputTNtuple(TNtuple *n, const char *command) {

    if (!AddCommand(command)) {
	return kFALSE;
    }

    fp_out[batch_pos-1] = n;

    //get the result
    key = makeDataBase()->GetEntry("batch_objects");
    //Examine NTuple file and get the branches

    TIter iter(fp_out[batch_pos-1]->GetListOfBranches());
    while(TBranch *br = (TBranch *)iter.Next()) {
	const char *name = br->GetName();

	//Each branch should be correlated to the batch key
	
	if (key_pos_out[batch_pos-1] == PROJECTOR_MAX_BRANCHES) {
	    Error("AddNTuple", "Too many branches in NTuple");
	    return kFALSE;
	}

	key_list_out[batch_pos-1][key_pos_out[batch_pos-1]] = makeDataBase()->GetEntry((char *) name);

	if (key_list_out[batch_pos-1][key_pos_out[batch_pos-1]] < 0) {
	    Warning("AddOutputNTuple", "Branch %s found in NTuple but not defined as a batch value", name);
	}

	key_pos_out[batch_pos-1]++;
    }
    return kTRUE;    
}

Bool_t PProjector::AddInputASCII(const char *filename, const char *command) {
    FILE * file = fopen(filename, "r");
    if (!file) {
	Error("AddInputASCII", "Could not open %s",filename);
	return kFALSE;
    }
    current_ascii_file = file;
    if (!AddCommand(command)) {
	return kFALSE;
    }
    batch[batch_pos-1]->SetToolObject(file);
    fPriority = FILEINPUT_PRIORITY;
    return kTRUE;   
}

Bool_t PProjector::AddOutputASCII(const char *filename, const char *command) {
    FILE * file = fopen(filename,"w");
    if (!file) {
	Error("AddOutputASCII", "Could not open %s", filename);
	return kFALSE;
    }
    current_ascii_file = file;
    if (!AddCommand(command)) {
	return kFALSE;
    }
    batch[batch_pos-1]->SetToolObject(file);
    return kTRUE;   
}

Bool_t PProjector::AddInputTNtuple(TNtuple *n,const  char *command) {

    if (!AddCommand(command)) {
	return kFALSE;
    }

    fp_in[batch_pos-1] = n;
    num_events_in[batch_pos-1] = n->GetEntries();
    num_events_in_c[batch_pos-1] = 0; //counted events

    //get the result
    key = makeDataBase()->GetEntry("batch_objects");
    //Examine NTuple file and get the branches

    TIter iter(fp_in[batch_pos-1]->GetListOfBranches());
    while(TBranch *br = (TBranch *)iter.Next()) {
	const char *name = br->GetName();

	//Each branch should be correlated to the batch key
	
	if (key_pos_in[batch_pos-1] == PROJECTOR_MAX_BRANCHES) {
	    Error("AddNTuple", "Too many branches in NTuple");
	    return kFALSE;
	}

	key_list_in[batch_pos-1][key_pos_in[batch_pos-1]] = makeDataBase()->GetEntry((char *) name);

	if (key_list_in[batch_pos-1][key_pos_in[batch_pos-1]] < 0) {
	    key_list_in[batch_pos-1][key_pos_in[batch_pos-1]] = makeStaticData()->
		MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, name);
        if (pluto_global::verbosity >= 3) {
            Info("AddInputTNtuple", "Created variable %s for the TNtuple branch", name);
        }
	}

	//Check if double is existing
	Double_t *val;
	if (!makeDataBase()->GetParamDouble(key_list_in[batch_pos-1][key_pos_in[batch_pos-1]], batch_value_param, &val)) {
	    Double_t *delme = new Double_t(0.);
	    makeDataBase()->SetParamDouble(key_list_in[batch_pos-1][key_pos_in[batch_pos-1]], "batch_value", delme);
	}

	//Set the branch adress to the floats
	n->SetBranchAddress(name, &(values_in[batch_pos-1][key_pos_in[batch_pos-1]]));

	key_pos_in[batch_pos-1]++;
    }

    fPriority = FILEINPUT_PRIORITY;

    return kTRUE;    
}

Int_t  PProjector::SetParticles(PParticle **mstack, int *, int *num, int, Int_t first_time) {
    //loop over batch object and see what I can do

    Int_t listkey=-1, *i_result, new_particles=0, counter_all=0, particle_key;

    if (first_time) {
 	for (int i=0; i<*num; i++) {
	    //cout << "FT:" << i << ":" << mstack[i]->ID() << endl;
 	    //Delete old default pos list
	    // Reset value on 1st call -> BUGBUG: can cause trouble when jumping over 2 projectors
 	    particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
 	    if (makeDataBase()->GetParamInt(particle_key, stream_default_pos_param, &i_result)) {
 		(*i_result) = 0;
 	    } 
 	}
	particle_key = makeStaticData()->GetParticleKey(0); //DUMMY=all particles
	if (makeDataBase()->GetParamInt(particle_key, stream_default_pos_param, &i_result)) {
	    (*i_result) = 0;
	} 
	for (int i=0; i<*num; i++) {
 	    //Delete old max pos list
 	    particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
 	    if (makeDataBase()->GetParamInt(particle_key, stream_max_pos_param, &i_result)) {
 		(*i_result) = 0;
 	    } 
 	}
	particle_key = makeStaticData()->GetParticleKey(0); //DUMMY=all particles
	if (makeDataBase()->GetParamInt(particle_key, stream_max_pos_param, &i_result)) {
	    (*i_result) = 0;
	} 
	for (int i=0; i<*num; i++) {
	    //Create max list
	    if (mstack[i]->IsActive()) { //count only active particles
		particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
		if (makeDataBase()->GetParamInt(particle_key, stream_max_pos_param, &i_result)) {
		    (*i_result)++;
		} else {
		    Int_t *dummy = new Int_t(1);
		    makeDataBase()->SetParamInt(particle_key, STREAM_MAX_POS, dummy);
		}
		particle_key = makeStaticData()->GetParticleKey(0); //DUMMY=all particles
		if (makeDataBase()->GetParamInt (particle_key, stream_max_pos_param, &i_result)) {
		    (*i_result)++;
		} else {
		    Int_t *dummy = new Int_t(1);
		    makeDataBase()->SetParamInt(particle_key, STREAM_MAX_POS, dummy);
		}
		counter_all++;
	    }
	    mstack[i]->SetScatterClone(kFALSE); //disable making clones when copy constructor is called -> memory leak
	} 
	//cout << "counted max " << counter_all << endl;
    }//END first_time

    while (makeDataBase()->MakeListIterator(key, NBATCH_NAME, LBATCH_NAME, &listkey)) {
        //loop over all particles
	//cout << "key=" << listkey << endl;
	if (makeDataBase()->GetParamInt (listkey, pid_param, &i_result)) {
	    Int_t pid = *i_result;

	    //fill object
	    Int_t pos = -1;
	    if (makeDataBase()->GetParamInt(listkey, link_param, &i_result)) {
		pos = *i_result-1;
	    }
	    //cout << "key=" << listkey << " pos:" << pos << endl;
	    //First clear the entry to avoid the use of old objects
	    if (pos >= 0) {
		makeDataBase()->SetParamTObj(listkey, batch_particle_param, NULL);
	    } else if ((pos == -112) && first_time && (*(proj_nr) == bulk_id)) { 
		//stumbled over "+"
		//cout << "CALLED + in "<< bulk_id << endl;
		new_particles++;
		if (new_particles == PROJECTOR_MAX_STACK) {
		    Warning("Modify", "PROJECTOR_MAX_STACK reached");
		    return kFALSE;
		}
		if (new_particles > stackpointer) {
		    //create new pparticle
		    stack[stackpointer] = new PParticle(0,0,0,0);
		    stackpointer++;
		    Info("SetParticles", "New particle created");
		}
		mstack[*num] = &(stack[new_particles-1]);
		stack[new_particles-1].SetID(pid);
		stack[new_particles-1].SetW(1.0);
		makeDataBase()->SetParamTObj(listkey, batch_particle_param, &(stack[new_particles-1]));
		(*num)++;
	    } else if (pos == -1) {
		//No pos, default 
		makeDataBase()->SetParamTObj(listkey, batch_particle_param, NULL);
		Int_t particle_key = makeStaticData()->GetParticleKey(pid);
		if (makeDataBase()->GetParamInt(particle_key, stream_default_pos_param, &i_result)) {
		    pos = (*i_result);
		} else
		    pos = 0;
	    } else if (pos < -1000) { //found link to variable
		//cout << "pos  is now:" << pos << edl;
		Double_t *res;
		if (makeDataBase()->GetParamDouble((-(pos+1))-1000, batch_value_param, &res)) {
		    //cout << "key: " << ((-(pos+1))-1000) << " res: " << *res << endl;
		    //makeDataBase()->ListEntries(-1,1,"*name,batch_value,*num_batch,*pid,*link");
		    pos = ((Int_t) *res)-1;
		} else {
		    pos = -1000;
		}
	    }

	    if (pos == -1000) {
		Error("SetParticles", "Unkown particle position for %s", makeDataBase()->GetName(listkey));
	    }

	    //cout << "pos  is now:" << pos << " for " << makeDataBase()->GetName(listkey) <<  endl;

	    for (int i=0; i<*num; i++) {
		//cout << mstack[i]->IsActive() << endl;
		if (mstack[i]->IsActive()) {
		    
		    //cout << "stack_pid:"<< mstack[i]->ID() << endl;
		    //mstack[i]->Print();
		    if ((mstack[i]->ID() == pid) || (pid==0)) { //0=DUMMY
			//cout << "match " << pid << " at " << pos << endl;
			if (pos == 0) {
			    TObject *delme =  (TObject *) mstack[i];
			    makeDataBase()->SetParamTObj(listkey, batch_particle_param, delme);
			    //makeDataBase()->ListEntries(listkey,1,"name,*pid,*batch_particle");
			} 
			pos--;
		    }
		}
	    }
	}	
    }
    return 0;
}

Bool_t PProjector::Modify(PParticle **mstack, int *decay_done, int * num, int stacksize) {
    //cout << "Modify " << endl;

    *w = current_weight;
    
    SetParticles(mstack, decay_done, num, stacksize, 1);
    //cout << "num:  " << *num << endl;
    //excuting batch
    Int_t startcommand = 0;
    Int_t retval = kFALSE;
    for (int i=0; i<batch_pos; i++) {

	//cout << "do now : " << i << endl;
	if (i<0) return -1;
	retval = batch[i]->Execute(startcommand, retval);
	//cout << "retval:" << retval << endl;
	startcommand = 0;

	if (retval & kTRUE) {
	    
	    current_weight = *w;

	    if (hist3[i] && fFillFlag[i]) {
		if (hist3[i]->GetSumw2()->GetSize() || fFillFlag[i]==1)
		    hist3[i]->Fill((*x), (*y), (*z), current_weight);
		else
		    hist3[i]->Fill((*x), (*y), (*z));
	    }
	    if (hist2[i] && fFillFlag[i]) {
		if (hist2[i]->GetSumw2()->GetSize() || fFillFlag[i]==1)
		    hist2[i]->Fill((*x), (*y), current_weight);
		else
		    hist2[i]->Fill((*x), (*y));
	    }
	    if (hist1[i] && fFillFlag[i]) {
		if (hist1[i]->GetSumw2()->GetSize() || fFillFlag[i]==1)
		    hist1[i]->Fill((*x), current_weight);
		else {
		    hist1[i]->Fill((*x));
		}
	    }
	    
	    if (fp_out[i]) {
		//fill the ntuple
		Double_t *val;

		for (int j=0; j<key_pos_out[i]; j++) {
		    if (key_list_out[i][j] > -1) {
			if (makeDataBase()->GetParamDouble(key_list_out[i][j], batch_value_param, &val))
			    values[j] = (Float_t)(*val);
			else values[j] = 0;
		    }
		    else values[j] = 0;
		}
		
		fp_out[i]->Fill(values);
	    }

	    if (fp_in[i]) {
		//read the ntuple
		if (num_events_in_c[i] == num_events_in[i]) {
            if (pluto_global::verbosity) {
                Info("Modify", "NTuple <%s>: number of events reached", fp_in[i]->GetTitle());
            }
		    return kFALSE;
		}

		fp_in[i]->GetEntry(num_events_in_c[i]);
		Double_t *val;

		for (int j=0; j<key_pos_in[i]; j++) {
		    if (key_list_in[i][j] > -1) {
			if (makeDataBase()->GetParamDouble(key_list_in[i][j], batch_value_param, &val))
			    *val = values_in[i][j];
		    }
		}
		num_events_in_c[i]++;
	    }
	} 
	
	int redo = 0;
	int setparticle = 0;
	if (retval & kGOTO) {
	    if (batch[i]->GetNewBulk() == bulk_id) {
		//stay in the same PProjector
		setparticle = 1;
		//SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles for "formore"
		Int_t new_batch = batch[i]->GetNewBatch();
		startcommand = batch[i]->GetNewCommand();
		i =  new_batch - 1;
	    } else {
		Error("Modify", "Jumping with a GOTO over different Projector not yet implemented");
	    }
	    retval &= ~kGOTO;
	} 
	if (retval & kPUSH) {
	    startcommand = batch[i]->GetOldCommand();
	    //cout << startcommand << endl;
	    if ((*num) >= stacksize) {
		Warning("Modify (kPUSH)", "Stack size too small, increase '_system_particle_stacksize'");
		return kTRUE;
	    }
	    if (batch[i]->GetCurrentParticle()) {
		*(mstack[*num]) = *(batch[i]->GetCurrentParticle());
		(*num)++;
		decay_done[*num]=0;
		setparticle = 2;
		//SetParticles(mstack, decay_done, num, stacksize, 1);  //reset particles
	    } else {
		Warning("Modify (kPUSH)", "Pushed particle not found");
	    }
	    //i -= 1; //stay in same batch
	    redo = 1;
	    retval &= ~kPUSH;
	} 
	if (retval & kPUSHBRANCH) {
	    startcommand = batch[i]->GetOldCommand();
	    if (batch[i]->GetBranch() > (*size_branches))  {
		Warning("Modify (kPUSHBRANCH)", "Branch number %i not known", *size_branches);
		return kTRUE;
	    }
	    if ( (*(current_size_branches[batch[i]->GetBranch()-1])) >= stacksize) {
		Warning("Modify (kPUSHBRANCH)", "Stack size too small, increase '_system_particle_stacksize'");
		return kTRUE;
	    }
	    
	    if (batch[i]->GetCurrentParticle()) {
		*((particle_array_branches[batch[i]->GetBranch()-1])
		  [*(current_size_branches[batch[i]->GetBranch()-1])] )
		    = *(batch[i]->GetCurrentParticle());
		(*(current_size_branches[batch[i]->GetBranch()-1]))++;
	    } else {
		Warning("Modify", "Pushed particle not found");
	    }
	    redo = 1;
	    retval &= ~kPUSHBRANCH;
	} 
	if (retval & kFOREACHEND) {
	    setparticle = 1;
	    //SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles like for "formore"
	    startcommand = 0;
	    retval &= ~kFOREACHEND;
	    retval &= ~kFOREACH;
	} 
	if (retval & kUPDATE) {
	    setparticle = 1;
	    //SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles like for "formore"
	    startcommand = batch[i]->GetOldCommand();
	    redo = 1;
	    retval &= ~kUPDATE;
	} 
	if (retval & kEOF) {
        if (pluto_global::verbosity) {
            Info("Modify", "EOF reached");
        }
	    return kFALSE;
	} 
	if (retval & kFOREACH) {
	    startcommand = batch[i]->GetOldCommand();
	    //cout << "sc foreach " << startcommand << endl;
	    redo = 1;
	    setparticle = 1;
	    //SetParticles(mstack, decay_done, num, stacksize, 0);  //reset particles like for "formore"
	    retval &= ~kFOREACHEND;
	} 
	if (retval & kELSE) {
	    if (batch[i]->GetElsePosition()>-1 &&  //only if else is provided
		batch[i]->GetCurrentPosition() < batch[i]->GetElsePosition() ) { //avoid deadlocks
		startcommand = batch[i]->GetElsePosition();
		redo = 1;
		//i -= 1; //stay in same batch		
	    }
	}
	if (setparticle == 1) SetParticles(mstack, decay_done, num, stacksize, 0);
	if (setparticle == 2) SetParticles(mstack, decay_done, num, stacksize, 1);
	if (redo) i -= 1; //stay in same batch	

	retval &= ~kTRUE; //clean true flag for next loop
    }

    Int_t *i_result;
    //Before leaving clean the max_list again, because next time we could have a different configuration of PIDS
    for (int i=0; i<*num; i++) {
	//Delete old max pos list
	Int_t particle_key = makeStaticData()->GetParticleKey(mstack[i]->ID());
	if (makeDataBase()->GetParamInt (particle_key, stream_max_pos_param, &i_result)) {
	    (*i_result) = 0;
	} 
    }

    return kTRUE;
}


void PProjector::Print(const Option_t*) const {
    for (int i=0; i<batch_pos; i++) 
	batch[i]->Print();
}


ClassImp(PProjector)
