/////////////////////////////////////////////////////////////////////
//
// Beam line transport
//
// Author:  T. Hennino (core), J. Biernat
//          Pluto implementation by IF
// Written: 14.05.2013
//
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "PBeamLineSimulation.h"
#include "PDefines.h"


ClassImp(PBeamLineSimulation)

PBeamLineSimulation::PBeamLineSimulation() {
};

PBeamLineSimulation::~PBeamLineSimulation() {
    if (mybeam)   delete (mybeam);
    if (mytarget) delete (mytarget);
};


PBeamLineSimulation::PBeamLineSimulation(Char_t *id, Char_t *de) :
    PBeamSmearing(id, de) {
    beam     = NULL;
    mybeam   = NULL;
    target   = NULL;
    mytarget = NULL;
    parent   = NULL;

    beam_counter     = 0;
    detector_counter = 0;
    num_target = -1;
    global_p = 0.;
    beamtube_size_x = beamtube_size_y = -1.;

    ii    = 0;
    ind1  = 0;
    ind23 = 0;
};

PDistribution* PBeamLineSimulation::Clone(const char *) const {
    return new PBeamLineSimulation((const PBeamLineSimulation &)* this);
};

void PBeamLineSimulation::MakeVars() {
    if (!projector) projector = new PProjector();

    in_xi =  makeStaticData()->GetBatchValue("_beam_x"); 
    in_yi =  makeStaticData()->GetBatchValue("_beam_y"); 
    in_dp =  makeStaticData()->GetBatchValue("_beam_dp"); 
    in_xpi = makeStaticData()->GetBatchValue("_beam_px"); 
    in_ypi = makeStaticData()->GetBatchValue("_beam_py"); 
    
    filter = makeStaticData()->GetBatchValue("#_beamline_filter");
};

Bool_t PBeamLineSimulation::AddEquation(char *command) {
    MakeVars();
    return projector->AddCommand(command);
};

Bool_t PBeamLineSimulation::AddHistogram(TH2 *histo, const char *command) {
    MakeVars();
    return  projector->AddHistogram(histo, command, 0);
};

Int_t PBeamLineSimulation::AddBeamParticle(char *name) {

    //create a clone of the string
    char * dummy = new char[strlen(name)+1];
    strcpy(dummy, name);

    //create key in data base:
    Int_t key1 = makeStaticData()->
	MakeDirectoryEntry("batch_objects", NBATCH_NAME, LBATCH_NAME, dummy);

    //check if we use this key already
    int used_id = -1;
    if (makeGlobalBulk()->GetSizeBranches() && makeGlobalBulk()->GetKeysBranches()) {
	for (int br=0; br<*(makeGlobalBulk()->GetSizeBranches()); br++) {
	    if ((makeGlobalBulk()->GetKeysBranches())[br] == key1) used_id = br;
	}
    } else {
	Error("AddBeamParticle","No access to branches (no PReaction?)");
    }
    if (used_id == -1) {
	if (*(makeGlobalBulk()->GetSizeBranches()) != MAX_NUM_BRANCHES) {
	    (makeGlobalBulk()->GetKeysBranches()[*(makeGlobalBulk()->GetSizeBranches())]) = key1;
	    used_id = *(makeGlobalBulk()->GetSizeBranches());
	    (*(makeGlobalBulk()->GetSizeBranches()))++;
	} else {
	    Error("AddDetector","MAX_NUM_BRANCHES reached");
	    return -1;
	}
    }
    
    //generate beam particles and put them to the stack

    if (beam_counter == MAX_NUM_BEAM_PARTICLES) {
	Error("AddBeamParticle", "MAX_NUM_BEAM_PARTICLES reached");
	return -1;
    }

    beam_particles[beam_counter] = new PParticle();
    beam_id[beam_counter] = used_id;
    beam_counter++;
    
    return beam_counter;
};

Int_t PBeamLineSimulation::AddDetector(char *name, Double_t distance) {
    //adds a detector "name" at place "distance" [in mm]
    //(distance must be negative if placed before the target)
    
    if (detector_counter == MAX_NUM_BEAM_PARTICLES) {
	Error("AddDetector","MAX_NUM_BEAM_PARTICLES reached");
	return -1;
    }

    detector_name[detector_counter]     = name;
    detector_distance[detector_counter] = distance;
    detector_counter++;

    return detector_counter - 1;
}

Bool_t PBeamLineSimulation::TargetIsElement(Int_t n) {
    if (n<1 || n>=nelem) return kFALSE;
    num_target = n;
    for (int i=1; i<nelem; i++)
	along[i] -= along[n];

    return kTRUE;
}

Bool_t PBeamLineSimulation::InitBeamLine(char *filename) {
    //opens the ASCII description of the beam line

    nelem = 1;

    std::ifstream plik1;
    plik1.open(filename);
    string str;

    Bool_t readfile = kTRUE;
    
    while (readfile) {
	
	str.clear();	
	int i = nelem;

	plik1>>along[i]>>txt2>>aa>>sig11>>txt5;
	plik1>>aa>>sig22>>txt5>>r12;
	plik1>>aa>>sig33>>txt5>>r13>>r23;
	plik1>>aa>>sig44>>txt5>>r14>>r24>>r34;
	plik1>>aa>>sig55>>txt5>>r15>>r25>>r35>>r45;
	plik1>>aa>>sig66>>txt5>>r16>>r26>>r36>>r46>>r56;
	
	getline(plik1,str,'\n'); 
	getline(plik1,str,'\n');
	
	plik1>>T11[i]>>T12[i]>>T13[i]>>T14[i]>>aa>>T16[i];
	plik1>>T21[i]>>T22[i]>>T23[i]>>T24[i]>>aa>>T26[i];
	plik1>>T31[i]>>T32[i]>>T33[i]>>T34[i]>>aa>>T36[i]>>std::fixed;
	// std::cout<<T31[i]<<endl;
	plik1>>T41[i]>>T42[i]>>T43[i]>>T44[i]>>aa>>T46[i];
	plik1>>aa>>aa>>aa>>aa>>aa>>aa;
	plik1>>aa>>aa>>aa>>aa>>aa>>aa;
	getline(plik1,str);
	
	//C
	//C     lecture second ordre T1ij
	//C
	getline(plik1,str);
	plik1>>ind1>>ind23>>T111[i];
	plik1>>ind1>>ind23>>T112[i]>>ind1>>ind23>>T122[i];
	plik1>>ind1>>ind23>>T113[i]>>ind1>>ind23>>T123[i]>>ind1>>ind23>>T133[i];
	plik1>>ind1>>ind23>>T114[i]>>ind1>>ind23>>T124[i]>>ind1>>ind23>>T134[i]>>ind1>>ind23>>T144[i];
	plik1>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa;
	plik1>>ind1>>ind23>>T116[i]>>ind1>>ind23>>T126[i]>>ind1>>ind23>>T136[i]>>ind1>>ind23>>T146[i]>>ind1>>ind23>>aa>>ind1>>ind23>>T166[i];
	getline(plik1,str);
	//plik2<<".."<<i<<".."<<T166[i]<<endl;
	
	//C
	//C     lecture second ordre T2ij
	//C
	plik1>>ii>>ii>>T211[i];
	plik1>>ii>>ii>>T212[i]>>ii>>ii>>T222[i];
	plik1>>ii>>ii>>T213[i]>>ii>>ii>>T223[i]>>ii>>ii>>T233[i];
	plik1>>ii>>ii>>T214[i]>>ii>>ii>>T224[i]>>ii>>ii>>T234[i]>>ii>>ii>>T244[i];
	plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
	plik1>>ii>>ii>>T216[i]>>ii>>ii>>T226[i]>>ii>>ii>>T236[i]>>ii>>ii>>T246[i]>>ii>>ii>>aa>>ii>>ii>>T266[i]; 
	getline(plik1,str);
	//plik2<<".."<<i<<".."<<T266[i]<<endl;
	//rT266 = T266[i];
	
	//C
	//C     lecture second ordre T3ij
	//C
	plik1>>ind1>>ind23>>T311[i];
	plik1>>ind1>>ind23>>T312[i]>>ind1>>ind23>>T322[i];
	plik1>>ind1>>ind23>>T313[i]>>ind1>>ind23>>T323[i]>>ind1>>ind23>>T333[i];
	plik1>>ind1>>ind23>>T314[i]>>ind1>>ind23>>T324[i]>>ind1>>ind23>>T334[i]>>ind1>>ind23>>T344[i];
	plik1>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa;
	plik1>>ind1>>ind23>>T316[i]>>ind1>>ind23>>T326[i]>>ind1>>ind23>>T336[i]>>ind1>>ind23>>T346[i]>>ind1>>ind23>>aa>>ind1>>ind23>>T366[i];
	getline(plik1,str);
	//plik2<<".."<<i<<".."<<T366[i]<<endl;
	//rT366 = T366[i];
	
	//C
	//C     lecture second ordre T4ij
	//C
	
	plik1>>ii>>ii>>T411[i];
	plik1>>ii>>ii>>T412[i]>>ii>>ii>>T422[i];
	plik1>>ii>>ii>>T413[i]>>ii>>ii>>T423[i]>>ii>>ii>>T433[i];
	plik1>>ii>>ii>>T414[i]>>ii>>ii>>T424[i]>>ii>>ii>>T434[i]>>ii>>ii>>T444[i];
	plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
	plik1>>ii>>ii>>T416[i]>>ii>>ii>>T426[i]>>ii>>ii>>T436[i]>>ii>>ii>>T446[i]>>ii>>ii>>aa>>ii>>ii>>T466[i];
	getline(plik1,str);
	//plik2<<".."<<i<<".."<<T466[i]<<endl;
	//rT466 = T466[i];
	//C
	//C     lecture second ordre T5ij
	//C
	
	plik1>>ii>>ii>>aa;
	plik1>>ii>>ii>>aa>>ii>>ii>>aa;
	plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
	plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
	plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
	plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
	getline(plik1,str);
	getline(plik1,str,'\n');
	
	if (!plik1.eof()) 
	    nelem++;
	else 
	    readfile = kFALSE;
	
    } //end nelem
  
    plik1.close();

    Info("InitBeamLine", "Read %i elements from file %s", nelem-1, filename);

    return kTRUE;
}


Bool_t PBeamLineSimulation::Init(void) {
    //Read beam, target and parent

    beam   = GetParticle("beam");
    target = GetParticle("target");
    parent = GetParticle("parent");

    if (!beam || !target) {
	Error("Init", "beam or target not found");
	return kFALSE;
    }

    if (!mybeam)
	mybeam = new PParticle(beam);
    else
	*mybeam = *beam;
    if (!mytarget)
	mytarget= new PParticle(target);
    else
	*mytarget = target;

    if (!parent){
	Error("Init","Parent not found");
	return kFALSE;
    }
  
    for (int i=0; i<detector_counter; i++) {
	AddBeamParticle(detector_name[i]);
    }

    if (!projector) {
	Error("Init", "PProjector not found, you have not specified the initial conditions");
	return kFALSE;
    }

    return kTRUE;
}

Bool_t PBeamLineSimulation::Prepare(void) {
    //We use the prepare function for smearing since it might affect
    //the mass and momentum sampling done in the next steps
    
    *filter = 1;
    
    //restore saved particles
    *beam   = *mybeam;
    *target = *mytarget;

    PParticle *real_beam = target;

    if (beam->P() > target->P()) {
	//change "beam"
	real_beam = beam;
    } 

    //the global event vertex cannot be set at this point.
    //this would require a major change in PReaction

    //change momentum of "real_beam" here.....

    //beam offset (in cm) at prod target in x
    //beam offset (in cm) at prod target in y
    Int_t num = 0;
    projector->Modify(NULL, NULL, &num, 0);
   
    xi  = *in_xi/10.;    //convert from Pluto standard (mm) to cm
    yi  = *in_xi/10.;    //convert from Pluto standard (mm) to cm
    xpi = *in_xpi*1000.; //convert from rad to mrad
    ypi = *in_ypi*1000.; //convert from rad to mrad
    dp  = *in_dp * 100.; //convert from dp/p to percent (as discussed with T.H. at CM Prague)
    //dp  = 0.08;
    
    //now loop over the number of elements
    for(Int_t k = 1; k < nelem; k++){

	z1 = T11[k]*xi;
	z2 = T12[k]*xpi;
	z3 = T13[k]*yi;
	z4 = T14[k]*ypi;
	z5 = T16[k]*dp;

	xf = z1+z2+z3+z4+z5;
	yf = T31[k]*xi+T32[k]*xpi+T33[k]*yi+T34[k]*ypi+T36[k]*dp;  //yf=T31(k)*xi+T32(k)*xpi+T33(k)*yi+T34(k)*ypi+T36(k)*dp
	//  std::cout<<"T31 = .. "<<T31[k]<<"..T32 = .. "<<T32[k]<<"..T33 = .. "<<T33[k]<<"..T34 = .. "<<T34[k]<<"..T36 = .. "<<T36[k]<<std::fixed<<endl;
	xpf  = T21[k]*xi+T22[k]*xpi+T23[k]*yi+T24[k]*ypi+T26[k]*dp;
	ypf  = T41[k]*xi+T42[k]*xpi+T43[k]*yi+T44[k]*ypi+T46[k]*dp;
	xf1  = xf;
	yf1  = yf;
	xpf1 = xpf;
	ypf1 = ypf;
       
	xf = xf1+T111[k]*xi*xi+T112[k]*xi*xpi+T113[k]*xi*yi;
	xf = xf+T114[k]*xi*ypi+T116[k]*xi*dp;
	xf = xf+T122[k]*xpi*xpi+T123[k]*xpi*yi+T124[k]*xpi*ypi;
	xf = xf+T126[k]*xpi*dp;
	xf = xf+T133[k]*yi*yi+T134[k]*yi*ypi+T136[k]*yi*dp;
	xf = xf+T144[k]*ypi*ypi+T146[k]*ypi*dp;
	xf = xf+T166[k]*dp*dp;
	
	// xf=xf1+T126[k]*xpi*dp+T166[k]*dp*dp;
	//  xf=xf1+T126[k]*xpi*dp;

	xpf = xpf1+T211[k]*xi*xi+T212[k]*xi*xpi+T213[k]*xi*yi;
	xpf = xpf+T214[k]*xi*ypi+T216[k]*xi*dp;
	xpf = xpf+T222[k]*xpi*xpi+T223[k]*xpi*yi+T224[k]*xpi*ypi;
	xpf = xpf+T226[k]*xpi*dp;
	xpf = xpf+T233[k]*yi*yi+T234[k]*yi*ypi+T236[k]*yi*dp;
	xpf = xpf+T244[k]*ypi*ypi+T246[k]*ypi*dp;
	xpf = xpf+T266[k]*dp*dp;
	
	yf = yf1+T311[k]*xi*xi+T312[k]*xi*xpi+T313[k]*xi*yi;
	yf = yf+T314[k]*xi*ypi+T316[k]*xi*dp;
	yf = yf+T322[k]*xpi*xpi+T323[k]*xpi*yi+T324[k]*xpi*ypi;
	yf = yf+T326[k]*xpi*dp;
	yf = yf+T333[k]*yi*yi+T334[k]*yi*ypi+T336[k]*yi*dp;
	yf = yf+T344[k]*ypi*ypi+T346[k]*ypi*dp;
	yf = yf+T366[k]*dp*dp;
	
	ypf = ypf1+T411[k]*xi*xi+T412[k]*xi*xpi+T413[k]*xi*yi;
	ypf = ypf+T414[k]*xi*ypi+T416[k]*xi*dp;
	ypf = ypf+T422[k]*xpi*xpi+T423[k]*xpi*yi+T424[k]*xpi*ypi;
	ypf = ypf+T426[k]*xpi*dp;
	ypf = ypf+T433[k]*yi*yi+T434[k]*yi*ypi+T436[k]*yi*dp;
	ypf = ypf+T444[k]*ypi*ypi+T446[k]*ypi*dp;
	ypf = ypf+T466[k]*dp*dp;
	//C
	//C      yf=yf1+T346[k]*ypi*dp+T366[k]*dp*dp
	//C      yf=yf1+T346[k]*ypi*dp
	// rxf = 
	xf2  = xf;
	xpf2 = xpf;
	yf2  = yf;
	ypf2 = ypf;
	
	out_xi[k]  = xf2*10.;    //convert back from cm. to mm
	out_yi[k]  = yf2*10.;    //convert back from cm. to mm
	out_dp[k]  = dp/100.;    //convert back from % to fraction
	out_xpi[k] = xpf2/1000.; //convert back from mrad to rad
	out_ypi[k] = ypf2/1000.; //convert back from mrad to rad
	
	//cout << out_xi[k] << ":" << out_yi[k] << ":" << out_xpi[k] << ":" << out_ypi[k] << endl;

    }

    Double_t p   = global_p * (1. + out_dp[num_target]);
    Double_t p_z = p / sqrt(1.+ tan(out_xpi[num_target])*tan(out_xpi[num_target]) + 
			    tan(out_ypi[num_target])*tan(out_ypi[num_target]));
    Double_t p_x = p_z * tan(out_xpi[num_target]);
    Double_t p_y = p_z * tan(out_ypi[num_target]);
    
    real_beam->SetXYZM(p_x, p_y, p_z, makeStaticData()->GetParticleMass(real_beam->ID()));
    
    //real_beam->Print();
    
    parent->Reconstruct();
    
    //parent->Print();
    
    //change particles at detector places
    for (int i=0; i<beam_counter; i++) {
	
	for (int k=1; k<nelem-1; k++) {
	    if (along[k]*1000. <= detector_distance[i] && detector_distance[i] < along[k+1]*1000.) {
		
		Double_t frac = (detector_distance[i] - along[k]*1000.) / (along[k+1]*1000. - along[k]*1000.);
		
		Double_t my_xi = out_xi[k] + frac*(out_xi[k+1] - out_xi[k]);
		Double_t my_yi = out_yi[k] + frac*(out_yi[k+1] - out_yi[k]);
		Double_t my_dp = out_dp[k] + frac*(out_dp[k+1] - out_dp[k]);
		Double_t my_xpi = out_xpi[k] + frac*(out_xpi[k+1] - out_xpi[k]);
		Double_t my_ypi = out_ypi[k] + frac*(out_ypi[k+1] - out_ypi[k]);
		p   = global_p * (1. + my_dp);
		
		//cout << p << ":" << my_xpi << ":" << my_ypi << endl;
		
		p_z = p / sqrt(1.+ tan(my_xpi)*tan(my_xpi) + tan(my_ypi)*tan(my_ypi));
		p_x = p_z * tan(my_xpi);
		p_y = p_z * tan(my_ypi);
		*(beam_particles[i]) = *(real_beam);
		beam_particles[i]->SetVertex(my_xi, my_yi, detector_distance[i]);
		beam_particles[i]->SetPx(p_x);
		beam_particles[i]->SetPy(p_y);
		beam_particles[i]->SetPz(p_z);
		
		//cout << p_x << ":" << p_y << ":" << p_z << endl;
		
	    }
	    
	}
    }
    
    for (int k=1; k<nelem-1; k++) {
    	if (beamtube_size_x > 0 && out_xi[k] > beamtube_size_x) {
	    *filter = 0;
	}		
	if (beamtube_size_y > 0 && out_yi[k] > beamtube_size_y) {
	    *filter = 0;
	}
    }
    
    for (int i=0; i<beam_counter; i++)
	(*(makeGlobalBulk()->GetBranchNum(beam_id[i]))) = 0;
    
    
    for (int i=0; i<beam_counter; i++) {
	*((makeGlobalBulk()->GetBranchArray(beam_id[i]))
	  [*(makeGlobalBulk()->GetBranchNum(beam_id[i]))] )
	    = *(beam_particles[i]);
	(*(makeGlobalBulk()->GetBranchNum(beam_id[i])))++;
    }

    return kTRUE;
}

