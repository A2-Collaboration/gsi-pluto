////////////////////////////////////////////////////////
// Fast conversion gamma -> e+e- in material
//
//
//                    Author:  
//                    Written: 12/01/2014
//                    Revised: 
//
////////////////////////////////////////////////////////

#include "PGammaConversion.h"
#include "PChannel.h"

#include "TFile.h"
#include "TH3F.h"
#include "TSystem.h"

PGammaConversion::PGammaConversion() {
    ConvProb = 0;
    Z = 41;   // Sep08 p+Nb run
    convpar[0] = 0.4319;
    convpar[1] = -0.02322;
    convpar[2] = 0.000359;
    ep = PParticle("e+");
    em = PParticle("e-");
    
    for (Int_t i=0; i<NEFUNCS; i++) 
	Elepton[i] = NULL;
    RGen = new TRandom3();
    ulep = new TF1("ulep","x*exp(-0.625*x)*(1+27*exp(-1.25*x))",0,20); // lepton angle w.r.t. photon
    ulep->SetNpx(300);
    ulep->SetLineWidth(1);

    nVertex   = 0;
    nbinTheta = 0;
    nbinR     = 0;
    nbinZ     = 0;
    hdimension  = 0;
    hvertexgeom = 0;
    x1 = y1 = z1 = 0;

}

PGammaConversion::PGammaConversion(Float_t z, Float_t p, Int_t run = 0) {
    ConvProb = p;
    Z = z;

    // following coeffs are fitted to results of a complete GEANT3 simulation:  ProbConv(Ephoton, Theta)
    if(run==0) {
	convpar[0] = 0.4319;     // Sep08 p+Nb run coeffs.
	convpar[1] = -0.02322;
	convpar[2] = 0.000359;
    } else if(run==1) {
	convpar[0] = 0.649;    // Sep05 Ar+KCl run coeffs.
	convpar[1] = 0.00740;
	convpar[2] = 3.97e-5;
    } else if(run==2) {
    convpar[0] =  0.380904   ; // 12 Au+Au run coeffs.
    convpar[1] = -0.0777771 ;
    convpar[2] =  0.00522763;
    convpar[3] = -1.33942e-4;
    convpar[4] =  1.49021e-6;
    convpar[5] = -6.11297e-9;
    } else {
	cout<< " no correct run number" <<endl;
        return;
    }
    ep = PParticle("e+");
    em = PParticle("e-");

    for (Int_t i=0; i<NEFUNCS; i++) 
	Elepton[i] = NULL;
    RGen = new TRandom3();
    ulep = new TF1("ulep", "x*exp(-0.625*x)*(1+27*exp(-1.25*x))", 0, 20); // lepton angle w.r.t. photon (GEANT3)
    ulep->SetNpx(300);
    ulep->SetLineWidth(1);

    nVertex   = 0;
    nbinTheta = 0;
    nbinR     = 0;
    nbinZ     = 0;
    hdimension = 0;
    hvertexgeom = 0;
    x1 = y1 = z1 = 0;

}


PGammaConversion::~PGammaConversion() {  // clean up
    for (Int_t i=0; i<NEFUNCS; i++) 
	if (Elepton[i] != NULL) 
	    Elepton[i]->Delete();
    ulep->Delete();
    RGen->Delete();
}


Float_t PGammaConversion::ConversionXS(Float_t E) { // return total conversion cross section in mb/atom (Geant3)

    const Float_t EMass = 5.109990615e-4; // electron mass in GeV
    if ( E<2.*EMass ) return 0;     // unphysical region
    Float_t Eph = (E<10 ? E : 10.);
    
    //  Coefficients from Geant3 routine GPRSGA
    const Float_t CPar[18] = {0.87842E-3, -0.19625E-2,  0.12949E-2, -0.20028E-3,  0.12575E-4, -0.28333E-6,
			      -0.10342E-4,  0.17692E-4, -0.82391E-5,  0.13063E-5, -0.90815E-7,  0.23586E-8,
			      -0.45263E-3,  0.11161E-2, -0.86749E-3,  0.21773E-3, -0.20467E-4,  0.65372E-6};
    Float_t CC[3] = {0, 0, 0};
    Float_t  X = log(Eph/EMass);
    
    // initialising parameters of CC1, CC2, CC3
    for (Int_t i=0; i<3; i++) {
	Float_t D = 1;
	for(int j=0; j<6; j++){
	    Int_t jx = 6*i + j;
	    CC[i] = CC[i] + D*CPar[jx];
	    D = D*X;
	}
    }
    
    //  total cross section for element with Z
    return (Z*(Z+1.0) * (CC[0] + CC[1]*Z + CC[2]/Z));
}


Float_t PGammaConversion::GetConversionProb(Float_t E, Float_t th, Int_t run = 0) {
    // return conversion prob for photon of energy E [GeV] and polar angle theta [deg]
    if (th<0. || th>90) 
	return 0;      // out of HADES acceptance
    if (ConvProb >= 1.0) 
	return 1;      // convert all
    if (ConvProb > 0) 
	return ConvProb;   // use constant conversion probability 
    else {
	if (run != 2)  
	    return 0.01 * convpar[0] * (1. + convpar[1]*th + convpar[2]*th*th) * ConversionXS(E);
	if (run == 2)  
	    return 0.01 * convpar[0] * (1. + convpar[1]*th + convpar[2]*th*th + 
					convpar[3]*th*th*th+convpar[4]*th*th*th*th) * ConversionXS(E);
    }
    return 0.;
}

Double_t PGammaConversion::dSde(Double_t *x, Double_t *par) {  // electron energy distribution (GEANT3)
    Double_t E = par[0];     // photon energy in GeV
    Double_t Z = par[1];     // average atomic charge of converter material

    Double_t eps = x[0]/E;   // electron fraction of photon energy
    if (eps<=0 || eps>=1 || Z<1) 
	return 0;

    const Double_t alpha = 1./137.036;   // fine structure constant
    const Double_t me = 0.0005109989;    // electron mass in GeV
    
    Double_t a   = pow(alpha*Z,2);
    Double_t fcZ = alpha*(1/(1+alpha) + 0.20206 - 0.0369*a + 0.0083*a*a - 0.002*a*a*a); // Coulomb correction
    Double_t FZ  = E<0.05 ? 4./3.*log(Z) : 4./3.*log(Z)+4.*fcZ;

    Double_t delta = 136.*me/(pow(Z,1./3.)*E*eps*(1-eps));
    Double_t phi1, phi2;  // screening functions
    if (delta<=1) {
       phi1 = 20.867 - 3.242*delta + 0.625*delta*delta;
       phi2 = 20.209 - 1.930*delta - 0.086*delta*delta;
    } else {
       phi1 = phi2 = 21.12 - 4.184*log(delta+0.952);
    }

    phi1 = max(phi1, FZ);  // to avoid negative result for eps close to 0 or 1
    phi2 = max(phi2, FZ);

    Double_t xiZ = log(1440./pow(Z,2./3.))/(log(183./pow(Z,1./3.)) - fcZ);

    Double_t val = alpha*Z*(1+xiZ)/(E*E) * ((eps*eps+(1-eps)*(1-eps))*(phi1-FZ)
                                           + 2./3.*eps*(1-eps)*(phi2-FZ));
    return val;
}

Int_t PGammaConversion::readHist(TString input) {
    
    if(gSystem->AccessPathName(input.Data())!=0) {
	Error("readHist()", "File %s not found!", input.Data());
        return -1;
    }

    TFile *file = new TFile(input, "READ");
    Int_t size = -1;

    hvertexgeom = (TH1F*) gDirectory->Get("hvertex");

    if(!hvertexgeom) {
	Error("readHist()", "hvertex hist for target geom not found! hopeless ....");
	exit(1);
    }

    nVertex = hvertexgeom->GetNbinsX();
    if(file) {
	for(Int_t binVertex=0; binVertex<nVertex; binVertex++) {
	    TH3F *h = (TH3F*)gDirectory->Get(Form("hVertex_full%i", binVertex+1));
	    if(h) {
		nbinTheta = h->GetNbinsZ();
		nbinR     = h->GetNbinsY();
		nbinZ     = h->GetNbinsX();
		
		if(size == -1) {
		    hdimension = h;

		    size = nVertex*nbinTheta;
		    for(Int_t i=0; i<size; i++){
			fhzr.push_back(0);
		    }
		    cout << "Histogram binning : theta " << nbinTheta << " r " << 
			nbinR << " z " << nbinZ << endl;
		}

		for(Int_t binTheta=0; binTheta<nbinTheta; binTheta++) {
		    h->GetZaxis()->SetRange(binTheta+1, binTheta+1);
		    fhzr[binVertex*nbinTheta+binTheta] = (TH2F*) h->Project3D("xy");
		    fhzr[binVertex*nbinTheta+binTheta]->SetName(Form("%s_%i_xy", h->GetName(), binTheta));
		    fhzr[binVertex*nbinTheta+binTheta]->SetDirectory(0);
		}

	    } else {
		Error("readHist()", "hvertex %i not found!", binVertex);
	    }

	}
    } else {
	Error("readHist()", "File point NULL!");
        return -1;
    }
    return 0;
}

Int_t PGammaConversion::findVertexBin(Float_t zVertexEvt) {
    Int_t bin = -1;
    if(hvertexgeom) {
	bin = hvertexgeom->GetXaxis()->FindBin(zVertexEvt);
	if(bin > hvertexgeom->GetNbinsX()) 
	    bin = -1;
	else 
	    bin -= 1;
        return bin;
    } else return bin;
}

Int_t PGammaConversion::getThetaBin(Double_t thetadeg) {

    if(hdimension) {
	Int_t bin = hdimension->GetZaxis()->FindBin(thetadeg);
	if(thetadeg==0||thetadeg>=nbinTheta) 
	    return -1;
        else 
	    return bin;
    }
    return -1;
}

TH2F *PGammaConversion::getHist(Double_t zVertexEvt, Double_t thetadeg)
{
    if(hdimension) {
	Int_t thbin = getThetaBin(thetadeg);
	if(thbin < 0) 
	    return 0;

        Int_t binVertex = findVertexBin(zVertexEvt);
        if(binVertex < 0) 
	    return 0;

	return fhzr[binVertex*nbinTheta+thbin] ;
    }
    return 0;
}

Bool_t PGammaConversion::Modify(PParticle **mstack, int *decay_done, int *num, int stacksize) {
    // Read the particles from the defined stack and copy this to the official
    // particle stream

    const Double_t R2D = 57.295780;    // rad to deg
    const Double_t me = 0.0005109989;  // electron mass in GeV
    
    int cur = *num;
    int gamma_id = makeStaticData()->GetParticleID("g");
    
    //    cout << cur << "   " << mstack[0]->ID() << endl;
    
    for (int i=0; i<cur; i++) {
	
	if (mstack[i]->HasID(gamma_id) && mstack[i]->IsActive()) {  // this is a photon
	    
	    if ((*num) == stacksize) {
		Warning("Modify", "stacksize reached");
		return kFALSE;
	    }

	    //make kinematics here....
	    PParticle *track = mstack[i];
	    TVector3 vertex  = track->GetVertex();

	    Double_t Ephot  = track->E();
	    Double_t thphot = R2D*track->Theta();

	    if ( thphot<0. || thphot>90 ) continue;
            /*
	    Double_t convprob = GetConversionProb(Ephot, thphot);  // get conversion probability
	    if ( RGen->Uniform()>convprob ) {  // do not convert this one
		decay_done[i] = 0;   // do not suppress unconverted photons from output
		continue;
	    }
            */
	    // Int_t mesonId = track->GetParentId();

	    if (Ephot < 2*me) 
		return kFALSE;
	    Int_t index = int(Ephot/0.01);          // select energy slice of 0.01 GeV/c
	    if (index < 0) 
		return kFALSE;
	    index = min(index, NEFUNCS-1);

	    if (Elepton[index] == NULL) {           // this slice does not yet exist
		Char_t funame[11] = "Elepton000";
		sprintf(funame, "Elepton%03d", index); // create function object
		Elepton[index] = new TF1(funame, PGammaConversion::dSde, me, double(index)*0.01+0.005, 2);
		Elepton[index]->SetNpx(200);
		Elepton[index]->SetParameters(Ephot, Z);  // and initialize it
	    }

	    // Do energy and angle sampling of the e+/e- pair

	    Double_t Eelec = Elepton[index]->GetRandom();     // lepton total energy
	    Eelec *= Ephot/(double(index)*0.01+0.005);        // correct for discretization of sampling
	    Double_t Pelec = sqrt(Eelec*Eelec - me*me);       // lepton momentum
	    Double_t Thelec = (me/Ephot) * ulep->GetRandom(); // lepton polar angle (w.r.t. photon)
	    Double_t Eposi = Ephot - Eelec;              // assume energy conservation
	    Double_t Pposi = sqrt(Eposi*Eposi - me*me);
	    Double_t arg = (Pelec/Pposi) * sin(Thelec);  // assume coplanar photon e+ e-
	    if (arg<-1.) arg = -1;
	    if (arg> 1.) arg =  1;
	    Double_t Thposi = TMath::ASin(arg);

	    Double_t Phelec = 6.28318531 * RGen->Uniform();  // sample uniform phi angle
	    Double_t Phposi = 6.28318531 - Phelec;           // assume that pair momentum is colinear
	    //       Double_t close = R2D * (Thelec + Thposi);        // opening angle of conversion pair
	    TVector3 Dirphot(track->Vect());  // photon direction
	    Dirphot = Dirphot.Unit();         // normalize photon vector

	    TVector3 lep3(Pelec*cos(Phelec)*sin(Thelec), Pelec*sin(Phelec)*sin(Thelec), Pelec*cos(Thelec));
	    TVector3 lep4(Pposi*cos(Phposi)*sin(Thposi), Pposi*sin(Phposi)*sin(Thposi), Pposi*cos(Thposi));
	    lep3.RotateUz(Dirphot);  // rotate lepton vectors into lab frame
	    lep4.RotateUz(Dirphot);

	    /*
	     Float_t ptot3 = lep3.Mag();
	     Float_t theta3 = R2D*lep3.Theta();
	     Float_t phi3 = R2D*TMath::ACos(lep3.X()/lep3.Perp());
	     if (lep3.Y() < 0.) phi3 = 360.-phi3;
	     Float_t ptot4 = lep4.Mag();
	     Float_t theta4 = R2D*lep4.Theta();
	     Float_t phi4 = R2D*TMath::ACos(lep4.X()/lep4.Perp());
	     if (lep4.Y() < 0.) phi4 = 360.-phi4;
	     */

	    Int_t id3, id4;
	    if (RGen->Uniform() < 0.5) { // decide which is e- and which is e+
		id3 = 3;
		id4 = 2;
	    } else {
		id3 = 2;
		id4 = 3;
	    }


	    ep.SetVectM(lep3,me);
	    ep.SetID(id3);
	    ep.SetSourceId(track->GetSourceId());
	    ep.SetParentId(1);
	    ep.SetParentIndex(track->GetIndex());
	    ep.SetW(track->W());
	    ep.SetVertex(0, 0, 0, 0);
	    ep.SetActive();

	    em.SetVectM(lep4,me);
	    em.SetID(id4);
	    em.SetSourceId(track->GetSourceId());
	    em.SetParentId(1);
	    em.SetParentIndex(track->GetIndex());
	    em.SetW(track->W());
	    em.SetVertex(0, 0, 0, 0);
	    em.SetActive();


	    if (nVertex > 0) {
		//----------------------------------
		// Convert vertex into target num
		TH2F *hzr = getHist(vertex.Z(), thphot);
		
		if(!hzr) {
		    Error("Modify()", "hzr at theta %f and vertex %f not found !", thphot, vertex.Z());
		    continue;
		}
		
		Double_t z,r;
		hzr->GetRandom2(z, r);
		
		if(!TMath::Finite(z) || !TMath::Finite(r)) {
		    ep.SetVertex(0, 0, 10000, 0);
		    em.SetVertex(0, 0, 10000, 0);
		} else {		    
		    Double_t phi = gRandom->Rndm()*2*TMath::Pi();
		    Double_t x   = cos(phi)*r;
		    Double_t y   = sin(phi)*r;;
		    ep.SetVertex(x, y, z, 0);
		    em.SetVertex(x, y, z, 0);		    
		    
		    cout<<" x " << x << " y " << y << " phi " << phi*TMath::RadToDeg() << 
			" r " << r << " th " << thphot << endl;		    
		}
	    } else {
		ep.SetVertex(x1, y1, z1, 0);
		em.SetVertex(x1, y1, z1, 0);	
	    }

 	    *(mstack[*num]) = ep;
	    (*num)++;
	    *(mstack[*num]) = em;
	    (*num)++;
	    decay_done[i] = 1;
	} 
    }
    
    return kTRUE;

}

ClassImp(PGammaConversion)
