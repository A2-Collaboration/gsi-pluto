////////////////////////////////////////////////////////
// dileptons (D0,D+,w,rho0,phi,eta,eta' and pi0)
// can be decayed via shining (translated for original UrQMD
// fortran routines, authors: Katharina Schmidt, Sascha Vogel,
// Christian Sturm 2007-2009) into e+/e-.
//
// USAGE: This classs is called via PHUrReader, see urqmd_f15_input.C
//
////////////////////////////////////////////////////////

#include "PHUrDilep.h"
#include "TMath.h"
#include "PParticle.h"
#include "PHUrReader.h"

ClassImp(PHUrDilep)


PHUrDilep::PHUrDilep() {
    omega    = 103;
    rho      = 104;
    phi      = 109;
    pion     = 101;
    etaprime = 107;
    eta      = 102;
    delta    = 17;
    omegadir = 1031;
    omegadal = 1032;

    alpha_em          =  (1./137.036);
    vacmass_pi0       = 0.1349764; // pi0
    mass_electron     = 0.00051099906; // e-
    mass_nucleon      = 0.93827231; // proton
    vacmass_eta       = 0.54745; // eta
    vacmass_etaprime  = 0.9577; // eta'
    vacmass_rho       = 0.7699; // rho0
    vacmass_omega     = 0.78194; // omega
    vacmass_phi       = 1.019413; // phi
    vacmass_delta     = 1.2320; // D+

    gev            = 0.2;
    g              = 5.44;     // delta

    lambda_omega   = 0.65;
    gamma_omega    = 0.04;
    gamma_photon   = 0.00075561;

    br_pi0         = 0.988;
    b_pi0          = 5.5;

    lambda_eta     = 0.72;
    br_eta         = 0.3933;

    lambda_etaprime = 0.76;
    gamma_etaprime  = 0.1;
    br_etaprime     = 0.0212;

    out           = NULL;
    outputLeptons = kTRUE;
    reader        = NULL;
}

PHUrDilep::~PHUrDilep() {
    if(out) fclose(out);
    out = 0;
}

void PHUrDilep::Output(TString infile,TString outfile) {
    if(outfile.CompareTo("") != 0) {
	out = fopen(outfile.Data(), "w");
	fprintf(out, "# InputFile  : %s\n", infile.Data());
	fprintf(out, "# OutputFile : %s\n", outfile.Data());
	fprintf(out, "# ityp weight          mres            p0_el_lab       px_el_lab       py_el_lab       pz_el_lab       p0_po_lab       px_po_lab       py_po_lab       pz_po_lab       tau             time_start      time_end       flag  dens_cre        dens_abs             noe  bev             acce            accp\n");
    }
}

void PHUrDilep::Dilep() {

    Double_t tau = (particleIn->t - particleOut->t)/particleIn->Gamma();

    // multi: multiplier for better statistics for rare resonances
    Int_t multi = 1;

    //     data input

    // tbeam_input is a char ... in GeV ... quick and dirty work around !!! Take care !!!
    // calculate beta from input tbeam

    Double_t ebeam = (evtheader->E_beam+mass_nucleon);
    Double_t pbeam = sqrt(pow(ebeam,2) - pow(mass_nucleon,2));

    Double_t beta  = pbeam/(ebeam+mass_nucleon);
    Double_t gamma = 1./sqrt(1.-beta*beta);

    Double_t weight = 0;                      // changed in dgamma_sum()
    Double_t smin, smax, tmax;                // changed in dgamma_sum()
    smin = smax = tmax = 0;

    Double_t s, t;                           // changed in gamma_star()
    Double_t tgamma;                         // changed in gamma_star()
    Double_t xgstar, ygstar, zgstar, tgstar; // changed in gamma_star()
    Double_t xgamma, ygamma, zgamma;         // changed in gamma_star()
    
    Double_t dgamma;       // changed in dalpi(),daleta(),daletaprime(),dalomega()
    Double_t mgstar;       // mgstar=s s from gamma_star()
    Double_t p0_gstar, p0_nucleon, p0_pion, p0_gamma; // local variables

    //     ccccc calculation of decay widths ccccc
    dgamma_sum(particleIn->id, tau, multi,
	       weight, smin, smax, tmax);  // output

    if(!outputLeptons) {
	if( particleIn->id == omega                            ||
	   (particleIn->id == rho   &&  particleIn->chrg == 0) ||
	    particleIn->id == phi                              ||
           (particleIn->id == pion  &&  particleIn->chrg == 0) ||
            particleIn->id == eta                              ||
            particleIn->id == etaprime                         ||
           (particleIn->id == delta &&  (particleIn->chrg == 0 || particleIn->chrg == 1) )
	  ) OutputROOT(particleIn->id,weight);
	return;
    }


    if (particleIn->id == omega)  {  //  decay of omega meson (direct and Dalitz)
        Double_t weight_omegadir = 0;
	for(Int_t n = 1; n <= multi; n++) {
	    diromega(tau, particleOut->mass, multi,
		     weight_omegadir);          // output

	    lobo_dir(beta, gamma, particleOut->mass); // em,ep output

	    Output(omegadir, weight_omegadir, particleOut->mass, multi,tau);

	} // for multi omega
	Int_t validm = 0;

	while(validm < multi){
	    t_omega(tau, multi,
		    smin, smax, tmax);  // output

	    gamma_star(smin, smax, tmax,
		       s, t,                             // output
		       xgstar, ygstar, zgstar, tgstar,   // output
		       xgamma, ygamma, zgamma, tgamma);  // output

	    dalomega(tau, s,
		     dgamma);                           // output

	    if (t <= dgamma) {
		mgstar = s;

		//    energy of the gamma* and the pion due to energy and
		//    momentum conservation

		p0_gstar = particleOut->mass/2. -(pow(vacmass_pi0,2) - pow(mgstar,2))/(2.*particleOut->mass);
		p0_pion  = particleOut->mass/2. +(pow(vacmass_pi0,2) - pow(mgstar,2))/(2.*particleOut->mass);

		//    check if energy is larger than mass (only for omega, for
		//    low mass omegas it might happen that p0_gstar < m_gstar

		if(p0_gstar > mgstar) {
		    validm++;
		    
		    lobo_dal(p0_gstar, p0_pion,
			     mgstar, beta, gamma,
			     xgstar, ygstar, zgstar, tgstar, xgamma, ygamma, zgamma, particleOut->mass);// em,ep output

		    Output(omegadal, weight, mgstar, multi, tau);

		} // if p0
	    }
	} // while vaidm < multi

    //----------------------------------------------------------------------
    } else if (particleIn->id == rho && particleIn->chrg == 0) {   // direct decay of the rho meson

	for(Int_t n = 1; n <= multi; n++){
	    dirrho(tau, particleOut->mass, multi,
		   weight);       // output

	    lobo_dir(beta, gamma, particleOut->mass);// em,ep output

	    Output(rho, weight, particleOut->mass, multi, tau);

	} //  for multi rho

    //----------------------------------------------------------------------
    } else if (particleIn->id == phi) {   // direct decay of the phi meson

	for(Int_t n = 1; n <= multi; n++){
	    dirphi(tau, particleOut->mass,multi,
		   weight);   // output

	    lobo_dir(beta, gamma, particleOut->mass);//  em,ep output

	    Output(phi, weight, particleOut->mass, multi, tau);

	} // for multi phi

    //----------------------------------------------------------------------
    } else if (particleIn->id == pion && particleIn->chrg == 0) { //  Dalitz decay of the pi0 meson

	Int_t validm = 0;

	while(validm < multi) {
	    gamma_star(smin, smax, tmax,
		       s, t,                           // output
		       xgstar, ygstar, zgstar, tgstar, // output
		       xgamma, ygamma, zgamma, tgamma);// output

	    dalpi(s,
		  dgamma);        // output

	    if (t <= dgamma) {
		Double_t mass_gamma = 0;
		mgstar              = s;

		p0_gstar = particleOut->mass/2 - (pow(mass_gamma,2) - pow(mgstar,2))/(2*particleOut->mass);
		p0_gamma = particleOut->mass/2 + (pow(mass_gamma,2) - pow(mgstar,2))/(2*particleOut->mass);

		validm = validm+1;

		lobo_dal(p0_gstar,p0_gamma,
			 mgstar,beta,gamma,
			 xgstar,ygstar,zgstar,tgstar,xgamma,ygamma,zgamma,particleOut->mass);// em,ep output

		Output (pion,weight,mgstar,multi,tau);
	    } // if t <= dgamma
	} // end while validm < multi

    //----------------------------------------------------------------------
    } else if (particleIn->id == etaprime) { //  Dalitz decay of the eta prime meson

	Int_t validm = 0;

	while (validm < multi) {
	    gamma_star(smin, smax, tmax,
		       s, t,                                 // output
		       xgstar, ygstar, zgstar, tgstar,       // output
		       xgamma, ygamma, zgamma, tgamma);      // output

	    daletaprime(s,
			dgamma);                            // output

	    if (t <= dgamma) {
		Double_t mass_gamma = 0;
		mgstar              = s;

		p0_gstar = particleOut->mass/2. - (pow(mass_gamma,2) - pow(mgstar,2))/(2.*particleOut->mass);
		p0_gamma = particleOut->mass/2. + (pow(mass_gamma,2) - pow(mgstar,2))/(2.*particleOut->mass);

		validm++;

		lobo_dal(p0_gstar, p0_gamma,
			 mgstar, beta, gamma,
			 xgstar, ygstar, zgstar, tgstar, xgamma, ygamma, zgamma, particleOut->mass);// em,ep output

		Output(etaprime, weight, mgstar, multi, tau);

	    } // if t <= dgamma
	} // end while validm < multi


    //----------------------------------------------------------------------
    } else if (particleIn->id == eta) { //  Dalitz decay of the eta meson
	Int_t validm = 0;

	while(validm < multi) {
	    gamma_star(smin, smax, tmax,
		       s, t,                            // output
		       xgstar, ygstar, zgstar, tgstar,  // output
		       xgamma, ygamma, zgamma, tgamma); // output

	    daleta(s,
		   dgamma);                            // output

	    if (t <= dgamma){
		Double_t mass_gamma = 0;
		mgstar              = s;

		p0_gstar = particleOut->mass/2. - (pow(mass_gamma,2) - pow(mgstar,2))/(2.*particleOut->mass);
		p0_gamma = particleOut->mass/2. + (pow(mass_gamma,2) - pow(mgstar,2))/(2.*particleOut->mass);

		validm++;

		lobo_dal(p0_gstar, p0_gamma, mgstar, beta, gamma,
			 xgstar, ygstar, zgstar, tgstar, xgamma, ygamma, zgamma, particleOut->mass);//  em,ep output

		Output (eta, weight, mgstar, multi, tau);

	    } // if t <= dgamma
	} // end validm < multi

    //----------------------------------------------------------------------
    } else if((particleIn->id == delta && particleIn->chrg == 0) ||
	      (particleIn->id == delta && particleIn->chrg == 1)) {  //     Dalitz decay of the Delta(1232) baryon

	Int_t validm = 0;

	dgamma_sum_delta(tau, particleOut->mass, multi,
			 weight,            // output
			 smin, smax, tmax); // output

	while(validm < multi) {

	    t_delta(tau, particleOut->mass, multi,
		    tmax);                  // output

	    gamma_star(smin, smax, tmax,
		       s, t,                                   // output
		       xgstar, ygstar, zgstar, tgstar,         // output
		       xgamma, ygamma, zgamma, tgamma);        // output

	    daldelta(tau, s, particleOut->mass,
		     dgamma);                                 // output


	    if (t <= dgamma) {
		
		mgstar = s ;

		p0_gstar   = particleOut->mass/2. - (pow(mass_nucleon,2) - pow(mgstar,2))/(2.*particleOut->mass);
		p0_nucleon = particleOut->mass/2. + (pow(mass_nucleon,2) - pow(mgstar,2))/(2.*particleOut->mass);

		validm++;

		lobo_dal(p0_gstar, p0_nucleon,
			 mgstar, beta, gamma,
			 xgstar, ygstar, zgstar, tgstar, xgamma, ygamma, zgamma, particleOut->mass);// em,ep output

		Output(delta, weight, mgstar, multi, tau);

	    } // if t <= dgamma
	} // end validm < multi
    } // end if case
}

void PHUrDilep::dalpi(Double_t mx, Double_t &dgamma) {
    // form factor from Landsberg 1985 (Phys.Rept.128,301-376)
   
    Double_t t1 = ((4.*alpha_em)/(3.*TMath::Pi()*mx))*
	sqrt(1.-(4.*pow(mass_electron,2))/(pow(mx,2)))*
	(1.+(2.*pow(mass_electron,2))/(pow(mx,2)));
    Double_t t2 = pow((1.-(pow(mx,2))/(pow(vacmass_pi0,2))),3);
    Double_t f1 = pow((1.+b_pi0*pow(mx,2)),2);
    dgamma      = t1*t2*f1*br_pi0;
}

void PHUrDilep::daleta(Double_t mx, Double_t &dgamma) {
                               
    //### integration dgamma, for more details please refer to (REFERENCE, MANUAL) ####

    Double_t t1 =((4.*alpha_em)/(3.*TMath::Pi()*mx))*
	sqrt(1.-(4.*pow(mass_electron,2))/(pow(mx,2)))*(1.+
							(2.*pow(mass_electron,2))/(pow(mx,2)));
    Double_t t2 = pow((1.-(pow(mx,2))/(pow(vacmass_eta,2))),3);

    //### form factor from Phys.Rept.308:65-233,1999 with a modification
    // from the diploma thesis of C.Ernst to ensure that F(0)=1  ###
    Double_t f1 = pow((1./(1.-(pow(mx,2))/(pow(lambda_eta,2)))),2);

    dgamma = t1*t2*f1*br_eta;
}

void PHUrDilep::daletaprime(Double_t mx, Double_t &dgamma) {

    //#### integration dgamma, for more details please refer to (REFERENCE, MANUAL) ####

    // ######### Landsberg 85
    Double_t t1 = ((4.*alpha_em)/(3.*TMath::Pi()*mx))
	*sqrt(1.-(4.*pow(mass_electron,2))/(pow(mx,2)))*(1.+
							 (2.*pow(mass_electron,2))/(pow(mx,2)));
    Double_t t2 = pow((1.-(pow(mx,2))/(pow(vacmass_etaprime,2))),3);

    //### form factor from Phys.Rept.308:65-233,1999 with a modification from the diploma thesis of C.Ernst to ensure that F(0)=1  ###

    Double_t  f1 = (pow(lambda_etaprime,2)*
		    (pow(lambda_etaprime,2)+pow(gamma_etaprime,2)))/
	(pow((pow(lambda_etaprime,2)-pow(mx,2)),2)+
	 (pow(lambda_etaprime,2)*pow(gamma_etaprime,2)));
    
    dgamma = t1*t2*f1*br_etaprime;
}

void PHUrDilep::daldelta(Double_t tau, Double_t mx, Double_t mres, Double_t &dgamma) {
    //
    // masses are in Gev/c2

    //    ######### Landsberg 85

    Double_t pf = sqrt((pow(mres,2)-pow((mass_nucleon+mx),2))*
		     (pow(mres,2)-pow((mass_nucleon-mx),2)))/2./mres;

    Double_t q0 = sqrt(pow(mx,2)+pow(pf,2));


    Double_t f   = (-3./2.)*(mres+mass_nucleon)/(mass_nucleon*
				    (pow((mres+mass_nucleon),2)-pow(mx,2)));

    Double_t e = sqrt(4.*TMath::Pi()*alpha_em);

    Double_t mt = pow(e*f*g,2)*pow(mres,2)/9./mass_nucleon *(q0*q0*(5.*mres-3.*
						   (q0+mass_nucleon))-pow(mx,2)*(mres+mass_nucleon+q0));

    Double_t ml = pow(e*f*g,2)*pow(mres,2)/9./mass_nucleon*pow(mx,2)*4.*
	(mres-mass_nucleon-q0);

    Double_t lambda = pow(mx,4)+pow(mass_nucleon,4)+pow(mres,4)-2.*
	(pow(mx,2)*pow(mass_nucleon,2)+pow(mx,2)*pow(mres,2)+pow(mass_nucleon,2)*pow(mres,2));

    Double_t gamma0 = sqrt(lambda)/(16.*TMath::Pi()*pow(mres,2))*mass_nucleon*(2.*mt+ml);
    dgamma = (2.*alpha_em)/(3.*TMath::Pi()*mx)*gamma0*(tau/gev);

}

void PHUrDilep::dalomega(Double_t tau, Double_t mx, Double_t &dgamma) {
    
    //######### Landsberg 85

    Double_t  t1 = ((2.*alpha_em)/(3.*TMath::Pi()*mx))*
	sqrt(1.-(4.*pow(mass_electron,2))/(pow(mx,2)))*(1.+
						    (2.*pow(mass_electron,2))/(pow(mx,2)));

    Double_t t2  = sqrt(pow((1.+(pow(mx,2))/(pow(vacmass_omega,2)-pow(vacmass_pi0,2))),2)-
			pow(((2.*vacmass_omega*mx)/(pow(vacmass_omega,2)-pow(vacmass_pi0,2))),2));

    Double_t f1 =(pow(lambda_omega,2)*(pow(lambda_omega,2)+pow(gamma_omega,2)))/
	(pow((pow(lambda_omega,2)-pow(mx,2)),2)+(pow(lambda_omega,2)*pow(gamma_omega,2)));

    dgamma = t1*pow(t2,3)*f1*gamma_photon*(tau/gev);
}

void  PHUrDilep::diromega(Double_t tau, Double_t mres, Int_t multi, Double_t &weight) {

    const Double_t cv = 0.000000778869;

    Double_t dwidth = cv/pow(mres,3)*pow(vacmass_omega,4);
    Double_t br     = dwidth*(tau/gev);

    weight = br/multi;
}

void PHUrDilep::dirphi(Double_t, Double_t mres, Int_t multi, Double_t &weight) {
            
    const Double_t gamma_tot = 0.00426;
    const Double_t cv        = 0.0000012411;

    Double_t dwidth = cv/pow(mres,3)*pow(vacmass_phi,4);
    Double_t br = dwidth/gamma_tot;

    weight      = br/multi;
}

void PHUrDilep::dirrho(Double_t tau, Double_t mres, Int_t multi, Double_t &weight) {
      
    const Double_t cv = 0.0000090545;

    Double_t dwidth = (cv/ pow(mres,3))*pow(vacmass_rho,4);
    Double_t br     =  dwidth*(tau/gev);

    weight = br/multi;  
}

void PHUrDilep::lobo_dal(Double_t p0_gstar, Double_t,
	      Double_t m_gstar, Double_t beta, Double_t gamma,
	      Double_t x_gstar, Double_t y_gstar, Double_t z_gstar, Double_t,
	      Double_t, Double_t, Double_t,
	      Double_t) {
    Double_t p0_res = particleOut->E();
    Double_t px_res = particleOut->Px();
    Double_t py_res = particleOut->Py();
    Double_t pz_res = particleOut->Pz();

// #### energies of the gamma* and the particle follow from energy and momentan conservation in the decay process ####

    Double_t p_gstar  = sqrt(p0_gstar*p0_gstar - m_gstar*m_gstar);
    Double_t px_gstar = x_gstar*p_gstar;
    Double_t py_gstar = y_gstar*p_gstar;
    Double_t pz_gstar = z_gstar*p_gstar;

    Double_t phi_el = 2.*TMath::Pi()*rndfunc();
    Double_t x      = rndfunc();
    Double_t sin_theta_el = (x-0.5)*2.;
    Double_t theta_el     = asin(sin_theta_el)+TMath::PiOver2();

    Double_t phi_po = phi_el + TMath::Pi() ;

    if(phi_po > 2*TMath::Pi()){
	phi_po = phi_po-2*TMath::Pi();
    }
    Double_t theta_po = 0;
    if(theta_el <= TMath::Pi()){
	theta_po = TMath::Pi()-theta_el;
    } else {
	theta_po = 0;
    }

    Double_t x_el = sin(theta_el)*cos(phi_el);
    Double_t y_el = sin(theta_el)*sin(phi_el);
    Double_t z_el = cos(theta_el);

    Double_t x_po = sin(theta_po)*cos(phi_po);
    Double_t y_po = sin(theta_po)*sin(phi_po);
    Double_t z_po = cos(theta_po);


    Double_t p0_el = m_gstar/2.;
    Double_t p0_po = m_gstar/2.;

    Double_t p_el = sqrt(p0_el*p0_el-mass_electron*mass_electron);
    Double_t p_po = sqrt(p0_po*p0_po-mass_electron*mass_electron);


    Double_t px_el = x_el*p_el;
    Double_t py_el = y_el*p_el;
    Double_t pz_el = z_el*p_el;

    Double_t px_po = x_po*p_po;
    Double_t py_po = y_po*p_po;
    Double_t pz_po = z_po*p_po;


    Double_t beta_x_gstar = -px_gstar/p0_gstar;
    Double_t beta_y_gstar = -py_gstar/p0_gstar;
    Double_t beta_z_gstar = -pz_gstar/p0_gstar;

    Double_t beta_gstar  = sqrt(beta_x_gstar*beta_x_gstar+beta_y_gstar*beta_y_gstar+beta_z_gstar*beta_z_gstar);
    Double_t gamma_gstar = 1./sqrt(1-beta_gstar*beta_gstar);

    Double_t p0_el_res = gamma_gstar*p0_el - beta_x_gstar*gamma_gstar*px_el -
	beta_y_gstar*gamma_gstar*py_el -beta_z_gstar*gamma_gstar*pz_el;

    Double_t px_el_res = (-beta_x_gstar*gamma_gstar)*p0_el +
	(1+(gamma_gstar-1)*((beta_x_gstar*beta_x_gstar)/
			    (beta_gstar*beta_gstar)))*px_el +
	(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
			 (beta_gstar*beta_gstar))*py_el +
	(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*pz_el;

    Double_t py_el_res = (-beta_y_gstar*gamma_gstar)*p0_el +
	(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
			 (beta_gstar*beta_gstar))*px_el +
	(1+(gamma_gstar-1)*((beta_y_gstar*beta_y_gstar)/
			    (beta_gstar*beta_gstar)))*py_el +
	(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*pz_el;

    Double_t pz_el_res = (-beta_z_gstar*gamma_gstar)*p0_el +
	(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*px_el +
	(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*py_el +
	(1+(gamma_gstar-1)*((beta_z_gstar*beta_z_gstar)/
			    (beta_gstar*beta_gstar)))*pz_el;

    Double_t p0_po_res = gamma_gstar*p0_po - beta_x_gstar*gamma_gstar*px_po -
	beta_y_gstar*gamma_gstar*py_po -beta_z_gstar*gamma_gstar*pz_po;

    Double_t px_po_res = (-beta_x_gstar*gamma_gstar)*p0_po +
	(1+(gamma_gstar-1)*((beta_x_gstar*beta_x_gstar)/
			    (beta_gstar*beta_gstar)))*px_po +
	(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
			 (beta_gstar*beta_gstar))*py_po +
	(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*pz_po;

    Double_t py_po_res = (-beta_y_gstar*gamma_gstar)*p0_po +
	(gamma_gstar-1)*((beta_x_gstar*beta_y_gstar)/
			 (beta_gstar*beta_gstar))*px_po +
	(1+(gamma_gstar-1)*((beta_y_gstar*beta_y_gstar)/
			    (beta_gstar*beta_gstar)))*py_po +
	(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*pz_po;

    Double_t pz_po_res = (-beta_z_gstar*gamma_gstar)*p0_po +
	(gamma_gstar-1)*((beta_x_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*px_po +
	(gamma_gstar-1)*((beta_y_gstar*beta_z_gstar)/
			 (beta_gstar*beta_gstar))*py_po +
	(1+(gamma_gstar-1)*((beta_z_gstar*beta_z_gstar)/
			    (beta_gstar*beta_gstar)))*pz_po;



    Double_t beta_x_res = -px_res/p0_res;
    Double_t beta_y_res = -py_res/p0_res;
    Double_t beta_z_res = -pz_res/p0_res;

    Double_t beta_res  = sqrt(beta_x_res*beta_x_res+beta_y_res*beta_y_res+beta_z_res*beta_z_res);

    Double_t gamma_res = 1./sqrt(1-beta_res*beta_res);

    Double_t p0_el_cm = gamma_res*p0_el_res-beta_x_res*gamma_res*px_el_res-
	beta_y_res*gamma_res*py_el_res-beta_z_res*gamma_res*pz_el_res;

    Double_t px_el_cm = (-beta_x_res*gamma_res)*p0_el_res+
	(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
			  (beta_res*beta_res)))*px_el_res +
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*py_el_res +
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*pz_el_res;

    Double_t py_el_cm = (-beta_y_res*gamma_res)*p0_el_res+
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*px_el_res+
	(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
			  (beta_res*beta_res)))*py_el_res+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*pz_el_res;

    Double_t pz_el_cm = (-beta_z_res*gamma_res)*p0_el_res+
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*px_el_res+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*py_el_res+
	(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
			  (beta_res*beta_res)))*pz_el_res;


    Double_t p0_po_cm = gamma_res*p0_po_res-beta_x_res*gamma_res*px_po_res-
	beta_y_res*gamma_res*py_po_res -beta_z_res*gamma_res*pz_po_res;

    Double_t px_po_cm = (-beta_x_res*gamma_res)*p0_po_res+
	(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
			  (beta_res*beta_res)))*px_po_res+
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*py_po_res+
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*pz_po_res;

    Double_t py_po_cm = (-beta_y_res*gamma_res)*p0_po_res+
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*px_po_res+
	(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
			  (beta_res*beta_res)))*py_po_res+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*pz_po_res;

    Double_t pz_po_cm = (-beta_z_res*gamma_res)*p0_po_res+
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*px_po_res+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*py_po_res+
	(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
			  (beta_res*beta_res)))*pz_po_res;


    Double_t beta_x_urqmd = 0;
    Double_t beta_y_urqmd = 0;
    Double_t beta_z_urqmd = -beta;


    Double_t beta_urqmd  = beta;
    Double_t gamma_urqmd = gamma;

    // output

    Double_t p0_el_lab = gamma_urqmd*p0_el_cm-beta_x_urqmd*gamma_urqmd*px_el_cm-
	beta_y_urqmd*gamma_urqmd*py_el_cm-beta_z_urqmd*gamma_urqmd*
	pz_el_cm;

    Double_t px_el_lab = (-beta_x_urqmd*gamma)*p0_el_cm +
	(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
			    (beta_urqmd*beta_urqmd)))*px_el_cm +
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_el_cm +
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_el_cm;

    Double_t py_el_lab = (-beta_y_urqmd*gamma)*p0_el_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_el_cm+
	(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
			    (beta_urqmd*beta_urqmd)))*py_el_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_el_cm;

    Double_t pz_el_lab = (-beta_z_urqmd*gamma)*p0_el_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_el_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_el_cm+
	(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
			    (beta_urqmd*beta_urqmd)))*pz_el_cm;

    Double_t p0_po_lab = gamma_urqmd*p0_po_cm-beta_x_urqmd*gamma_urqmd*px_po_cm-
	beta_y_urqmd*gamma_urqmd*py_po_cm -beta_z_urqmd*gamma_urqmd*
	pz_po_cm;

    Double_t px_po_lab = (-beta_x_urqmd*gamma)*p0_po_cm+
	(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
			    (beta_urqmd*beta_urqmd)))*px_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_po_cm;

    Double_t py_po_lab = (-beta_y_urqmd*gamma)*p0_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_po_cm+
	(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
			    (beta_urqmd*beta_urqmd)))*py_po_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_po_cm;

    Double_t pz_po_lab = (-beta_z_urqmd*gamma)*p0_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_po_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_po_cm+
	(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
			    (beta_urqmd*beta_urqmd)))*pz_po_cm;

       em.SetPxPyPzE(px_el_lab,py_el_lab,pz_el_lab,p0_el_lab);
       ep.SetPxPyPzE(px_po_lab,py_po_lab,pz_po_lab,p0_po_lab);
}

void PHUrDilep::lobo_dir(Double_t beta, Double_t gamma, Double_t m_res) {

    Double_t p0_res = particleOut->E();
    Double_t px_res = particleOut->Px();
    Double_t py_res = particleOut->Py();
    Double_t pz_res = particleOut->Pz();

    Double_t phi_el       = 2.*TMath::Pi()*rndfunc();
    Double_t x            = rndfunc();
    Double_t sin_theta_el = (x-0.5)*2.;
    Double_t theta_el     = asin(sin_theta_el)+TMath::PiOver2();

    Double_t phi_po = phi_el + TMath::Pi();
    if(phi_po > 2*TMath::Pi()) {
	phi_po = phi_po-2*TMath::Pi();
    }
    Double_t theta_po = 0;
    if(theta_el <= TMath::Pi()) {
	theta_po = TMath::Pi()-theta_el;
    } else {
	theta_po = 0;
    }

    Double_t x_el = sin(theta_el)*cos(phi_el);
    Double_t y_el = sin(theta_el)*sin(phi_el);
    Double_t z_el = cos(theta_el);

    Double_t x_po = sin(theta_po)*cos(phi_po);
    Double_t y_po = sin(theta_po)*sin(phi_po);
    Double_t z_po = cos(theta_po);


    Double_t p0_el = m_res/2.;
    Double_t p0_po = m_res/2.;

    Double_t p_el = sqrt(pow(p0_el,2)-pow(mass_electron,2));
    Double_t p_po = sqrt(pow(p0_po,2)-pow(mass_electron,2));

    Double_t px_el = x_el*p_el;
    Double_t py_el = y_el*p_el;
    Double_t pz_el = z_el*p_el;

    Double_t px_po = x_po*p_po;
    Double_t py_po = y_po*p_po;
    Double_t pz_po = z_po*p_po;


    Double_t beta_x_res = -px_res/p0_res;
    Double_t beta_y_res = -py_res/p0_res;
    Double_t beta_z_res = -pz_res/p0_res;

    Double_t beta_res  = sqrt(beta_x_res*beta_x_res+beta_y_res*beta_y_res+beta_z_res*beta_z_res);

    Double_t gamma_res = 1./sqrt(1-beta_res*beta_res);

    Double_t p0_el_cm = gamma_res*p0_el-beta_x_res*gamma_res*px_el-
	beta_y_res*gamma_res*py_el-beta_z_res*gamma_res*pz_el;

    Double_t  px_el_cm = (-beta_x_res*gamma_res)*p0_el+
	(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
			  (beta_res*beta_res)))*px_el+
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*py_el+
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*pz_el;

    Double_t py_el_cm = (-beta_y_res*gamma_res)*p0_el+
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*px_el+
	(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
			  (beta_res*beta_res)))*py_el+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*pz_el;

    Double_t pz_el_cm = (-beta_z_res*gamma_res)*p0_el+
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*px_el+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*py_el+
	(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
			  (beta_res*beta_res)))*pz_el;

    Double_t p0_po_cm = gamma_res*p0_po-beta_x_res*gamma_res*px_po-
	beta_y_res*gamma_res*py_po -beta_z_res*gamma_res*pz_po;

    Double_t px_po_cm = (-beta_x_res*gamma_res)*p0_po+
	(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
			  (beta_res*beta_res)))*px_po+
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*py_po+
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*pz_po;

    Double_t py_po_cm = (-beta_y_res*gamma_res)*p0_po+
	(gamma_res-1)*((beta_x_res*beta_y_res)/
		       (beta_res*beta_res))*px_po+
	(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
			  (beta_res*beta_res)))*py_po+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*pz_po;

    Double_t pz_po_cm = (-beta_z_res*gamma_res)*p0_po+
	(gamma_res-1)*((beta_x_res*beta_z_res)/
		       (beta_res*beta_res))*px_po+
	(gamma_res-1)*((beta_y_res*beta_z_res)/
		       (beta_res*beta_res))*py_po+
	(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
			  (beta_res*beta_res)))*pz_po;


    Double_t beta_x_urqmd = 0;
    Double_t beta_y_urqmd = 0;
    Double_t beta_z_urqmd = -beta;

    Double_t beta_urqmd  = beta;
    Double_t gamma_urqmd = gamma;

    // output

    Double_t p0_el_lab = gamma_urqmd*p0_el_cm-beta_x_urqmd*gamma_urqmd*px_el_cm-
	beta_y_urqmd*gamma_urqmd*py_el_cm-beta_z_urqmd*gamma_urqmd*
	pz_el_cm;

    Double_t px_el_lab = (-beta_x_urqmd*gamma)*p0_el_cm +
	(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
			    (beta_urqmd*beta_urqmd)))*px_el_cm +
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_el_cm +
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_el_cm;

    Double_t py_el_lab = (-beta_y_urqmd*gamma)*p0_el_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_el_cm+
	(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
			    (beta_urqmd*beta_urqmd)))*py_el_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_el_cm;

    Double_t pz_el_lab = (-beta_z_urqmd*gamma)*p0_el_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_el_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_el_cm+
	(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
			    (beta_urqmd*beta_urqmd)))*pz_el_cm;

    Double_t p0_po_lab = gamma_urqmd*p0_po_cm-beta_x_urqmd*gamma_urqmd*px_po_cm-
	beta_y_urqmd*gamma_urqmd*py_po_cm -beta_z_urqmd*gamma_urqmd*
	pz_po_cm;

    Double_t px_po_lab = (-beta_x_urqmd*gamma)*p0_po_cm+
	(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
			    (beta_urqmd*beta_urqmd)))*px_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_po_cm;

    Double_t py_po_lab = (-beta_y_urqmd*gamma)*p0_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_po_cm+
	(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
			    (beta_urqmd*beta_urqmd)))*py_po_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*pz_po_cm;

    Double_t pz_po_lab = (-beta_z_urqmd*gamma)*p0_po_cm+
	(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*px_po_cm+
	(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
			 (beta_urqmd*beta_urqmd))*py_po_cm+
	(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
			    (beta_urqmd*beta_urqmd)))*pz_po_cm;

     em.SetPxPyPzE(px_el_lab,py_el_lab,pz_el_lab,p0_el_lab);
     ep.SetPxPyPzE(px_po_lab,py_po_lab,pz_po_lab,p0_po_lab);

}

void  PHUrDilep::dgamma_sum(Int_t ityp, Double_t tau, Int_t multi,
			    Double_t &weight,
			    Double_t &smin, Double_t &smax, Double_t &tmax) {
    const Double_t dx  = 0.001;

    Double_t dgamma = 0, dgamma_sum = 0;
    Double_t mx = 0, t1 = 0, t2 = 0, f1 = 0;
    Int_t z = 0;
    if(ityp == pion) { //     pi0
 
	smin = 0.0013;
	smax = 0.1348;
	tmax = 1.9065;

	z = (Int_t)((smax-smin)/dx);

	for(Int_t i = 1; i <= z; i++) {
	    mx = smin+dx*i;

	    t1 = ((4*alpha_em)/(3.*TMath::Pi()*mx)) *
		sqrt(1-(4*pow(mass_electron,2))/(pow(mx,2)))*(1+(2*pow(mass_electron,2))/(pow(mx,2)));

	    t2 = pow((1-(pow(mx,2))/(pow(vacmass_pi0,2))),3);
	    f1 = pow((1+b_pi0*pow(mx,2)),2);

	    dgamma     = t1*t2*f1*br_pi0;
	    dgamma_sum = dgamma_sum + dgamma;
	}

	weight = dgamma_sum*dx/multi;
	
    } else if(ityp == omega) {  //     omega

	smin = 0.0013;
	smax = 0.6446;
	tmax = 0.086086;

	z = (Int_t)((smax-smin)/dx);

	for(Int_t i=1; i <= z; i++) {
	    mx = smin+dx*i;

	    t1 = ((2*alpha_em)/(3*TMath::Pi()*mx))*
		sqrt(1-(4*pow(mass_electron,2))/
		     (pow(mx,2)))*(1+(2*pow(mass_electron,2))/(pow(mx,2)));

	    t2 = sqrt(pow((1+(pow(mx,2))/(pow(vacmass_omega,2)-
					  pow(vacmass_pi0,2))),2)
		      -pow(((2*vacmass_omega*mx)/(pow(vacmass_omega,2)-pow(vacmass_pi0,2))),2));

	    f1 = (pow(lambda_omega,2)*(pow(lambda_omega,2)+pow(gamma_omega,2)))/
		(pow((pow(lambda_omega,2)-pow(mx,2)),2)+
		 (pow(lambda_omega,2)*pow(gamma_omega,2)));

	    dgamma     = t1*pow(t2,3)*f1*gamma_photon*(tau/gev);
	    dgamma_sum = dgamma_sum + dgamma;
	}

	weight = dgamma_sum*dx/multi;

    } else if(ityp == etaprime) { //    eta'

	smin = 0.0013;
	smax = 0.958;
	tmax = 0.04092;

        z = (Int_t)((smax-smin)/dx);

	for(Int_t i=1; i <= z; i++) {
	    mx = smin+dx*i;

	    t1 = ((4*alpha_em)/(3*TMath::Pi()*mx))*
		sqrt(1-(4*pow(mass_electron,2))/
		     (pow(mx,2)))*(1+(2*pow(mass_electron,2))/(pow(mx,2)));

	    t2 = pow((1-(pow(mx,2))/(pow(vacmass_etaprime,2))),3);

	    f1 = (pow(lambda_etaprime,2)*
		  (pow(lambda_etaprime,2)+pow(gamma_etaprime,2)))/
		(pow((pow(lambda_etaprime,2)-pow(mx,2)),2)+
		 (pow(lambda_etaprime,2)*pow(gamma_etaprime,2)));

	    dgamma     = t1*t2*f1*br_etaprime;
	    dgamma_sum = dgamma_sum + dgamma;

	}

	weight=dgamma_sum*dx/multi;

    } else if(ityp == eta) { //    eta

	smin = 0.0013;
	smax = 0.547;
	tmax = 0.757592;

	z = (Int_t)((smax-smin)/dx);

	for(Int_t i=1; i <= z;i++) {
	    mx = smin+dx*i;

	    t1 = ((4.*alpha_em)/(3.*TMath::Pi()*mx))*
		sqrt(1.-(4.*pow(mass_electron,2))/
		     (pow(mx,2))) * (1.+(2.*pow(mass_electron,2))/(pow(mx,2)));

	    t2 = pow((1.-(pow(mx,2))/(pow(vacmass_eta,2))),3);
	    f1 = pow((1./(1.-(pow(mx,2))/(pow(lambda_eta,2)))),2);

	    dgamma     = t1*t2*f1*br_eta;
	    dgamma_sum = dgamma_sum + dgamma;

	}

	weight = dgamma_sum*dx/multi;
    }
}

void PHUrDilep::dgamma_sum_delta(Double_t tau,Double_t mres,Int_t multi,
				 Double_t& weight_delta,
				 Double_t& smin_del,Double_t& smax_del,Double_t& tmax_del) {

      const Double_t dx = 0.000001;

      //  delta

      Double_t mx, pf, q0, e, f, mt, ml, lambda, gamma0, dgamma, dgamma_sum=0;
      //Double_t q,qr;
      //     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
  
      smin_del = 0.0013;
      tmax_del = 0.00903;
      smax_del = mres-mass_nucleon;


      Int_t z_del = (Int_t)((smax_del-smin_del)/dx);

      for(Int_t i=1; i <= z_del; i++) {
      
	  mx = smin_del+dx*i;
      
	  pf = sqrt((pow(mres,2)-pow((mass_nucleon+mx),2))*
		  (pow(mres,2)-pow((mass_nucleon-mx),2)))/2./mres;

	  q0 = sqrt(mx*mx+pf*pf);

/*	  q = sqrt((pow(mres,2)-pow((mass_nucleon+vacmass_pi0),2))*
		   (pow(mres,2)-pow((mass_nucleon-vacmass_pi0),2))/(4.*pow(mres,2)));

	  qr = sqrt((pow(vacmass_delta,2)-pow((mass_nucleon+vacmass_pi0),2))*
		    (pow(vacmass_delta,2)-pow((mass_nucleon-vacmass_pi0),2))/
		    (4.*pow(vacmass_delta,2)));
*/
	  f = (-3./2.)*(mres+mass_nucleon)/
	      (mass_nucleon*(pow((mres+mass_nucleon),2)-mx*mx));
      
	  e = sqrt(4.*TMath::Pi()*alpha_em);
      
	  mt = pow((e*f*g),2)*pow(mres,2)/9./mass_nucleon*
	      (pow(q0,2)*(5.*mres-3.*(q0+mass_nucleon))-
	       pow(mx,2)*(mres+mass_nucleon+q0));
      
	  ml = pow((e*f*g),2)*pow(mres,2)/9./mass_nucleon*mx*mx*4.*
	      (mres-mass_nucleon-q0);
      
	  lambda = pow(mx,4)+pow(mass_nucleon,4)+pow(mres,4)-2.*
	      (mx*mx*pow(mass_nucleon,2)+mx*mx*pow(mres,2)+pow(mass_nucleon,2)*mres*mres);


	  if(fabs(lambda) <= 1e-15) {
	      lambda = 0;
	  }
            
	  if(mx <= q0) {
	      gamma0     = sqrt(lambda)/(16.*TMath::Pi()*mres*mres)*mass_nucleon*(2.*mt+ml);
	      dgamma     = (2.*alpha_em)/(3.*TMath::Pi()*mx)*gamma0*(tau/gev);
	      dgamma_sum = dgamma_sum + dgamma;
	  } else {
	      dgamma_sum = 0; // added
	  }
      }
 
      weight_delta = dgamma_sum*dx/multi;
}

void  PHUrDilep::gamma_star(Double_t smin, Double_t smax, Double_t tmax,
			    Double_t& s, Double_t& t,
			    Double_t& xgstar, Double_t& ygstar, Double_t& zgstar, Double_t& tgstar,
			    Double_t& xparticle, Double_t& yparticle, Double_t& zparticle, Double_t& tparticle) {

    Double_t phi_gstar        = 2.*TMath::Pi()*rndfunc();
    Double_t y                = rndfunc();
    Double_t sin_theta_gstar  = (y-0.5)*2.;
    Double_t theta_gstar      = asin(sin_theta_gstar)+TMath::PiOver2();

    Double_t phi_particle     = phi_gstar + TMath::Pi();
    if (phi_particle > 2*TMath::Pi()) {
	phi_particle = phi_particle-2*TMath::Pi();
    }
    Double_t theta_particle = 0;
    if (theta_gstar <= TMath::Pi()) {
	theta_particle = TMath::Pi()-theta_gstar;
    } else {
	theta_particle = 0;
    }

    xgstar = sin(theta_gstar)*cos(phi_gstar);
    ygstar = sin(theta_gstar)*sin(phi_gstar);
    zgstar = cos(theta_gstar);
    tgstar = sin(theta_gstar);

    xparticle = sin(theta_particle)*cos(phi_particle);
    yparticle = sin(theta_particle)*sin(phi_particle);
    zparticle = cos(theta_particle);
    tparticle = sin(theta_particle);


    t = rndfunc()*tmax;
    s = smin + rndfunc()*(smax-smin);

}

void PHUrDilep::t_delta(Double_t tau, Double_t mres, Int_t, Double_t &tmax_del) {
    const Double_t smin_del = 0.0013;

    //     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    Double_t mx = smin_del;

    Double_t pf = sqrt((pow(mres,2)-pow((mass_nucleon+mx),2))*
		       (pow(mres,2)-pow((mass_nucleon-mx),2)))/2./mres;

    Double_t q0 = sqrt(mx*mx+pf*pf);


    Double_t f = (-3./2.)*(mres+mass_nucleon)/
	(mass_nucleon*(pow((mres+mass_nucleon),2)-pow(mx,2)));

    Double_t e = sqrt(4.*TMath::Pi()*alpha_em);

    Double_t mt = pow(e*f*g,2)*pow(mres,2)/9./mass_nucleon*
	(q0*q0*(5.*mres-3.*(q0+mass_nucleon))-
	 pow(mx,2)*(mres+mass_nucleon+q0));

    Double_t ml = pow(e*f*g,2)*pow(mres,2)/9./mass_nucleon*pow(mx,2)*4.*
	(mres-mass_nucleon-q0);

    Double_t lambda = pow(mx,4)+pow(mass_nucleon,4)+pow(mres,4)-2.*
	(pow(mx,2)*pow(mass_nucleon,2)+pow(mx,2)*pow(mres,2)+pow(mass_nucleon,2)*pow(mres,2));

    if(fabs(lambda) <= 1e-15) {
	lambda = 0;
    }

    if(mx <= q0) {
	Double_t gamma0 = sqrt(lambda)/(16.*TMath::Pi()*pow(mres,2))*mass_nucleon*(2.*mt+ml);
	Double_t dgamma = (2.*alpha_em)/(3.*TMath::Pi()*mx)*gamma0*(tau/gev);
	tmax_del = dgamma;
    }
}

void PHUrDilep::t_omega(Double_t tau,Int_t,
			Double_t& smin_omega,Double_t& smax_omega,Double_t& tmax_omega) {

    //     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    smin_omega = 0.0013;
    smax_omega = 0.6446;


    Double_t mx_omega = smin_omega;
    Double_t t1_omega = ((2.*alpha_em)/(3.*TMath::Pi()*mx_omega))*
	sqrt(1-(4.*pow(mass_electron,2))/
	     (pow(mx_omega,2))) * (1+(2.*pow(mass_electron,2))/(pow(mx_omega,2)));

    Double_t t2_omega = sqrt(pow((1+pow(mx_omega,2)/(pow(vacmass_omega,2)-pow(vacmass_pi0,2))),2)
			     -pow(((2.*vacmass_omega*mx_omega)/(pow(vacmass_omega,2)-pow(vacmass_pi0,2))),2));

    Double_t f1_omega = (pow(lambda_omega,2)*(pow(lambda_omega,2)+pow(gamma_omega,2)))/
	(pow((pow(lambda_omega,2)-pow(mx_omega,2)),2)+(pow(lambda_omega,2)*pow(gamma_omega,2)));

    Double_t dgamma_omega = t1_omega*pow(t2_omega,3)*f1_omega*gamma_photon*(tau/gev);

    tmax_omega = dgamma_omega;

}

void PHUrDilep::Output(Int_t ityp, Double_t weight, Double_t mres, Int_t, Double_t tau) {

    /*          .   ,
      write(*,'(I5,1X,13(E16.8E3,1X),I4,1X,2(E16.8E3,1X),I9,1X,F12.8,1X,2(E16.8,1x))')
      ityp,weight,mres,
      p0_el_lab,px_el_lab,py_el_lab,pz_el_lab,
      p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,
      tau,time_start,time_end,
      flag,
      dens_cre,dens_abs,
      noe,bev,acce,accp
    */
    //# ityp    weight           mres             p0_el_lab        px_el_lab        py_el_lab        pz_el_lab        p0_po_lab        px_po_lab        py_po_lab        pz_po_lab        tau              time_start       time_end     flag     dens_cre         dens_abs           noe bev              acce             accp
    //   17     9.377040e-06     1.940621e-03     3.300701e-01    -1.145014e-01    -1.798632e-01     2.519616e-01     1.025126e-01    -3.523874e-02    -5.656054e-02     7.789554e-02     3.753789e-01     3.410117e+00     3.792198e+00    2     3.075000e-01     3.360000e-01         1 3.941000e+00     1.000000e+00     1.000000e+00

    if(out) {
	fprintf(out, "%5i %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e ",
		ityp, weight, mres,
		em.E(), em.Px(), em.Py(), em.Pz(),
		ep.E(), ep.Px(), ep.Py(), ep.Pz(),
		tau, particleOut->t, particleIn->t
	       );
	fprintf(out, "%4i %15.8e %15.8e %9i %15.8e %15.8e %15.8e \n",
		collisionIn->n_in, collisionOut->baryon_density, collisionIn->baryon_density,
		evtheader->evtnum, evtheader->impact_par, 1., 1.
	       );
    }

    OutputROOT(ityp, weight);
}

void PHUrDilep::OutputROOT(Int_t, Double_t weight) {

    Int_t nRequest             = 1;
    if(outputLeptons) nRequest = 3;

    Int_t nfree = reader->GetNFreeSlots();
    if(nfree < nRequest) {

	Error("output()", "event reached max stacksize ! Will skip ....");
        return;
    }
    Int_t index = 0;
    PParticle *pMoth    = reader->CreateParticle(index);
    PHUrAddon *pMothAdd = reader->CreateAddon(index);

    Int_t pdg = particleIn->pdg;

    pMoth->Reset(pdg, particleIn->Px(), particleIn->Py(), particleIn->Pz(), particleIn->E(), weight);

    Int_t stable = 0;
    if(particleOut->t == particleIn ->t) stable =1;

    pMothAdd->t_cre    = particleOut ->t;
    pMothAdd->t_abs    = particleIn  ->t;
    pMothAdd->dens_cre = collisionOut->baryon_density;
    pMothAdd->dens_abs = collisionIn ->baryon_density;
    pMothAdd->fr_cre   = particleOut ->fr;
    pMothAdd->fr_abs   = particleIn  ->fr;
    pMothAdd->stable   = stable;
    pMothAdd->n_in        = collisionIn->n_in;
    pMothAdd->n_out       = collisionIn->n_out;
    pMothAdd->process_id  = collisionIn->process_id;
    pMothAdd->n_collision = collisionIn->n_collision;
    pMothAdd->instance    = particleIn->instance;
    pMothAdd->first       = particleIn->first;

    if(outputLeptons) {
        Int_t indexlep1 = 0;
        Int_t indexlep2 = 0;
	PParticle *pep    = reader->CreateParticle(indexlep1);
	PParticle *pem    = reader->CreateParticle(indexlep2);
	PHUrAddon *pepAdd = reader->CreateAddon(indexlep1);
	PHUrAddon *pemAdd = reader->CreateAddon(indexlep2);

	pep->Reset(-11, ep.Px(), ep.Py(), ep.Pz(), ep.E(), weight);       // ep
	pem->Reset( 11, em.Px(), em.Py(), em.Pz(), em.E(), weight);       // em
	pem->SetParentId(pdg);
	pem->SetParentIndex(index);
	pep->SetParentId(pdg);
	pep->SetParentIndex(index);

	*pepAdd = *pMothAdd;
	*pemAdd = *pMothAdd;
    }
}

Double_t PHUrDilep::rndfunc() {

    static long int init = 0;

    if(init == 0) {
	srand(1);
    }

    // 0 to RAND_MAX
    double val = ((double)rand())/RAND_MAX;

    init++;
    return val;
}
