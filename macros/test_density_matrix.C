{

    //Creates the matrix object:
    PDensityMatrix *matrix = new PDensityMatrix();
    
    //Open the HSD-file:
    matrix->ReadDensityMatrix("fort_0.5.925", 3, kTRUE, 0.49, 0.51);
    //Options:
    // * Filename
    // * Number of dimension-columns (Important!)
    // * Correct integral for the bin width
    // * Lower bound for section selection
    // * Upper bound for section selection
    
    //Create a seed particle:
    PParticle *seed = new PParticle("dilepton", 0., 0., 0.);
    //matrix->AddParticle(seed); //Do not use that when in already taken in PReaction as seed!!!

    //Use matrix 12:
    matrix->SetMatrix(12);
    //(counting from 0, the columns for the dimensions must be substracted)
    
    //Convert _x, _y and _z into physical variables:
    matrix->Do("mt  = sqrt(_z*_z + _x*_x)");
    matrix->Do("p3  = mt * sinh(_y)");
    matrix->Do("phi = sampleFlat() * TMath::Pi() * 2.0");
    matrix->Do("p1  = _z * sin(phi)");
    matrix->Do("p2  = _z * cos(phi)");
    matrix->Do("[dilepton]->SetPx(p1)");
    matrix->Do("[dilepton]->SetPy(p2)");
    matrix->Do("[dilepton]->SetPz(p3)");
    matrix->Do("[dilepton]->SetM(_x)");
    
    //Start the reaction:
    PParticle *ep = new PParticle("e+");
    PParticle *em = new PParticle("e-");
    PParticle *p[]={seed,em,ep};
    PChannel dilepton_decay(p,2);
    dilepton_decay.AddPrologueBulk(matrix);

    PChannel *c[]={&dilepton_decay};
    PReaction my_reaction(c, "hsd", 1);
    
    //Just for debugging...
    TH1F * histo1 = new TH1F ("histo1","ee mass",100,0.01,1.5);
    my_reaction.Do(histo1, "_x = ([e+]+[e-])->M()");
    TH2F * histo2 = new TH2F ("histo2","Y vs. Pt", 100, -5.0, 5.0, 300, 0.,1.);
    my_reaction.Do(histo2, "_x = ([e+]+[e-])->Rapidity(); _y = ([e+]+[e-])->Pt()");

    my_reaction.Loop(1000000);

    histo1->Draw();

}
