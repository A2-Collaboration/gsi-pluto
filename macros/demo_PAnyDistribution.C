{

    //First, we create the general purpose distribution
    //model:
    PAnyDistribution *mydecay = 
        new PAnyDistribution("pt_rho",
			     "A function to define a new pt-distribution");
    mydecay->Add("q,     parent");
    mydecay->Add("p,     daughter");
    mydecay->Add("p,     daughter");
    mydecay->Add("rho0,  daughter");

    //This is the cache for the undistorted data
    //It is needed because the mandelstam t is a non-uniform
    //distribution. The size and binning of the cache must
    //be chosen such, that during runtime (or even better: during Preheating) the statistic
    //is sufficiently filled
    TH1F *cache  = new TH1F ("cache", "Rho0 pt cache", 400, 0, 1.0);

    //For the following it could be important to know that the daughters are in 
    //their rest frame (i.e. the parent).
    //But the parent itself is in lab frame.
    //The parent is indicated by "_parent"
    //Moreover, all particles of the decay (daughters and parent)
    //can be accessed via the usual '[parname]' identifier.
    //N.B.: "x,y,z,t" are reserved in TFormula, do NOT use it
    mydecay->AddEquation(cache, "_x = [rho0]->Pt();");
    
    //This is the final equation. The distribution (the probability function)
    //must be stored in "_f"
    mydecay->AddEquation("_f = 1; if ([rho0]->Pt() > 0.4) _f = 0.5;");
    //This function is of course completely arbitrary!

    //Remember, AnyDistribution is a rejection method. Therefore
    //it can happen, that parts of the phase space, where _f has a 
    //large probability, is not well populated by the generated events.
    //In this case, the event loop will run forever, as Pluto tries
    //to match the shape defined by _f.
    //The following factor is the maximum enhancement factor to avoid such
    //deadlocks.
    //N.B.: It directly scales with the computing time!!!
    mydecay->SetMaxEnhancementFactor(10);

    //Add this model to the Pluto data base:
    makeDistributionManager()->Add(mydecay);

    //Construct the reaction, as usual:
    PReaction my_reaction("_T1 = 2.2", "p", "p", "p p rho0 [pi+ pi-]");

    TH1F *histo1 = new TH1F("histo1", "rho0 pt", 40, 0, 0.7);
    TH1F *histo3 = new TH1F("histo3", "cos theta of rho0", 50, -1., 1.);

    my_reaction.Do(histo1, "_x = [rho0]->Pt();");
    my_reaction.Do(histo3, "_rho=[rho0]; _rho->Boost([p+p]); _x= cos(_rho->Theta())");

    my_reaction.Print();

    //Make a dummy loop to fill the AnyDistribution with some statistics:
    my_reaction.Preheating(500);
    //The event loop prints warnings - that's normal
    //The MaxFactor should not changed very much after the PreHeating finished
    
    my_reaction.Loop(10000);

    histo1->Draw();
}
