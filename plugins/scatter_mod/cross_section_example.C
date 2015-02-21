//TITLE Example for a 2-body scattering with total and differential cross sections

{

    //Add a model for the reaction "g+p -> p+eta":
    PScatterCrossSection * model = new PScatterCrossSection("mymodel","My cross section");
    model->Add("g,grandparent,beam");
    model->Add("p,grandparent,target");
    model->Add("q,parent");
    model->Add("p,daughter");
    model->Add("eta,daughter,primary");

    //Define the range of the c.m. sampling:
    model->SetRange(1.4,2.0); 

    //Now start to fill the helper histogram
    //The axis *can* be like in the PScatterCrossSection model (see below)
    //but this is not to be required. You can use the batch equations 
    //to make a conversion
    TH1F *distribution = new TH1F("distribution", "Angular distribution", 10, -1 , 1 );
    distribution->SetBinContent(1,20.);
    distribution->SetBinContent(2,16.);
    distribution->SetBinContent(3,11.);
    distribution->SetBinContent(4,8.);
    distribution->SetBinContent(5,5.);
    distribution->SetBinContent(6,4.);
    distribution->SetBinContent(7,3.);
    distribution->SetBinContent(8,2.5);
    distribution->SetBinContent(9,2.);
    distribution->SetBinContent(10,1.);
    
    //Now add the histogram to the model class, and define an equation
    //Input:  _x is cos(theta), _y is the c.m. energy
    //Output: _f: cross section
    model->AddHistogram(distribution,"value = Eval(_x); _f = _y * value");
    //This equation only a placeholder. If the histogram has a
    //different axis, it is possible to convert the axis, e.g.: 
    //
    //my_x = f(_x); value = Eval(my_x); 
    //
    //with f() as any user-defined function which is possible (see
    //script manual for details). A TH2 can also be used, in this case
    //one should the 2-dimensional version of Eval(x,y)
    //
    //The function can also be staged:
    //model->AddHistogram(distribution,"value = Eval(_x);");
    //model->AddEquation(" _f = _y * value");
    //produces the same result as above
    //
    //The function _f can be folded with the beam profile. But one has
    //to keep in mind that the beam profile must be converted into a
    //function based on the c.m. energy (e.g. by calculating the T_lab
    //based on q)
    //
    //The following example demonstate this. Let's say that a
    //histogram "profile" contains the beam profile and the x-axis is
    //T_lab. Therefore, one has to do a small calculation like:
    //
    //model->AddEquation("t_lab = (_y*_y - g.mass*g.mass - p.mass*p.mass)/(2*p.mass*p.mass) - g.mass;");
    //model->AddHistogram(profile,"_f = _f * Eval(t_lab);");
    //(the gamma mass "g.mass" is only added for completeness)
    //
    //if one would like to use an analytic function instead of a
    //histogram is also possible. 


    //Add the model to the manager:
    makeDistributionManager()->Add(model);
    
    
    PReaction my_reaction("_P1 = 2.2","g","p","p eta [dilepton [e+ e-] g]");
    TH2F * histo2 = new TH2F ("histo2","c.m. vs cos_theta",100,1.0,3.0,20,-1,1);
    my_reaction.Do(histo2,"_x = [g+p]->M(); myeta = [eta]; myeta->Boost([g+p]); _y = myeta->CosTheta()");
    my_reaction.Print();


    my_reaction.Loop(100000);


}
