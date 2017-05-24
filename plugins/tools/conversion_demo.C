{
    

    PReaction my_reaction("_T1 = 2.2","p","p","p p eta [g g]", "eta_conv",0,0,1,1);
    
    Float_t Zmean = 79;
    Float_t pconv = 20.;   // 0 => do lookup, >0 => use fixed value
    PGammaConversion *conv = new PGammaConversion(Zmean, pconv,2);
    conv->SetVertex(1,2,3);

    //if(conv->readHist("HGeantSim_white_gamma_vertex.root")<0){
    //	cout<<"COULD NOT READ VERTEXHIST!"<<endl;	
    //} else {
    //	cout<<"VERTEXHIST INITIALIZED!"<<endl;
    //}

    conv->SetPriority(51);   // needed to get conversion done last
    
    my_reaction.AddBulk(conv);
    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(5);
}
