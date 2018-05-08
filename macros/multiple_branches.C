//TITLE Using more than one ROOT TBranch to store undistorted sim particles

{
    
    PUtils::SetSeed(123); 
    PReaction my_reaction("_T1 = 2.2", "p", "p", "p p w [e+ e-]", "w_multiple_branches");
    //for syntax see parse_eta.C

    //save particles on a new branch, called "particles_stage1":
    my_reaction.Do("foreach(*); [*]->Push(Branch(particles_stage1))");

    //apply a filter (keeping the w out):
    my_reaction.Do("loop: theta_ep = ([*]->Theta() * 180.)/TMath::Pi()");
    my_reaction.Do("theta_em = ([*]->Theta() * 180.)/TMath::Pi()");
    my_reaction.Do("if ((theta_ep<18 || theta_ep>85 || theta_em<18 || theta_em>85) && ![*]->HasID(w.pid)); [*]->SetInActive(); ");
    my_reaction.Do("formore(*); goto loop;");

    //reconstruct w, which "survived" the filter
    my_reaction.Do("[w]->SetM(([e+]+[e-])->M()); else; [w]->SetInActive();");

    //save particles on a new branch, called "particles_stage2", before doing the smearing:
    my_reaction.Do("foreach(*); [*]->Push(Branch(particles_stage2))");

    //10% momentum resolution:
    my_reaction.Do("foreach(*); mom = [*]->P(); newmom = sampleGaus(mom,0.1*mom);[*]->SetMom(newmom);");

    //reconstruct w, after smearing:
    my_reaction.Do("[w]->SetM( ([e+]+[e-])->M() );");

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(1000);

#if 0
    //root -l w_multiple_branches.root
    data->Draw("particles_stage1.M() >> htemp1","particles_stage1.ID() == 52");
    htemp1->SetLineColor(1);
    data->Draw("particles_stage2.M() >> htemp2","particles_stage2.ID() == 52","same");
    htemp2->SetLineColor(2);
    data->Draw("Particles.M() >> htemp3","Particles.ID() == 52","same");
#endif

}
