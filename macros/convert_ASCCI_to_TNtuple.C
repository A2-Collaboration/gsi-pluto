
{

    TFile *f = new TFile("ntuple.root", "RECREATE");    //Create a file for the NTuple
    TNtuple *ntuple = new TNtuple("ntuple", 
                                  "data from ASCII file", 
                                  "p1px:p1py:p1pz:p1e:p1id");   

    PReaction my_reaction;
    
    PProjector *input = new PProjector();
    input->AddInputASCII("ascii.txt", 
			 "readline{@delme}; ");
    input->Do("readline{@p1px @p1py @p1pz @p1e @p1id}; ");
    input->Do("readline{@p2px @p2py @p2pz @p2e @p2id}; ");
    input->Do("readline{@p3px @p3py @p3pz @p3e @p3id}; ");
    input->Do("readline{@p4px @p4py @p4pz @p4e @p4id}; ");
    my_reaction.AddPrologueBulk(input); //The "prolog" is done before the decay
    

    //my_reaction.Do("echo $p1id, $p2id, $p3id, $p4id"); //only for debugging
    my_reaction.Output(ntuple);     


    my_reaction.Print();

    cout << my_reaction.Loop() << " events converted" << endl;
    
    
    f->Write();

}
