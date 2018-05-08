//TITLE Reading an UniGen file

{

    *(makeStaticData()->GetBatchValue("_system_particle_stacksize"))=1000;

    PUniGenInput *input = new PUniGenInput();
    input->Input("out-auau-6gev-10k-minbias.root");

    PReaction my_reaction("out-auau-6gev-10k-minbias_pluto");
 
    my_reaction.AddBulk(input);
    

    my_reaction.Do("foreach(*); id = [*]->ID(); echo $id");

    my_reaction.Print();
    cout << my_reaction.Loop() << " events converted" << endl;


}
