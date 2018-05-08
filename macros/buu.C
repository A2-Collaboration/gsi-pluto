//this is a short demo macro the use 
//Pluto as a simple filter, without generating particles

{

    //Uncomment this for debug info:
    //makeDistributionManager()->Startup("_filter_debug=1"); 

    //Compile the input/output classes
    gSystem->CompileMacro( "PGiBuuInput.C");
    gSystem->CompileMacro( "PGiBuuOutput.C");

    //Unpack and attach the HADES filter
    //makeDistributionManager()->Unpack("pluto_ee_filter_may07.root"); //for d+p @ 1.25AGeV
    makeDistributionManager()->Unpack("pluto_ee_filter_apr06.root"); //for p+p @ 1.25GeV

    //Create the dummy "reaction"
    PReaction dummy;
    //Add the IO objects
    dummy.AddBulk(new PGiBuuInput("buu_dil.dat"));
    PGiBuuOutput *out = new PGiBuuOutput("buu_dil2.dat");
    dummy.AddBulk(out);
    dummy.Print();
    
    //Loop, and see what happens:
    cout << dummy.Loop(100) << " events read from file" << endl;
    out->CloseFile();

}
