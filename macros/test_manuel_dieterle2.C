{

    makeDistributionManager()->Exec("nucleus_fermi:gamma");

    //PReaction *Reac = new PReaction (1.5,"g","12C","(g n) n pi0 [g g] (11C)","delme");
    PReaction *Reac = new PReaction (1.5,"g","12C","(g p) p pi0 [g g] (11B)","delme");

    Reac->Print();        // Write to .root file      
    Reac->loop(10);  // Number of events

}
