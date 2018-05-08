//TITLE Using the pdg-code

{

    //Analyze example:
     PParticle *array[70];

    for (int i=0; i<70; i++) {
	array[i] = new PParticle(i);
    }

    makeDistributionManager()->Exec("pdg:init");

    for (int i=0; i<70; i++) {
	cout << array[i]->Name() << ", old id: " << array[i]->ID() << ", pdg code: " << array[i]->GetDBInt("pdg") << endl;
    }

}
