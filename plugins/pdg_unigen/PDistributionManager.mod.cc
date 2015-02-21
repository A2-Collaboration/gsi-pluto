//this line is added to the PDistributionManager.cc
PPDGPlugin *pdg = new PPDGPlugin("pdg","PDG extension");
AddPlugin(pdg);
Enable("pdg");
PluginInfo("PDG/UNIGEN classes available");

