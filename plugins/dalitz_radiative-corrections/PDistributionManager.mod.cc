//this line is added to the PDistributionManager.cc
AddPlugin(new PDalitzCorrectionsPlugin("dalitz_corrections","Plugin for radiative corrections in Dalitz decays"));
//Enable("dalitz_corrections"); //Auto-enabled

PluginInfo("Plugin to apply radiative corrections to Dalitz decays of pi^0, eta, and eta'");
