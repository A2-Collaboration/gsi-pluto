//this line is added to the PDistributionManager.cc
PEtaDecaysPlugin *eta_decays = new PEtaDecaysPlugin("eta_decays", "Plugin for (rare) eta decays");
AddPlugin(eta_decays);
Enable("eta_decays"); //Auto-enabled
PluginInfo("Rare eta decays are enabled");

