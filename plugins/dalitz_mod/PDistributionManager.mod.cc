//this line is added to the PDistributionManager.cc
PDalitzModPlugin * plugin = new PDalitzModPlugin("dalitz_mod","Plugin to modify Dalitz decays");
AddPlugin(plugin);
Enable("dalitz_mod"); //Auto-enabled

PluginInfo("Plugin for Dalitz decay (generator & new D Dalitz) available");
PluginInfo("Dalitz decays of N* activated");

