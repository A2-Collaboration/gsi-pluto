# Makefile for PLUTO 5.30
#


######################
#Directories and Names
######################
SRC = src
INSTALL = lib
DEPENDFILE := .depend

OPTIMIZE = -g
#OPTIMIZE = -O3


###########################
#Compiler flags and options
###########################
CXX           = g++

CXXFLAGS      = $(OPTIMIZE) -rdynamic -c -Wall -fno-exceptions -fPIC -DROOTVER=$(ROOTVER) \
		-I$(shell root-config --incdir) -I\. -I$(SRC) $(PLUGIN:%=-I%)

ARFLAGS       = r

LD            = g++
SOFLAGS       = $(OPTIMIZE) -rdynamic -shared -Wl,-soname,$@

###########################
#include basics and plugins
###########################

include Makefile.plugins

#check for additional user plugins:
PLUGIN := $(STD_PLUGIN)
ifeq ($(origin PLUTO_PLUGINS), environment)
	PLUGIN := $(PLUTO_PLUGINS)
endif
ifeq ($(origin PLUTO_USER_PLUGINS), environment)
	PLUGIN := $(PLUGIN) $(PLUTO_USER_PLUGINS)
endif


include Makefile.base
include $(PLUGIN:%=%/Makefile)


ROOTVERSION       := $(shell root-config --version)
ROOTLINEARVERSION := $(subst /,,$(subst .,,$(ROOTVERSION)))
ROOTVER           := $(word 1,$(subst /, ,$(subst ., ,$(ROOTVERSION))))
ROOTSUBVER        := $(word 2,$(subst /, ,$(subst ., ,$(ROOTVERSION))))
ROOTSUBSUBVER     := $(word 3,$(subst /, ,$(subst ., ,$(ROOTVERSION))))

PLUTOVERSION      = $(shell sed 's/^[^"]*"//' Version.h | sed  's/".*//')



ifdef PYTHIA6
  CXXFLAGS += -DUSE_PYTHIA6
endif

# Linker:


# Root libs
RLIBS = $(shell root-config --glibs) -lProof -lTreePlayer -lEG \
                                     -lm -ldl
ifdef PYTHIA6
  RLIBS += -lEGPythia6 -L$(PYTHIA6) -lPythia6
endif

SNAMES = $(sort $(BASE_CLASSES))

#ANAMES         =  PHADES PTrackH PTrack
ANAMES         =  

HDRS	      = $(addsuffix .h, $(addprefix $(SRC)/, $(SNAMES) ) ) 
#
PHDRS	      = $(addsuffix .h, $(PLUGIN_CLASSES) )
PSRC	      = $(addsuffix .cc, $(PLUGIN_CLASSES) )
SOBJS         = $(addsuffix .o, $(SNAMES) ) \
	        $(addsuffix .o, $(PLUGIN_CLASSES) )
TSOBJS        = $(addsuffix .o, $(addprefix $(INSTALL)/, $(SNAMES) ) ) 
PSOBJS        =	$(addsuffix .o, $(addprefix $(INSTALL)/, $(PLUGIN_CLASSES_NAMES) ) )

PRULES1	      = $(join $(addsuffix : , $(PSOBJS) ) , $(PSRC))
PRULES	      = $(join $(addsuffix + , $(PRULES1) ) , $(PHDRS))

#PLUGIN_COLLECTION = $(filter %Plugin, $(PLUGIN_CLASSES_NAMES))






# Output files:
LIBA          = libPluto.a
LIBSO         = libPluto.so

# Root:
ROOTCINT      = rootcint

# Default action:
all : $(INSTALL) $(LIBSO) $(LIBA)
	@echo Libraries done.

config :
	echo Pluto-Version: $(PLUTOVERSION)
	@echo Makefiles:
	@echo $(MAKEFILE_LIST)
	@echo \*\*\* Plugin classes to be compiled in:
	@echo $(PLUGIN_CLASSES_NAMES)
	@echo \*\*\* Plugin classes with fullpath:
	@echo $(PLUGIN_CLASSES)
	@echo \*\*\* Plugin target objects:
	@echo $(PSOBJS)
	@echo \*\*\* Plugin header files:
	@echo $(PHDRS)
	@echo \*\*\* Plugin source files:
	@echo $(PSRC)
	@echo \*\*\* Plugin rules:
	@echo $(PRULES)
	@echo $(MAKEFILE_LIST)

depend :
	rm -f $(DEPENDFILE)
	make $(DEPENDFILE)


$(DEPENDFILE) : $(PLUGIN:%=%/Makefile) Plugins.h
	@echo $(PRULES) | sed 's/+/\ /g'  | sed 's/\.h/\.h\n/g' > $(DEPENDFILE)

Plugins.cc : $(PLUGIN:%=%/Makefile)
	@echo $(addprefix \#include \", $(addsuffix \"+ , $(PLUGIN_COLLECTION)))| sed 's/+/\n/g'  > Plugins.cc

Plugins.h : $(PLUGIN:%=%/Makefile)
	@echo > Plugins.h
	@echo $(addprefix \#include \", $(addsuffix \"+, $(PHDRS))) | sed 's/+/\n/g'  > Plugins.h
#	@echo $(addprefix \#include \", $(addsuffix \"+, $(PHDRS))) > Plugins.h


$(INSTALL)/PDistributionManager.o : $(PLUGIN:%=%/Makefile) Plugins.h Plugins.cc

# Static library:
$(LIBA) : $(TSOBJS) $(INSTALL)/PlutoCint.o
	@rm -f $@
	@$(AR) $(ARFLAGS) $@ $^
# Shared library:
$(LIBSO) : $(TSOBJS) $(PSOBJS) $(INSTALL)/PlutoCint.o
	@echo Linking libPluto.so
	@$(LD) $(SOFLAGS) -o $@ $^

$(INSTALL) :
	@mkdir $(INSTALL)

# Rules for PlutoCint:
$(INSTALL)/PlutoCint.o : $(INSTALL)/PlutoCint.cc
	@echo Compiling dictionary
	@$(CXX) $(CXXFLAGS) $< -o $@
$(INSTALL)/PlutoCint.cc : $(HDRS) $(PHDRS) PlutoLinkdef.h 
	@echo Creating dictionary
	@echo $(addprefix "#pragma link C++ class ", $(addsuffix ";", $(PLUGIN_CLASSES_NAMES))) | sed 's/;/;\n/g' > PluginLinkdef.h
	@$(ROOTCINT) -f $@ -c  -I\. -I$(SRC) $(PLUGIN:%=-I%) $^


# Static pattern rule for object file dependency on sources:
$(TSOBJS): $(INSTALL)/%.o : $(SRC)/%.cc $(SRC)/%.h
	 @echo Compiling Base Class $*
	 @echo char *date_string = "(char*)"\"$$(date +"%e %B %Y")\"\; > Compiled.h
	 @$(CXX) $(CXXFLAGS) $< -o $@

#$(PSOBJS): $(INSTALL)/%.o : $(shell grep $* $(PLUGIN:%=%/config))
$(PSOBJS): $(INSTALL)/%.o : $(DEPENDFILE)
	 @echo Compiling Plugin     $*
	 @$(CXX) $(CXXFLAGS) $(filter %$*.cc, $(PSRC)) -o $@

Pluto:	$(LIBSO) pluto.cc 
	@$(CXX) $(CXXFLAGS) pluto.cc -o pluto.o
	@$(LD) -g pluto.o -L. -L$(ROOTSYS)/lib $(RLIBS)  -lPluto -o Pluto

docs:   macros/makeClassDocs.C Pluto.so
	@echo Generating the documentation...
	@root -q -l makeClassDocs.C

tar:    
	tar cvf pluto.tar src/*.cc src/*.h Makefile PlutoLinkdef.h README AUTHORS REFERENCES macros/*.C; gzip pluto.tar

# cleanup:
clean:
	rm -f $(INSTALL)/*.o
	rm -f $(INSTALL)/*Cint.*
	rm -f $(INSTALL)/$(LIBSO) $(INSTALL)/$(LIBA)
	rm -f Plugins.h Plugins.cc

pluginclean:
	rm -f $(PSOBJS)

fairroot:
	rm -fr $(SIMPATH)/generators/Pluto.$(PLUTOVERSION)*
	rm -fr $(SIMPATH)/generators/pluto
	mkdir -p $(SIMPATH)/generators/pluto
	cp -r src $(SIMPATH)/generators/pluto
	cp -r plugins/ $(SIMPATH)/generators/pluto
	cp $(PHDRS) $(SIMPATH)/generators/pluto/src
	cp -r macros $(SIMPATH)/generators/pluto
	cp plugins/fairroot/PFairGenerator.h $(SIMPATH)/generators/pluto/src
	cp  Makefile Makefile.base Makefile.plugins Makefile.fairsoft newline.awk  PlutoLinkdef.h README AUTHORS REFERENCES Version.h $(SIMPATH)/generators/pluto/
	cd $(SIMPATH)/generators/;  tar cvf pluto.tar pluto/* --exclude-vcs; gzip pluto.tar; mv pluto.tar.gz Pluto.$(PLUTOVERSION).tar.gz
	rm -fr $(SIMPATH)/generators/pluto
	rm -fr $(SIMPATH)/generators/libPluto.*


sinclude $(DEPENDFILE)




