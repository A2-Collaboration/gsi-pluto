////////////////////////////////////////////////////////
//  
//
//                    Author: I. Froehlich
//                    Written: 6.10.2007
//
////////////////////////////////////////////////////////


#include "PCommandList.h"
#include "PUtils.h"
#include <iostream>
#include "TROOT.h"


using namespace std;

PCommandList::PCommandList() {
    fNextP = NULL;
    fLevel = 0;
}

PCommandList::PCommandList(TString name) : TObjString(name) {
    fNextP = NULL;
    fLevel = 0;
};

PCommandList::PCommandList(TObject *obj, TString name, TString cmd) : TObjString(name) {
    fCmd.SetString(cmd);
    fNextP = NULL;
    fLevel = 0;
    if (obj) {
	Info("PCommandList", "Added object %s", obj->GetName());
	fTool.SetString(obj->GetName());
    }
};

PCommandList::PCommandList(TString name, TString cmd) : TObjString(name) {
    fCmd.SetString(cmd);
    fNextP = NULL;
    fLevel = 0;
};

PCommandList::~PCommandList() {
};

Bool_t PCommandList::AddCommand(const char *cmd, const char *basename, int level) {
    return AddCommand(NULL, cmd, basename, level);
};

Bool_t PCommandList::AddCommand(TObject *obj, const char *cmd, const char *basename, int level) {
    //if (basename && !fNextP) cout << "PCommandList::AddCommand " << basename << ":" << level << endl;
    //else if (!fNextP) cout << "PCommandList::AddCommand " << GetName() << ":" << level << endl;
    if (fNextP) {
	if (basename && level) return fNextP->AddCommand(obj,cmd,basename,level+1);
	else return fNextP->AddCommand(obj,cmd,GetName(),1);
    } else {
	if (basename && level) {
	    TString tmp(basename);
	    tmp += "_";
	    tmp += (level+1);
	    fNextP = new PCommandList(obj, tmp, cmd);    
	    fNext.SetString(PUtils::NewString(tmp));
	    return kTRUE;
	} else {
	    TString tmp(GetName());
	    tmp += "_1";
	    fNextP = new PCommandList(obj, tmp, cmd);  
	    fNext.SetString(PUtils::NewString(tmp));
	    return kTRUE;
	}
    }
    return kFALSE;
};

Bool_t PCommandList::GetCommand(char **cmd, int level, TObject **obj) {
    if (level == 0) {
	(*cmd) = PUtils::NewString(fCmd.GetString().Data());
	if (obj) {
	    if (fTool.GetString().Data() && strlen(fTool.GetString().Data())) {
		//Try to recover tool object
		*obj = gROOT->FindObject(fTool.GetString().Data());
		if (!(*obj)) {
		    Error("GetCommand", "Object %s could not be recovered", fTool.GetString().Data());
		    return kFALSE;
		}
		return kTRUE;
	    } else (*obj)=NULL;
	}
	return kTRUE;
    }

    if (fNextP) return fNextP->GetCommand(cmd, level-1, obj);
 
    if (fNext.GetString().Data() && strlen(fNext.GetString().Data())) {
	fNextP = (PCommandList*)gROOT->FindObject(fNext.GetString().Data());
	if (!fNextP) return kFALSE;
    }

    if (fNextP) return fNextP->GetCommand(cmd, level-1, obj);
    return kFALSE;
}

void PCommandList::Print(const Option_t *) const {
    
    if (fTool.GetString().Data() && strlen(fTool.GetString().Data())) {
	cout<< "<" <<  fTool.GetString().Data() << "> " << fCmd.GetString() << endl;
    } else {
	cout << "> " << fCmd.GetString() << endl;
    }
    //print following:
    if (fNext.GetString().Data() && strlen(fNext.GetString().Data())) {
	PCommandList *fNextP = (PCommandList*)gROOT->FindObject(fNext.GetString().Data());
	if (!fNextP) {
	    Error("Print", "Object %s not found", fNext.GetString().Data());
	    return;
	}
	fNextP->Print();
    }
};


Int_t PCommandList::Write(const char *name, Int_t option, Int_t bufsize) const {
    
    if (fNextP) 
	fNextP->Write(name, option, bufsize);
    return TObject::Write(name, option, bufsize);
}


ClassImp(PCommandList)
