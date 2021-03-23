#!/usr/bin/env python
# Analyzer for WWG Analysis based on nanoAOD tools

import os, sys
import math
import ROOT
from math import sin, cos, sqrt
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import countHistogramsProducer

class WWG_Producer(Module):
    def __init__(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.out = wrappedOutputTree

        self.out.branch("event",  "i")
        self.out.branch("run",  "i")
        self.out.branch("lumi",  "i")
	self.out.branch("pass_selection",  "B");
	self.out.branch("channel",  "I");

        self.out.branch("lep1_pid",  "I")
        self.out.branch("lep2_pid",  "I")
        self.out.branch("lep1pt",  "F")
        self.out.branch("lep2pt",  "F")
        self.out.branch("lep1eta",  "F")
        self.out.branch("lep2eta",  "F")
        self.out.branch("lep1phi",  "F")
        self.out.branch("lep2phi",  "F")
        self.out.branch("lepton1_isprompt", "I")
        self.out.branch("lepton2_isprompt", "I")
        self.out.branch("photonet",  "F")
        self.out.branch("photoneta",  "F")
        self.out.branch("photonphi",  "F")
        self.out.branch("photon_isprompt", "I")
        self.out.branch("photon_gen_matching", "I")
        self.out.branch("mll",  "F")
        self.out.branch("mllg",  "F")
        self.out.branch("ptll",  "F")
        self.out.branch("met",  "F")
        self.out.branch("metup",  "F")
        self.out.branch("puppimet","F")
        self.out.branch("puppimetphi","F")
        self.out.branch("rawmet","F")
        self.out.branch("rawmetphi","F")
        self.out.branch("metphi","F")
        self.out.branch("gen_weight","F")
        self.out.branch("npu",  "I");
        self.out.branch("ntruepu",  "F");
        self.out.branch("n_pos", "I")
        self.out.branch("n_minus", "I")
        self.out.branch("n_num", "I")
        self.out.branch("MET_pass","I")
        self.out.branch("npvs","I")
        self.out.branch("njets","I")
        self.out.branch("njets50","I")
        self.out.branch("njets40","I")
        self.out.branch("njets30","I")
        self.out.branch("njets20","I")
        self.out.branch("njets15","I")
        self.out.branch("pass_Ele23Ele12","I")
        self.out.branch("pass_Ele23Ele12_DZ","I")
        self.out.branch("pass_Ele35","I")
        self.out.branch("pass_Ele32","I")
        self.out.branch("pass_Mu17Mu8_3p8","I")
        self.out.branch("pass_Mu17Mu8","I")
        self.out.branch("pass_Mu27","I")
        self.out.branch("pass_Mu12Ele23","I")
        self.out.branch("pass_Mu23Ele12","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.out.fillBranch("event",event.event)
        self.out.fillBranch("lumi",event.luminosityBlock)
        self.out.fillBranch("run",event.run)
#        print event.event,event.luminosityBlock,event.run
        if hasattr(event,'Generator_weight'):
            if event.Generator_weight > 0 :
                n_pos=1
                n_minus=0
            else:
                n_minus=1
                n_pos=0
            self.out.fillBranch("gen_weight",event.Generator_weight)
            self.out.fillBranch("n_pos",n_pos)
            self.out.fillBranch("n_minus",n_minus)
        else:    
            self.out.fillBranch("gen_weight",0)
            self.out.fillBranch("n_pos",0)
            self.out.fillBranch("n_minus",0)

        pass_Ele23Ele12 = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
        pass_Ele23Ele12_DZ = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
        pass_Ele35 = event.HLT_Ele35_WPTight_Gsf
        pass_Ele32 = event.HLT_Ele32_WPTight_Gsf_L1DoubleEG

        pass_Mu17Mu8_3p8 = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
        pass_Mu17Mu8 = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8
        pass_Mu27 = event.HLT_IsoMu27

        pass_Mu12Ele23 = event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
        pass_Mu23Ele12 = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL

        if not (pass_Ele23Ele12 or pass_Ele23Ele12_DZ or pass_Ele35 or pass_Ele32 or pass_Mu17Mu8_3p8 or pass_Mu17Mu8 or pass_Mu27 or pass_Mu12Ele23 or pass_Mu23Ele12):
           self.out.fillBranch("pass_selection",0)
           return True
        self.out.fillBranch("pass_Ele23Ele12",pass_Ele23Ele12)
        self.out.fillBranch("pass_Ele23Ele12_DZ",pass_Ele23Ele12_DZ)
        self.out.fillBranch("pass_Ele35",pass_Ele35)
        self.out.fillBranch("pass_Ele32",pass_Ele32)
        self.out.fillBranch("pass_Mu17Mu8_3p8",pass_Mu17Mu8_3p8)
        self.out.fillBranch("pass_Mu17Mu8",pass_Mu17Mu8)
        self.out.fillBranch("pass_Mu27",pass_Mu27)
        self.out.fillBranch("pass_Mu12Ele23",pass_Mu12Ele23)
        self.out.fillBranch("pass_Mu23Ele12",pass_Mu23Ele12)

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")
        jets = Collection(event, "Jet")
	if hasattr(event, 'nGenPart'):
           genparts = Collection(event, "GenPart")

        jet_select = [] 
        dileptonp4 = ROOT.TLorentzVector()
        photons_select = []
        electrons_select = []
        muons_select = [] 
        jets_select = []
        leptons_select=[]

        #selection on muons
        muon_pass =0
        for i in range(0,len(muons)):
            if muons[i].pt < 20:
                continue
            if abs(muons[i].eta) > 2.5:
                continue
            if muons[i].pfRelIso04_all > 0.25:
                continue   
            if muons[i].mediumId == True:
                muons_select.append(i)
                muon_pass += 1
                leptons_select.append(i)

        # selection on electrons
        electron_pass=0
        for i in range(0,len(electrons)):
            if electrons[i].pt < 20:
                continue
            if abs(electrons[i].eta + electrons[i].deltaEtaSC) > 2.5:
                continue
            if (abs(electrons[i].eta + electrons[i].deltaEtaSC) < 1.479 and abs(electrons[i].dz) < 0.1 and abs(electrons[i].dxy) < 0.05) or (abs(electrons[i].eta + electrons[i].deltaEtaSC) > 1.479 and abs(electrons[i].dz) < 0.2 and abs(electrons[i].dxy) < 0.1):
                if electrons[i].cutBased >= 3:
                    electrons_select.append(i)
                    electron_pass += 1
                    leptons_select.append(i)

#        print 'the number of leptons: ',len(electrons_select)+len(muons_select)
        if len(electrons_select)+len(muons_select) != 2:      #reject event if there are not exactly two leptons
	   self.out.fillBranch("pass_selection",0)
	   return True

        # selection on photons
	photon_pass=0
        for i in range(0,len(photons)):
            if photons[i].pt < 20:
                continue
            if abs(photons[i].eta) > 2.5:
                continue
            if not (photons[i].isScEtaEE or photons[i].isScEtaEB):
                continue
            if photons[i].pixelSeed:
                continue
            pass_lepton_dr_cut = True
            for j in range(0,len(muons_select)):
                if deltaR(muons[muons_select[j]].eta,muons[muons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False
            for j in range(0,len(electrons_select)):
                if deltaR(electrons[electrons_select[j]].eta,electrons[electrons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False
            if not pass_lepton_dr_cut:
                continue
            if photons[i].cutBased >=2:
                photons_select.append(i)
                photon_pass += 1

#        print 'the number of photons: ',len(photons_select)
        if  len(photons_select)<1:
            self.out.fillBranch("pass_selection",0)
            return True

        pass_lepton_dr_cut = True
        njets = 0
        njets50 = 0
        njets40 = 0
        njets30 = 0
        njets20 = 0
        njets15 = 0
        for i in range(0,len(jets)):
            if jets[i].btagDeepB > 0.4184 and i<=6 :  # DeepCSVM
               self.out.fillBranch("pass_selection",0)
               return True
            if abs(jets[i].eta) > 4.7:
               continue
            if jets[i].pt<30:
               continue
	    if deltaR(jets[i].eta,jets[i].phi,photons[photons_select[0]].eta,photons[photons_select[0]].phi) < 0.5:
	       continue;
            for j in range(0,len(electrons_select)):
                if deltaR(jets[i].eta,jets[i].phi,electrons[electrons_select[j]].eta,electrons[electrons_select[j]].phi) < 0.5:
                   pass_lepton_dr_cut = False
            for j in range(0,len(muons_select)):
                if deltaR(jets[i].eta,jets[i].phi,muons[muons_select[j]].eta,muons[muons_select[j]].phi) < 0.5:
                   pass_lepton_dr_cut = False

            if  not pass_lepton_dr_cut == True:
	        continue
            if jets[i].jetId >> 1 & 1:
               jets_select.append(i)
               njets += 1
            if jets[i].pt > 50:
                njets50+=1
            if jets[i].pt > 40:
                njets40+=1
            if jets[i].pt > 30:
                njets30+=1
            if jets[i].pt > 20:
                njets20+=1
            if jets[i].pt > 15:
                njets15+=1
#        print len(jets),("njets",njets)
#        if njets >=2 :
#            self.out.fillBranch("pass_selection",0)
#            return True

        isprompt_mask = (1 << 0) #isPrompt
        isdirectprompttaudecayproduct_mask = (1 << 5) #isDirectPromptTauDecayProduct
        isdirecttaudecayproduct_mask = (1 << 4) #isDirectTauDecayProduct
        isprompttaudecayproduct = (1 << 3) #isPromptTauDecayProduct
        isfromhardprocess_mask = (1 << 8) #isPrompt

        channel = 0 
        # emu:     1
        # ee:      2
        # mumu:    3

        # emu
        mll = -10
        ptll = -10
        lepton1_isprompt=-10
        lepton2_isprompt=-10
        if len(muons_select)==1 and len(electrons_select)==1:  # emu channel 
            if deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) < 0.5:
	       self.out.fillBranch("pass_selection",0)
               return True 
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,muons[muons_select[0]].eta,muons[muons_select[0]].phi) < 0.5:
                self.out.fillBranch("pass_selection",0)
                return True
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) < 0.5:
                self.out.fillBranch("pass_selection",0)
                return True
            if muons[muons_select[0]].charge * (electrons[electrons_select[0]].charge) >= 0:
                self.out.fillBranch("pass_selection",0)
                return True
#           print 'test emu channel',len(genparts)
            if hasattr(event, 'nGenPart'):
                print 'calculate the lepton flag in channel emu'
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1
                       break 
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1 
                       break 
            channel = 1
            self.out.fillBranch("channel",channel)
	    self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
	    self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lep1_pid",13)
            self.out.fillBranch("lep2_pid",11)
            self.out.fillBranch("lep1pt",muons[muons_select[0]].pt)
            self.out.fillBranch("lep1eta",muons[muons_select[0]].eta)
            self.out.fillBranch("lep1phi",muons[muons_select[0]].phi)
            self.out.fillBranch("lep2pt",electrons[electrons_select[0]].pt)
            self.out.fillBranch("lep2eta",electrons[electrons_select[0]].eta)
            self.out.fillBranch("lep2phi",electrons[electrons_select[0]].phi)
            self.out.fillBranch("mll",(muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()).M())
            self.out.fillBranch("ptll",(muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()).Pt())
            self.out.fillBranch("photonet",photons[photons_select[0]].pt)
            self.out.fillBranch("photoneta",photons[photons_select[0]].eta)
            self.out.fillBranch("photonphi",photons[photons_select[0]].phi)
            self.out.fillBranch("mllg",(muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()+photons[photons_select[0]].p4()).M())
        # ee
        elif len(muons_select)==0 and len(electrons_select)==2:
            if deltaR(electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi,electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi)<0.5:
	       self.out.fillBranch("pass_selection",0)
               return True 
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) < 0.5:
                self.out.fillBranch("pass_selection",0)
                return True
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi) < 0.5:
                self.out.fillBranch("pass_selection",0)
                return True
            if electrons[electrons_select[0]].charge * electrons[electrons_select[1]].charge >=0:
	        self.out.fillBranch("pass_selection",0)
                return True 
#           print 'test',len(genparts)
	    if hasattr(event, 'nGenPart'):
                print 'calculate the lepton flag in channel ee'
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1 
                       break 
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1 
                       break 
            channel = 2
            self.out.fillBranch("channel",channel)
	    self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
	    self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lep1_pid",11)
            self.out.fillBranch("lep2_pid",11)
            self.out.fillBranch("lep1pt",electrons[electrons_select[0]].pt)
            self.out.fillBranch("lep1eta",electrons[electrons_select[0]].eta)
            self.out.fillBranch("lep1phi",electrons[electrons_select[0]].phi)
            self.out.fillBranch("lep2pt",electrons[electrons_select[1]].pt)
            self.out.fillBranch("lep2eta",electrons[electrons_select[1]].eta)
            self.out.fillBranch("lep2phi",electrons[electrons_select[1]].phi)
            self.out.fillBranch("mll",(electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()).M())
            self.out.fillBranch("ptll",(electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()).Pt())
            self.out.fillBranch("photonet",photons[photons_select[0]].pt)
            self.out.fillBranch("photoneta",photons[photons_select[0]].eta)
            self.out.fillBranch("photonphi",photons[photons_select[0]].phi)
            self.out.fillBranch("mllg",(electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()+photons[photons_select[0]].p4()).M())


        # mumu 
        elif len(electrons_select)==0 and len(muons_select)==2:
            if deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,muons[muons_select[1]].eta,muons[muons_select[1]].phi)<0.5:
	       self.out.fillBranch("pass_selection",0)
               return True 
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,muons[muons_select[0]].eta,muons[muons_select[0]].phi) < 0.5:
                self.out.fillBranch("pass_selection",0)
                return True
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,muons[muons_select[1]].eta,muons[muons_select[1]].phi) < 0.5:
                self.out.fillBranch("pass_selection",0)
                return True
            if muons[muons_select[0]].charge * muons[muons_select[1]].charge >= 0:
	       self.out.fillBranch("pass_selection",0)
               return True 
	    if hasattr(event, 'nGenPart'):
                print 'calculate the lepton flag in channel mumu'
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1
                       break 
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muons_select[1]].eta,muons[muons_select[1]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1
                       break 
            channel = 3
            self.out.fillBranch("channel",channel)
	    self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
	    self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lep1_pid",13)
            self.out.fillBranch("lep2_pid",13)
            self.out.fillBranch("lep1pt",muons[muons_select[0]].pt)
            self.out.fillBranch("lep1eta",muons[muons_select[0]].eta)
            self.out.fillBranch("lep1phi",muons[muons_select[0]].phi)
            self.out.fillBranch("lep2pt",muons[muons_select[1]].pt)
            self.out.fillBranch("lep2eta",muons[muons_select[1]].eta)
            self.out.fillBranch("lep2phi",muons[muons_select[1]].phi)
            self.out.fillBranch("mll",(muons[muons_select[0]].p4() + muons[muons_select[1]].p4()).M())
            self.out.fillBranch("ptll",(muons[muons_select[0]].p4() + muons[muons_select[1]].p4()).Pt())
            self.out.fillBranch("photonet",photons[photons_select[0]].pt)
            self.out.fillBranch("photoneta",photons[photons_select[0]].eta)
            self.out.fillBranch("photonphi",photons[photons_select[0]].phi)
            self.out.fillBranch("mllg",(muons[muons_select[0]].p4() + muons[muons_select[1]].p4()+photons[photons_select[0]].p4()).M())

        else:
            self.out.fillBranch("pass_selection",0)
            return True
        photon_gen_matching=-10
        photon_isprompt =-10
        if hasattr(photons[photons_select[0]],'genPartIdx') :
            print 'calculate the photon flag'
            if photons[photons_select[0]].genPartIdx >= 0 and genparts[photons[photons_select[0]].genPartIdx].pdgId  == 22: 
                if ((genparts[photons[photons_select[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photons_select[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)) and (genparts[photons[photons_select[0]].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask):
                    photon_gen_matching = 6
                elif ((genparts[photons[photons_select[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photons_select[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):       
                    if (genparts[photons[photons_select[0]].genPartIdx].genPartIdxMother >= 0 and (abs(genparts[genparts[photons[photons_select[0]].genPartIdx].genPartIdxMother].pdgId) == 11 or abs(genparts[genparts[photons[photons_select[0]].genPartIdx].genPartIdxMother].pdgId) == 13 or abs(genparts[genparts[photons[photons_select[0]].genPartIdx].genPartIdxMother].pdgId) == 15)):
                        photon_gen_matching = 4
                    else:    
                        photon_gen_matching = 5
                else:
                    photon_gen_matching = 3
            elif photons[photons_select[0]].genPartIdx >= 0 and abs(genparts[photons[photons_select[0]].genPartIdx].pdgId) == 11:     
                if ((genparts[photons[photons_select[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photons_select[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):  
                    photon_gen_matching = 1
                else:
                    photon_gen_matching = 2
                    
            else:
                assert(photons[photons_select[0]].genPartFlav == 0)
                photon_gen_matching = 0
        if hasattr(event, 'nGenPart') :
            for j, genpart in enumerate(genparts):
	        if photons[photons_select[0]].genPartIdx >=0 and genpart.pt > 5 and abs(genpart.pdgId) == 22 and ((genparts[photons[photons_select[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photons_select[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) or (genparts[photons[photons_select[0]].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) and deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,genpart.eta,genpart.phi) < 0.3:
                   photon_isprompt =1
                   break

        self.out.fillBranch("photon_gen_matching",photon_gen_matching)
        self.out.fillBranch("photon_isprompt",photon_isprompt)

        if hasattr(event,'Pileup_nPU'):    
            self.out.fillBranch("npu",event.Pileup_nPU)
        else:
            self.out.fillBranch("npu",0)
    
        if hasattr(event,'Pileup_nTrueInt'):    
            self.out.fillBranch("ntruepu",event.Pileup_nTrueInt)
        else:
            self.out.fillBranch("ntruepu",0)

        print 'channel', channel,'mu_pass:',muon_pass,' ele_pass:',electron_pass,' photon_pass:',photon_pass,' is lepton1 real ',lepton1_isprompt,' is lepton2 real ',lepton2_isprompt,' is photon real ',photon_isprompt,' or ',photon_gen_matching

        self.out.fillBranch("njets50",njets50)
        self.out.fillBranch("njets40",njets40)
        self.out.fillBranch("njets30",njets30)
        self.out.fillBranch("njets20",njets20)
        self.out.fillBranch("njets15",njets15)
        self.out.fillBranch("npvs",event.PV_npvs)
        self.out.fillBranch("met",event.MET_pt)
        self.out.fillBranch("metup",sqrt(pow(event.MET_pt*cos(event.MET_phi) + event.MET_MetUnclustEnUpDeltaX,2) + pow(event.MET_pt*sin(event.MET_phi) + event.MET_MetUnclustEnUpDeltaY,2)))
        self.out.fillBranch("puppimet",event.PuppiMET_pt)
        self.out.fillBranch("puppimetphi",event.PuppiMET_phi)
        self.out.fillBranch("rawmet",event.RawMET_pt)
        self.out.fillBranch("rawmetphi",event.RawMET_phi)
        self.out.fillBranch("metphi",event.MET_phi)
        self.out.fillBranch("pass_selection",1)
        return True

WWG_Module = lambda: WWG_Producer()

