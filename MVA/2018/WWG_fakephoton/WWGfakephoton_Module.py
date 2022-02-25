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

        self.out.branch("event",  "F")
        self.out.branch("run",  "F")
        self.out.branch("lumi",  "F")
	self.out.branch("photon_selection",  "I");
	self.out.branch("pass_selection",  "B");

        self.out.branch("channel",  "I")
        self.out.branch("lepton1_pid",  "I")
        self.out.branch("lepton2_pid",  "I")
        self.out.branch("lepton1_pt",  "F")
        self.out.branch("lepton2_pt",  "F")
        self.out.branch("lepton1_eta",  "F")
        self.out.branch("lepton2_eta",  "F")
        self.out.branch("lepton1_phi",  "F")
        self.out.branch("lepton2_phi",  "F")
        self.out.branch("lepton1_mvaTTH",  "F")
        self.out.branch("lepton2_mvaTTH",  "F")
        self.out.branch("lepton1_miniISO",  "F")
        self.out.branch("lepton2_miniISO",  "F")
        self.out.branch("is_lepton_tight", "I")
	self.out.branch("lepton1_isprompt", "I")
	self.out.branch("lepton2_isprompt", "I")
        self.out.branch("n_bjets","I")
        self.out.branch("n_photon", "I")
        self.out.branch("photonet",  "F")
        self.out.branch("photoneta",  "F")
        self.out.branch("photonphi",  "F")
        self.out.branch("photonchiso",  "F")
        self.out.branch("photonsieie",  "F")
        self.out.branch("drl1a",  "F")
        self.out.branch("drl2a",  "F")
        self.out.branch("photon_isprompt", "I")
        self.out.branch("photon_gen_matching", "I")
        self.out.branch("mll",  "F")
        self.out.branch("mllg",  "F")
        self.out.branch("ptll",  "F")
        #self.out.branch("mt",  "F")
        #self.out.branch("puppimt",  "F")
        self.out.branch("met",  "F")
        self.out.branch("puppimet","F")
        self.out.branch("gen_weight","F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        self.out.fillBranch("event",event.event)
        self.out.fillBranch("lumi",event.luminosityBlock)
        self.out.fillBranch("run",event.run)
#        print event.event,event.luminosityBlock,event.run
        if hasattr(event,'Generator_weight'):
            self.out.fillBranch("gen_weight",event.Generator_weight)
        else:    
            self.out.fillBranch("gen_weight",0)

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")

	isprompt_mask = (1 << 0) #isPrompt
	isdirectprompttaudecayproduct_mask = (1 << 5) #isDirectPromptTauDecayProduct
	isdirecttaudecayproduct_mask = (1 << 4) #isDirectTauDecayProduct
	isprompttaudecayproduct = (1 << 3) #isPromptTauDecayProduct
	isfromhardprocess_mask = (1 << 8) #isPrompt
	   
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
        loose_but_not_tight_muons = []
        loose_but_not_tight_electrons = []

        #selection on muons
        muon_pass =0
        for i in range(0,len(muons)):
            if muons[i].pt < 10:
                continue
            if abs(muons[i].eta) > 2.5:
                continue
            if muons[i].mediumId == True and muons[i].mvaTTH > -0.2 and muons[i].miniPFRelIso_all < 0.4:
                muons_select.append(i)
                muon_pass += 1
                leptons_select.append(i)
            elif muons[i].mediumId == True and muons[i].mvaTTH < -0.2 and muons[i].miniPFRelIso_all < 0.4:
                 loose_but_not_tight_muons.append(i)

        # selection on electrons
        electron_pass=0
        loose_electron_pass=0
        for i in range(0,len(electrons)):
            if electrons[i].pt < 10:
                continue
            if abs(electrons[i].eta + electrons[i].deltaEtaSC) > 2.5:
                continue
            if (abs(electrons[i].eta + electrons[i].deltaEtaSC) < 1.479 and abs(electrons[i].dz) < 0.1 and abs(electrons[i].dxy) < 0.05) or (abs(electrons[i].eta + electrons[i].deltaEtaSC) > 1.479 and abs(electrons[i].dz) < 0.2 and abs(electrons[i].dxy) < 0.1):
                if electrons[i].mvaFall17V2Iso_WP80==True:
                    electrons_select.append(i)
                    electron_pass += 1
                    leptons_select.append(i)
                elif electrons[i].cutBased >= 1:
                    loose_but_not_tight_electrons.append(i)

        # selection on photons for tight photon
	photon_pass=0
	for i in range(0,len(photons)):
           # if photons[i].pt < 20:
           #     continue
           # if abs(photons[i].eta) > 2.5:
           #     continue
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

            #| pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho |
            mask1 = 0b11111111111111 # full medium ID
            mask2 = 0b00111111111111 # fail Isopho
            mask3 = 0b11001111111111 # fail IsoNeu
            mask4 = 0b11110011111111 # fail Isoch
            mask5 = 0b11111100111111 # fail sigma ieie

            bitmap = photons[i].vidNestedWPBitmap & mask1

            #the photon pass the full ID or only fail Isoch/sieie
            #if not ((bitmap == mask1) or (bitmap == mask2) or (bitmap == mask3) or (bitmap == mask4) or (bitmap == mask5)):
            if not ((bitmap == mask1) or (bitmap == mask4) or (bitmap == mask5)):
                continue

            photons_select.append(i)
            photon_pass += 1

        self.out.fillBranch("n_photon",photon_pass) #the number of medium photons and control photons


	lepton1_isprompt=-10
	lepton2_isprompt=-10
        channel=-10
        photon_selection=-10
        if len(muons_select) == 2 and len(electrons_select)==0:
	    channel=1
            if len(muons_select) == 2:
                muon_index1 = muons_select[0]
                muon_index2 = muons_select[1]
                self.out.fillBranch("is_lepton_tight",1)
            if hasattr(event, 'nGenPart'):
	       for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muon_index1].eta,muons[muon_index1].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1
                       break

	       for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muon_index2].eta,muons[muon_index2].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1
                       break
            njets=0
            n_bjets=0
            pass_lepton_dr_cut = True
            for i in range(0,len(jets)):
                if jets[i].btagDeepB > 0.4168 and i<=6 :  # DeepCSVM
                   n_bjets += 1
                if abs(jets[i].eta) > 4.7:
                   continue
                if jets[i].pt<20:
                   continue
		if deltaR(jets[i].eta,jets[i].phi,muons[muon_index1].eta,muons[muon_index1].phi) < 0.3 or  deltaR(jets[i].eta,jets[i].phi,muons[muon_index2].eta,muons[muon_index2].phi) < 0.3:
                       pass_lepton_dr_cut = False

                if  not pass_lepton_dr_cut == True:
	            continue
                if jets[i].jetId >> 1 & 1:
                   jets_select.append(i)
                   njets += 1
#           print 'fake muon, the number of jets ',njets
            #if njets <1 :
            #   return False
            self.out.fillBranch("n_bjets",n_bjets)
            self.out.fillBranch("channel",channel) # muon-muon channel
            self.out.fillBranch("lepton1_pt",muons[muon_index1].pt)
            self.out.fillBranch("lepton2_pt",muons[muon_index2].pt)
            self.out.fillBranch("lepton1_eta",muons[muon_index1].eta)
            self.out.fillBranch("lepton2_eta",muons[muon_index2].eta)
            self.out.fillBranch("lepton1_phi",muons[muon_index1].phi)
            self.out.fillBranch("lepton2_phi",muons[muon_index2].phi)
            self.out.fillBranch("lepton1_pid",muons[muon_index1].pdgId)
            self.out.fillBranch("lepton2_pid",muons[muon_index2].pdgId)
            self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
            self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lepton1_mvaTTH",muons[muon_index1].mvaTTH)
            self.out.fillBranch("lepton2_mvaTTH",muons[muon_index2].mvaTTH)
            self.out.fillBranch("lepton1_miniISO",muons[muon_index1].miniPFRelIso_all)
            self.out.fillBranch("lepton2_miniISO",muons[muon_index2].miniPFRelIso_all)
            self.out.fillBranch("mll",(muons[muon_index1].p4() + muons[muon_index2].p4()).M())
            self.out.fillBranch("ptll",(muons[muon_index1].p4() + muons[muon_index2].p4()).Pt())

            self.out.fillBranch("met",event.MET_pt)
            self.out.fillBranch("puppimet",event.PuppiMET_pt)
	elif (len(electrons_select)  == 2) and (len(muons_select)  == 0):

	    channel=2
            if len(electrons_select) == 2:
                electron_index1 = electrons_select[0]
                electron_index2 = electrons_select[1]
                self.out.fillBranch("is_lepton_tight",1)


            if hasattr(event, 'nGenPart'):
	       for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electron_index1].eta,electrons[electron_index1].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1
                       break
	       for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electron_index2].eta,electrons[electron_index2].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1
                       break

            njets=0
            n_bjets=0
            pass_lepton_dr_cut = True
            for i in range(0,len(jets)):
                if jets[i].btagDeepB > 0.4184 and i<=6 :  # DeepCSVM
                   n_bjets +=1       
                if abs(jets[i].eta) > 4.7:
                   continue
                if jets[i].pt<20:
                   continue
		if deltaR(jets[i].eta,jets[i].phi,electrons[electron_index1].eta,electrons[electron_index1].phi) < 0.3 or deltaR(jets[i].eta,jets[i].phi,electrons[electron_index2].eta,electrons[electron_index2].phi) < 0.3 :
                       pass_lepton_dr_cut = False

                if  not pass_lepton_dr_cut == True:
	            continue
                if jets[i].jetId >> 1 & 1:
                   jets_select.append(i)
                   njets += 1
#           print 'fake ele, the number of jets ',njets
            #if njets <1 :
            #   return False
            self.out.fillBranch("n_bjets",n_bjets)
            self.out.fillBranch("channel",channel) #electron-electron channel
            self.out.fillBranch("lepton1_pt",electrons[electron_index1].pt)
            self.out.fillBranch("lepton2_pt",electrons[electron_index2].pt)
            self.out.fillBranch("lepton1_eta",electrons[electron_index1].eta)
            self.out.fillBranch("lepton2_eta",electrons[electron_index2].eta)
            self.out.fillBranch("lepton1_phi",electrons[electron_index1].phi)
            self.out.fillBranch("lepton2_phi",electrons[electron_index2].phi)
            self.out.fillBranch("lepton1_pid",electrons[electron_index1].pdgId)
            self.out.fillBranch("lepton2_pid",electrons[electron_index2].pdgId)
            self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
            self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lepton1_mvaTTH",electrons[electron_index1].mvaTTH)
            self.out.fillBranch("lepton2_mvaTTH",electrons[electron_index2].mvaTTH)
            self.out.fillBranch("lepton1_miniISO",electrons[electron_index1].miniPFRelIso_all)
            self.out.fillBranch("lepton2_miniISO",electrons[electron_index2].miniPFRelIso_all)
            self.out.fillBranch("mll",(electrons[electron_index1].p4() + electrons[electron_index2].p4()).M())
            self.out.fillBranch("ptll",(electrons[electron_index1].p4() + electrons[electron_index2].p4()).Pt())
            self.out.fillBranch("met",event.MET_pt)
            self.out.fillBranch("puppimet",event.PuppiMET_pt)
	else:
	    return False

        self.out.fillBranch("met",event.MET_pt)
        self.out.fillBranch("puppimet",event.PuppiMET_pt)
        #print 'lepton1 is prompt', lepton1_isprompt,'lepton2 is prompt', lepton2_isprompt,' lepton1_pid',lepton1_pid,' the number of jets ',njets,'-> this event is saved'
        photon_gen_matching=-10
        photon_isprompt =-10
        if photon_pass>0:
           for j in range(0,len(photons_select)):
               if channel==1 and deltaR(photons[photons_select[j]].eta,photons[photons_select[j]].phi,muons[muons_select[0]].eta,muons[muons_select[0]].phi) > 0.5 and deltaR(photons[photons_select[j]].eta,photons[photons_select[j]].phi,muons[muons_select[1]].eta,muons[muons_select[1]].phi) > 0.5:
                  photon_index=photons_select[j]
                  break
               elif channel==2 and deltaR(photons[photons_select[j]].eta,photons[photons_select[j]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) > 0.5 and deltaR(photons[photons_select[j]].eta,photons[photons_select[j]].phi,electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi) > 0.5:
                  photon_index=photons_select[j] 
                  break
               else:
                  photon_pass=0
                  break
        if photon_pass>0:
           if hasattr(photons[photon_index],'genPartIdx') :
               print 'calculate the photon flag'
               if photons[photon_index].genPartIdx >= 0 and genparts[photons[photon_index].genPartIdx].pdgId  == 22: 
                   if ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)) and (genparts[photons[photon_index].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask):
                       photon_gen_matching = 6
                   elif ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):       
                       if (genparts[photons[photon_index].genPartIdx].genPartIdxMother >= 0 and (abs(genparts[genparts[photons[photon_index].genPartIdx].genPartIdxMother].pdgId) == 11 or abs(genparts[genparts[photons[photon_index].genPartIdx].genPartIdxMother].pdgId) == 13 or abs(genparts[genparts[photons[photon_index].genPartIdx].genPartIdxMother].pdgId) == 15)):
                           photon_gen_matching = 4
                       else:    
                           photon_gen_matching = 5
                   else:
                       photon_gen_matching = 3
               elif photons[photon_index].genPartIdx >= 0 and abs(genparts[photons[photon_index].genPartIdx].pdgId) == 11:     
                   if ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):  
                       photon_gen_matching = 1
                   else:
                       photon_gen_matching = 2
                       
               else:
                   assert(photons[photon_index].genPartFlav == 0)
                   photon_gen_matching = 0
           if hasattr(event, 'nGenPart') and len(photons_select)>0 :
               for j, genpart in enumerate(genparts):
	           if photons[photon_index].genPartIdx >=0 and genpart.pt > 5 and abs(genpart.pdgId) == 22 and ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) and deltaR(photons[photon_index].eta,photons[photon_index].phi,genpart.eta,genpart.phi) < 0.3:
                      photon_isprompt =1
                      break
           mask1 = 0b11111111111111 # full medium ID
           mask2 = 0b00111111111111 # fail Isopho
           mask3 = 0b11001111111111 # fail IsoNeu
           mask4 = 0b11110011111111 # fail Isoch
           mask5 = 0b11111100111111 # fail sigma ieie
        
	   bitmap = photons[photon_index].vidNestedWPBitmap & mask1   
           if (bitmap == mask1):
               photon_selection=1
               self.out.fillBranch("photon_selection",1) #all cuts applied
               print ' channel',channel,'photon_selection ',photon_selection,' the number of jets ',njets,' chiso ',photons[photon_index].pfRelIso03_chg,'-> this event is saved'
           #elif (bitmap == mask2):
           #    self.out.fillBranch("photon_selection",2) # fail Isopho
           #elif (bitmap == mask3):
           #    self.out.fillBranch("photon_selection",3) # fail IsoNeu
           elif (bitmap == mask4):
               photon_selection=4
               self.out.fillBranch("photon_selection",4) # fail Isoch
           elif (bitmap == mask5):
               photon_selection=5
               self.out.fillBranch("photon_selection",5) # fail sigma ieie
	   else:
               assert(0)
           #(photon_selection==2 || photon_selection==3 || photon_selection==4 || photon_selection ==5 )->build fake photon enriched sample 
	   if channel ==1:
               self.out.fillBranch("drl1a",deltaR(photons[photon_index].eta,photons[photon_index].phi,muons[muons_select[0]].eta,muons[muons_select[0]].phi))
               self.out.fillBranch("drl2a",deltaR(photons[photon_index].eta,photons[photon_index].phi,muons[muons_select[1]].eta,muons[muons_select[1]].phi))
               self.out.fillBranch("mllg",(muons[muons_select[0]].p4() + muons[muons_select[1]].p4()+photons[photon_index].p4()).M())
               self.out.fillBranch("photonchiso",photons[photon_index].pfRelIso03_chg)
	       self.out.fillBranch("photonet",photons[photon_index].pt)
               self.out.fillBranch("photoneta",photons[photon_index].eta)
               self.out.fillBranch("photonphi",photons[photon_index].phi)
               self.out.fillBranch("photonsieie",photons[photon_index].sieie)
               print ' channel',channel,'photon_selection ',photon_selection,' the number of jets ',njets,' chiso ',photons[photon_index].pfRelIso03_chg,'-> this event is saved'
           elif channel ==2:
	       self.out.fillBranch("drl1a",deltaR(photons[photon_index].eta,photons[photon_index].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi))
	       self.out.fillBranch("drl2a",deltaR(photons[photon_index].eta,photons[photon_index].phi,electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi))
               self.out.fillBranch("mllg",(electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()+photons[photon_index].p4()).M())
               self.out.fillBranch("photonchiso",photons[photon_index].pfRelIso03_chg)
	       self.out.fillBranch("photonet",photons[photon_index].pt)
               self.out.fillBranch("photoneta",photons[photon_index].eta)
               self.out.fillBranch("photonphi",photons[photon_index].phi)
               self.out.fillBranch("photonsieie",photons[photon_index].sieie)
          # self.out.fillBranch("mllg",(muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()+photons[photon_index].p4()).M())
               print ' channel',channel,'photon_selection ',photon_selection,' the number of jets ',njets,' chiso ',photons[photon_index].pfRelIso03_chg,'-> this event is saved'
        else:
           self.out.fillBranch("photon_selection",0) #if there is no photons selected
           self.out.fillBranch("drl1a",-10)
           self.out.fillBranch("drl2a",-10)
           self.out.fillBranch("photonet",-10)
           self.out.fillBranch("photoneta",-10)
           self.out.fillBranch("photonphi",-10)
           self.out.fillBranch("mllg",-10)
           self.out.fillBranch("photonsieie",-10)
           self.out.fillBranch("photonchiso",-10)
           
        self.out.fillBranch("photon_gen_matching",photon_gen_matching)
        self.out.fillBranch("photon_isprompt",photon_isprompt)
        #print ' channel',channel,'photon_selection ',photon_selection,' the number of jets ',njets,' chiso ',photons[photon_index].pfRelIso03_chg,'-> this event is saved'

        return True
WWGfakephoton_Module = lambda: WWG_Producer()
