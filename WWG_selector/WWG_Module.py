#!/usr/bin/env python
# Analyzer for WWG Analysis based on nanoAOD tools

import os, sys
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import countHistogramsProducer
Np=0
Nm=0
test=0
MET_pass = 0
photon_pass = 0
electron_pass = 0
muon_pass = 0
none_lepton_reject = 0
none_2lepton_reject = 0
dilepton_pass = 0
emu_pass = 0
ee_pass = 0
mumu_pass = 0
btagjet_reject = 0
deltar_reject = 0
mll_reject = 0
pt_reject = 0
none_photon_reject = 0
different_charge_reject = 0
minus_mll=0
n_posi=0
n_minus=0
njet_reject = 0
n_num = 0
class WWG_Producer(Module):
    def __init__(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("event",  "i")
        self.out.branch("MET",  "F")
        self.out.branch("MET_emu",  "F")
        self.out.branch("MET_mumu",  "F")
        self.out.branch("MET_ee",  "F")
        self.out.branch("photon_pt",  "F")
        self.out.branch("photon_pt_emu",  "F")
        self.out.branch("photon_pt_ee",  "F")
        self.out.branch("photon_pt_mumu",  "F")
        self.out.branch("photon_eta",  "F")
        self.out.branch("photon_eta_ee",  "F")
        self.out.branch("photon_eta_emu",  "F")
        self.out.branch("photon_eta_mumu",  "F")
        self.out.branch("photon_phi",  "F")
        self.out.branch("photon_phi_ee",  "F")
        self.out.branch("photon_phi_emu",  "F")
        self.out.branch("photon_phi_mumu",  "F")

        self.out.branch("photon_sieie",  "F")
        self.out.branch("lepton1_pt_ee",  "F")
        self.out.branch("lepton2_pt_ee",  "F")
        self.out.branch("lepton1_pt_emu",  "F")
        self.out.branch("lepton2_pt_emu",  "F")
        self.out.branch("lepton1_pt_mumu",  "F")
        self.out.branch("lepton2_pt_mumu",  "F")


        self.out.branch("Njets",  "i")
        self.out.branch("Njets_emu",  "i")
        self.out.branch("Njets_mumu",  "i")
        self.out.branch("Njets_ee",  "i")
        

        self.out.branch("dilepton_mass_mumu",  "F")
        self.out.branch("dilepton_g_mass_mumu",  "F")
        self.out.branch("dilepton_pt_mumu",  "F")
        self.out.branch("dilepton_mass_emu",  "F")
        self.out.branch("dilepton_g_mass_emu",  "F")
        self.out.branch("dilepton_pt_emu",  "F")
        self.out.branch("dilepton_mass_ee",  "F")
        self.out.branch("dilepton_g_mass_ee",  "F")
        self.out.branch("dilepton_pt_ee",  "F")
        self.out.branch("dilepton_mass",  "F")
        self.out.branch("dilepton_g_mass",  "F")
        self.out.branch("dilepton_pt",  "F")
        self.out.branch("Generator_weight","F")
        self.out.branch("lepton_motherid_mumu",  "F")
        self.out.branch("lepton_motherid_emu",  "F")
        self.out.branch("electrons_is_real", "I")
        self.out.branch("muons_is_real", "I")
        self.out.branch("photons_is_real", "I")

        self.out.branch("channel_mark","i")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")
        jets = Collection(event, "Jet")

        jet_select = [] 
        dileptonp4 = ROOT.TLorentzVector()
        photons_select = []
        electrons_select = []
        muons_select = [] 
        jets_select = []
        leptons_select=[]
        if hasattr(event,'nGenPart'):
            genparts = Collection(event, "GenPart")

	# Record the pass numbers for each cut. Noticed that for efficiency, those who can't pass the MET cut may not be counted because it will pass to next event directly.
        
        global MET_pass
        global photon_pass
        global electron_pass
        global muon_pass
        global none_lepton_reject
        global none_2lepton_reject
        global dilepton_pass
        global ee_pass 
        global emu_pass 
        global mumu_pass 
        global btagjet_reject
        global mll_reject
        global njet_reject
        global pt_reject
        global deltar_reject
        global test
        global none_photon_reject
        global different_charge_reject
        global minus_mll 
        global n_posi 
        global n_minus 
        # selection on MET. Pass to next event directly if fail.
        global n_num
        n_num +=1
        #if event.Generator_weight > 0 :
        #    n_posi +=1
        #else:
        #    n_minus +=1

    #if event.Generator_weight > 0 :
    #    n_posi +=1
    #else:
    #    n_minus +=1

        if  event.PuppiMET_pt>20:
            MET_pass += 1
        else:
            return False  

        #selection on muons
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

        if len(electrons_select)+len(muons_select) != 2:      #reject event if there are not exactly two leptons
            none_2lepton_reject += 1
            return False

        # selection on photons
        for i in range(0,len(photons)):
            #if not photons[i].electronVeto:
            #    continue
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
            if photons[i].cutBasedBitmap >=1: #1==medium
            #if photons[i].cutBased >=2:
                photons_select.append(i)
                photon_pass += 1

        if  len(photons_select)<1:
            none_photon_reject +=1 
            return False                        #reject event if there is not exact one photon in the event 

        njets = 0
        #print len(jets)
        for i in range(0,len(jets)):
            btag_cut = False
            #if jets[i].btagCMVA > -0.5884:  # cMVAv2L
            # if jets[i].btagCMVA > 0.4432:  # cMVAv2M
            # if jets[i].btagCSVV2 > 0.5426:  # CSVv2L
            # if jets[i].btagCSVV2 > 0.8484:  # CSVv2M
            #if jets[i].btagDeepB > 0.2219:  # DeepCSVL
            if jets[i].btagDeepB > 0.4184:  # DeepCSVM
            # if jets[i].btagDeepB > 0.6324:  # DeepCSVL
                btag_cut = True      #initialize
            if abs(jets[i].eta) > 4.7:
                continue
            if jets[i].pt<30:
                continue
            for j in range(0,len(photons_select)):          # delta R cut, if all deltaR(lep,jet) and deltaR(gamma,jet)>0.3, consider jet as a b jet
                if deltaR(jets[i].eta,jets[i].phi,photons[photons_select[j]].eta,photons[photons_select[j]].phi) < 0.5:
                    btag_cut = True
                    #jets_select.append(i)
                    #njets +=1
            for j in range(0,len(electrons_select)):
                if deltaR(jets[i].eta,jets[i].phi,electrons[electrons_select[j]].eta,electrons[electrons_select[j]].phi) < 0.5:
                    btag_cut = True
                    #jets_select.append(i)
                    #njets +=1
            for j in range(0,len(muons_select)):
                if deltaR(jets[i].eta,jets[i].phi,muons[muons_select[j]].eta,muons[muons_select[j]].phi) < 0.5  :
                    btag_cut = True
            if btag_cut == True:
                btagjet_reject += 1
                #return False
            else:
                if jets[i].jetId >> 1 & 1:
                    jets_select.append(i)
                    njets += 1
        #njets = njets + 1
        print ("njets",njets)
        if njets >=2 :
            njet_reject +=1
            return False

        electrons_is_real=-99
        photons_is_real=-99
        muons_is_real = -99

        isprompt_mask = (1 << 0) #isPrompt
        isdirectprompttaudecayproduct_mask = (1 << 5) #isDirectPromptTauDecayProduct
        isfromhardprocess_mask = (1 << 8) #isPrompt

        if hasattr(event, 'nGenPart'):
            genparts = Collection(event, "GenPart")
            is_real_flag = 0
            for i,lep in enumerate(electrons_select):
                for j,genpart in enumerate(genparts):
                    if genpart.pt>5 and abs(genpart.pdgId)==11 and deltaR(electrons[electrons_select[i]].eta, electrons[electrons_select[i]].phi, genpart.eta, genpart.phi) < 0.3:
                        is_real_flag=1
                        break
            if is_real_flag == 1: electrons_is_real=1
            is_real_flag = 0
            for i, mu in enumerate(muons_select):
                for j, genpart in enumerate(genparts):
                    if genpart.pt > 5 and abs(genpart.pdgId) == 13 and deltaR(muons[muons_select[i]].eta, muons[muons_select[i]].phi,genpart.eta,genpart.phi) < 0.3:
                            is_real_flag = 1
                            break
            if is_real_flag == 1:  muons_is_real=1
            is_real_flag = 0
            for i, pho in enumerate(photons_select):
                for j, genpart in enumerate(genparts):
                    if photons[photons_select[i]].genPartIdx >=0:
                        if genpart.pt > 5 and abs(genpart.pdgId) == 22 and ((genparts[photons[photons_select[i]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photons_select[i]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) or (genparts[photons[photons_select[i]].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) and deltaR(photons[photons_select[i]].eta,photons[photons_select[i]].phi,genpart.eta,genpart.phi) < 0.3:
                            is_real_flag =1
                            break
            if is_real_flag == 1: photons_is_real=1
        #initial value
        dileptonmass_emu = -99
        dileptonpt_emu = -99
        met_emu = -99
        dileptonmass_ee = -99
        dileptonpt_ee = -99
        met_ee = -99
        dileptonmass_mumu = -99
        dileptonpt_mumu = -99
        met_mumu = -99
        #dilepton mass selection and channel selection
        channel = 0 
        # emu:     1
        # ee:      2
        # mumu:    3

        # emu
        if len(muons_select)==1 and len(electrons_select)==1:  # emu channel 
            #print (muons[muons_select[0]].pdgId,electrons[electrons_select[0]].pdgId)
            if deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) < 0.5:
                deltar_reject +=1
                return False
            if muons[muons_select[0]].charge * (electrons[electrons_select[0]].charge) > 0:
                different_charge_reject += 1
                return False
            dileptonmass_emu = (muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()).M()
            dileptongmass_emu = (muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()+photons[photons_select[0]].p4()).M()
            dileptonpt_emu = (muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()).Pt()
            met_emu = event.PuppiMET_pt
            if dileptonmass_emu >= 50 and dileptonmass_emu <= 100 :
                mll_reject +=1
                return False
            if dileptonpt_emu <= 40:
                pt_reject +=1
                return False
            if met_emu <60:
                return False
                # if dileptonmass >= 60 and dileptonmass <= 120:
                # print "a=",photons_select, "e=",electrons_select, "mu=",muons_select
            if dileptonmass_emu > 0:
                channel = 1
                emu_pass += 1
            else :
                minus_mll +=1
                return False

        # ee
        if len(muons_select)==0 and len(electrons_select)==2:
            if deltaR(electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi,electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi)<0.5:
                deltar_reject +=1
                return False
            if electrons[electrons_select[0]].charge * electrons[electrons_select[1]].charge >0:
                different_charge_reject += 1
                return False

            dileptonmass_ee = (electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()).M()
            dileptongmass_ee = (electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()+photons[photons_select[0]].p4()).M()
            dileptonpt_ee = (electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()).Pt()
            met_ee = event.PuppiMET_pt
            if dileptonmass_ee >= 50 and dileptonmass_ee <= 100 :
                mll_reject +=1
                return False
            if dileptonpt_ee <= 40:
                pt_reject +=1
                return False
            if met_ee <60:
                return False
            if dileptonmass_ee > 0:
                channel = 2
                ee_pass += 1
            else:
                minus_mll +=1
                return False


        # mumu 
        if len(electrons_select)==0 and len(muons_select)==2:
            if deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,muons[muons_select[1]].eta,muons[muons_select[1]].phi)<0.5:
                deltar_reject +=1
                return False
            if muons[muons_select[0]].charge * muons[muons_select[1]].charge >0:
                different_charge_reject +=1
                return False
            dileptonmass_mumu = (muons[muons_select[0]].p4() + muons[muons_select[1]].p4()).M()
            dileptongmass_mumu = (muons[muons_select[0]].p4() + muons[muons_select[1]].p4()+photons[photons_select[0]].p4()).M()
            met_mumu = event.PuppiMET_pt
            dileptonpt_mumu = (muons[muons_select[0]].p4() + muons[muons_select[1]].p4()).Pt()
                # print "a=",photons_select, "e=",electrons_select, "mu=",muons_select
            if dileptonmass_mumu >= 50 and dileptonmass_mumu <= 100 :
                mll_reject +=1
                return False
            if dileptonpt_mumu <= 40:
                pt_reject +=1
                return False
            if met_mumu <60:
                return False
            if dileptonmass_mumu > 0:
                channel = 3
                mumu_pass += 1
            else :
                minus_mll +=1
                return False

 #    check channel to fill branch
        if channel == 0:
            print len(electrons_select),len(muons_select)
            return False




        #self.out.fillBranch("MET",event.PuppiMET_pt)
        #self.out.fillBranch("photon_pt",photons[photons_select[0]].pt)
        #self.out.fillBranch("photon_eta",photons[photons_select[0]].eta)
        #self.out.fillBranch("photon_phi",photons[photons_select[0]].phi)
        #self.out.fillBranch("photon_sieie",photons[photons_select[0]].sieie)
        #self.out.fillBranch("Njets",njets)
        if channel == 1:
            pirnt("dilepton_mass_emu", dileptonmass_emu)
            self.out.fillBranch("electrons_is_real",electrons_is_real)
            self.out.fillBranch("photons_is_real", photons_is_real)
            self.out.fillBranch("muons_is_real", muons_is_real)
            self.out.fillBranch("dilepton_mass_emu",dileptonmass_emu)
            self.out.fillBranch("dilepton_g_mass_emu",dileptongmass_emu)
            self.out.fillBranch("dilepton_pt_emu",dileptonpt_emu)
            self.out.fillBranch("Njets_emu",njets)
            self.out.fillBranch("MET_emu",met_emu)
            self.out.fillBranch("photon_pt_emu",photons[photons_select[0]].pt)
            self.out.fillBranch("photon_eta_emu",photons[photons_select[0]].eta)
            self.out.fillBranch("photon_phi_emu",photons[photons_select[0]].phi)
            #lep1_pt lep2_pt
            if muons[muons_select[0]].pt > electrons[electrons_select[0]].pt:
                self.out.fillBranch("lepton1_pt_emu",muons[muons_select[0]].pt)
                self.out.fillBranch("lepton2_pt_emu",electrons[electrons_select[0]].pt)
            else:
                self.out.fillBranch("lepton2_pt_emu",muons[muons_select[0]].pt)
                self.out.fillBranch("lepton1_pt_emu",electrons[electrons_select[0]].pt)
            self.out.fillBranch("channel_mark", channel)
            return True
            #self.out.fillBranch("photon_sieie",photons[photons_select[0]].sieie)

        elif channel == 2:
            #self.out.fillBranch("Njets",njets)
            pirnt ("dilepton_mass_ee",dileptonmass_ee)
            self.out.fillBranch("electrons_is_real",electrons_is_real)
            self.out.fillBranch("photons_is_real", photons_is_real)
            self.out.fillBranch("muons_is_real", muons_is_real)
            self.out.fillBranch("dilepton_mass_ee", dileptonmass_ee)
            self.out.fillBranch("dilepton_g_mass_ee",dileptongmass_ee)
            self.out.fillBranch("dilepton_pt_ee",dileptonpt_ee)
            self.out.fillBranch("Njets_ee",njets)
            self.out.fillBranch("MET_ee",met_ee)
            self.out.fillBranch("photon_pt_ee",photons[photons_select[0]].pt)
            self.out.fillBranch("photon_eta_ee",photons[photons_select[0]].eta)
            self.out.fillBranch("photon_phi_ee",photons[photons_select[0]].phi)
            #lep1_pt lep2_pt
            if electrons[electrons_select[0]].pt > electrons[electrons_select[1]].pt:
                self.out.fillBranch("lepton1_pt_ee",electrons[electrons_select[0]].pt)
                self.out.fillBranch("lepton2_pt_ee",electrons[electrons_select[1]].pt)
            else:
                self.out.fillBranch("lepton2_pt_ee",electrons[electrons_select[0]].pt)
                self.out.fillBranch("lepton1_pt_ee",electrons[electrons_select[1]].pt)
            self.out.fillBranch("channel_mark", channel)
            return True
        elif channel == 3:
            #self.out.fillBranch("lepton_motherid_mumu",genparts[genparts[muons[muons_select[0]].genPartIdx].genPartIdxMother].pdgId)
            #self.out.fillBranch("lepton_motherid_mumu",genparts[genparts[muons[muons_select[1]].genPartIdx].genPartIdxMother].pdgId)
            #self.out.fillBranch("Njets",njets)
            pirnt("dilepton_mass_mumu", dileptonmass_mumu)
            self.out.fillBranch("electrons_is_real",electrons_is_real)
            self.out.fillBranch("photons_is_real", photons_is_real)
            self.out.fillBranch("muons_is_real", muons_is_real)
            self.out.fillBranch("dilepton_mass_mumu",dileptonmass_mumu)
            self.out.fillBranch("dilepton_g_mass_mumu",dileptongmass_mumu)
            self.out.fillBranch("dilepton_pt_mumu",dileptonpt_mumu)
            self.out.fillBranch("Njets_mumu",njets)
            self.out.fillBranch("MET_mumu",met_mumu)
            self.out.fillBranch("photon_pt_mumu",photons[photons_select[0]].pt)
            self.out.fillBranch("photon_eta_mumu",photons[photons_select[0]].eta)
            self.out.fillBranch("photon_phi_mumu",photons[photons_select[0]].phi)
            #lep1_pt lep2_pt
            if muons[muons_select[0]].pt > muons[muons_select[1]].pt:
                self.out.fillBranch("lepton1_pt_mumu",muons[muons_select[0]].pt)
                self.out.fillBranch("lepton2_pt_mumu",muons[muons_select[1]].pt)
            else:
                self.out.fillBranch("lepton2_pt_mumu",muons[muons_select[0]].pt)
                self.out.fillBranch("lepton1_pt_mumu",muons[muons_select[1]].pt)
            self.out.fillBranch("channel_mark", channel)
            return True

        #self.out.fillBranch("event",event.event)
        #self.out.fillBranch("dilepton_mass",dileptonmass)
        #self.out.fillBranch("dilepton_g_mass",dileptongmass)
        #self.out.fillBranch("dilepton_pt",dileptonpt)
        #self.out.fillBranch("Generator_weight",event.Generator_weight)
        else:
            print "Unkown channel",channel
            return False


print "MET_pass","\t","=","\t",MET_pass
print "muon_pass","\t","=","\t",muon_pass 
print "electron_pass","\t","=","\t",electron_pass
print "photon_pass","\t","=","\t",photon_pass
print
print "none_photon_reject","\t","=","\t",none_photon_reject
print "none_lepton_reject","\t","=","\t",none_lepton_reject
print "none_2lepton_reject","\t","=","\t",none_2lepton_reject
print "different_charge_reject","\t","=","\t",different_charge_reject
print "minus_mll ","\t","=","\t",minus_mll 
print "all reject",none_photon_reject+none_2lepton_reject+different_charge_reject+minus_mll


print "emu_pass","\t","=","\t",emu_pass
print "ee_pass","\t","=","\t",ee_pass
print "mumu_pass","\t","=","\t",mumu_pass
print "btagjet_reject","\t","=","\t",btagjet_reject
print "deltar_reject","\t","=","\t",deltar_reject
print "n_posi","\t","=","\t",n_posi
print "n_minus","\t","=","\t",n_minus
print "n_num","\t","=","\t",n_num


