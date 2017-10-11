#!/usr/bin/env python
import sys
import ROOT, math
import numpy as np
from NtupleDataFormat import HGCalNtuple
from helperTools import *
# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

max_events = 1000

z_half = +1
minE = .5
min_pt = 5

min_fbrem = 0.9
#max_fbrem = 0.1
max_fbrem = 10.3
#max_dR = 0.3
max_dR = 0.3

def main(fname = "hgcalNtuple-El15-100_noReClust.root"):
    ntuple = HGCalNtuple(fname)

    foutname = fname.replace(".root","_ele.root")
    tfile = ROOT.TFile(foutname,"recreate")

    tot_nevents = 0
    tot_genpart = 0
    tot_rechit = 0
    tot_rechit_raw = 0
    tot_cluster2d = 0
    tot_ele = 0
    tot_multiclus = 0
    tot_simcluster = 0
    tot_pfcluster = 0
    tot_calopart = 0
    tot_track = 0

    hist_data = {}

    # store rechits
    good_cluster_rechits = {}
    bad_cluster_rechits = {}

    for event in ntuple:
        #if stop_run: break
        if tot_nevents >= max_events: break

        if tot_nevents % 100 == 0: print("Event %i" % tot_nevents)
        # print "Event", event.entry()
        tot_nevents += 1


        #####
        # Select good genparticles

        genParts = event.genParticles()
        #tot_genpart += len(genParts)

        good_genparts = []
        for part in genParts:

            if part.gen() < 1: continue
            if part.reachedEE() < 2: continue
            if part.fbrem() > max_fbrem: continue

            good_genparts.append(part)

        #print("Event %i" % event.event())
        #print("Found %i good gen particles" %len(good_genparts))
        tot_genpart += len(good_genparts)

        #####
        # Select good electrons

        electrons = event.electrons()
        #pfmultiClusters = event.pfClustersFromMultiCl()

        if len(electrons) == 0: continue
        #print("Number of electrons: %i" %len(electrons))

        good_electrons = []
        for i_ele, ele in enumerate(electrons):

            if ele.isEB(): continue

            #n_multi =  len([mc for mc in ele.clustersFromMultiCl()])
            #print ("# Ele %i of pt %0.2f has %i multiclusters" %(i_ele,ele.pt(),n_multi ))
            #print ("# Ele %i of pt %0.2f has %i multiclusters" %(i_ele,ele.pt(),ele.numClinSC() ))

            ## select SCs with 1 cluster (seed=electron)
            #if ele.numClinSC() != 1: continue

            #print ele.pt()

            '''
            for mc in ele.clustersFromMultiCl():
                #print mc.eta(), mc.energy(), mc.pt()
                if abs(mc.pt() - ele.pt())/ele.pt() < 0.5:
                    print mc.eta(), ele.eta(), mc.pt(), ele.pt()

            '''
            #ele.seedeta()
            good_electrons.append(ele)

        if len(good_electrons) == 0: continue
        #tot_ele += len(good_electrons)

        ## Select good multiclusters
        multiClusters = event.multiClusters()
        good_mclust = []

        for i_mcl, mclust in enumerate(multiClusters):
            if mclust.NLay() < 10: continue
            if mclust.pt() < 10: continue

            good_mclust.append(mclust)

        tot_multiclus += len(good_mclust)

        ##########################################
        ############## MATCHING ##################
        ##########################################

        for part in good_genparts:

            part_tlv = ROOT.TLorentzVector()
            part_tlv.SetPtEtaPhiE(part.pt(), part.eta(), part.phi(), part.energy())

            for ele in good_electrons:

                if ele.eta() * part.eta() < 0: continue

                ele_tlv = ROOT.TLorentzVector()
                ele_tlv.SetPtEtaPhiE(ele.pt(), ele.eta(), ele.phi(), ele.energy())
                ## use seed instead of ele!
                seed_tlv = ROOT.TLorentzVector()
                seedpt = ele.seedenergy() * ele.pt()/ele.energy()
                seed_tlv.SetPtEtaPhiE(seedpt, ele.seedeta(), ele.seedphi(), ele.energy())

                dR = part_tlv.DeltaR(ele_tlv)
                addDataPoint(hist_data,"part_ele_dR",dR)

                if dR < 0.3: tot_ele += 1

                dR = part_tlv.DeltaR(seed_tlv)
                addDataPoint(hist_data,"part_seed_dR",dR)

                for mclust in good_mclust:

                    if mclust.eta() * part.eta() < 0: continue

                    mclust_tlv = ROOT.TLorentzVector()
                    mclust_tlv.SetPtEtaPhiE(mclust.pt(), mclust.eta(), mclust.phi(), mclust.energy())

                    dR = part_tlv.DeltaR(mclust_tlv)
                    addDataPoint(hist_data,"part_mclust_dR",dR)
                    #addDataPoint(hist_data,"part_mclust_dR_pt",(dR,mclust.pt()))

                    dR = seed_tlv.DeltaR(mclust_tlv)
                    addDataPoint(hist_data,"seed_mclust_dR",dR)
                    #addDataPoint(hist_data,"seed_mclust_dR_pt",(dR,mclust.pt()))

                    if dR > 0.3: continue

                    tot_multiclus += 1

                    # compare variabs
                    delta_NLay = ele.ele_nlay() - mclust.help_NLay()
                    addDataPoint(hist_data,"seed_mclust_NLay",(ele.ele_nlay(),mclust.help_NLay()))
                    addDataPoint(hist_data,"seed_mclust_eta_delta_NLay",(abs(ele.eta()),delta_NLay))

                    if delta_NLay != 0: continue
                    delta_sigvv = ele.ele_sigvv() - mclust.help_sigvv()
                    addDataPoint(hist_data,"seed_mclust_sigvv",(ele.ele_sigvv(),mclust.help_sigvv()))
                    addDataPoint(hist_data,"seed_mclust_eta_delta_sigvv",(abs(ele.eta()),delta_sigvv))
                    addDataPoint(hist_data,"seed_sigvv",ele.ele_sigvv())
                    addDataPoint(hist_data,"mclust_sigvv",mclust.help_sigvv())
                    addDataPoint(hist_data,"seed_sigvv_vs_eta",(abs(ele.eta()),ele.ele_sigvv()))
                    addDataPoint(hist_data,"mclust_sigvv_vs_eta",(abs(mclust.eta()),mclust.help_sigvv()))


        #print("Found good: %i gen particles and %i ele and %i multi" %(len(good_genparts),len(good_electrons),len(good_mclust)))

    #print("Found %i gen particles and %i multicl" %(tot_genpart,tot_multiclus))
    #print("Event %i" % event.event())
    #print("Found %i good gen particles" %len(good_genparts))
    #print("Found good: %i gen particles and %i ele and %i multi" %(len(good_genparts),len(good_electrons),len(good_mclust)))
    print("Found %i gen particles and %i ele and %i multicl" %(tot_genpart,tot_ele, tot_multiclus))


    hists = []
    for data_name in hist_data:
        print("Plotting hist for data: %s" %data_name)
        canv = ROOT.TCanvas("canv_" + data_name,data_name,800,600)
        hist = getHisto(hist_data[data_name],data_name,data_name)
        hist.Draw("colz")
        canv.Update()
        #hists.append(hist)
        ROOT.SetOwnership(canv,0)

        canv.Write()
        hist.Write()
    q = raw_input("exit")

    tfile.Close()

if __name__ == "__main__":

    #if '-b' in sys.argv: sys.argv = [sys.argv[0]]

    if len(sys.argv) > 1:
        if '-b' in sys.argv:
            fname = sys.argv[1]
        else:
            fname = sys.argv[1]
        print '# Input file is', fname
        main(fname)
    else:
        print("No input files given!")
        main()
