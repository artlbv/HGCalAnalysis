#!/usr/bin/env python
import os,sys
import ROOT, math
import numpy as np
from NtupleDataFormat import HGCalNtuple
from helperTools import *
# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

max_events = 10000

z_half = +1
minE = .5
min_pt = 5

min_fbrem = 0.9
#max_fbrem = 0.1
max_fbrem = 10.3
#max_dR = 0.3
max_dR = 0.3

#mydir = "/grid_mnt/data__data.polcms/cms/lobanov/hgcal/clustering/tuples/myplots/"
mydir = "/eos/user/a/alobanov/www/HGCAL/reco/eleID/compare_ele_ID_vars/dRcomp/dRcuts/myplots/"

def main(fname = "hgcalNtuple-El15-100_noReClust.root"):
    ntuple = HGCalNtuple(fname)

    foutname = fname.replace(".root","_ele_ArbMatch.root")

    if "lobanov" not in foutname:
        print "Saving in another file",
        foutname = mydir + os.path.basename(foutname)
        print foutname

    tfile = ROOT.TFile(foutname,"recreate")

    ## sample type flag
    is_signal = True
    if "qcd" in fname.lower(): is_signal = False
    if "pi" in fname.lower(): is_signal = False

    if is_signal:
        print "# This sample is signal"
    else:
        print "# This sample is background"

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

    ## Gen particle collection prefix
    gen_prefix = "gen"
    #gen_prefix = "genpart"

    ele_cnt = 0
    gen_cnt = 0

    n_duplicates = 0
    n_uniq_ele = 0

    for event in ntuple:
        #if stop_run: break
        if tot_nevents >= max_events: break

        if tot_nevents % 100 == 0: print("Event %i" % tot_nevents)
        # print "Event", event.entry()
        tot_nevents += 1

        #####
        good_genparts = []
        # Select good genparticles
        #genParts = event.genParticles()
        if is_signal:
            genParts = event.genParticles(gen_prefix)
            #tot_genpart += len(genParts)
            #print "Contains %i genparts" % len(genParts)

            for i_part, part in enumerate(genParts):
                if not is_signal: break

                if abs(part.eta()) < 1.5: continue
                if abs(part.eta()) > 2.9: continue
                #if i_part > 10: continue
                #print part.status()
                #if part.pt() < 10: continue

                if gen_prefix == "genpart":
                    if part.gen() < 1: continue
                    if part.reachedEE() < 2: continue
                elif gen_prefix == "gen":
                    #if part.status() != 23 and part.status() != 1: continue
                    #if part.status() != 23: continue
                    if part.status() != 1: continue
                    if abs(part.pdgid()) != 11: continue

                #if part.fbrem() > max_fbrem: continue
                #if abs(part.pid()) != 11: continue

                good_genparts.append(part)

                if len(good_genparts) > 1: break

            #print("Event %i" % event.event())
            #print("Found %i good gen particles" %len(good_genparts))
            tot_genpart += len(good_genparts)

        if is_signal and tot_genpart == 0: continue
        #####
        # Select good electrons

        electrons = event.electrons()
        #pfmultiClusters = event.pfClustersFromMultiCl()

        if len(electrons) == 0: continue
        #print("Number of electrons: %i" %len(electrons))

        good_electrons = []
        for i_ele, ele in enumerate(electrons):

            if ele.isEB(): continue
            if ele.ele_siguu() == -1:
                #print("Problem!")
                continue

            if ele.pt() < 5: continue

            if hasattr(ele,"energyEE"):
                et = ele.energyEE()/math.cosh(ele.eta())
                #if et < 5: continue
                #if ele.energyEE() < 10: continue
            #if ele.seedenergy() < 20: continue

            ## select SCs with 1 cluster (seed=electron)
            #if ele.numClinSC() != 1: continue
            good_electrons.append(ele)

        if len(good_electrons) == 0: continue
        #tot_ele += len(good_electrons)

        ##########################################
        ############## MATCHING ##################
        ##########################################

        ele_match_ids = { i:[] for i in range(len(good_genparts)) }
        matched_electrons = []

        ## Match gsfElectrons to gen electrons
        if is_signal:
            found_match = False

            for i_gen, part in enumerate(good_genparts):

                part_tlv = ROOT.TLorentzVector()
                part_tlv.SetPtEtaPhiE(part.pt(), part.eta(), part.phi(), part.energy())

                for i_ele,ele in enumerate(good_electrons):
                    ele_tlv = ROOT.TLorentzVector()
                    ele_tlv.SetPtEtaPhiE(ele.pt(), ele.eta(), ele.phi(), ele.energy())

                    if ele.eta() * part.eta() < 0: continue

                    dR = part_tlv.DeltaR(ele_tlv)
                    addDataPoint(hist_data,"part_ele_dR",dR)

                    if dR > 0.1: continue

                    # found match
                    found_match = True
                    ele_match_ids[i_gen].append(i_ele)
                    matched_electrons.append(ele)

                ## end ele loop
            ## end gen loop
            if not found_match: continue

        else:
            # if background: save all electrons
            matched_electrons = good_electrons

        # Arbitrate in case of multiple matches per gen particle / overlapping gsfElectons - remove duplicates
        clean_electrons = []

        for i1, ele1 in enumerate(matched_electrons):
            is_best = True

            ele1_tlv = ROOT.TLorentzVector(); ele1_tlv.SetPtEtaPhiE(ele1.pt(), ele1.eta(), ele1.phi(), ele1.energy())

            for i2, ele2 in enumerate(matched_electrons):
                if i1 == i2: continue
                if ele1.eta() * ele2.eta() < 0: continue
                if abs(ele1.seedenergy() - ele2.seedenergy()) > 0.1: continue

                ele2_tlv = ROOT.TLorentzVector(); ele2_tlv.SetPtEtaPhiE(ele2.pt(), ele2.eta(), ele2.phi(), ele2.energy())

                dR = ele1_tlv.DeltaR(ele2_tlv)

                addDataPoint(hist_data,"gsf_ele_dR",dR)

                if dR > 0.3: continue

                #n_duplicates += 1

                ## choose better ele
                if hasattr(ele1,"energyEE") and hasattr(ele2,"energyEE"):
                    if ele1.energyEE() < ele2.energyEE():
                        is_best = False
                        #continue
                        break
                #else: is_best = False
                if abs(ele1.eEleClusterOverPout() - 1) > abs(ele2.eEleClusterOverPout() - 1):
                    is_best = False
                    break

            if is_best:
                #n_uniq_ele += 1
                clean_electrons.append(ele1)
                break

        n_duplicates += len(matched_electrons) - len(clean_electrons)
        n_uniq_ele += len(clean_electrons)

        if len(clean_electrons) > 2:
            #print clean_electrons
            print
            for ele in clean_electrons:
                #print ele.eta(), ele.phi(), ele.pt(), ele.seedenergy(), ele.energyEE(), ele.seedeta(), ele.seedphi()
                print ele.eta(), ele.phi(), ele.pt(), ele.seedenergy(), ele.seedeta(), ele.seedphi()

        '''
        clean_electrons = []
        if is_signal:
            for gen in ele_match_ids:
                n_match = len(ele_match_ids[gen])
                gen_cnt += 1
                ele_cnt += n_match

                addDataPoint(hist_data,"gen_nmatch_ele",n_match)

                if n_match == 1:
                    # clean_electrons.append(good_electrons[ele_match_ids[gen][0]))
                    continue

                ele_idxs = ele_match_ids[i_gen]

                if len(ele_idxs) < 2: continue

                # Arbitrate in case of multiple matches per gen particle
                best_dEoP = -1
                best_idx = ele_idxs[0]

                for i_ele in ele_idxs:

                    if i_ele == best_idx: continue

                    # previous best ele
                    ele1 = good_electrons[best_idx]
                    # previous best ele
                    ele2 = good_electrons[i_ele]



        else:
            clean_electrons = good_electrons

        for ele in clean_electrons:

            tot_ele += 1

            if ele.ele_realDepth() > 500: continue

            et = ele.seedenergy()/math.cosh(ele.eta()

            addDataPoint(hist_data,"ele_seedET_vs_PT",(ele.pt(),ele.seedenergy()/math.cosh(ele.eta())))
            if hasattr(ele,"energyEE"):
                addDataPoint(hist_data,"ele_energy",ele.energyEE())
                addDataPoint(hist_data,"ele_ET",et)
                addDataPoint(hist_data,"ele_seedenergy_vs_eneEE",(ele.seedenergy(),ele.energyEE()))

            if is_signal:
                addDataPoint(hist_data,"part_ele_PT",(part.pt(),ele.pt()))

                if hasattr(ele,"energyEE"):
                    addDataPoint(hist_data,"ele_ET",et)
                    addDataPoint(hist_data,"ele_ET_vs_PT",(ele.pt(),et))

                    addDataPoint(hist_data,"part_PT_vs_ele_ET",(part.pt(),et))
                    if ele.energyEE() < 100:
                        addDataPoint(hist_data,"part_PT_vs_ele_EE",(part.pt(),ele.energyEE()))

            addDataPoint(hist_data,"ele_siguu",ele.ele_siguu())
            addDataPoint(hist_data,"ele_sigvv",ele.ele_sigvv())
            addDataPoint(hist_data,"ele_NLay",ele.ele_nlay())
            addDataPoint(hist_data,"ele_firstLay",ele.ele_firstlay())
            addDataPoint(hist_data,"ele_lastLay",ele.ele_lastlay())

            addDataPoint(hist_data,"ele_posZ",abs(ele.ele_pcaPosZ()))

            if ele.fbrem() > 0 and ele.fbrem() < 1:
                addDataPoint(hist_data,"ele_sigvv_vs_fbrem",(ele.fbrem(),ele.ele_sigvv()))

            addDataPoint(hist_data,"ele_varZ",ele.ele_pcaEigSig3())
            addDataPoint(hist_data,"ele_pcaSig1",ele.ele_pcaEigSig1())
            addDataPoint(hist_data,"ele_pcaSig2",ele.ele_pcaEigSig2())
            addDataPoint(hist_data,"ele_pcaSig3",ele.ele_pcaEigSig3())
            addDataPoint(hist_data,"ele_pcaEig1",ele.ele_pcaEigVal1())
            addDataPoint(hist_data,"ele_pcaEig2",ele.ele_pcaEigVal2())
            addDataPoint(hist_data,"ele_pcaEig3",ele.ele_pcaEigVal3())

            addDataPoint(hist_data,"ele_layEfrac10",ele.ele_layEfrac10())
            addDataPoint(hist_data,"ele_layEfrac90",ele.ele_layEfrac90())
            addDataPoint(hist_data,"ele_layEfrac10_wrt1",ele.ele_layEfrac10() - ele.ele_firstlay())
            addDataPoint(hist_data,"ele_layEfrac90_wrt1",ele.ele_layEfrac90() - ele.ele_firstlay())
            #addDataPoint(hist_data,"ele_outEnergyFrac",ele.ele_outEnergy()/ele.seedenergy())
            #addDataPoint(hist_data,"ele_outEnergyFrac",ele.ele_outEnergy()/ele.energyEE())

            addDataPoint(hist_data,"ele_depthCompat",ele.ele_depthCompat())
            addDataPoint(hist_data,"ele_realDepth",ele.ele_realDepth())
            addDataPoint(hist_data,"ele_predDepth",ele.ele_predDepth())
            addDataPoint(hist_data,"ele_predDepthSigma",ele.ele_predDepthSigma())

            addDataPoint(hist_data,"ele_EE4overE",ele.ele_EE4overEE())
            addDataPoint(hist_data,"ele_FHoverE",ele.ele_FHoverEE())
            addDataPoint(hist_data,"ele_HoverE",ele.ele_HoverEE())

            addDataPoint(hist_data,"ele_numCl",ele.numClinSC())
            ## track vars
            addDataPoint(hist_data,"ele_EoverPout",ele.eEleClusterOverPout())
            addDataPoint(hist_data,"ele_dEtaEle",ele.deltaEtaEleClusterTrackAtCalo())
            addDataPoint(hist_data,"ele_dPhiEle",ele.deltaPhiEleClusterTrackAtCalo())

            #if "qcd" in fname.lower(): break

        #print ele_match_ids

        '''
        for gen in ele_match_ids:
            n_match = len(ele_match_ids[gen])
            gen_cnt += 1
            ele_cnt += n_match

            addDataPoint(hist_data,"gen_nmatch_ele",n_match)

            if n_match > 1:
                #print "# Event", event.event()
                #print "(ele.charge(),ele.eta(),ele.phi(),ele.pt(),ele.sceta()),ele.sceta(),ele.seedenergy(),ele.seedlayer())"
                #for vals in ele_match_ids[gen]:
                #    print vals
                pass

    print("Found %i gen particles and %i ele and %i multicl" %(tot_genpart,tot_ele, tot_multiclus))
    print("Match %i gen particles and %i ele" %(gen_cnt,ele_cnt))
    print(" also %i gsfEle duplicates and %i unique" %(n_duplicates, n_uniq_ele))

    #Set nbins,xmin,xmax
    ranges = {}
    ranges["part_ele_dR"] = 100,0,0.5
    ranges["part_seed_dR"] = 100,0,0.5
    ranges["ele_siguu"] = 100,0.5,1.5
    ranges["ele_sigvv"] = 100,0.5,1.5
    ranges["ele_NLay"] = 50,0,50
    ranges["ele_firstLay"] = 25,1,26
    ranges["ele_lastLay"] = 52,1,53
    ranges["ele_posZ"] = 30,320,350
    ranges["ele_varZ"] = 100,0,10

    ranges["ele_pcaEig1"] = 50,0.9,1
    ranges["ele_pcaEig2"] = 50,0,0.2
    ranges["ele_pcaEig3"] = 50,0,0.2

    ranges["ele_pcaSig1"] = 100,0,3
    ranges["ele_pcaSig2"] = 100,0,3
    ranges["ele_pcaSig3"] = 100,0,10

    ranges["ele_EE4overE"] = 100,0,0.4
    ranges["ele_FHoverE"] = 100,0,0.05
    ranges["ele_HoverE"] = 100,0,0.05

    ranges["ele_EoverPout"] = 100,0,4
    ranges["ele_dEtaEle"] = 100,-0.05,0.05
    ranges["ele_dPhiEle"] = 100,-0.05,0.05

    ranges["ele_numCl"] = 20,0,20

    ranges["ele_layEfrac10"] = 30,1,31
    ranges["ele_layEfrac90"] = 50,1,51

    ranges["ele_depthCompat"] = 100,-10,10
    ranges["ele_realDepth"] = 100,0,20
    ranges["ele_predDepth"] = 100,5,15
    ranges["ele_predDepthSigma"] = 100,0.85,1.2


    hists = []
    for data_name in hist_data:
        print("Plotting hist for data: %s" %data_name)
        canv = ROOT.TCanvas("canv_" + data_name,data_name,800,600)

        if data_name in ranges:
            nbins,xmin,xmax = ranges[data_name]
            hist = make1DHist(hist_data[data_name],data_name,data_name, nbins,xmin,xmax)
        else:
            hist = getHisto(hist_data[data_name],data_name,data_name)

        hist.Draw("colz")
        canv.Update()
        #hists.append(hist)
        ROOT.SetOwnership(canv,0)

        #canv.Write()
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
