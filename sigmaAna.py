#!/usr/bin/env python
import sys
import ROOT, math
import numpy as np
from NtupleDataFormat import HGCalNtuple
from helperTools import *
# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

max_events = 500

z_half = +1
minE = .5
min_pt = 5

min_fbrem = 0.9
#max_fbrem = 0.1
max_fbrem = 1.0
#max_dR = 0.3
max_dR = 0.3

def main(fname = "hgcalNtuple-El15-100_noReClust.root"):
    ntuple = HGCalNtuple(fname)

    foutname = fname.replace(".root","_sigma.root")
    tfile = ROOT.TFile(foutname,"recreate")

    tot_nevents = 0
    tot_genpart = 0
    tot_rechit = 0
    tot_rechit_raw = 0
    tot_cluster2d = 0
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

        genParts = event.genParticles()
        #tot_genpart += len(genParts)

        multiClusters = event.multiClusters()
        #tot_multiclus += len(multiClusters)

        # 2d clusters in layers
        #layerClusters = event.layerClusters()
        #tot_cluster2d += len(layerClusters)

        # rechits
        #recHits = event.recHits()
        #tot_rechit += len(recHits)

        for i_mcl, multicl in enumerate(multiClusters):

            found_part = False
            for part in genParts:

                #if part.eta() * z_half < 0: continue
                if part.gen() < 1: continue
                if not part.reachedEE(): continue
                #if part.fbrem() < min_fbrem: continue
                if part.fbrem() > max_fbrem: continue

                if multicl.z() * part.eta() < 0: continue

                found_part = True
                if i_mcl == 0: tot_genpart += 1

                break

            if not found_part: continue

            part_tlv = ROOT.TLorentzVector()
            part_tlv.SetPtEtaPhiE(part.pt(), part.eta(), part.phi(), part.energy())

            #if multicl.energy() < minE: continue
            #if multicl.pt() < min_pt: continue
            #if multicl.pt() > 10: continue
            #if multicl.z() * part.eta() < 0: continue
            if len(multicl.cluster2d()) < 3: continue

            mcl_tlv = ROOT.TLorentzVector()
            mcl_tlv.SetPtEtaPhiE(multicl.pt(), multicl.eta(), multicl.phi(), multicl.energy())

            dR = float(part_tlv.DeltaR(mcl_tlv))
            #dR2 =  math.hypot(multicl.eta()-part.eta(), multicl.phi()-part.phi())
            #print dR, dR2

            #if dR < 5:
            addDataPoint(hist_data,"part_mcl_dR",dR)
            #addDataPoint(hist_data,"part_mcl_dR_fbrem",(dR,part.fbrem()))

            #if dR > 2 and dR < 10:
            if dR < 0.5:
                tot_multiclus += 1
                '''
                #continue
                addDataPoint(hist_data,"mcl_pt_sige",(multicl.sigee(),multicl.pt()))
                addDataPoint(hist_data,"mcl_pt_sigp",(multicl.sigpp(),multicl.pt()))
                addDataPoint(hist_data,"mcl_pt_sigu",(multicl.siguu(),multicl.pt()))
                addDataPoint(hist_data,"mcl_pt_sigv",(multicl.sigvv(),multicl.pt()))
                '''

                if multicl.sigee() > 0.05 or multicl.siguu() > 4 or multicl.sigpp() > 0.05: continue

                if multicl.pt() > 10:
                    addDataPoint(hist_data,"mcl_dR_sige",(multicl.sigee(),dR))
                    addDataPoint(hist_data,"mcl_dR_sigp",(multicl.sigpp(),dR))
                    addDataPoint(hist_data,"mcl_dR_sigu",(multicl.siguu(),dR))
                    addDataPoint(hist_data,"mcl_dR_sigv",(multicl.sigvv(),dR))
                else:
                    addDataPoint(hist_data,"mcl_lowpt_dR_sige",(multicl.sigee(),dR))
                    addDataPoint(hist_data,"mcl_lowpt_dR_sigp",(multicl.sigpp(),dR))
                    addDataPoint(hist_data,"mcl_lowpt_dR_sigu",(multicl.siguu(),dR))
                    addDataPoint(hist_data,"mcl_lowpt_dR_sigv",(multicl.sigvv(),dR))
                    '''
                if dR < 0: print dR
                continue
                '''

            '''
                #addDataPoint(hist_data,"mcl_axisZ",multicl.pcaAxisZ())

                if part.fbrem() > max_fbrem:
                    addDataPoint(hist_data,"high_brem_sigu",multicl.siguu())
                    addDataPoint(hist_data,"high_brem_sigv",multicl.sigvv())
                    addDataPoint(hist_data,"high_brem_sigp",multicl.sigpp())
                    addDataPoint(hist_data,"high_brem_sige",multicl.sigee())
                elif part.fbrem() < 0.1:
                    addDataPoint(hist_data,"low_brem_sigu",multicl.siguu())
                    addDataPoint(hist_data,"low_brem_sigv",multicl.sigvv())
                    addDataPoint(hist_data,"low_brem_sigp",multicl.sigpp())
                    addDataPoint(hist_data,"low_brem_sige",multicl.sigee())
            elif dR > 1:
                if part.fbrem() > max_fbrem:
                    addDataPoint(hist_data,"far_high_brem_sigu",multicl.siguu())
                    addDataPoint(hist_data,"far_high_brem_sigv",multicl.sigvv())
                    addDataPoint(hist_data,"far_high_brem_sigp",multicl.sigpp())
                    addDataPoint(hist_data,"far_high_brem_sige",multicl.sigee())
                elif part.fbrem() < 0.1:
                    addDataPoint(hist_data,"far_low_brem_sigu",multicl.siguu())
                    addDataPoint(hist_data,"far_low_brem_sigv",multicl.sigvv())
                    addDataPoint(hist_data,"far_low_brem_sigp",multicl.sigpp())
                    addDataPoint(hist_data,"far_low_brem_sige",multicl.sigee())
            '''
        #break

    print("Found %i gen particles and %i multicl" %(tot_genpart,tot_multiclus))

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
