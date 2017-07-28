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
        for part in genParts:

            #if part.eta() * z_half < 0: continue
            if part.gen() < 1: continue
            #if part.gen() > 0: continue
            if not part.reachedEE(): continue
            #if part.fbrem() < min_fbrem: continue
            #if part.fbrem() > max_fbrem: continue

            part_vect = ROOT.TVector3(
                part.posx()[1]-part.posx()[0],
                part.posy()[1]-part.posy()[0],
                part.posz()[1]-part.posz()[0],
            )

            track_vect = ROOT.TVector3(
                part.posx()[0],
                part.posy()[0],
                part.posz()[0],
            )

            part_angle =  part_vect.Angle(track_vect)
            addDataPoint(hist_data,"part_angle",part_angle)
            #addDataPoint(hist_data,"part_angle_fbrem",(part_angle,part.fbrem()))

        #continue

        multiClusters = event.multiClusters()

        for i_mcl, multicl in enumerate(multiClusters):

            #if len(multicl.cluster2d()) < 3: continue
            if multicl.NLay() < 3: continue

            found_part = False
            for part in genParts:

                #if part.eta() * z_half < 0: continue
                if part.gen() < 1: continue
                #if part.gen() > 0: continue
                if not part.reachedEE(): continue
                #if part.fbrem() < min_fbrem: continue
                #if part.fbrem() > max_fbrem: continue

                if multicl.z() * part.eta() < 0: continue

                found_part = True
                if i_mcl == 0: tot_genpart += 1

                break

            if not found_part: continue


            #######
            # IP-vector based stuff
            #######
            part_tlv = ROOT.TLorentzVector()
            part_tlv.SetPtEtaPhiE(part.pt(), part.eta(), part.phi(), part.energy())

            if multicl.energy() < 1: continue
            #if multicl.pt() < min_pt: continue
            #if multicl.pt() < 10: continue
            #if multicl.pt() < 1: continue
            #if multicl.z() * part.eta() < 0: continue

            mcl_tlv = ROOT.TLorentzVector()
            mcl_tlv.SetPtEtaPhiE(multicl.pt(), multicl.eta(), multicl.phi(), multicl.energy())

            dR = float(part_tlv.DeltaR(mcl_tlv))
            #dR2 =  math.hypot(multicl.eta()-part.eta(), multicl.phi()-part.phi())
            #print dR, dR2

            #if dR < 5:
            addDataPoint(hist_data,"part_mcl_dR",dR)
            #addDataPoint(hist_data,"part_mcl_dR_fbrem",(dR,part.fbrem()))


            #######
            # real vector based stuff
            #######
            if abs(multicl.pcaAxisZ() < 0.1): continue

            ## 1. get particle propagation vector
            part_vect = ROOT.TVector3(
                part.posx()[1]-part.posx()[0],
                part.posy()[1]-part.posy()[0],
                part.posz()[1]-part.posz()[0],
            )

            ## 2. get mulcutluster vector
            mclut_vect = ROOT.TVector3(
                multicl.pcaAxisX(),
                multicl.pcaAxisY(),
                multicl.pcaAxisZ(),
            )

            #if abs(mclut_vect.Mag()-1) > 0.1: continue
            angle = part_vect.Angle(mclut_vect)
            #print part_vect.Mag(), mclut_vect.Mag()
            #print part_vect.Angle(mclut_vect), mclut_vect.Angle(part_vect)
            #print multicl.pcaAxisX(),multicl.pcaAxisY(),multicl.pcaAxisZ()

            addDataPoint(hist_data,"part_mcl_angle",angle)
            addDataPoint(hist_data,"part_mcl_angle_pt",(angle,multicl.pt()))
            addDataPoint(hist_data,"part_mcl_angle_pt_br",(angle,multicl.pt(),part.fbrem()))

            addDataPoint(hist_data,"part_mcl_angle_pt_br",(angle,multicl.pt(),part.fbrem()))
            addDataPoint(hist_data,"part_mcl_angle_br",(angle,part.fbrem()))
            addDataPoint(hist_data,"part_mcl_dR_br",(dR,part.fbrem()))

            addDataPoint(hist_data,"dR_vs_angle",(angle,dR))
            #addDataPoint(hist_data,"dR_vs_angle_vs_pt",(angle,dR,multicl.pt()))
            addDataPoint(hist_data,"dR_vs_angle_vs_ene",(angle,dR,multicl.energy()))


    print("Found %i gen particles and %i multicl" %(tot_genpart,tot_multiclus))

    hists = []
    for data_name in hist_data:
        print("Plotting hist for data: %s" %data_name)
        canv = ROOT.TCanvas("canv_" + data_name,data_name,800,600)
        hist = getHisto(hist_data[data_name],data_name,data_name)
        #hist.SetDirectory(0)
        ROOT.SetOwnership(hist,0)

        if "TH" in hist.ClassName():
            hist.Draw("colz")
        elif "TGraph" in hist.ClassName():
            hist.Draw("pcolz")
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
