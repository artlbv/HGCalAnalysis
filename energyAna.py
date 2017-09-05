#!/usr/bin/env python
import sys
import ROOT, math
import numpy as np
from NtupleDataFormat import HGCalNtuple
from helperTools import *
# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

#ROOT.gROOT.SetBatch(1)

max_events = 1000

z_half = +1
minE = .5
min_pt = 1

min_fbrem = 0.9
#max_dR = 0.3
max_dR = 0.3

def main(fname = "hgcalNtuple-El15-100_noReClust.root"):
    ntuple = HGCalNtuple(fname)

    foutname = fname.replace(".root","_ene.root")
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
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)

        # rechits
        recHits = event.recHits()
        tot_rechit += len(recHits)

        ## clustered rechits
        clustered_rechits_idx = []
        for mcl in multiClusters:
            for cl_idx in mcl.cluster2d():
                cluster = layerClusters[cl_idx]
                clustered_rechits_idx += [indx for indx in cluster.rechits()]

        #print len(clustered_rechits_idx), len(recHits)
        clustered_rechits = []
        #for rh in recHits:
        for rh_indx in clustered_rechits_idx:
            rh = recHits[rh_indx]
            if rh.flags() > 2: continue
            if rh.layer() > 28: continue
            clustered_rechits.append(rh)


        for i_mcl, multicl in enumerate(multiClusters):
            #if multicl.energy() < minE: continue
            if multicl.pt() < min_pt: continue
            #if len(multicl.cluster2d()) < 3: continue
            if multicl.NLay() < 10: continue

            #if abs(multicl.eta()) > 2.5: continue

            found_part = False
            for part in genParts:

                #if part.eta() * z_half < 0: continue
                if part.gen() < 1: continue
                #if part.gen() > 0: continue
                if part.reachedEE() < 2: continue
                #if part.fbrem() < min_fbrem: continue
                #if part.fbrem() > max_fbrem: continue

                if multicl.z() * part.eta() < 0: continue

                found_part = True
                if i_mcl < 2: tot_genpart += 1

                break

            if not found_part: continue

            if multicl.eta() * part.eta() < 0: continue

            # make lorentz vector for particle
            part_tlv = ROOT.TLorentzVector()
            part_tlv.SetPtEtaPhiE(part.pt(), part.eta(), part.phi(), part.energy())

            # make vector for cluster
            mcl_tlv = ROOT.TLorentzVector()
            mcl_tlv.SetPtEtaPhiE(multicl.pt(), multicl.eta(), multicl.phi(), multicl.energy())
            dR = part_tlv.DeltaR(mcl_tlv)

            #dR =  math.hypot(multicl.eta()-part.eta(), multicl.phi()-part.phi())

            if dR > max_dR: continue

            #a_p,a_v = get_angles(multicl,part)
            #if abs(a_v) < 0.03: continue
            #addDataPoint(hist_data,"part_mcl_a_v",a_v)
            #if abs(a_v) < 0.03: continue

            tot_multiclus += 1


            addDataPoint(hist_data,"part_mcl_dR",dR)
            #addDataPoint(hist_data,"mcl_axisZ",multicl.pcaAxisZ())

            ##### ENERGY
            mcl_ene_core = 0 # energy sum of all core rechits
            mcl_ene_halo = 0 # energy sum of all halo rechits
            mcl_ene_candh = 0 # energy sum of all core and halo rechits
            mcl_ene_cyl = 0 # energy sum of rechits within cylinder
            mcl_ene_cyl_all = 0 # energy sum of all rechits within cylinder

            ## play with cluster layers
            clusters = [layerClusters[clust_idx] for clust_idx in multicl.cluster2d()]

            #mis_layers = get_missing_layers([cl.layer() for cl in clusters])

            rh_cyl_cnt = 0
            #for rh in recHits:
            for rh in clustered_rechits:
                if rh.flags() > 2: continue
                if rh.layer() > 28: continue
                if rh.z() * multicl.pcaPosZ() < 0: continue
                ## filter by distance to cluster axis
                dist = math.hypot(rh.x() - multicl.pcaPosX(), rh.y() - multicl.pcaPosY())
                if dist < 10:
                    mcl_ene_cyl_all += rh.energy()
                    rh_cyl_cnt += 1

            ## Rechit plots
            rhcl_cyl_cnt = 0

            rechits = []
            for i_cl,cluster in enumerate(clusters):

                if cluster.layer() > 28: continue

                #flagsum = sum([1 for rh_idx in cluster.rechits() if recHits[rh_idx].flags() == 0])
                e_sum = sum([recHits[rh_idx].energy() for rh_idx in cluster.rechits()])

                #if flagsum == 0 and e_sum :
                if cluster.energy() == 0 and e_sum > 1:
                    print "Halo only cluster", event.event(), cluster.z(), cluster.layer(), e_sum, multicl.pt(), multicl.NLay()

                rhits = [recHits[rh_idx] for rh_idx in cluster.rechits()]# if recHits[rh_idx].flags() == 0]
                rechits += rhits

                ## analyze halo hits
                cl_seed = recHits[cluster.rechitSeed()]
                #print cl_seed.x(), cl_seed.y(), cl_seed.energy()
                nclu = sum([cl.layer()==cluster.layer() for cl in clusters])

                for rh in rhits:

                    if rh.flags() == 0:
                        mcl_ene_core += rh.energy()
                        mcl_ene_candh += rh.energy()
                    elif rh.flags() == 2:
                        mcl_ene_halo += rh.energy()
                        mcl_ene_candh += rh.energy()

                    ## filter by distance to cluster axis
                    dist = math.hypot(rh.x() - multicl.pcaPosX(), rh.y() - multicl.pcaPosY())
                    if dist < 10:
                        mcl_ene_cyl += rh.energy()
                        rhcl_cyl_cnt += 1

            print rh_cyl_cnt, rhcl_cyl_cnt

            ### fill energy
            #addDataPoint(hist_data,"mcl_ene_orig",multicl.energy())
            addDataPoint(hist_data,"mcl_pt_orig", multicl.pt())

            eta_factor = multicl.pt() / multicl.energy()

            #addDataPoint(hist_data,"mcl_ene_core", mcl_ene_core)
            addDataPoint(hist_data,"mcl_pt_core", mcl_ene_core * eta_factor)

            #addDataPoint(hist_data,"mcl_ene_halo", mcl_ene_halo)
            addDataPoint(hist_data,"mcl_pt_halo", mcl_ene_halo * eta_factor)

            #addDataPoint(hist_data,"mcl_ene_candh", mcl_ene_candh)
            addDataPoint(hist_data,"mcl_pt_candh", mcl_ene_candh * eta_factor)

            #addDataPoint(hist_data,"mcl_ene_cyl", mcl_ene_cyl)
            addDataPoint(hist_data,"mcl_pt_cyl", mcl_ene_cyl * eta_factor)

            #addDataPoint(hist_data,"mcl_ene_cyl", mcl_ene_cyl_all)
            addDataPoint(hist_data,"mcl_pt_cyl_all", mcl_ene_cyl_all * eta_factor)

    print("Found %i gen particles and %i multicl" %(tot_genpart,tot_multiclus))

    hists = []
    for data_name in hist_data:
        print("Plotting hist for data: %s" %data_name)
        canv = ROOT.TCanvas("canv_" + data_name,data_name,800,600)
        hist = getHisto(hist_data[data_name],data_name,data_name)
        hist.Draw("colz")
        canv.Update()
        #hists.append(hist)
        hist.Write()
        ROOT.SetOwnership(canv,0)
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
