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
min_pt = 5

min_fbrem = 0.9
max_fbrem = 0.3
#max_dR = 0.3
max_dR = 0.3

def main(fname = "hgcalNtuple-El15-100_noReClust.root"):
    ntuple = HGCalNtuple(fname)

    foutname = fname.replace(".root","_profiles.root")
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

        multiClusters = event.multiClusters()
        #tot_multiclus += len(multiClusters)

        # 2d clusters in layers
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)

        # rechits
        recHits = event.recHits()
        tot_rechit += len(recHits)

        for i_mcl, multicl in enumerate(multiClusters):
            #if multicl.energy() < minE: continue
            if multicl.pt() < min_pt: continue
            #if len(multicl.cluster2d()) < 3: continue
            if multicl.NLay() < 3: continue

            #if abs(multicl.eta()) > 2.0: continue

            found_part = False
            for part in genParts:

                #if part.eta() * z_half < 0: continue
                if part.gen() < 1: continue
                #if part.gen() > 0: continue
                if part.reachedEE() < 2: continue
                #if part.fbrem() < min_fbrem: continue
                if part.fbrem() > max_fbrem: continue

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

            tot_multiclus += 1

            #addDataPoint(hist_data,"part_mcl_a_v",a_v)

            #if abs(a_v) < 0.03: continue

            #addDataPoint(hist_data,"part_mcl_dR",dR)
            #addDataPoint(hist_data,"mcl_axisZ",multicl.pcaAxisZ())

            ## play with cluster layers
            clusters = [layerClusters[clust_idx] for clust_idx in multicl.cluster2d()]

            #gr_rh_XYZ = ROOT.TGraph2D(); gr_rh_XYZ.SetMarkerStyle(20)

            ## Rechit plots
            rechits = []
            for cluster in clusters:
                rechits += [recHits[rh_idx] for rh_idx in cluster.rechits()]# if recHits[rh_idx].flags() == 0]

            ## sort indeces by z
            rh_ind_sort_z = sorted(xrange(len(rechits)), key=lambda x: abs(rechits[x].z()))
            z_positions = sorted(set([abs(int(rh.z()*10)/10.) for rh in rechits]))

            '''
            print "first and last z:" , rechits[rh_ind_sort_z[0]].z(), rechits[rh_ind_sort_z[-1]].z()
            print "z positions", z_positions
            '''

            lay_energies = np.array([0.]*52)
            #z_energies = [0]*len(z_positions)
            z_energies = np.array([0.]*len(z_positions))
            #z_energies = {z_pos:0 for z_pos in z_positions}

            for irh in rh_ind_sort_z:
                rh = rechits[irh]

                #if rh.layer() > 28: continue
                lay_energies[rh.layer()] += rh.energy()

                z_pos = abs(int(rh.z()*10)/10.)
                z_indx = z_positions.index(z_pos)
                z_energies[z_indx] += rh.energy()

            # normalize by multiclust energy and plot
            lay_energies[:] = lay_energies/multicl.energy()
            lay_cumene = np.cumsum(lay_energies)

            z_energies[:] = z_energies/multicl.energy()
            z_cumene = np.cumsum(z_energies)

            #print z_energies
            #print lay_energies
            #print lay_cumene

            ## calculate inclination of multicluster to z_axis
            ## 0. vector normal for XY plane
            norm_vect = ROOT.TVector3(0,0,abs(part.eta())/part.eta())
            ## 1. get mulcutluster vector
            mclut_vect = ROOT.TVector3(
                multicl.pcaAxisX(),
                multicl.pcaAxisY(),
                multicl.pcaAxisZ(),
            )
            # 1. multicluster center
            mcl_cent = (multicl.pcaPosX(), multicl.pcaPosY(), multicl.pcaPosZ())
            a_p,a_v = calc_angles(mcl_cent,mclut_vect,norm_vect)
            a_comb = math.hypot(a_p,a_v)
            a_norm = norm_vect.Angle(mclut_vect)

            cosThets_mcl = math.cos(a_norm)
            #print cosThets_mcl, abs(mcl_tlv.CosTheta())
            addDataPoint(hist_data,"cos_vs",(cosThets_mcl, abs(mcl_tlv.CosTheta())))

            first_z = abs(rechits[rh_ind_sort_z[0]].z())
            ## plot cumulative 50% :
            z_cumul50 = -1

            # add points to plot
            for lay,ene in enumerate(lay_energies):
                addDataPoint(hist_data,"ene_prof_lay",(lay,ene))
                addDataPoint(hist_data,"ene_cumul_lay",(lay,lay_cumene[lay]))

            for z_indx,z_pos in enumerate(z_positions):
                #z_pos = z_positions[z_indx]
                # correct z position for cosTheta
                #z_pos /= abs(mcl_tlv.CosTheta())
                #z_pos /= abs(part_tlv.CosTheta())
                #z_pos = (z_pos - first_z) /abs(mcl_tlv.CosTheta()) # / cosThets_mcl

                z_pos = (z_pos - 320.75500) / cosThets_mcl

                ene = z_energies[z_indx]
                cumul = z_cumene[z_indx]
                addDataPoint(hist_data,"ene_prof_z",(z_pos,ene))
                addDataPoint(hist_data,"ene_cumul_z",(z_pos,cumul))

                if z_cumul50 == -1 and cumul > 0.5: z_cumul50 = z_pos

            addDataPoint(hist_data,"z_cumul50",z_cumul50)

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
