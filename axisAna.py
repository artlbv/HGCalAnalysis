#!/usr/bin/env python
import sys
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
max_fbrem = 0.3
#max_fbrem = 1.0
#max_dR = 0.3
max_dR = 0.3

def main(fname = "hgcalNtuple-El15-100_noReClust.root"):
    ntuple = HGCalNtuple(fname)

    foutname = fname.replace(".root","_axis.root")
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

        for i_mcl, multicl in enumerate(multiClusters):

            if multicl.energy() < 1: continue
            #if multicl.pt() < min_pt: continue
            #if multicl.pt() < 10: continue
            if multicl.pt() < 5: continue
            #if multicl.z() * part.eta() < 0: continue

            #if len(multicl.cluster2d()) < 3: continue
            if multicl.NLay() < 3: continue

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

            #if abs(multicl.pcaAxisZ() < 0.1): continue

            #######
            # IP-vector based stuff
            #######
            part_tlv = ROOT.TLorentzVector()
            part_tlv.SetPtEtaPhiE(part.pt(), part.eta(), part.phi(), part.energy())

            mcl_tlv = ROOT.TLorentzVector()
            mcl_tlv.SetPtEtaPhiE(multicl.pt(), multicl.eta(), multicl.phi(), multicl.energy())

            dR = part_tlv.DeltaR(mcl_tlv)

            #if dR < 5:
            if dR > 1: continue

            addDataPoint(hist_data,"part_mcl_dR",dR)
            #addDataPoint(hist_data,"part_mcl_dR_fbrem",(dR,part.fbrem()))

            #######
            # local vector based stuff
            #######

            ## 0. vector normal for XY plane
            norm_vect = ROOT.TVector3(0,0,abs(part.eta())/part.eta())

            ## 1. get particle propagation vector
            part_vect = ROOT.TVector3(
                part.posx()[10]-part.posx()[0],
                part.posy()[10]-part.posy()[0],
                part.posz()[10]-part.posz()[0],
            )

            ## 2. get mulcutluster vector
            mclut_vect = ROOT.TVector3(
                multicl.pcaAxisX(),
                multicl.pcaAxisY(),
                multicl.pcaAxisZ(),
            )

            #if abs(mclut_vect.Mag()-1) > 0.1: continue
            angle = part_vect.Angle(mclut_vect)

            ## Calculate veritcal/perp angles
            # 1. multicluster center
            mcl_cent = (multicl.pcaPosX(), multicl.pcaPosY(), multicl.pcaPosZ())
            a_p,a_v = calc_angles(mcl_cent,mclut_vect,part_vect)

            an_p,an_v = calc_angles(mcl_cent,mclut_vect,norm_vect)

            if abs(a_p) > 0.1 or abs(a_v) > 0.1: continue

            tot_multiclus += 1

            '''
            addDataPoint(hist_data,"part_mcl_angle",angle)

            addDataPoint(hist_data,"part_mcl_a_perp",a_p)
            addDataPoint(hist_data,"part_mcl_a_vert",a_v)
            addDataPoint(hist_data,"part_mcl_avert_vs_aperp",(a_p,a_v))

            addDataPoint(hist_data,"mcl_an_perp",an_p)
            addDataPoint(hist_data,"mcl_an_vert",an_v)
            addDataPoint(hist_data,"mcl_anvert_vs_anperp",(an_p,an_v))
            '''

            addDataPoint(hist_data,"mcl_an_perp_vs_eta",(abs(part.eta()),an_p))
            addDataPoint(hist_data,"mcl_an_vert_vs_eta",(abs(part.eta()),an_v))

            ## entry position
            z_front = 320.75500 * part.eta()/abs(part.eta())
            x_entry, y_entry = get_entry_point(mcl_cent,mclut_vect,z_front)

            #addDataPoint(hist_data,"mcl_entry_x",x_entry)
            #addDataPoint(hist_data,"mcl_entry_y",y_entry)
            addDataPoint(hist_data,"mcl_entry_x_vs_y",(x_entry,y_entry))

            entry_vect = ROOT.TVector3(x_entry,y_entry,z_front)
            ## entry_vector Eta/Phi,dR
            dR_entr_part = part_vect.DeltaR(entry_vect)
            dR_entr_mcl = mclut_vect.DeltaR(entry_vect)

            addDataPoint(hist_data,"part_entr_dR", dR_entr_part)
            addDataPoint(hist_data,"mcl_entr_dR", dR_entr_mcl)

            addDataPoint(hist_data,"part_mcl_dR_vs_eta", (abs(part.eta()),dR))
            addDataPoint(hist_data,"part_entr_dR_vs_eta", (abs(part.eta()),dR_entr_part))
            addDataPoint(hist_data,"mcl_entr_dR_vs_eta", (abs(part.eta()),dR_entr_mcl))

            addDataPoint(hist_data,"mcl_entr_dR_vs_part_mcl_dR", (dR,dR_entr_mcl))

            ## propagate particle to cluster centre
            #x_centr, y_centr = get_entry_point(mcl_cent,part_vect,mcl_cent[2])
            prop_entr = [
                part.posx()[0],
                part.posy()[0],
                part.posz()[0],
            ]

            x_centr, y_centr = get_entry_point(prop_entr,part_vect,mcl_cent[2])
            #addDataPoint(hist_data,"part_centr_x_vs_y",(x_centr,y_centr))

            drho_cent = math.hypot(x_centr-mcl_cent[0], y_centr-mcl_cent[1])
            addDataPoint(hist_data,"centr_drho",drho_cent)

            drho_entr = math.hypot(x_entry-prop_entr[0], y_entry-prop_entr[1])
            addDataPoint(hist_data,"entr_drho",drho_entr)

            addDataPoint(hist_data,"centr_vs_entr_drho",(drho_cent,drho_entr))

            centr_vect = ROOT.TVector3(x_centr,y_centr,mcl_cent[2])
            dR_cent_mcl = mclut_vect.DeltaR(centr_vect)

            addDataPoint(hist_data,"mcl_cent_dR", dR_cent_mcl)

            '''
            a_comb = math.hypot(a_p,a_v)
            addDataPoint(hist_data,"part_mcl_a_combined",a_comb)

            addDataPoint(hist_data,"part_mcl_a_perp_pcaZ",(a_p,abs(multicl.pcaAxisZ())))
            addDataPoint(hist_data,"part_mcl_a_vert_pcaZ",(a_v,abs(multicl.pcaAxisZ())))

            addDataPoint(hist_data,"part_mcl_pt_eta",(multicl.pt(),abs(multicl.eta())))
            addDataPoint(hist_data,"part_mcl_pt_peta",(multicl.pt(),abs(part.eta())))

            addDataPoint(hist_data,"part_mcl_a_perp_eta",(a_p,abs(multicl.eta())))
            addDataPoint(hist_data,"part_mcl_a_vert_eta",(a_v,abs(multicl.eta())))

            addDataPoint(hist_data,"part_mcl_a_perp_peta",(a_p,abs(part.eta())))
            addDataPoint(hist_data,"part_mcl_a_vert_peta",(a_v,abs(part.eta())))

            addDataPoint(hist_data,"part_mcl_a_perp_ene",(a_p,multicl.energy()))
            addDataPoint(hist_data,"part_mcl_a_vert_ene",(a_v,multicl.energy()))

            addDataPoint(hist_data,"part_mcl_a_p_vs_a_norm",(a_p,norm_vect.Angle(part_vect)))
            addDataPoint(hist_data,"part_mcl_a_v_vs_a_norm",(a_v,norm_vect.Angle(part_vect)))

            addDataPoint(hist_data,"part_mcl_a_perp_fbrem",(a_p,part.fbrem()))
            addDataPoint(hist_data,"part_mcl_a_vert_fbrem",(a_v,part.fbrem()))

            addDataPoint(hist_data,"part_mcl_a_perp_sigu",(a_p,multicl.siguu()))
            addDataPoint(hist_data,"part_mcl_a_perp_sigv",(a_p,multicl.sigvv()))

            addDataPoint(hist_data,"part_mcl_a_vert_sigu",(a_v,multicl.siguu()))
            addDataPoint(hist_data,"part_mcl_a_vert_sigv",(a_v,multicl.sigvv()))
            '''


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
        elif "TGraph2D" in hist.ClassName():
            hist.Draw("pcolz")
        elif "TGraph" in hist.ClassName():
            hist.Draw("ap")
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
