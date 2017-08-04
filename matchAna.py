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

def calc_angles_0(point1 = (1,2,3), vect1 = (2,3,4), vect2 = (4,5,6)):
    v1 = ROOT.TVector3(vect1[0],vect1[1],vect1[2])
    v2 = ROOT.TVector3(vect2[0],vect2[1],vect2[2])

    ## 1. get radial vector from point
    v_rad = ROOT.TVector3(point1[0],point1[1],0)

    ## 2. perpendicular vector to r and vect1
    v_perp = v1.Cross(v_rad)

    ## 3. vertical vector to perp and vect1
    v_vert = v_perp.Cross(v1)

    ## get angles
    a_perp = v_perp.Angle(v2)
    a_vert = v_vert.Angle(v2)

    print a_perp,a_vert

def calc_angles(point,v1,v2):

    ## 1. get radial vector from point
    v_rad = ROOT.TVector3(point[0],point[1],0)

    ## 2. perpendicular vector to r and vect1
    v_perp = v1.Cross(v_rad)

    ## 3. vertical vector to perp and vect1
    v_vert = v_perp.Cross(v1)

    ## get angles
    a_perp = ROOT.TMath.Pi()/2 - v_perp.Angle(v2)
    a_vert = ROOT.TMath.Pi()/2 - v_vert.Angle(v2)

    return a_perp,a_vert

def get_angles(multicl,part):

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


    ## Calculate veritcal/perp angles
    # 1. multicluster center
    mcl_cent = (multicl.pcaPosX(), multicl.pcaPosY(), multicl.pcaPosZ())
    a_p,a_v = calc_angles(mcl_cent,mclut_vect,part_vect)

    return a_p,a_v


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

        '''
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
        '''
        #continue

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
                #if part.fbrem() > max_fbrem: continue

                if multicl.z() * part.eta() < 0: continue

                found_part = True
                if i_mcl < 2: tot_genpart += 1

                break

            if not found_part: continue


            #######
            # IP-vector based stuff
            #######
            part_tlv = ROOT.TLorentzVector()
            part_tlv.SetPtEtaPhiE(part.pt(), part.eta(), part.phi(), part.energy())

            mcl_tlv = ROOT.TLorentzVector()
            mcl_tlv.SetPtEtaPhiE(multicl.pt(), multicl.eta(), multicl.phi(), multicl.energy())

            dR = float(part_tlv.DeltaR(mcl_tlv))
            #dR2 =  math.hypot(multicl.eta()-part.eta(), multicl.phi()-part.phi())
            #print dR, dR2

            #if dR < 5:
            if dR > 1: continue

            addDataPoint(hist_data,"part_mcl_dR",dR)
            #addDataPoint(hist_data,"part_mcl_dR_fbrem",(dR,part.fbrem()))

            #######
            # real vector based stuff
            #######
            #if abs(multicl.pcaAxisZ() < 0.1): continue

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

            if abs(a_p) > 0.1 or abs(a_v) > 0.1: continue

            tot_multiclus += 1

            addDataPoint(hist_data,"part_mcl_angle",angle)

            addDataPoint(hist_data,"part_mcl_a_perp",a_p)
            addDataPoint(hist_data,"part_mcl_a_vert",a_v)
            addDataPoint(hist_data,"part_mcl_avert_vs_aperp",(a_p,a_v))

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

            continue

            ## construct "dR" from dR and angle
            dR2 = math.hypot(angle,dR)
            addDataPoint(hist_data,"part_mcl_dR2",dR2)

            if angle > 0.3: continue

            addDataPoint(hist_data,"part_mcl_angle_vs_dR",(dR,angle))

            if dR > 0.06: continue

            addDataPoint(hist_data,"part_mcl_angle",angle)

            #print part_vect.Mag(), mclut_vect.Mag()
            #print part_vect.Angle(mclut_vect), mclut_vect.Angle(part_vect)
            #print multicl.pcaAxisX(),multicl.pcaAxisY(),multicl.pcaAxisZ()

            addDataPoint(hist_data,"part_norm_angle",norm_vect.Angle(part_vect))
            addDataPoint(hist_data,"mcl_norm_angle",norm_vect.Angle(mclut_vect))

            addDataPoint(hist_data,"part_mcl_angle",angle)

            '''
            addDataPoint(hist_data,"part_norm_angle_eta",(part.eta(),norm_vect.Angle(part_vect)))
            addDataPoint(hist_data,"mcl_norm_angle_eta",(multicl.eta(),norm_vect.Angle(mclut_vect)))
            addDataPoint(hist_data,"angle_vs_eta",(part.eta(),angle))

            addDataPoint(hist_data,"part_mcl_angle_pt",(angle,multicl.pt()))
            addDataPoint(hist_data,"part_mcl_angle_pt_br",(angle,multicl.pt(),part.fbrem()))

            addDataPoint(hist_data,"part_mcl_angle_pt_br",(angle,multicl.pt(),part.fbrem()))
            addDataPoint(hist_data,"part_mcl_angle_br",(angle,part.fbrem()))
            addDataPoint(hist_data,"part_mcl_dR_br",(dR,part.fbrem()))

            addDataPoint(hist_data,"dR_vs_angle",(angle,dR))
            #addDataPoint(hist_data,"dR_vs_angle_vs_pt",(angle,dR,multicl.pt()))
            #addDataPoint(hist_data,"dR_vs_angle_vs_ene",(angle,dR,multicl.energy()))

            addDataPoint(hist_data,"siguu_vs_eta",(part.eta(),multicl.siguu()))
            addDataPoint(hist_data,"sigvv_vs_eta",(part.eta(),multicl.sigvv()))

            addDataPoint(hist_data,"siguu_vs_angle",(angle,multicl.siguu()))
            addDataPoint(hist_data,"sigvv_vs_angle",(angle,multicl.sigvv()))
            '''
            addDataPoint(hist_data,"siguu_vs_angle",(angle,multicl.siguu()))
            addDataPoint(hist_data,"sigvv_vs_angle",(angle,multicl.sigvv()))

            #addDataPoint(hist_data,"hbr_siguu_vs_angle_vs_fbrem",(angle,multicl.siguu(),part.fbrem()))
            #addDataPoint(hist_data,"hbr_sigvv_vs_angle_vs_fbrem",(angle,multicl.sigvv(),part.fbrem()))
            addDataPoint(hist_data,"siguu_vs_fbrem",(multicl.siguu(),part.fbrem()))
            addDataPoint(hist_data,"sigvv_vs_fbrem",(multicl.sigvv(),part.fbrem()))
            addDataPoint(hist_data,"angle_vs_fbrem",(angle,part.fbrem()))

            addDataPoint(hist_data,"sigvv_vs_fbrem_vs_pt",(multicl.sigvv(),part.fbrem(),multicl.pt()))

            addDataPoint(hist_data,"sigpp_vs_fbrem",(multicl.sigpp(),part.fbrem()))
            addDataPoint(hist_data,"sigee_vs_fbrem",(multicl.sigee(),part.fbrem()))

            continue

            if part.fbrem() > max_fbrem:
                addDataPoint(hist_data,"hbr_siguu_vs_angle",(angle,multicl.siguu()))
                addDataPoint(hist_data,"hbr_sigvv_vs_angle",(angle,multicl.sigvv()))
            else:
                addDataPoint(hist_data,"lbr_siguu_vs_angle",(angle,multicl.siguu()))
                addDataPoint(hist_data,"lbr_sigvv_vs_angle",(angle,multicl.sigvv()))


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
