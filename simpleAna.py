#!/usr/bin/env python
import ROOT, math
import numpy as np
from NtupleDataFormat import HGCalNtuple
from helperTools import *
# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

max_events = 10

z_half = +1
minE = .5
min_pt = 2

min_fbrem = 0.9
#max_dR = 0.3
max_dR = 0.3

def main():
    #ntuple = HGCalNtuple("/Users/clange/CERNBox/partGun_PDGid211_x120_E80.0To80.0_NTUP_9.root")
    ntuple = HGCalNtuple("hgcalNtuple-pca-1000.root")
    #ntuple = HGCalNtuple("../hgcalNtuple-pca.root")

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

    h_mclust_dR = ROOT.TH1F("h_mclust_dR","multi dR; dR(p,mclust)",100,0,3.2)

    h_clust_dEta = ROOT.TH1F("h_clust_dEta","h_clust_dEta;dEta(mclust,clust)",100,-0.1,0.1)
    h_clust_dPhi = ROOT.TH1F("h_clust_dPhi","h_clust_dPhi;dPhi(mclust,clust)",100,-0.1,0.1)
    h_clust_dR = ROOT.TH1F("h_clust_dR","h_clust_dR;dR(mclust,clust)",100,0,0.2)

    h_mclust_dR_dZ = ROOT.TH2F("h_mclust_dR_dZ","multi dR; dR(p,mclust);dZ (clusts)",100,0,3.2,100,0,20)
    #h_mclust_dXdY = ROOT.TH2F("h_mclust_dXdY","dX/dY; layer; dX/dY",28,0,28,100,0,30)
    h_mclust_dX = ROOT.TH2F("h_mclust_dX","dX; layer; dX",28,1,29,100,0,10)
    h_mclust_dY = ROOT.TH2F("h_mclust_dY","dY; layer; dY",28,1,29,100,0,10)

    h_mclust_str_E = ROOT.TH2F("h_mclust_str_E","start position; layer; E thr",28,1,29,20,0,2)

    #h_Event = ROOT.TH3F("h_Event","event; x; y",100,-80,80,100,-80,80,30,300,330)
    ## data storages for hists
    mgr_Event_X = ROOT.TMultiGraph("ev_x","X good")
    mgr_Event_Y = ROOT.TMultiGraph("ev_y","Y good")
    mgr_Event2_X = ROOT.TMultiGraph("ev2_x","X bad")
    mgr_Event2_Y = ROOT.TMultiGraph("ev2_y","Y bad")

    gr_Event = ROOT.TGraph2D(); gr_Event.SetTitle("good axisZ"); gr_Event.SetMarkerStyle(20)
    gr_Event2 = ROOT.TGraph2D(); gr_Event2.SetTitle("bad axisZ"); gr_Event2.SetMarkerStyle(20)
    gr_cnt = 0
    gr_cnt2 = 0

    gr2ds = []
    gr2d_2s = []

    hist_data = {}

    col_cnt = 0

    stop_run = False

    for event in ntuple:
        #if stop_run: break
        if tot_nevents >= max_events: break

        if tot_nevents % 100 == 0: print("Event %i" % tot_nevents)
        # print "Event", event.entry()
        tot_nevents += 1

        genParts = event.genParticles()
        #tot_genpart += len(genParts)

        found_part = False
        for part in genParts:

            if part.eta() * z_half < 0: continue
            if part.gen() < 1: continue
            if not part.reachedEE(): continue
            if part.fbrem() < min_fbrem: continue

            found_part = True
            tot_genpart += 1
            break

        if not found_part: continue

        multiClusters = event.multiClusters()
        #tot_multiclus += len(multiClusters)

        # 2d clusters in layers
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)

        # rechits
        recHits = event.recHits()
        tot_rechit += len(recHits)

        #print part.posx()[0],part.posy()[0],part.posz()[0]

        for i_mcl, multicl in enumerate(multiClusters):
            if multicl.energy() < minE: continue
            if multicl.pt() < min_pt: continue
            if multicl.eta() * part.eta() < 0: continue

            #if abs(multicl.pcaAxisZ()) > 0.4: continue

            if len(multicl.cluster2d()) < 3: continue

            #if abs(multicl.pcaAxisZ()) < 0.4: continue

            dR =  math.hypot(multicl.eta()-part.eta(), multicl.phi()-part.phi())

            '''
            h_mclust_dR.Fill(dR)
            addDataPoint(hist_data,"clust_dR",dR)

            #drho = math.hypot(multicl.slopeX()-part.posx()[0],multicl.slopeY()-part.posy()[0])
            drho = calcDeltaRho(part,multicl)
            #print drho,drho2

            addDataPoint(hist_data,"clust_drho",drho)
            addDataPoint(hist_data,"clust_dR_drho",(dR,drho))

            addDataPoint(hist_data,"multi_dR_sigV",(dR,multicl.sigvv()))
            addDataPoint(hist_data,"multi_dR_sigU",(dR,multicl.siguu()))
            '''

            if dR > max_dR: continue
            tot_multiclus += 1

            addDataPoint(hist_data,"mcl_axisZ",multicl.pcaAxisZ())
            #if abs(multicl.pcaAxisZ()) > 0.4: continue

            ## play with cluster layers
            clusters = [layerClusters[clust_idx] for clust_idx in multicl.cluster2d()]
            #print len(clusters), " clusters in multi"

            # plot event display
            #gr_Event = ROOT.TGraph2D()
            gr_Event_X = ROOT.TGraph()
            gr_Event_Y = ROOT.TGraph()
            gr2d = ROOT.TGraph2D()
            gr2d_2 = ROOT.TGraph2D()

            #col_frac = multicl.pt()/part.pt()
            #col_frac = multicl.pt()/multiClusters[0].pt()
            #col_frac = i_mcl/len(multiClusters)
            col_frac = col_cnt/60.
            col = ROOT.gStyle.GetColorPalette(int(col_frac * 255))
            col_cnt += 10
            gr_Event_X.SetMarkerColor(col)
            gr_Event_X.SetMarkerStyle(20)
            gr_Event_Y.SetMarkerColor(col)
            gr_Event_Y.SetMarkerStyle(20)

            gr2d.SetMarkerColor(col)
            gr2d.SetMarkerStyle(20)
            gr2d_2.SetMarkerColor(col)
            gr2d_2.SetMarkerStyle(20)

            rh_cnt = 0

            for i,cluster in enumerate(clusters):
                for j,rh_idx in enumerate(cluster.rechits()):
                    rechit = recHits[rh_idx]

                    #if rechit.energy() < 0.1: continue

                    #ind = i * len(cluster.rechits()) + j
                    gr_Event_Y.SetPoint(rh_cnt,rechit.x(),rechit.z())
                    gr_Event_X.SetPoint(rh_cnt,rechit.x(),rechit.y())

                    if abs(multicl.pcaAxisZ()) > 0.5:
                        gr_Event.SetPoint(gr_cnt,rechit.z(),rechit.x(),rechit.y())
                        gr2d.SetPoint(rh_cnt,rechit.z(),rechit.x(),rechit.y())
                        gr_cnt+=1
                    else:
                        gr_Event2.SetPoint(gr_cnt2,rechit.z(),rechit.x(),rechit.y())
                        gr2d_2.SetPoint(rh_cnt,rechit.z(),rechit.x(),rechit.y())
                        gr_cnt2+=1

                    rh_cnt += 1

            if abs(multicl.pcaAxisZ()) > 0.5:
                mgr_Event_X.Add(gr_Event_X)
                mgr_Event_Y.Add(gr_Event_Y)
                gr2ds.append(gr2d)
            else:
                print event.event(), event.run(), i_mcl
                mgr_Event2_X.Add(gr_Event_X)
                mgr_Event2_Y.Add(gr_Event_Y)
                gr2d_2s.append(gr2d_2)


            stop_run = True
            '''
            for i,cluster in enumerate(clusters):
                #gr_Event.SetPoint(i+1,cluster.x(),cluster.y(),cluster.z())
                gr_Event_Y.SetPoint(i,cluster.y(),cluster.z())
                gr_Event_X.SetPoint(i,cluster.x(),cluster.z())

                if abs(multicl.pcaAxisZ()) > 0.5:
                    gr_Event.SetPoint(gr_cnt,cluster.z(),cluster.x(),cluster.y())
                    gr_cnt+=1
                else:
                    gr_Event2.SetPoint(gr_cnt2,cluster.z(),cluster.x(),cluster.y())
                    gr_cnt2+=1

            if abs(multicl.pcaAxisZ()) > 0.5:
                mgr_Event_X.Add(gr_Event_X)
                mgr_Event_Y.Add(gr_Event_Y)
            else:
                mgr_Event2_X.Add(gr_Event_X)
                mgr_Event2_Y.Add(gr_Event_Y)
            '''
            break

            #print compDeltaVar(clusters,'z',"mean"),
            dZ = compDeltaVar(clusters,'z',"std")
            addDataPoint(hist_data,"clust_dR_dZ",(dR,dZ))
            h_mclust_dR_dZ.Fill(dR,dZ)

            dX = compDeltaVar(clusters,'x',"std")
            addDataPoint(hist_data,"clust_dR_dX",(dR,dX))
            dY = compDeltaVar(clusters,'y',"std")
            addDataPoint(hist_data,"clust_dR_dY",(dR,dY))
            addDataPoint(hist_data,"clust_dX_dY",(dX,dY))


            # sigma EtaEta
            dEta = compDeltaVar(clusters,'eta',"std")
            addDataPoint(hist_data,"clust_dEta_sigU",(dEta,multicl.siguu()))
            addDataPoint(hist_data,"clust_dEta_sigV",(dEta,multicl.sigvv()))

            sigmaEtaEta = 0
            #for rh in [recHits[rh_idx] for rh_idx in cluster.rechits()

            #continue
            #if len(clusters) < 25: continue

            cluster_layers = [[] for layer in range(28)]
            rechits_layers = [[] for layer in range(28)]
            #cluster_layers = [28*[]]
            #rechits_layers = [28*[]]

            # count number of layers
            n_lay = len(set([cluster.layer() for cluster in clusters if cluster.layer() < 29]))
            #print n_lay
            addDataPoint(hist_data,"clust_Zax_Nlay",(multicl.pcaAxisZ(),n_lay))

            ##
            # Fill layer-wise cluster/rechit collections
            ##
            for cluster in clusters:
                if cluster.layer() > 28: continue
                #if cluster.energy() < minE/100.: continue
                cluster_layers[cluster.layer()-1].append(cluster)
                # rechits in layer
                rechits_layers[cluster.layer()-1]+=[recHits[rh_idx] for rh_idx in cluster.rechits()]

            ### determine start position
            for e_thr_step in range(21):
                e_thr = e_thr_step * 0.1

                for lay,clusts in enumerate(cluster_layers):
                    if len(clusts) == 0: continue

                    layer_e = sum([clust.energy() for clust in clusts])
                    if layer_e > e_thr:
                        #addDataPoint(hist_data,"shower_start_thrE",(e_thr,lay))
                        addDataPoint(hist_data,"shower_start_thrE",(lay+1,e_thr))
                        h_mclust_str_E.Fill(lay+1,e_thr)
                        break


            for lay,rhits in enumerate(rechits_layers):
                if len(rhits) == 0: continue
                dX = compDeltaVar(rhits,'x',"std")
                dY = compDeltaVar(rhits,'y',"std")

                h_mclust_dX.Fill(lay+1,dX)
                h_mclust_dY.Fill(lay+1,dY)
            '''
            for lay,clusts in enumerate(cluster_layers):
                if len(clusts) == 0: continue
                dX = compDeltaVar(clusts,'x',"std")
                dY = compDeltaVar(clusts,'y',"std")

                h_mclust_dX.Fill(lay,dX)
                h_mclust_dY.Fill(lay,dY)


                print("%i - %i , dX = %f" %(lay, len(clusts), dX))

                #for prop in ["std","mean"]:
                #    print prop
                #    for var in ["x","y","z"]:
                #        print var, compDeltaVar(clusts, var, prop),

            #print
            '''

        # for genPart in genParts:
        #     print tot_nevents, "genPart pt:", genPart.pt()

    print "Processed %d events" % tot_nevents
    print "On average %f generator particles" % (float(tot_genpart) / tot_nevents)
    print "On average %f reconstructed hits" % (float(tot_rechit) / tot_nevents)
    print "On average %f raw reconstructed hits" % (float(tot_rechit_raw) / tot_nevents)
    print "On average %f layer clusters" % (float(tot_cluster2d) / tot_nevents)
    print "On average %f multi-clusters" % (float(tot_multiclus) / tot_nevents)
    print "On average %f sim-clusters" % (float(tot_simcluster) / tot_nevents)
    print "On average %f PF clusters" % (float(tot_pfcluster) / tot_nevents)
    print "On average %f calo particles" % (float(tot_calopart) / tot_nevents)
    print "On average %f tracks" % (float(tot_track) / tot_nevents)

    canv_Event = ROOT.TCanvas("canv_Event","Event",1000,800)
    canv_Event.Divide(2,2)
    canv_Event.cd(1);    mgr_Event_X.Draw("Ap")
    canv_Event.cd(2);    mgr_Event_Y.Draw("Ap")
    canv_Event.cd(3);    mgr_Event2_X.Draw("Ap")
    canv_Event.cd(4);    mgr_Event2_Y.Draw("Ap")
    canv_Event.Update()

    canv_Event2 = ROOT.TCanvas("canv_Event2","Event",1000,800)
    canv_Event2.Divide(2,1)
    canv_Event2.cd(1);
    gr_Event.Draw("p")
    #gr2ds[0].Draw("p")
    for gr in gr2ds:
        gr.Draw("psame")
    canv_Event2.cd(2);
    gr_Event2.Draw("pcol")
    for gr in gr2d_2s:
        gr.Draw("psame")
    canv_Event2.Update()

    '''
    canv_dR = ROOT.TCanvas("canv_dR","canv",800,600)
    h_mclust_dR_dZ.Draw("colz")
    canv_dR.Update()

    canv_dXY = ROOT.TCanvas("canv_dXY","canv",1000,600)
    canv_dXY.Divide(2,1)

    canv_dXY.cd(1)
    h_mclust_dX.Draw("colz")
    canv_dXY.cd(2)
    h_mclust_dY.Draw("colz")

    canv_dXY.Update()

    canv_LayE = ROOT.TCanvas("canv_LayE","canv",1000,600)
    h_mclust_str_E.Draw("colz")
    canv_LayE.Update()
    '''

    #print hist_data
    hists = []
    for data_name in hist_data:
        print("Plotting hist for data: %s" %data_name)
        canv = ROOT.TCanvas("canv_" + data_name,data_name,800,600)
        hist = getHisto(hist_data[data_name],data_name,data_name)
        hist.Draw("colz")
        canv.Update()
        #hists.append(hist)
        ROOT.SetOwnership(canv,0)

    q = raw_input("exit")

if __name__ == "__main__":
    main()
