#!/usr/bin/env python
import os, sys
import ROOT
from array import array
import phase2tdrStyle
#import CMS_lumi, tdrstyle
import numpy as np

#ROOT.gROOT.LoadMacro("plot_style/HttStyles.cc")
#ROOT.gROOT.LoadMacro("plot_style/CMS_lumi.C")
#ROOT.setTDRStyle()

labels = {}

labels["ele_pT"] = "p_{T} (GeV)"
labels["ele_eta"] = "|#eta|"
labels["ele_scleta"] = "|#eta|"
labels["Nvtx"] = "Number of reco. vertices"
#labels["NPU"] = "Number of gen. vertices"
labels["NPU"] = "Number of PU interactions"
labels["PU_density"] = "Density (events / mm)"

def getScores(fname):

    treename = "TestTree"
    tfile = ROOT.TFile(fname)
    tree = tfile.Get(treename)
    print "Found tree " + treename + " with %i events" % tree.GetEntries()

    print "Reading scores...",

    sig_scores = []
    bkg_scores = []

    for i_entry,event in enumerate(tree):

        #if i_entry > 100: continue
        #if abs(event.ele_deltaetaele) > 0.02: continue
        #if (event.ele_sigmauu > 0.7 and event.ele_sigmavv > 0.8) : continue
        #if event.ele_sigmavv > 0.8 : continue
        #if event.ele_hgc_nlay < 20: continue
        #if event.ele_ET < 5: continue

        #if event.classID == 0: sig_scores.append(event.BDT)
        if event.classID == 0:
            #if abs(event.ele_deltaetaele) > 0.02: continue
            #if event.ele_ET < 5: continue
            sig_scores.append(event.BDT)
        elif event.classID == 1:
            bkg_scores.append(event.BDT)

    print "done"

    tfile.Close()

    #print "Read scores"
    print "Recorded events:", len(sig_scores), len(bkg_scores)

    sig_scores = np.array(sig_scores)
    bkg_scores = np.array(bkg_scores)

    print "Signal mean BDT score: ", np.mean(sig_scores)
    print "Signal WP80% BDT score: ", np.percentile(sig_scores,100-80)
    print "Signal WP90% BDT score: ", np.percentile(sig_scores,100-90)
    print "Background mean BDT score: ", np.mean(bkg_scores)

    return sig_scores, bkg_scores

def getWPscore(scores, eff = 50):

    wp_score =  np.percentile(scores,100-eff)
    print "Score for %0.2f efficiency is %0.3f" %(eff,wp_score)

    return wp_score

def getROC(sig_scores, bkg_scores):

    roc_data = []

    min_perc = 50
    max_perc = 100
    for i_eff in range(min_perc,max_perc,1):

        sig_eff = i_eff / 1.
        sig_rej = 100 - sig_eff
        #print i_eff, sig_eff, sig_rej
        bdt_score = np.percentile(sig_scores,sig_rej)

        #bkg_eff = len([score for score in ])/len(bkg_scores)
        #bkg_eff = len([score for score in bkg_scores if score > bdt_score])/float(len(bkg_scores))
        bkg_eff = sum(bkg_scores > bdt_score)/float(len(bkg_scores))

        #print i_eff, bdt_score, sig_eff, bkg_eff
        roc_data.append((bdt_score, sig_eff, bkg_eff))

    return roc_data

def plotROC(roc_data, grname = "gr"):

    #gr = ROOT.TGraph(len(roc_data))
    gr = ROOT.TGraph()
    gr.SetName(grname)

    print "Plotting ROC"
    for i,point in enumerate(roc_data):

        score,sig_eff,bkg_eff = point
        #print point
        gr.SetPoint(i,sig_eff,1-bkg_eff)

    print "Filled ROC graph with %i points" %(gr.GetN())

    cname = "c" + grname
    canv = phase2tdrStyle.setCanvas()
    #canv = getCanvas(cname)
    gr.Draw("apl")
    #ROOT.SetDirectory(gr,0)
    #gr.SetOwnership(0)

    phase2tdrStyle.drawCMS(True)
    phase2tdrStyle.drawEnPu()

    canv.Update()
    canv.Draw()
    ROOT.SetOwnership(canv,0)

    canv.SaveAs("ROC.pdf")

    q = raw_input("q")

    return gr#canv

def getHistFromTree(tree, var = "ele_pT", cuts = "", hname = ""):

    hname_pref = 'h_' + var + '_' + hname
    #histList = [] # going to store sig_ref, sig_sel, bkg_ref, bkg_sel

    if "ele_pT" in var:
        #hRef = ROOT.TH1F(hname_pref + "_ref", "", 20, 10, 20)
        #hRef = ROOT.TH1F(hname_pref + "_ref", "", 20, 20, 80)
        #hRef = ROOT.TH1F(hname_pref, "", 40, 10, 90)
        #pt_bins = range(10,60,2) + range(60,80,5)
        pt_bins = range(10,50,2) + range(50,80,5)
        hRef = ROOT.TH1F(hname_pref, "", len(pt_bins)-1, array('f',pt_bins))
    elif "ele_eta" in var:
        var = "abs(" + var + ")"
        hRef = ROOT.TH1F(hname_pref, "", 27, 1.55, 2.9)
        #hRef = ROOT.TH1F(hname_pref, "", 40, 0, 1.5)
        #hRef = ROOT.TH1F(hname_pref, "", 60, 0, 3.0)
        #bins = range(15,29,1)/10 + [2.9,3.0]
        #hRef = ROOT.TH1F(hname_pref, "", len(bins)-1, array('f',bins))
    elif "ele_scleta" in var:
        var = "abs(" + var + ")"
        hRef = ROOT.TH1F(hname_pref, "", 24, 1.6, 2.8)
    elif "ele_ET" in var:
        hRef = ROOT.TH1F(hname_pref, "", 20, 0, 60)
    elif "Nvtx" in var:
        hRef = ROOT.TH1F(hname_pref, "", 16, 105, 185)
    elif "PU_density" in var:
        hRef = ROOT.TH1F(hname_pref, "", 10, 0, 1.3)
    elif "NPU" in var:
        #hRef = ROOT.TH1F(hname_pref, "", 18, 155, 245)
        bins = range(155,175,10) + range(175,245,10) + range(245,255,10)
        hRef = ROOT.TH1F(hname_pref, "", len(bins)-1, array('f',bins))
    elif "BDT" in var:
        hRef = ROOT.TH1F(hname_pref, "", 100, -1, 1)
    elif "ele_deltaetaele" in var:
        hRef = ROOT.TH1F(hname_pref, "", 100, -0.1, 0.1)
    else:
        hRef = ROOT.TH1F(hname_pref, "", 100, -100, 100)

    #hRef.Sumw2()

    tree.Draw(var + '>>' + hRef.GetName(),cuts,"e1 goff")
    print "Got %i events in hist %s" %(hRef.GetEntries(), hname_pref)
    return hRef


def plotEff(fname, var = "ele_pT", score = "0."):

    treename = "TestTree"
    tfile = ROOT.TFile(fname)
    tree = tfile.Get(treename)
    print "Found tree " + treename + " with %i events" % tree.GetEntries()

    outdir = os.path.basename(fname).replace(".root","_plotScale5")
    if not os.path.exists(outdir): os.makedirs(outdir)
    print "Storing output files in ", outdir

    otfile = ROOT.TFile(outdir + "/plots_" + var + ".root","recreate")

    # plot option
    plotOpt = 'e1'

    #cuts = "abs(ele_deltaetaele) < 0.02 &&"
    cuts = ""
    #cuts = "ele_ET > 5 &&"
    #cuts = "!(ele_sigmavv > 0.8 && ele_sigmauu > 0.7) && ele_hgc_nlay > 20 && "
    #cuts = "abs(ele_deltaetaele) < 0.02 && !(ele_sigmavv > 0.8 && ele_sigmauu > 0.7) && ele_hgc_nlay > 20 && "

    hRef = getHistFromTree(tree,var,cuts + "classID == 0","refSig")
    hSel = getHistFromTree(tree,var,cuts + "classID == 0 && BDT > %f " % score ,"selSig")

    #cuts = ""
    #cuts = "!(ele_sigmavv > 0.8 && ele_sigmauu > 0.7) && ele_hgc_nlay > 20 && "

    hRef2 = getHistFromTree(tree,var,cuts + "classID == 1","refBkg")
    hSel2 = getHistFromTree(tree,var,cuts + "classID == 1 && BDT > %f " % score,"selBkg")

    tfile.Close()

    otfile.cd()

    cname = "c"
    canv = phase2tdrStyle.setCanvas()
    #phase2tdrStyle.drawCMS(True)
    #phase2tdrStyle.drawEnPu()
    #canv = getCanvas(cname)
    #CMS_lumi.CMS_lumi(canv, iPeriod, iPos)
    #canv = ROOT.MakeCanvas("canv", "histograms", 800, 600)

    hRefEff = hRef.Clone(hRef.GetName()+'Eff')
    hRefEff.Divide(hRef)
    label = labels[var] if var in labels else var
    hRefEff.GetXaxis().SetTitle(label)
    hRefEff.GetYaxis().SetTitle("Selection efficiency")
    hRefEff.GetYaxis().SetRangeUser(0.,1.25)
    hRefEff.SetLineWidth(0)
    hRefEff.Draw()

    tEff = ROOT.TEfficiency(hSel,hRef)
    tEff.SetName(hRef.GetName() + "_Eff")
    tEff.SetLineColor(ROOT.kBlue+2)
    tEff.SetMarkerColor(ROOT.kBlue+2)
    tEff.SetMarkerStyle(7)
    tEff.Draw("same")

    #print hSel2.GetNbinsX(), hRef2.GetNbinsX()
    #hSel2.Scale(5)
    hSel2.Scale(5)
    #print hSel2.GetNbinsX(), hRef2.GetNbinsX()
    #hRef2.Scale(1/10.)
    tEff2 = ROOT.TEfficiency(hSel2,hRef2)
    tEff2.SetName(hRef2.GetName() + "_Eff")
    tEff2.SetLineColor(ROOT.kRed+2)
    tEff2.SetMarkerColor(ROOT.kRed+2)
    tEff2.SetMarkerStyle(7)
    tEff2.Draw("same")

    #ROOT.CMS_lumi( canv, 4, 10 )

    canv.Update()
    canv.Draw()
    ROOT.SetOwnership(canv,0)

    phase2tdrStyle.drawCMS(True)
    phase2tdrStyle.drawEnPu(pileup = 200)

    ## Legend
    #leg = ROOT.TLegend(0.45,0.45,0.9,0.55)
    #leg = ROOT.TLegend(0.65,0.4,1.,0.65)
    leg = ROOT.TLegend(0.6,0.4,0.95,0.7)
    ROOT.SetOwnership(leg,0)

    leg.SetBorderSize(0)
    #leg.SetTextFont(62)
    #leg.SetTextSize(0.05)
    leg.SetFillStyle(0)

    #leg.AddEntry(tEff,"DY: Z #rightarrow ee, <PU> = 200","pl")
    #leg.AddEntry(tEff2,"QCD, <PU> = 200","pl")
    #leg.SetHeader("<PU> = 200")
    #leg.SetHeader("10 < p_{T} < 20 GeV, 1.6 < |#eta| < 2.8")
    leg.SetHeader("p_{T} > 20 GeV, 1.6 < |#eta| < 2.8")
    leg.AddEntry(tEff,"Z #rightarrow ee","ple")
    leg.AddEntry(tEff2,"QCD multijets x5","ple")

    leg.Draw()

    hRefEff.Write()

    gr1 = tEff.GetPaintedGraph()
    gr1.SetName(tEff.GetName())
    gr1.Write()

    gr2 = tEff2.GetPaintedGraph()
    gr2.SetName(tEff2.GetName())
    gr2.Write()

    canv.SaveAs(outdir + "/eff_" + var + ".pdf")
    #canv.SaveAs(outdir + "/eff_" + var + ".root")

    #q = raw_input("q")

    otfile.Close()

def plotHist(fname, var = "BDT"):

    treename = "TestTree"
    tfile = ROOT.TFile(fname)
    tree = tfile.Get(treename)
    print "Found tree " + treename + " with %i events" % tree.GetEntries()

    outdir = os.path.basename(fname).replace(".root","")
    if not os.path.exists(outdir): os.makedirs(outdir)
    print "Storing output files in ", outdir

    otfile = ROOT.TFile(outdir + "/plots_hist.root","recreate")

    # plot option
    plotOpt = 'e1'

    #var = "ele_pT"
    #var = "ele_eta"
    #var = "Nvtx"
    #var = "ele_hgc_nlay"
    #var = "BDT"
    #var = "ele_deltaetaele"

    #cuts = "abs(ele_deltaetaele) < 0.02 &&"
    cuts = ""
    hRef = getHistFromTree(tree,var,cuts + "classID == 0","refSig")
    hRef2 = getHistFromTree(tree,var,cuts + "classID == 1","refBkg")

    hRef.SetLineColor(ROOT.kBlue+2)
    hRef2.SetLineColor(ROOT.kRed+2)

    #cuts = "abs(ele_deltaetaele) < 0.01 &&"
    #cuts = "!(ele_sigmauu > 0.6 && ele_sigmavv > 0.8) && "
    #cuts = "!(ele_sigmauu > 0.6 && ele_sigmavv > 0.8) && ele_hgc_nlay > 20 && "
    #cuts = "!(ele_sigmavv > 0.8) && ele_hgc_nlay > 20 && "
    hRef3 = getHistFromTree(tree,var,cuts + "classID == 0","refSig2")
    hRef4 = getHistFromTree(tree,var,cuts + "classID == 1","refBkg2")

    hRef3.SetLineColor(ROOT.kBlue+2)
    hRef4.SetLineColor(ROOT.kRed+2)

    tfile.Close()


    otfile.cd()

    cname = "c"
    canv = phase2tdrStyle.setCanvas()

    hRef.GetXaxis().SetTitle(var)
    #hRef.GetYaxis().SetTitle("Efficiency")
    hRef.DrawNormalized()

    hRef2.DrawNormalized("same")
    hRef3.DrawNormalized("same")
    hRef4.DrawNormalized("same")

    phase2tdrStyle.drawCMS(True)
    phase2tdrStyle.drawEnPu()

    canv.Update()
    canv.Draw()
    ROOT.SetOwnership(canv,0)

    hRef.Write()

    canv.SaveAs(outdir + "/" + var + ".pdf")

    q = raw_input("q")

    otfile.Close()

def plotROCs(rocs, names):

    cname = "c"
    canv = phase2tdrStyle.setCanvas()

    for i, roc in enumerate(rocs):

        if i == 0:
            roc.Draw("apl")
        else:
            roc.Draw("pl")

        roc.SetName(names[i])

    phase2tdrStyle.drawCMS(True)
    phase2tdrStyle.drawEnPu()

    canv.Update()
    canv.Draw()
    ROOT.SetOwnership(canv,0)

def main(fname = "/Users/artur/cernbox/www/HGCAL/reco/eleID/MVA/HGCTDRTMVA_1020_trackepshowerlonghgcaltdrV2DR01presel.root"):

    print "Reading file", fname
    '''
    outdir = os.path.basename(fname).replace(".root","_plots/")
    if not os.path.exists(outdir): os.makedirs(outdir)
    print "Storing output files in ", outdir
    '''

    ## Calculate ROC and define WP
    sig_scores, bkg_scores = getScores(fname)
    score_wp95 = getWPscore(sig_scores, 95)
    #score_wp95 = getWPscore(sig_scores, 80)
    print "Rejecting bck:", sum(bkg_scores > score_wp95)/float(len(bkg_scores))

    #var = "ele_pT"
    #var = "ele_eta"
    #var = "ele_ET"
    #var = "Nvtx"
    #var = "PU_density"
    #var = "NPU"
    #plotEff(fname, var, score_wp95)

    variabs = ["ele_pT","ele_eta","ele_scleta","Nvtx","PU_density","NPU","ele_ET"]

    for var in variabs:
        plotEff(fname, var, score_wp95)

    #plotHist(fname, var)

    '''
    roc_data = getROC(sig_scores, bkg_scores)
    #canv = plotROC(roc_data)
    gr = plotROC(roc_data)

    otfile = ROOT.TFile(outdir + "/" + "ROC.root","recreate")
    print "Writing to ", otfile.GetName()
    gr.Write()
    otfile.Close()
    #canv.SaveAs(outdir + "/" + "ROC.pdf")
    '''

if __name__ == "__main__":

    if len(sys.argv) > 1:
        if '-b' in sys.argv:
            fname = sys.argv[1]
        else:
            fname = sys.argv[1]
        print '# Input file is', fname
        main(fname)
    else:
        main()
