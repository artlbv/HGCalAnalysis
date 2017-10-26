#!/usr/bin/env python
import os, sys
import ROOT
import phase2tdrStyle

def addOverflow(hist):

    #for ibin in range(1,hist.GetNbinsX()+1):
    nbins = hist.GetNbinsX()

    hist.SetBinContent(nbins,hist.GetBinContent(nbins) + hist.GetBinContent(nbins+1))
    hist.SetBinContent(1,hist.GetBinContent(1) + hist.GetBinContent(0))


indir = "/home/llr/cms/lobanov/HGCAL/reco/clustering/flat_ntuple/CMSSW_9_3_2/src/RecoNtuples/HGCalAnalysis/test/ntuples/"

samples = {}

#samples["1_sig"] = ("eleIDntuple_ele15_noPU_n10000_ele.root", "Ele, 15 GeV, no PU",ROOT.kRed+2)
#samples["1_sig"] = ("eleIDntuple_ele15_ele.root", "Ele, 15 GeV, no PU",ROOT.kRed+2)
samples["1_sig"] = ("eleIDntuple_ele15_ele.root", "Ele, 15 GeV, no PU",ROOT.kRed+2)
samples["2_sig"] = ("eleIDntuple_ele15_PU140_n1000_ele.root", "Ele, 15 GeV, PU140",ROOT.kMagenta+2)
samples["3_bkg"] = ("eleIDntuple_qcd_n10000_ele.root", "QCD, 15-7000 GeV, no PU",ROOT.kCyan+2)
samples["4_bkg"] = ("eleIDntuple_pi25_n10000_ele.root", "Pi, 25 GeV, no PU",ROOT.kYellow+2)
samples["5_bkg"] = ("eleIDntuple_pi25_pu140_n10000_ele.root", "Pi, 25 GeV, 140 PU",ROOT.kYellow-7)
'''
samples["1_sig"] = ("eleIDntuple_ele15_noPU_n1000_ele.root", "Ele, 15 GeV, Multicl",ROOT.kRed+2)
samples["2_sig"] = ("eleIDntuple_ele15_noPU_n1000_SIMCL_ele.root", "Ele, 15 GeV, SIMCL",ROOT.kCyan+2)
'''

hists = {}

for sample in sorted(samples):

    fname, title, color = samples[sample]

    print fname, title

    tfile = ROOT.TFile(indir + fname)

    for key in tfile.GetListOfKeys():
        obj = tfile.Get(key.GetName())
        if obj.ClassName() != "TH1F": continue

        ## Copy hist to storage
        hname = sample + "_" + obj.GetName()
        hist = obj.Clone(hname)
        hist.SetDirectory(0)

        hist.SetTitle(title)
        hist.GetXaxis().SetTitle(obj.GetName().replace("ele_",""))

        hist.SetLineWidth(2)
        hist.SetLineColor(color)
        hist.SetFillColorAlpha(color,0.5)
        hist.SetMarkerColor(color)

        ## add overflow bins
        addOverflow(hist)

        if obj.GetName() in hists:
            hists[obj.GetName()].append(hist)
        else:
            hists[obj.GetName()] = [hist]

    tfile.Close()

#outdir = "/eos/user/a/alobanov/www/HGCAL/reco/eleID/compare_ele_ID_vars/"
outdir = "compare_vars_wOverflow/"
#outdir = "compare_vars/"
#outdir = "compare_vars_simcl/"

if not os.path.exists(outdir): os.makedirs(outdir)

for hname in hists:

    cname = "c" + hname
    canv = phase2tdrStyle.setCanvas()
    canv.SetName(cname)
    canv.SetTitle(hname)

    for i, hist in enumerate(hists[hname]):

        if i == 0: hist.DrawNormalized()
        else: hist.DrawNormalized("same")

    leg = canv.BuildLegend()
    leg.Draw()
    leg.SetX1(0.65)
    leg.SetX2(0.9)
    leg.SetY1(0.7)
    leg.SetY2(0.9)

    phase2tdrStyle.drawCMS(True)
    phase2tdrStyle.drawEnPu()

    canv.Update()
    canv.Draw()
    #canv.SetLogy()

    #q = raw_input("q")

    canv.SaveAs(outdir + cname + ".pdf")
    canv.SaveAs(outdir + cname + ".png")
