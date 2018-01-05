#!/usr/bin/env python
import os, sys
import ROOT
#import phase2tdrStyle

ROOT.gROOT.LoadMacro("plot_style/HttStyles.cc")
ROOT.gROOT.LoadMacro("plot_style/CMS_lumi.C")
ROOT.setTDRStyle()


if __name__ == "__main__":

    fname = "merge_v6_ele_pT.root"

    tfile = ROOT.TFile(fname)

    cname = "c"
    #canv = phase2tdrStyle.setCanvas()
    #canv.SetName(cname)
    #canv.SetTitle("")
    canv = ROOT.MakeCanvas("canv", "histograms", 800, 600)

    h_ele_pT_refSigEff = tfile.Get("h_ele_pT_refSigEff")
    h_ele_pT_refSig_Eff = tfile.Get("h_ele_pT_refSig_Eff")
    h_ele_pT_refBkg_Eff = tfile.Get("h_ele_pT_refBkg_Eff")

    h_ele_pT_refSigEff.Draw()
    h_ele_pT_refSigEff.GetXaxis().SetRangeUser(10,60)
    h_ele_pT_refSigEff.GetYaxis().SetRangeUser(0,1.2)

    h_ele_pT_refSig_Eff.Draw("pe1same")
    h_ele_pT_refBkg_Eff.SetTitle("DY")
    h_ele_pT_refBkg_Eff.Draw("pe1same")
    h_ele_pT_refBkg_Eff.SetTitle("QCD")

    #phase2tdrStyle.drawCMS(True)
    #phase2tdrStyle.drawEnPu(pileup = 200)
    ROOT.CMS_lumi( canv, 4, 0 )

    ## Legend
    #leg = ROOT.TLegend(0.45,0.45,0.9,0.55)
    #leg = ROOT.TLegend(0.65,0.4,1.,0.65)
    leg = ROOT.TLegend(0.6,0.4,0.9,0.65)

    leg.SetFillStyle(0)


    ROOT.SetOwnership(leg,0)

    leg.SetBorderSize(0)
    #leg.SetTextFont(62)
    #leg.SetTextSize(0.05)

    #leg.AddEntry(h_ele_pT_refSig_Eff,"DY: Z #rightarrow ee, <PU> = 200","pl")
    #leg.AddEntry(h_ele_pT_refBkg_Eff,"QCD, <PU> = 200","pl")
    #leg.SetHeader("<PU> = 200")
    leg.SetHeader("1.6 < |#eta| < 2.8")
    leg.AddEntry(h_ele_pT_refSig_Eff,"Z #rightarrow ee","pl")
    leg.AddEntry(h_ele_pT_refBkg_Eff,"QCD multijets x10","pl")

    leg.Draw()

    canv.Update()
    canv.Draw()
    #canv.SetLogy()

    #q = raw_input("q")

    canv.SaveAs("pt_eff_merge.pdf")
