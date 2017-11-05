#!/usr/bin/env python
import os, sys
import ROOT
import phase2tdrStyle

if __name__ == "__main__":

    fname = "merge_v3_ele_pT.root"

    tfile = ROOT.TFile(fname)

    cname = "c"
    canv = phase2tdrStyle.setCanvas()
    canv.SetName(cname)
    canv.SetTitle("")

    h_ele_pT_refSigEff = tfile.Get("h_ele_pT_refSigEff")
    h_ele_pT_refSig_Eff = tfile.Get("h_ele_pT_refSig_Eff")
    h_ele_pT_refBkg_Eff = tfile.Get("h_ele_pT_refBkg_Eff")

    h_ele_pT_refSigEff.Draw()
    h_ele_pT_refSigEff.GetXaxis().SetRangeUser(10,70)
    h_ele_pT_refSigEff.GetYaxis().SetRangeUser(0,1.2)

    h_ele_pT_refSig_Eff.Draw("same")
    h_ele_pT_refBkg_Eff.SetTitle("DY")
    h_ele_pT_refBkg_Eff.Draw("same")
    h_ele_pT_refBkg_Eff.SetTitle("QCD")

    '''
    leg = canv.BuildLegend()
    leg.Draw()
    leg.SetX1(0.65)
    leg.SetX2(0.9)
    leg.SetY1(0.7)
    leg.SetY2(0.9)
    '''

    phase2tdrStyle.drawCMS(True)
    phase2tdrStyle.drawEnPu()

    canv.Update()
    canv.Draw()
    #canv.SetLogy()

    #q = raw_input("q")

    canv.SaveAs("pt_eff_merge.pdf")
