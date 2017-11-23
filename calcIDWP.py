#!/usr/bin/env python

from plotIDeff import *

def main(fname = "/Users/artur/cernbox/www/HGCAL/reco/eleID/MVA/HGCTDRTMVA_1020_trackepshowerlonghgcaltdrV2DR01presel.root"):

    print "Reading file", fname

    ## Calculate ROC and define WP
    sig_scores, bkg_scores = getScores(fname)

    bkg_effs = [0.1,1.0,10]

    for bkg_eff in bkg_effs:

        score = getWPscore(bkg_scores, bkg_eff)
        #print "Rejecting bck:", sum(bkg_scores > score_wp95)/float(len(bkg_scores))
        #print "Signal eff:", sum(sig_scores > score)/float(len(sig_scores))
        sig_eff = sum(sig_scores > score)/float(len(sig_scores))

        #print bkg_eff, score, sig_eff
        print "Bkg eff %0.2f, sig eff %0.2f, score %0.3f" %(bkg_eff, sig_eff, score)

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
