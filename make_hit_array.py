#!/usr/bin/env python
import sys
import ROOT, math
import numpy as np
#import cPickle as pickle

from NtupleDataFormat import HGCalNtuple
from helperTools import *
# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

ROOT.gROOT.SetBatch(1)

max_events = 100

def fill_event_data(event_list, particles, hits):
    ### Skim particle/hit data

    if not ( len(particles) > 0 and len(hits) > 0) : return 0

    skim_parts = []
    skim_hits = []

    for particle in particles:

        skim_parts.append(
            [ float(particle.posx()[10]),
              float(particle.posy()[10]),
              float(abs(particle.posz()[10])),
              float(particle.energy()),
              float(abs(particle.eta())),
              float(particle.phi()),
              ]
        )

    for hit in hits:

        skim_hits.append(
            [ float(hit.x()),
              float(hit.y()),
              float(abs(hit.z())),
              float(hit.energy()),
              float(abs(hit.eta())),
              float(hit.phi()),
            ]
        )

    event_list.append( [skim_parts, skim_hits] )

    return 1

def load_tree():

    print "loading tree"
    a = np.load("test.npy")

    print len(a)
    print a.shape
    print a.ndim


def main(fname = "../ntuples/hgcalNtuple_ele15_n100_testhelper.root"):
    ntuple = HGCalNtuple(fname)

    #foutname = fname.replace(".root","_hittree.root")
    #tfile = ROOT.TFile(foutname,"recreate")

    tot_nevents = 0
    tot_genpart = 0
    tot_rechit = 0

    hist_data = {}

    part_hit_data = []

    for event in ntuple:
        #if stop_run: break
        if tot_nevents >= max_events: break

        if tot_nevents % 100 == 0:
            print("Event %i" % tot_nevents)

        tot_nevents += 1

        genParts = event.genParticles()
        #tot_genpart += len(genParts)

        good_genparts = []
        for part in genParts:

            if part.gen() < 1: continue
            if part.pt() < 2: continue
            if part.reachedEE() != 2: continue

            good_genparts.append(part)

        if len(good_genparts) < 1: continue
        tot_genpart += len(good_genparts)

        ##### CLEAN CLOSEBY GEN PARTICLES (cannot resolve in HGCAL)
        ##
        # TODO
        ##

        ## filter flag 3 rechits (below 3sigma noise)
        good_hits = []#rh for rh in event.recHits() if rh.flags() < 3]
        for i_hit, hit in enumerate(event.recHits()):

            if i_hit > 1000: continue
            if hit.flags() > 2: continue
            if hit.layer() > 28: continue

            good_hits.append(hit)

        tot_rechit += len(good_hits)

        # split z+ and z- into two separate events
        for z_sign in [-1,+1]:

            particles = [part for part in good_genparts if part.eta() * z_sign > 0]
            if len(particles) == 0: continue

            hits = [hit for hit in good_hits if hit.z() * z_sign > 0] #[:100] ## take just 10 hits
            if len(hits) == 0: continue

            fill_event_data(part_hit_data, particles, hits)

    print "Finished event loop / event filling"

    #print np.array(part_hit_data)
    event_data = np.array(part_hit_data)

    np.save("test",event_data)


    ## save to pickle

    print("Found %i gen particles and %i rechits in %i events" %(tot_genpart,tot_rechit, tot_nevents))

    #tfile.Close()

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

    # load tree
    #load_tree()
