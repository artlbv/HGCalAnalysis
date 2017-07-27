#!/usr/bin/env python
# import ROOT
from NtupleDataFormat import HGCalNtuple

# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple


def main():
    #ntuple = HGCalNtuple("/Users/clange/CERNBox/partGun_PDGid211_x120_E80.0To80.0_NTUP_9.root")
    ntuple = HGCalNtuple("../hgcalNtuple-pca.root")

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

    for event in ntuple:
        if tot_nevents >= 10: break
        # print "Event", event.entry()
        tot_nevents += 1
        '''
        genParts = event.genParticles()
        tot_genpart += len(genParts)
        recHits = event.recHits()
        tot_rechit += len(recHits)
        if (ntuple.hasRawRecHits()):
            recHitsRaw = event.recHits("rechit_raw")
            tot_rechit_raw += len(recHitsRaw)
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)
        '''
        multiClusters = event.multiClusters()
        tot_multiclus += len(multiClusters)
        '''
        simClusters = event.simClusters()
        tot_simcluster += len(simClusters)
        pfClusters = event.pfClusters()
        tot_pfcluster += len(pfClusters)
        pfClusters = event.pfClusters()
        tot_pfcluster += len(pfClusters)
        caloParts = event.caloParticles()
        tot_calopart += len(caloParts)
        tracks = event.tracks()
        tot_track += len(tracks)
        '''

        for multicl in multiClusters:
            for item in multicl.pcaAxisX: #pxaPosX
                print item

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


if __name__ == "__main__":
    main()
