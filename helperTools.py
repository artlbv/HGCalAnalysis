#!/usr/bin/env python
import ROOT, math
import numpy as np

def addDataPoint(data,key,point):
    if key in data: data[key].append(point)
    else: data[key] = [point]

def getHisto(values, hname = "hist", htitle = "hist"):

    if len(values) == 0:
        print("Error! no values to build histo from!")
        return 0

    nbins = 100
    #nbins = len(values)*10/10

    # detect value type (1D/2D)
    #if "tuple" in type(values[0]):
    if isinstance(values[0], tuple):

        if len(values[0]) == 2:
            #hist_type = "2d"
            nbins /= 2

            # define histo
            xmin = min([val[0] for val in values])
            xmax = max([val[0] for val in values])
            ymin = min([val[1] for val in values])
            ymax = max([val[1] for val in values])

            if not xmin > 0:
                print hname,htitle,nbins,xmin,xmax,nbins,ymin,ymax

            hist = ROOT.TH2F(hname,htitle,nbins,xmin,xmax,nbins,ymin,ymax)

            # fill
            for val in values: hist.Fill(val[0],val[1])
        else:
            gr = ROOT.TGraph2D()
            gr.SetName(hname)
            gr.SetTitle(hname)
            gr.SetMarkerStyle(20)

            for i,val in enumerate(values):

                gr.SetPoint(i,val[0],val[1],val[2])

            return gr

    else:
        #hist_type = "1d"

        # define histo
        xmin,xmax  = min(values), max(values)
        hist = ROOT.TH1F(hname,htitle,nbins,xmin,xmax)

        # fill
        for val in values: hist.Fill(val)

    #if "int" in type(value[0]) or "float" in type(value[0]):
    #else:
    #    print("Unknown type %s" % type(value[0]))
    ROOT.SetOwnership(hist,0)
    return hist

def compDeltaVar(items, var = "x", prop = "std"):

    values = [getattr(item,var)() for item in items]
    #values = [item.x() for item in items]
    #print values

    res = getattr(np,prop)(values)
    return res

def calcDeltaRho(particle,cluster):

    drho = 999

    for layer in range(1,len(particle.posz())):
        if abs( particle.posz()[layer] - cluster.pcaPosZ() ) < 2:
            drho = math.hypot(cluster.slopeX()-particle.posx()[layer],cluster.slopeY()-particle.posy()[layer])

    return drho

    ### LATER
    ## cluster center: pcaPosZ
