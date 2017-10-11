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

        elif len(values[0]) == 2:
            gr = ROOT.TGraph()
            gr.SetName(hname)
            gr.SetTitle(hname)
            gr.SetMarkerStyle(20)

            for i,val in enumerate(values):
                gr.SetPoint(i,val[0],val[1])

            return gr

        elif len(values[0]) == 3:
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

def make1DHist(values, hname = "hist", htitle = "htitle", nbins = None, xmin = None, xmax = None):

    if nbins == None: return getHisto(values, hname, htitile)

    hist = ROOT.TH1F(hname,htitle,nbins,xmin,xmax)
    for val in values: hist.Fill(val)

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

    ## Calculate veritcal/perp angles
    # 1. multicluster center
    mcl_cent = (multicl.pcaPosX(), multicl.pcaPosY(), multicl.pcaPosZ())
    a_p,a_v = calc_angles(mcl_cent,mclut_vect,part_vect)

    return a_p,a_v

#### axis line equation:
"""
point: x0,y0,z0
vector: a,b,c

eq:
x = x0 + ta
y = y0 + tb
z = z0 + tc
"""

def get_entry_point(point,axis_vector,plane_z):
    """
    z0 = plane_z, solve eq

    t = (z-z0)/c
    x = x0 + a*t
    y = y0 + b*t
    """

    if axis_vector.Z() == 0: return 0

    t = (plane_z - point[2]) / axis_vector.Z()
    x = point[0] + axis_vector.X() * t
    y = point[1] + axis_vector.Y() * t

    return x,y

### layer manipulation

def get_missing_layers(layers):

    cnt = 0
    for lay in range(min(layers), max(layers)+1):
        #if lay not in range(5,22): continue
        if lay not in layers: cnt +=1

    return cnt
