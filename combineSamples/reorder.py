#!/bin/env python2


setpoints = []

values = []

POINTS = False
FIELD = False

with open("postProcessing/sets/0.5/point_p.vtk", "rb") as f:
    for row in f:
        row.strip()
        if row.startswith("POINTS"):
            print "p"
            POINTS = True
            
        elif row.startswith("DATA"):
            print "-p"
            POINTS = False
        
        elif row.startswith("Temp"):
            print "t"
            FIELD = True
        
        elif POINTS:
#            setpoints.append([round(float(x), 6) for x in row.strip().split()])
            setpoints.append(row)
        
        elif FIELD:
            print "yes"
            values = [round(float(x), 6) for x in row.strip().split()]
        
#print setpoints

POINTS = False
FIELD = False

surfpoints = []

with open("postProcessing/surfaces/0.5/p_surfaces.vtk", "rb") as f:
    for row in f:
        if row.startswith("POINTS"):
            POINTS = True
            
        elif row.startswith("POLYGONS"):
            POINTS = False
            break
        
        elif POINTS:
            surfpoints.append([round(float(x), 6) for x in row.strip().split()])

    
#for point in surfpoints:
#    if point in setpoints:
#        print values[setpoints.index(point)]


