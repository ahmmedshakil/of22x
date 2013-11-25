/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    cutHexMesh

Description
    Utility to mesh STL-Surfaces using a cut-cell approach.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "cellCuts.H"
#include "cellSet.H"
#include "meshCutter.H"
#include "triSurfaceSearch.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "surfaceFeatures.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "triangle.H"
#include "polyTopoChange.H"
#include "OFstream.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label findNextIntersection
(
   const triSurface& surf,
   const label startTri,
   const labelHashSet& nextTris
)
{
    label nextIntersection = -1;
    labelHashSet visited(1000);
    
    labelList faceEdges = surf.faceEdges()[startTri];
    label thisEdge =  faceEdges[1];
    label lastEdge1 = faceEdges[0];
    label lastEdge2 = faceEdges[2];
    
    label lastTri = startTri;
    
    for(int i = 0; i < 1000; i++)
    {
        label tryTri1 = triSurfaceTools::otherFace(surf, lastTri, lastEdge1);
        label tryTri2 = triSurfaceTools::otherFace(surf, lastTri, lastEdge2);
        
        if(visited[tryTri1] && tryTri2 != -1)
        {
            thisEdge = lastEdge2;
            lastTri = tryTri2;
        }
        else if (tryTri1 != -1)
        {
            thisEdge = lastEdge1;
            lastTri = tryTri1;
        }
        
        triSurfaceTools::otherEdges(
            surf, 
            lastTri, 
            thisEdge, 
            lastEdge1, 
            lastEdge2
        );
        
        if(nextTris[lastTri])
        {
            nextIntersection = lastTri;
            break;
        }
        
        visited.insert(lastTri);
    }
    
    return nextIntersection;
}
    



int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    argList::noParallel();
    argList::validArgs.append("input surfaceFile");
//    argList::validArgs.append("feature angle");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");
    const fileName surfName = args[1];
    triSurface surf(runTime.constantPath()/"triSurface"/surfName);  
    

    pointField points = mesh.points();
    labelList edgeLabels(mesh.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
    const edgeList edges = mesh.edges();
    
    
    triSurfaceSearch querySurf(surf);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
    
    
    
//    DynamicList<label> allCutPoints(mesh.nPoints());
    List<labelHashSet> pointTris(mesh.nPoints());
//    DynamicList<label> allCutEdges(mesh.nEdges()); 
    List<DynamicList<scalar> > edgeWeights(mesh.nEdges()); 
    List<DynamicList<label> > edgeTris(mesh.nEdges()); 
    List<bool> pointCuts(mesh.nPoints(), false); 
    List<bool> edgeCuts(mesh.nEdges(), false);
    scalar weight;
    
    Info << "find hits" << nl;
    
        
//    OFstream pointStream("cutPoints.obj");
//    meshTools::writeOBJ(pointStream, points[i]);

            

    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];
        const edge e = edges[edgeI];
        const point pStart = points[e.start()] ;
        const point pEnd = points[e.end()] ;
        
        point p0 = pStart;
        const point p1 = pEnd;
        
        const vector eVec(pEnd - pStart);
        const scalar eMag = mag(eVec);
        const vector n(eVec/(eMag + VSMALL));
        const point tolVec = 1e-6*eVec;
              
        while(true)
        {
            pointIndexHit pHit = tree.findLine(p0, p1);
            
            if(pHit.hit())
            {
                
                if (mag(pHit.hitPoint() - pStart) < 0.01 * eMag)
                {
                    const label startPoint = e.start();
//                    allCutPoints.append(startPoint);
                    pointTris[startPoint].insert(pHit.index());
                    pointCuts[startPoint] = true;
                }
                else if (mag(pHit.hitPoint() - pEnd) < 0.01 * eMag)
                {
                    const label endPoint = e.end();
//                    allCutPoints.append(endPoint);
                    pointTris[endPoint].insert(pHit.index());
                    pointCuts[endPoint] = true;
                    break;
                }
                else
                {
                    const vector eVec(pEnd - pStart);
                    const vector pVec(pHit.hitPoint() - pStart);
                    
                    weight = mag(pVec)/mag(eVec);
//                    allCutEdges.append(edgeI);
                    edgeWeights[edgeI].append(weight);
                    edgeTris[edgeI].append(pHit.index());
                    edgeCuts[edgeI] = true;
                }
                p0 = pHit.hitPoint() + tolVec;
            }
            else
            {
                // No hit.
                break;
            }
        }
    }
    
    DynamicList<label> cutPoints(mesh.nPoints()); 
    forAll(pointCuts, pointCutI)
    {
        bool cutPoint = pointCuts[pointCutI];
        if(cutPoint)
        {
            cutPoints.append(pointCutI);
        }
    }
    
    DynamicList<label> cutEdges(mesh.nEdges()); 
    forAll(edgeCuts, edgeCutI)
    {
        bool cutEdge = edgeCuts[edgeCutI];
        if(cutEdge)
        {
            cutEdges.append(edgeCutI);
        }
    }
    
   
    Info << "ready" << nl;
    
    List<label> nFaceCuts(mesh.nFaces(), 0);
    
    labelListList cutFacePoints(mesh.nFaces());
    forAll(cutPoints, cutPointI)
    {
        label cutPoint = cutPoints[cutPointI];
        label nCuts = pointTris[cutPoint].size();
        
        labelList cutFaces = mesh.pointFaces()[cutPoint];
        
        forAll(cutFaces, cutFaceI)
        {
            label cutFace = cutFaces[cutFaceI];
            nFaceCuts[cutFace] += nCuts;
            cutFacePoints[cutFace].append(cutPoint);
        }   
    }
    
    Info << "cutFacePoints" << cutFacePoints << nl;
    
    labelListList cutFaceEdges(mesh.nFaces());   
    forAll(cutEdges, cutEdgeI)
    {
        label cutEdge = cutEdges[cutEdgeI];
        labelList cutFaces = mesh.edgeFaces()[cutEdge];
        label nCuts = edgeTris[cutEdge].size();
        
        forAll(cutFaces, cutFaceI)
        {
            label cutFace = cutFaces[cutFaceI];
            nFaceCuts[cutFace] += nCuts;
            cutFaceEdges[cutFace].append(cutEdge);
        }   
    }
    
    Info << "cutFaceEdges" << cutFaceEdges << nl;
    
    Info << "cutPoints" << cutPoints << nl;
    Info << "pointTris" << pointTris << nl;
    Info << "cutEdges" << cutEdges << nl;
    Info << "edgeTris" << edgeTris << nl;
    Info << "edgeWeights" << edgeWeights << nl;
    
    
    
    Info << "nFaceCuts" << nFaceCuts << nl;
    
    
    labelList n(20, 0);
    
    forAll(nFaceCuts, faceI)
    {
        label nCuts = nFaceCuts[faceI];
        n[nCuts]++;
    }
    
    Info << "nCuts" << n << nl;
    
    DynamicList<label> allCutPoints(mesh.nPoints());
    DynamicList<label> allCutEdges(mesh.nEdges());
    DynamicList<scalar> cutEdgeWeights(mesh.nEdges());
    
    List<bool> addedEdges(mesh.nEdges(), false);
    List<bool> addedPoints(mesh.nPoints(), false);
    List<bool> removedEdges(mesh.nEdges(), false);
    List<bool> removedPoints(mesh.nPoints(), false);
    
    forAll(nFaceCuts, cutFace)
    {
        label faceCuts = nFaceCuts[cutFace];
        if(faceCuts == 2)
        {
            labelList cutEdges = cutFaceEdges[cutFace];
            forAll(cutEdges, cutEdgeI)
            {
                label edge = cutEdges[cutEdgeI];
                scalarList weights = edgeWeights[edge];
                if(weights.size() == 1 && !addedEdges[edge])
                {
                    allCutEdges.append(edge);
                    scalar weight = weights[0];
                    cutEdgeWeights.append(weight);
                    addedEdges[edge] = true;
                }
            }
            labelList cutPoints = cutFacePoints[cutFace];
            forAll(cutPoints, cutPointI)
            {
                label point = cutPoints[cutPointI];
                if(pointTris[point].size() == 1 && !addedPoints[point])
                {
                    allCutPoints.append(point);
                    addedPoints[point] = true;
                }
            }
        }
        else if(faceCuts > 2)
        {
            labelList cutEdges = cutFaceEdges[cutFace];
            labelList cutPoints = cutFacePoints[cutFace];
            label startEdge = -1;
            label startPoint = -1;
            
            if(cutEdges.size() > 0)
            {
                forAll(cutEdges, edgeI)
                {
                    label edge = cutEdges[edgeI];
                    if(!removedEdges[edge])
                    {
                        startEdge = edge;
                        break;
                    }
                }
            }
            else if(cutPoints.size() > 0)
            {
                forAll(cutPoints, pointI)
                {
                    label point = cutPoints[pointI];
                    if(!removedPoints[point])
                    {
                        startPoint = point;
                        break;
                    }
                }
            }
                
            if(startEdge != -1)
            {
                labelHashSet nextTris;
                
                forAll(cutEdges, edgeI)
                {
                    label edge = cutEdges[edgeI];
                    if(!removedEdges[edgeI] && edge != startEdge)
                    {
                        nextTris.insert(edgeTris[edge]);
                    }
                }
                forAll(cutPoints, pointI)
                {
                    label point = cutPoints[pointI];
                    if(!removedPoints[pointI] && point != startPoint)
                    {
                        nextTris.insert(pointTris[point].toc());
                    }
                }
                
                label nextTri = findNextIntersection(
                    surf,
                    startEdge,
                    nextTris
                );
                
                forAll(cutEdges, edgeI)
                {
                    label edge = cutEdges[edgeI];
                    labelList tris = edgeTris[edge];
                    scalarList weights = edgeWeights[edge];
                    forAll(tris, triI)
                    {
                        label tri = tris[triI];
                        
                        if(tri == nextTri)
                        {
                            if(!addedEdges[edge])
                            {
                                allCutEdges.append(edge);
                                scalar weight = weights[triI];
                                cutEdgeWeights.append(weight);
                                Info << "start: " << startEdge << nl 
                                     << "next: " << edge << nl << nl;
                                addedEdges[edge] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    
    Info  << "allCutEdges" << allCutEdges << nl;
    Info  << "cutEdgeWeights" << cutEdgeWeights << nl;
    Info  << "allCutPoints" << allCutPoints << nl;

    scalarField allCutEdgeWeights;
    allCutEdgeWeights.transfer(cutEdgeWeights);
    cutEdgeWeights.clear();
//   
    
//    forAll(cutPoints, i)
//    {
//        Info << mesh.pointPoints()[cutPoints[i]];
//    }


//    Finde die nachbarn und entferne Punkt wenn er 2 nachbarn hat
    
    
    cellCuts cuts
    (
        mesh,
        allCutPoints,
        allCutEdges,
        allCutEdgeWeights
    );
    
    polyTopoChange meshMod(mesh);

    // Cutting engine
    meshCutter cutter(mesh);

    // Insert mesh refinement into polyTopoChange.
    cutter.setRefinement(cuts, meshMod);

    // Do all changes
    Info<< "Morphing ..." << endl;


    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    // Update stored labels on meshCutter.
    cutter.updateMesh(morphMap());
       
    if (!overwrite)
    {
        runTime++;
    }
    
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
//    Info<< "Writing refined morphMesh to time " << runTime.timeName() << endl;

    mesh.write(); 
    
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
