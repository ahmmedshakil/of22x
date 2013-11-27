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
#include "geometryCut.H"
#include "unitConversion.H"


using namespace Foam;


namespace Foam
{
    scalar alignedCos_ = cos(degToRad(89.0));
}

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
    const vectorField& normals = surf.faceNormals();
    

    pointField points = mesh.points();
    labelList edgeLabels(mesh.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
    const edgeList edges = mesh.edges();
    
    
    triSurfaceSearch querySurf(surf);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
    
    DynamicList<GeometryCut> cuts(mesh.nPoints()*4);
    
    Info << "find hits" << nl;
    
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
                    const GeometryCut newCut(startPoint, pHit.index());
                    cuts.append(newCut);
                }
                else if (mag(pHit.hitPoint() - pEnd) < 0.01 * eMag)
                {
                    const label endPoint = e.end();
                    GeometryCut newCut(endPoint, pHit.index());
                    cuts.append(newCut);
                    break;
                }
                else if (mag(n & normals[pHit.index()]) < alignedCos_)
                {
                    break;
                }
                else
                {
                    const vector eVec(pEnd - pStart);
                    const vector pVec(pHit.hitPoint() - pStart);
                    const scalar weight = mag(pVec)/mag(eVec);
                    const GeometryCut newCut(edgeI, pHit.index(), weight);
                    cuts.append(newCut);
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
    
//    Info << "inter " << cuts << nl;

    Info << "ready" << nl;
    
    cuts.shrink();
    
    
    
//  find faceCuts

    labelListList facesCuts(mesh.nFaces());
    forAll(cuts, cutI)
    {
        GeometryCut cut = cuts[cutI];
        labelList faces;
        
        if(cut.isEdge())
        {
            label edge = cut.geometry();
            faces = mesh.edgeFaces()[edge];
        }
        else
        {
            label point = cut.geometry();
            faces = mesh.pointFaces()[point];
        }
        forAll(faces, faceI)
        {
            label face = faces[faceI];
            facesCuts[face].append(cutI);
        }
    }
//    Info << "facesCuts" << facesCuts << nl;
    
    
//    List<bool> addEdges(mesh.nEdges(), false);
//    List<bool> addPoints(mesh.nPoints(), false);
//    List<bool> rmEdges(mesh.nEdges(), false);
//    List<bool> rmPoints(mesh.nPoints(), false);
    List<bool> addedCuts(cuts.size(), false);
    List<bool> removedCuts(cuts.size(), false);

    DynamicList<label> allCutPoints(mesh.nPoints());
    DynamicList<label> allCutEdges(mesh.nEdges());
    DynamicList<scalar> cutEdgeWeights(mesh.nEdges());
    

    DynamicList<label> cutFaces(mesh.nFaces());
    forAll(facesCuts, faceI)
    {
        labelList faceCuts = facesCuts[faceI];
        label nCuts = faceCuts.size();
        
        Info << nl << "nCuts " << faceI << ": " << nCuts << nl;
        
        if(nCuts > 1)
        {
//            if(nCuts == 2)
//            {
//                forAll(faceCuts, faceCutI)
//                {
//                    label cutI = faceCuts[faceCutI];
//                    GeometryCut cut = cuts[cutI];
//                    
////                    if(cut.isEdge())
////                    {
////                        addedCuts[cutI] = true;
////                    }
//                    if(!removedCuts[cutI])
//                    {
//                        addedCuts[cutI] = true;
//                    }
//                }
//            }
//            else
//            {
                label startCutI = -1;
                label nextCutI = -1;
                DynamicList<label> otherCutsI(faceCuts.size());
                
                forAll(faceCuts, faceCutI)
                {
                    label cutI = faceCuts[faceCutI];
                    GeometryCut cut = cuts[cutI];
                    Info << "cut " << faceCutI << ": (" << cutI << ") " << cut << nl;
                    
                    if(!removedCuts[cutI])
                    {
                        if(cut.isEdge() && startCutI == -1)
                        {
                            addedCuts[cutI] = true;
                            startCutI = cutI;
                        }
                        else
                        {
//                            Info << "append this" << nl;
                            otherCutsI.append(cutI);
                        }
                    }
                }
                
                if(startCutI != -1)
                {
                    Info << "startCut " << cuts[startCutI] << nl;
//                    Info << "nOther " << otherCutsI.size() << nl;
                    forAll(otherCutsI, otherCutI)
                    {
                        label cutI = otherCutsI[otherCutI];
                        Info << cuts[cutI] << nl;
                    }
//                    Info << "otherCutsI " << otherCutsI << nl;
                    nextCutI = cuts[startCutI].findNext(surf, cuts, otherCutsI);
                }
                
                
                if(nextCutI != -1)
                {
                    Info << "nextCut " << cuts[nextCutI] << nl;
                    addedCuts[nextCutI] = true;
                    
                    forAll(otherCutsI, otherCutI)
                    {
                        label cutI = otherCutsI[otherCutI];
                        if(cutI != startCutI && cutI != nextCutI)
                        {
                            removedCuts[cutI] = true;
                        }
                    }
                }
                
//            }
        }   
    }
    
//    Info << "addedCuts   " << addedCuts << nl;
//    Info << "removedCuts " << removedCuts << nl;
    
    
    forAll(addedCuts, cutI)
    {
        bool addCut = addedCuts[cutI];
        if(addCut)
        {
            GeometryCut cut = cuts[cutI];
            if(cut.isEdge())
            {
                allCutEdges.append(cut.geometry());
                cutEdgeWeights.append(cut.weight());
            }
            else
            {
                allCutPoints.append(cut.geometry());
            }
        }
    }
    
    
    
    
    
    allCutPoints.shrink();
    allCutEdges.shrink();
    
    scalarField allCutEdgeWeights;
    allCutEdgeWeights.transfer(cutEdgeWeights);
    cutEdgeWeights.clear();
    
//    Info << "allCutPoints" << allCutPoints << nl;
//    Info << "allCutEdges" << allCutEdges << nl;
//    Info << "allCutEdgeWeights" << allCutEdgeWeights << nl;
    
    
    OFstream pointStream("cutPoints.obj");
    
    forAll(allCutPoints, pointI)
    {
        label point = allCutPoints[pointI];
        Info << points[point] << nl;
        meshTools::writeOBJ(pointStream, points[point]);
    }
    
    forAll(allCutEdges, edgeI)
    {
        scalar weight = allCutEdgeWeights[edgeI];
        label edgeL = allCutEdges[edgeI];
        const edge e = edges[edgeL];
        const point pStart = points[e.start()];
        const point pEnd = points[e.end()];
        const point cutPoint = pStart + (pEnd - pStart) * weight;
        
        meshTools::writeOBJ(pointStream, cutPoint);
    }
    
//    forAll(cuts, cutI)
//    {
//        GeometryCut cut = cuts[cutI];
//        
//        if(cut.isEdge())
//        {
//            label edge = cut.elem();
//            labelList faces = mesh.edgeFaces(edge);
//            forAll(faces, faceI)
//            {
//                label face = faces[faceI];
//                
//            }
//        }
//        
//        
//    }
    
//    
    
    cellCuts cut
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
    cutter.setRefinement(cut, meshMod);

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
