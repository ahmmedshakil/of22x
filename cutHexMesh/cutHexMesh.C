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
/*
TODO

Finde für jede geschnittene Zelle die innenliegenden triSurface Flächen.

Zelle ist geschnitten und beinhaltet Flächenpunkte.
    -> Schnittflächen bestehen aus mehreren triangles.
    -> Finde alle triangles, welche mit den innenliegenden Punkten verbunden sind.
    -> speichere triangles als zusammengehörig, wenn sie punkte teilen, welche innerhalb der zelle liegen.
    -> getrennte flächen haben keine gemeinsamen punkte innerhalb der Zelle.
    -> Ergebnis sind zusammenhängende flächen, pro schnitt durch eine Zelle.
    
    FINDE AUCH FLÄCHEN SCHNEIDENDE TRISURFACE KANTEN
    
Methode:
    Überprüfe in welchen Zellen sich die Flächenpunkte befinden. 
    nehme die zugehörigen triangles und füge sie einer liste hinzu.
    Überprüfe, welche Flächen die kanten schneiden. 
    nehme die zugehörigen triangles und füge sie ebenfalls den jeweiligen zellenlisten hinzu.



Falls geschnittene Punkte mehrmals auftreten entferne redundante cuts.
entweder punkte und triangles sind identisch -> identischer cut
Falls punkte identisch aber triangles unterschiedlich:
finde nachbar triangles (mit eckpunkten)
falls nachbartriangle identisch -> identischer cut
ansonsten: cut durch andere fläche -> mehrfacher cut


PackedBoolList elemToRemove(l.size());

forAll(l, i)
{
    if (checkCondition)
    {
        elemToRemove[i] = true;
    }
}

inplaceSubset(elemToRemove, l);
*/


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
#include "SortableList.H"


using namespace Foam;


namespace Foam
{
    scalar alignedCos_ = cos(degToRad(89.0));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


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
    DynamicList<label> cutPoints(mesh.nPoints());
    DynamicList<label> cutPointsTriangle(mesh.nPoints());
    DynamicList<label> cutEdges(mesh.nEdges());
    DynamicList<label> cutEdgesTriangle(mesh.nEdges());
    DynamicList<label> cutEdgesWeight(mesh.nEdges());
    
    
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
//                    pointCuts.append(startPoint);
//                    pointCutsTriangle.append(pHit.index());
                    const GeometryCut newCut(startPoint, pHit.index());
                    cuts.append(newCut);
                }
                else if (mag(pHit.hitPoint() - pEnd) < 0.01 * eMag)
                {
                    const label endPoint = e.end();
//                    pointCuts.append(endPoint);
//                    pointCutsTriangle.append(pHit.index());
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
//                    edgeCuts.append(edgeI);
//                    edgeCutsTriangle.append(pHit.index());
//                    edgeCutsWeight.append(weight);
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

    Info << "ready" << nl;
    
    cuts.shrink();
    sort(cuts);
    
    // Remove duplicate cuts
    labelList keepCuts;
    uniqueOrder(cuts, keepCuts);
    
    
    
        
    
    for (int i = 0; i < keepCuts.size() - 1; i++)
    {
        bool keep = true;
        label keepCut = keepCuts[i];
        label nextKeepCut = keepCuts[i+1];
        
        if (cuts[keepCut].isPoint() && cuts[nextKeepCut].isPoint())
        {
            if (cuts[keepCut].geometry() == cuts[nextKeepCut].geometry())
            {
                Info << cuts[keepCut] << nl;
                
                label triangle1 = cuts[keepCut].triangle();
                labelList pointList1(surf[triangle1]);
                labelHashSet points1;
                points1.insert(pointList1);
                
                label triangle2 = cuts[nextKeepCut].triangle();
                labelList pointList2(surf[triangle2]);
                labelHashSet points2;
                points2.insert(pointList2);
                
                if( (points1 & points2).size() > 0)
                {
                    Info << "remove " << keepCut << nl;
                }
                
                
//                points1 &= points2;
                
//                Info << result << nl;
                
//                Info << tri << nl;
                
            }
        }
        
    }
    
    
    
    PackedBoolList cutsToKeep(cuts.size() , false);
    forAll(keepCuts, keepCutsI)
    {
        label cutI = keepCuts[keepCutsI];
        cutsToKeep[cutI] = true;
    }
    inplaceSubset(cutsToKeep, cuts);
    
    
    
    Info << cuts << nl;

    OFstream cutPointStream("cuts.obj");
    
    forAll(cuts, cutI)
    {
        point addPoint;
        GeometryCut cut = cuts[cutI];
        if(cut.isPoint())
        {
            label p = cut.geometry();
            addPoint = points[p];
        }
        else
        {
            label cutEdge = cut.geometry();
            scalar weight = cut.weight();
            const edge e = edges[cutEdge];
            const point pStart = points[e.start()];
            const point pEnd = points[e.end()];
            addPoint = pStart + (pEnd - pStart) * weight;
        }
//        Info << points[point] << nl;
        meshTools::writeOBJ(cutPointStream, addPoint);
    }
    
    
    
//    List<bool> addEdges(mesh.nEdges(), false);
//    List<bool> addPoints(mesh.nPoints(), false);
//    List<bool> rmEdges(mesh.nEdges(), false);
//    List<bool> rmPoints(mesh.nPoints(), false);
//    
//    List<bool> removedCuts(cuts.size(), false);
//    
//    List<bool> addedCuts(cuts.size(), false);
//    List<bool> visitedFaces(mesh.nFaces(), false);
//    
    
//    //  find faceCuts

//    labelListList facesCuts(mesh.nFaces());
//    forAll(cuts, cutI)
//    {
//        GeometryCut cut = cuts[cutI];
//        labelList faces;
//        
//        if(cut.isEdge())
//        {
//            label edge = cut.geometry();
//            faces = mesh.edgeFaces()[edge];
//        }
//        else
//        {
//            label point = cut.geometry();
//            faces = mesh.pointFaces()[point];
//        }
//        forAll(faces, faceI)
//        {
//            label face = faces[faceI];
//            facesCuts[face].append(cutI);
//        }
//    }
////    Info << "facesCuts" << facesCuts << nl;
//    
    
    
    
//    DynamicList<label> checkCuts;
//    
//    forAll(cuts, cutI)
//    {
//        if(cuts[cutI].isEdge())
//        {
//            checkCuts.append(cutI);
//            break;
//        }
//    }
    
    
    
//    DynamicList<label> newCheckCuts = checkCuts;
//    
//    while(newCheckCuts.size() > 0)
//    {
//        checkCuts.clear();
//        checkCuts.transfer(newCheckCuts);
////        checkCuts = newCheckCuts;
//        newCheckCuts.clear();
//        
//        
////        Info << nl << "checkCuts " << checkCuts << nl;
//        
//        forAll(checkCuts, checkCutI)
//        {
//            label cutI = checkCuts[checkCutI];
//            
////            if(!addedCuts[cutI] && !removedCuts[cutI])
////            {
//                addedCuts[cutI] = true;
//                GeometryCut cut = cuts[cutI];
//                labelList cutFaces;
//    //            Info << cut << nl;
//                if(cut.isEdge())
//                {
//                    label cutEdge = cut.geometry();
//                    cutFaces.append(mesh.edgeFaces()[cutEdge]);
//                }
//                else
//                {
//                    label cutPoint = cut.geometry();
//                    cutFaces.append(mesh.pointFaces()[cutPoint]);
//                }
//                
////                Info << "cutFaces " << cutFaces << nl;
//                
//                forAll(cutFaces, faceI)
//                {
//                    label cutFace = cutFaces[faceI];
//                    if(!visitedFaces[cutFace])
//                    {
//                        labelList faceCuts = facesCuts[cutFace];
//                        label nCuts = faceCuts.size();
//                        
////                        Info << nl << "nCuts " << faceI << ": " << nCuts << nl;
//                        
//                        if(nCuts > 1)
//                        {
//                            label startCutI = cutI;
//                            label nextCutI = -1;
//                            DynamicList<label> otherCutsI(faceCuts.size());
//                            
//                            forAll(faceCuts, i)
//                            {
//                                label faceCutI = faceCuts[i];
//                                GeometryCut faceCut = cuts[faceCutI];
////                                Info << "cut " << i << ": (" << faceCutI << ") " << faceCut << nl;
//                                
//                                if(!removedCuts[faceCutI] && faceCutI != startCutI)
//                                {
//                                    otherCutsI.append(faceCutI);
//                                }
//                            }
//    //                        
//    //                        if(startCutI != -1)
//    //                        {
////                            Info << "startCut " << cuts[startCutI] << nl;
//            //                    Info << "nOther " << otherCutsI.size() << nl;
//                            forAll(otherCutsI, otherCutI)
//                            {
//                                label cutI = otherCutsI[otherCutI];
////                                Info << "nextCuts " << cuts[cutI] << nl;
//                            }
//            //                    Info << "otherCutsI " << otherCutsI << nl;
//                            nextCutI = cuts[startCutI].findNext(surf, cuts, otherCutsI);
//    //                        }
//                            
//                            
//                            if(nextCutI != -1)
//                            {
////                                Info << "nextCut  " << cuts[nextCutI] << nl;
//                                addedCuts[nextCutI] = true;
//                                newCheckCuts.append(nextCutI);
//                                
//                                forAll(otherCutsI, i)
//                                {
//                                    label otherCutI = otherCutsI[i];
//                                    if(otherCutI != startCutI && otherCutI != nextCutI)
//                                    {
//                                        removedCuts[cutI] = true;
//                                    }
//                                }
//                            }
//                        }  
//                    }
//                    
//                    visitedFaces[cutFace] = true;
//                }
////            }
//        }
//    }



//    DynamicList<label> allCutPoints(mesh.nPoints());
//    DynamicList<label> allCutEdges(mesh.nEdges());
//    DynamicList<scalar> cutEdgeWeights(mesh.nEdges());

//    forAll(addedCuts, cutI)
//    {
//        bool addCut = addedCuts[cutI];
//        if(addCut)
//        {
//            GeometryCut cut = cuts[cutI];
//            if(cut.isEdge())
//            {
//                allCutEdges.append(cut.geometry());
//                cutEdgeWeights.append(cut.weight());
//            }
//            else
//            {
//                allCutPoints.append(cut.geometry());
//            }
//        }
//    }
//    
//    
    
    
//    
//    allCutPoints.shrink();
//    allCutEdges.shrink();
//    
//    scalarField allCutEdgeWeights;
//    allCutEdgeWeights.transfer(cutEdgeWeights);
//    cutEdgeWeights.clear();
//    
//    
//    OFstream pointStream("cutPoints.obj");
//    
//    forAll(allCutPoints, pointI)
//    {
//        label point = allCutPoints[pointI];
//        meshTools::writeOBJ(pointStream, points[point]);
//    }
//    
//    forAll(allCutEdges, edgeI)
//    {
//        scalar weight = allCutEdgeWeights[edgeI];
//        label edgeL = allCutEdges[edgeI];
//        const edge e = edges[edgeL];
//        const point pStart = points[e.start()];
//        const point pEnd = points[e.end()];
//        const point cutPoint = pStart + (pEnd - pStart) * weight;
//        
//        meshTools::writeOBJ(pointStream, cutPoint);
//    }
//    

    
//    cellCuts cut
//    (
//        mesh,
//        allCutPoints,
//        allCutEdges,
//        allCutEdgeWeights
//    );
////    
//    polyTopoChange meshMod(mesh);

//    // Cutting engine
//    meshCutter cutter(mesh);

//    // Insert mesh refinement into polyTopoChange.
//    cutter.setRefinement(cut, meshMod);

//    // Do all changes
//    Info<< "Morphing ..." << endl;


//    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

//    if (morphMap().hasMotionPoints())
//    {
//        mesh.movePoints(morphMap().preMotionPoints());
//    }

//    // Update stored labels on meshCutter.
//    cutter.updateMesh(morphMap());
       
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
