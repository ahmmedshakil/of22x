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
#include "meshSearch.H"
#include "triSurfaceSearch.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "surfaceFeatures.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "treeDataCell.H"
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

void computeCuts
(
    DynamicList<GeometryCut>& cuts, 
    const polyMesh& mesh, 
    const triSurface& surf
)
{
    GeometryCut::setTriSurface(surf);
    const vectorField& normals = surf.faceNormals();
    

    labelList edgeLabels(mesh.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
        
    triSurfaceSearch querySurf(surf);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

//  Find all cuts
    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];
        const edge e = mesh.edges()[edgeI];
        const point pStart = mesh.points()[e.start()] ;
        const point pEnd = mesh.points()[e.end()] ;
        
        const vector eVec(pEnd - pStart);
        const scalar eMag = mag(eVec);
        const vector n(eVec/(eMag + VSMALL));
        const point tolVec = 1e-6*eVec;
        
        point p0 = pStart - tolVec;
        const point p1 = pEnd + tolVec;
              
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
    cuts.shrink();
    
//  Remove duplicate cuts
    labelList keepCuts;
    uniqueOrder(cuts, keepCuts);
    PackedBoolList cutsToKeep(cuts.size() , false);
    forAll(keepCuts, keepCutsI)
    {
        bool keep = true;
        label keepCut = keepCuts[keepCutsI];
        if(keepCutsI < keepCuts.size() - 1)
        {
            label nextKeepCut = keepCuts[keepCutsI + 1];
            
            if (cuts[keepCut].isEqual(cuts[nextKeepCut]))
                keep = false;
        }
        cutsToKeep[keepCut] = keep;
    }
    inplaceSubset(cutsToKeep, cuts);
}



void writeCuts
(
    const DynamicList<GeometryCut>& cuts, 
    const polyMesh& mesh
)
{
    OFstream cutPointStream("cuts.obj");
    forAll(cuts, cutI)
    {
        point addPoint;
        GeometryCut cut = cuts[cutI];
        if(cut.isPoint())
        {
            label p = cut.geometry();
            addPoint = mesh.points()[p];
        }
        else
        {
            label cutEdge = cut.geometry();
            scalar weight = cut.weight();
            const edge e = mesh.edges()[cutEdge];
            const point pStart = mesh.points()[e.start()];
            const point pEnd = mesh.points()[e.end()];
            addPoint = pStart + (pEnd - pStart) * weight;
        }
        meshTools::writeOBJ(cutPointStream, addPoint);
    }
}


void computeTrianglesPerCell
(
    List<labelHashSet>& cutsPerCell, 
    const polyMesh& mesh, 
    const triSurface& surf
)
{
    labelList edgeLabels(surf.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
    Map<label> map(surf.meshPointMap());
    
    treeBoundBox allBb(mesh.points());
    // Extend domain slightly (also makes it 3D if was 2D)
    scalar bbTol = 1e-6 * allBb.avgDim();

    point& bbMin = allBb.min();
    bbMin.x() -= bbTol;
    bbMin.y() -= bbTol;
    bbMin.z() -= bbTol;

    point& bbMax = allBb.max();
    bbMax.x() += 2*bbTol;
    bbMax.y() += 2*bbTol;
    bbMax.z() += 2*bbTol;

    indexedOctree<treeDataFace> faceTree
    (
        treeDataFace(false, mesh),
        allBb, // overall search domain
        8, // maxLevel
        10, // leafsize
        3.0 // duplicity
    );

    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];
        const edge e = surf.edges()[edgeI];
        const point pStart = surf.localPoints()[e.start()];
        const point pEnd = surf.localPoints()[e.end()];
        const vector eVec(pEnd - pStart);
        const scalar eMag = mag(eVec);
        const vector n(eVec/(eMag + VSMALL));
        const vector tolVec = 1e-6*eVec;
        
        point p0 = pStart - tolVec;
        const point p1 = pEnd + tolVec;
        
        pointIndexHit pHit;
        do
        {
            pHit = faceTree.findLine(p0, p1);
            
            if (pHit.hit())
            {
                const label faceI = pHit.index();
                
                labelList triangles = surf.edgeFaces()[edgeI];
                
                forAll(triangles, i)
                {
                    label triangle = triangles[i];
                    cutsPerCell[mesh.faceOwner()[faceI]].insert(triangle);
                    if (mesh.isInternalFace(faceI))
                    {
                        cutsPerCell[mesh.faceNeighbour()[faceI]].insert(triangle);
                    }
                } 
            
                const vector& area = mesh.faceAreas()[pHit.index()];
                scalar typDim = Foam::sqrt(mag(area));
                if ((mag(pHit.hitPoint() - pEnd)/typDim) < SMALL)
                {
                    break;
                }
                p0 = pHit.hitPoint() + tolVec;
            }
            
        } while (pHit.hit());
    }
    
    
    labelList pointLabels(surf.nPoints());
    forAll(pointLabels, i)
    {
        pointLabels[i] = i;
    }
    
    indexedOctree<treeDataCell> cellTree
    (
        treeDataCell(false, mesh, polyMesh::FACEDIAGTETS),
        allBb, // overall search domain
        8, // maxLevel
        10, // leafsize
        3.0 // duplicity
    ); 
    
    forAll(pointLabels, pointI)
    {
        const point searchPoint = surf.points()[pointI];
        const label cellI = cellTree.findInside(searchPoint);
        labelList triangles = surf.pointFaces()[map[pointI]];
        if (cellI != -1)
        {
            forAll(triangles, i)
            {
                const label triangle = triangles[i];
                cutsPerCell[cellI].insert(triangle);
            }
        }
    }
}

void agglomerateTriangles(
    const labelList& triangles, 
    const triSurface& surf, 
    labelList& agglomeration
)
{
    List<labelHashSet> surfaces(triangles.size());
    List<labelHashSet> points(triangles.size());
    
    forAll(triangles, triangleI)
    {
        labelList trianglePoints(surf[triangles[triangleI]]);
        
        surfaces[triangleI].insert(triangleI);
        points[triangleI].insert(trianglePoints);
    }
    
    bool foundConnected;
    
    do
    {
        foundConnected = false;
        for (int i = 0; i < surfaces.size(); i++)
        {
            if (surfaces[i].size() > 0)
            {
                for (int j = i + 1; j < surfaces.size(); j++)
                {
                    if ((points[i] & points[j]).size() > 0)
                    {
                        points[i] += points[j];
                        points[j].clear();
                        surfaces[i] += surfaces[j];
                        surfaces[j].clear();
                        foundConnected = true;
                    }
                }
            }
        }
    }
    while (foundConnected);
    
    
    label group = 0;
    forAll(surfaces, surfaceI)
    {
        labelList triangles = surfaces[surfaceI].toc();
        if (triangles.size() > 0)
        {
            forAll(triangles, triangleI)
            {
                label triangle = triangles[triangleI];
                agglomeration[triangle] = group;
            }
            group++;
        }
    }
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
    
    DynamicList<GeometryCut> cuts(mesh.nPoints()*4);
    computeCuts(cuts, mesh, surf);
    writeCuts(cuts, mesh);

    
    
    
//    //  find facesCuts
    List<DynamicList<label> > facesCuts(mesh.nFaces());
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
    
    DynamicList<label> cutFaces(mesh.nFaces());
    forAll(facesCuts, faceI)
    {
        label nFaceCuts = facesCuts[faceI].size();
        if (nFaceCuts > 1)
        {
            cutFaces.append(faceI);
        }
    }
    cutFaces.shrink();
//    Info << "cutFaces " << cutFaces << nl;

    List<DynamicList<label> > cellsCutFaces(mesh.nCells());
    forAll(cutFaces, i)
    {
        label faceI = cutFaces[i];
        label owner = mesh.faceOwner()[faceI];
        cellsCutFaces[owner].append(faceI);
        
        if (mesh.isInternalFace(faceI))
        {
            label neighbour = mesh.faceNeighbour()[faceI];
            cellsCutFaces[neighbour].append(faceI);
        }
    }
//    Info << "cellsCutFaces " << cellsCutFaces << nl;
    
    DynamicList<label> cutCells(mesh.nCells());
    forAll(cellsCutFaces, cellI)
    {
        label nCellCuts = cellsCutFaces[cellI].size();
        
        
        if (nCellCuts > 2)
        {
            cutCells.append(cellI);
        }
    }
    cutCells.shrink();
//    Info << "cutCells " << cutCells << nl;
    
    
    
    List<labelHashSet> trianglesPerCell(mesh.nCells());
    computeTrianglesPerCell(trianglesPerCell, mesh, surf);
    
    List<DynamicList<DynamicList<label> > > cutsPerCell(mesh.nCells());
    
    forAll(trianglesPerCell, cellI)
    {
        labelList triangles = trianglesPerCell[cellI].toc();
        if (triangles.size() > 0)
        {
            labelList agglomeration(triangles.size());
            agglomerateTriangles(triangles, surf, agglomeration);
            Info << agglomeration << nl << nl;
        }
    }
    
    
//    forAll(cutsPerCell, i)
//    {
//        Info << cutsPerCell[i].toc() << nl;
//    }
    
    
//    Map<label> map(surf.meshPointMap());

//    forAll(cutsPerCell, cellI)
//    {
//        
//        labelList faces = cutsPerCell[cellI].toc();
//        if (faces.size() > 0)
//        {
//            Info << name(cellI) << " " << faces << nl;
//            
//            OFstream cutPointStream("cell" + name(cellI) + ".obj");
//            
//            
//            forAll(faces, i)
//            {
//                labelList points(surf[faces[i]]);
//                forAll(points, i)
//                {
//                    meshTools::writeOBJ(
//                        cutPointStream, 
//                        surf.points()[points[i]]
//                    );
//                }
//            }
//        }
//    }
    
//    forAll(cutsPerCell, cellI)
//    {
//        
//        labelList points = cutsPerCell[cellI].toc();
//        if (points.size() > 0)
//        {
//            Info << points << nl;
//            
//            OFstream cutPointStream("cell" + name(cellI) + ".obj");
//            
//            forAll(points, i)
//            {
//                meshTools::writeOBJ(
//                    cutPointStream, 
//                    surf.points()[points[i]]
//                );
//            }
//        }
//    }
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
