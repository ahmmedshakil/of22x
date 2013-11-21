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
#include "surfaceFeatures.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "triangle.H"
#include "polyTopoChange.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
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
        
        triSurfaceTools::otherEdges(surf, lastTri, thisEdge, lastEdge1, lastEdge2);
        
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
//    const scalar featureAngle = args.argRead<scalar>(2);
    
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
    triSurfaceSearch searchSurf(surf);
    
    
    
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
    
    label nHits = 0;
    
    
    DynamicList<label> allCutEdges(mesh.nEdges()); 
//    DynamicList<label> allCutEdges(mesh.nEdges()); 
    List<bool> pointCuts(mesh.nPoints(), false); 
    List<List<pointIndexHit> > edgeIntersections(mesh.nEdges());
    List<List<pointIndexHit> > pointIntersections(mesh.nPoints());
    
    
    List<List<label> > intersectionTypes(mesh.nEdges());
    
    Info << "find hits" << nl;
    
    
    
    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];
        const edge e = edges[edgeI];
        const point pStart = points[e.start()] ;
        const point pEnd = points[e.end()] ;
        
        bool isCut = false;
        
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
                
                if (mag(pHit.hitPoint() - pStart) < 0.1 * eMag)
                {
                    const label startPoint = e.start();
                    pointIntersections[startPoint].append(pHit);
//                    allCutPoints.insert(startPoint);
                    pointCuts[startPoint] = true;
                }
                else if (mag(pHit.hitPoint() - pEnd) < 0.1 * eMag)
                {
                    const label endPoint = e.end();
                    pointIntersections[endPoint].append(pHit);
//                    allCutPoints.insert(endPoint);
                    pointCuts[endPoint] = true;
                    break;
                }
//                else if (mag(n & normals[pHit.index()]) < alignedCos_)
//                {
//                    edgeEnd = 2;
//                }
                else
                {
                    edgeIntersections[edgeI].append(pHit);
                    isCut = true;
                }
                p0 = pHit.hitPoint() + tolVec;
            }
            else
            {
                // No hit.
                break;
            }
        }
        
        if(isCut)
        {
            allCutEdges.append(edgeI);
        }
    }
    
    DynamicList<label> cutPoints(mesh.nPoints()); 
    
    forAll(pointCuts, i)
    {
        bool cutPointI = pointCuts[i];
        if(cutPointI)
        {
            cutPoints.append(i);
        }
    }
    
    
    
    cutPoints.shrink();
    List<label> allCutPoints;
    allCutPoints.transfer(cutPoints);
    cutPoints.clear();
    
    
    Info << "ready" << nl;
    




    List<label> nFaceCuts(mesh.nFaces(), 0);
    
    
//    List<bool> cutFaces(mesh.nFaces(), false);
    
    labelListList cutFaceEdges(mesh.nFaces());    
    forAll(allCutEdges, i)
    {
        const label cutEdgeI = allCutEdges[i];
        const labelList cutEdgeFaces = mesh.edgeFaces()[cutEdgeI];
        
        forAll(cutEdgeFaces, i)
        {
            const label cutEdgeFaceI = cutEdgeFaces[i];
            nFaceCuts[cutEdgeFaceI]++;
            
            cutFaceEdges[cutEdgeFaceI].append(cutEdgeI);
        }
    }
    
    
    
    
    
    labelListList cutFacePoints(mesh.nFaces());
    forAll(allCutPoints, i)
    {
        const label cutPointI = allCutPoints[i];
        const labelList cutPointFaces = mesh.pointFaces()[cutPointI];
        
        forAll(cutPointFaces, i)
        {
            const label cutPointFaceI = cutPointFaces[i];
            nFaceCuts[cutPointFaceI]++;
            
            cutFacePoints[cutPointFaceI].append(cutPointI);
        }
    }
    
    
    DynamicList<label> cutFaces(mesh.nFaces());    
    forAll(nFaceCuts, i)
    {
        label nFaceCutsI = nFaceCuts[i];
        
        if(nFaceCutsI >= 2)
        {
            cutFaces.append(i);
        }
    }
    
    cutFaces.shrink();
    List<label> allCutFaces;
    allCutFaces.transfer(cutFaces);
    cutFaces.clear();
    
    
    
    const labelList owner = mesh.faceOwner();
    const labelList neighbour = mesh.faceNeighbour();
    const label nNeighbour = neighbour.size();
    
    labelListList cutCellFaces(mesh.nCells());
    List<label> nCellCuts(mesh.nCells(), 0);
    forAll(allCutFaces, i)
    {
        const label cutFaceI = allCutFaces[i];
        const label cutFaceOwner = owner[cutFaceI];
        
        nCellCuts[cutFaceOwner]++;
        cutCellFaces[cutFaceOwner].append(cutFaceI);
        
        if(mesh.isInternalFace(cutFaceI))
        {
            const label cutFaceNeighbour = neighbour[cutFaceI];
            nCellCuts[cutFaceNeighbour]++;
            cutCellFaces[cutFaceNeighbour].append(cutFaceI);
        }
    }
    
    DynamicList<label> cutCells(mesh.nCells());    
    forAll(nCellCuts, i)
    {
        label nCellCutsI = nCellCuts[i];
        
        if(nCellCutsI >= 3)
        {
            cutCells.append(i);
        }
    }
    
    cutCells.shrink();
    
    List<label> allCutCells;
    allCutCells.transfer(cutCells);
    cutCells.clear();
    
//    scalarField cutEdgeWeights;
//    cutEdgeWeights.transfer(allCutEdgeWeights);
//    allCutEdgeWeights.clear();
    
    
    Info << "allCutEdges " << allCutEdges << nl;
    Info << "edgeIntersections " << edgeIntersections << nl;
    
    Info << "allCutPoints " << allCutPoints << nl;
    Info << "pointIntersections " << pointIntersections << nl;
    
    Info << "allCutFaces " << allCutFaces << nl;
    Info << "cutFacePoints " << cutFacePoints << nl;
    Info << "cutFaceEdges " << cutFaceEdges << nl;
    
    Info << "allCutCells " << allCutCells << nl;
    Info << "cutCellFaces " << cutCellFaces << nl;
    
    
    
    
//    const labelList owner = mesh.faceOwner();
//    const labelList neighbour = mesh.faceNeighbour();
//    const label nNeighbour = neighbour.size();
//    
//    List<label> nCellsCutFaces(mesh.nCells(), 0);
//    labelHashSet allCutCells(mesh.nCells());
//    labelListList cellsCutFaces(mesh.nCells());
    
//    labelListList alll = allCutFaces.toc();
//    
//    labelList cutFaces = allCutFaces.toc();

    
//    allCutFaces.clearStorage();
//    
//    
//    Info << "allCutFaces : " << allCutFaces << nl;
//    
//    forAll(allCutFaces, i)
//    {
//        const label cutFaceI = allCutFaces.toc()[i];
//        const label cutFaceOwner = owner[cutFaceI];
//        
//        nCellsCutFaces[cutFaceOwner]++;
//        cellsCutFaces[cutFaceOwner].append(cutFaceI);
//        if(nCellsCutFaces[cutFaceOwner] >= 3)
//        {
//            allCutCells.insert(cutFaceOwner);
//        }
//        
//        
//        if(cutFaceI < nNeighbour)
//        {
//            const label cutFaceNeighbour = neighbour[cutFaceI];
//            nCellsCutFaces[cutFaceNeighbour]++;
//            cellsCutFaces[cutFaceNeighbour].append(cutFaceI);
//            if(nCellsCutFaces[cutFaceNeighbour] >= 3)
//            {
//                allCutCells.insert(cutFaceNeighbour);
//            }
//        }
//    }
//    
//    allCutCells.shrink();
    
//    Info  << cellsCutFaces << nl;
    
    
    
    
//    label startTriFace = 9;
//    labelHashSet nextTriFaces(3);
//    
//    nextTriFaces.insert(1052);
//    nextTriFaces.insert(1086);
//    nextTriFaces.insert(805);
//    
//    label next = findNextIntersection(surf, startTriFace, nextTriFaces);
//    
//    Info << "Returned " << next << nl;

    
    
    forAll(allCutCells, i)
    {
        label cutCellI = allCutCells[i];
        labelList cutFaces = cutCellFaces[cutCellI];
        
        forAll(cutFaces, i)
        {
            label firstFace = cutFaces[i];
            labelList cutEdges = cutFaceEdges[firstFace];
            label firstEdge = cutEdges[0];
            List<pointIndexHit> cutEdgePoints = edgeIntersections[firstEdge];
        
            Info <<  nl;
                
            forAll(cutEdgePoints, i)
            {
                pointIndexHit firstPoint = cutEdgePoints[i];
                Info << "this : " << firstPoint << nl;
            }
            
            if(cutEdges.size() >= 2)
            {
                for(int k = 1; k < cutEdges.size(); k++)
                {
                    label nextEdge = cutEdges[k];
                    
                    List<pointIndexHit> nextPoints = edgeIntersections[nextEdge];
                    
                    forAll(nextPoints, i)
                    {
                        pointIndexHit nextPoint = nextPoints[i];
                        Info << "next : " << nextPoint << nl;
                    }
                }
            }
            else
            {
                Info << "only one Edge" << nl;
            }
        }
    }
    
    
    forAll(allCutCells, i)
    {
        label cutCellI = allCutCells[i];
        labelList cutFaces = cutCellFaces[cutCellI];
        label thisFace = cutFaces[0];
        label thisEdge = cutFaceEdges[thisFace][0];
        
        forAll(cutFaces, i)
        {
            label thisFace = cutFaces[i];
            labelList cutEdges = cutFaceEdges[thisFace];
            labelList cutPoints = cutFacePoints[thisFace];
            DynamicList<pointIndexHit> intersections;
            DynamicList<label> nextEdges;
            
            forAll(cutEdges, i)
            {
                label cutEdgeI = cutEdges[i];
                intersections.append(edgeIntersections[cutEdgeI]);
            }
            forAll(cutPoints, i)
            {
                label cutPointI = cutPoints[i];
                intersections.append(pointIntersections[cutPointI]);
            }
            
            
            
            
            Info << intersections << nl;
            
            
            
            if(intersections.size() == 2)
            {
                
            }
        }
    }
    

    


    
    
//    Info << "allCutEdges : " << allCutEdges << nl;
//    Info << "allCutCells : " << allCutCells << nl;
    
//    const faceList& faces = mesh.faces();
    
//    labelList nCellsCutFaces(mesh.nCells(), 0);
//    labelListList cellsCutFaces(mesh.nCells());
    
    
    
    
//    Info << "owner : " << owner << nl;
//    Info << neighbour << nl;
    
//    for(int i = 0; i < mesh.nFaces(); i++)
//    {
//    
////        Info << owner[i] << neighbour[i] << nl;

//        if(nFacesCutEdges[i] >= 2)
//            {
//            nCellsCutFaces[owner[i]]++;
//            if(i < nNeighbour)
//            { 
//                nCellsCutFaces[neighbour[i]]++;
//            }
//        }
//    
//    }
    
//    forAll(
    
    
    
    
//    Info << "nCellsCutFaces : " << nCellsCutFaces << nl;
//    Info << "nCellsCutFaces : " << nCellsCutFaces << nl;
    
    
    
    
    
//    Info << "nFacesCutEdges : " << nFacesCutEdges << nl;
//    Info << "FacesCutEdges : " << FacesCutEdges << nl;
    
    
//    labelListList* faceCutsPtr_;
//    const faceList& faces = mesh.faces();
//    faceCutsPtr_ = new labelListList(mesh.nFaces());
//    labelListList& faceCuts = *faceCutsPtr_;

    
    
//    polyTopoChange meshMod(mesh);
//    
//    point pc1(0,    0,    0.05);
//    point pc2(0.1,  0,    0.05);
//    point pc3(0.1,  0.1,  0.05);
//    point pc4(0,    0.1,  0.05);
//    
//    
//    label p1 = meshMod.addPoint(pc1, -1, -1, -1);
//    label p2 = meshMod.addPoint(pc2, -1, -1, -1);
//    label p3 = meshMod.addPoint(pc3, -1, -1, -1);
//    label p4 = meshMod.addPoint(pc4, -1, -1, -1);
//    
//    
//    
////    label fc1[4] = {p1,p2,p3,p4};
//    
//    labelList fc1(4);
//    
//    fc1[0] = p1;
//    fc1[1] = p2;
//    fc1[2] = p3;
//    fc1[3] = p4;
//    
//    face fcc1(fc1);
//    
//    label newCell = meshMod.addCell(-1,-1,-1,-1,-1);
//    
//    label f1 = meshMod.addFace(fcc1,0,newCell,-1,-1,-1,false,-1,-1,false);
//    
//    
//    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);


    
    
//    meshMod.addPoint
    
    
    
    cellCuts cuts
    (
        mesh,
        surf,
        allCuts,
        cutTris,
        cutEdges,
        cutPoints,
        

        allCutPoints,       // cut vertices
        allCutEdges,        // cut edges
        cutEdgeWeights,     // weight on cut edges
        cutTris
    );
    
    
    

    // Transfer DynamicLists to straight ones.
//    scalarField cutEdgeWeights;
//    cutEdgeWeights.transfer(allCutEdgeWeights);
//    allCutEdgeWeights.clear();

//     Gets cuts across cells from cuts through edges.
//    DynamicList<scalar> allcutEdgeWeights(3);
//    
//    allcutEdgeWeights.append(0.5);
//    allcutEdgeWeights.append(0.5);
//    allcutEdgeWeights.append(0.5);
//    allcutEdgeWeights.append(0.5);
//    
//    scalarField cutEdgeWeights;    
//    cutEdgeWeights.transfer(allcutEdgeWeights);


//    cellCuts cuts
//    (
//        mesh,
//        allCutPoints.toc(),       // cut vertices
//        allCutEdges,        // cut edges
//        cutEdgeWeights      // weight on cut edges
//    );

//    polyTopoChange meshMod(mesh);

//    // Cutting engine
//    meshCutter cutter(mesh);

//    // Insert mesh refinement into polyTopoChange.
//    cutter.setRefinement(cuts, meshMod);

//    // Do all changes
//    Info<< "Morphing ..." << endl;


//    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

//    if (morphMap().hasMotionPoints())
//    {
//        mesh.movePoints(morphMap().preMotionPoints());
//    }

//    // Update stored labels on meshCutter.
//    cutter.updateMesh(morphMap());
//   
       
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
