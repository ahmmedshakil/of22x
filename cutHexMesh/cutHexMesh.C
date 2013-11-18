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
#include "triSurface.H"
#include "unitConversion.H"


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
    
    
//    DynamicList<scalar> allCutEdgeWeights;
    DynamicList<label> allCutEdges(mesh.nEdges()); 
    labelHashSet allCutPoints(mesh.nPoints()); 
    List<List<pointIndexHit> > intersections(mesh.nEdges());
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
    

        DynamicList<pointIndexHit> currentIntersections(100);        
        DynamicList<label> currentIntersectionTypes(100);
                
        while(true)
        {
            pointIndexHit pHit = tree.findLine(p0, p1);
            
//            Info << pHit << nl;
            
            if(pHit.hit())
            {
                
                if (mag(pHit.hitPoint() - pStart) < 0.1 * eMag)
                {
                    const label startPoint = e.start();
                    pointIntersections[startPoint].append(pHit);
                    allCutPoints.insert(startPoint);
                }
                else if (mag(pHit.hitPoint() - pEnd) < 0.1 * eMag)
                {
                    const label endPoint = e.end();
                    pointIntersections[endPoint].append(pHit);
                    allCutPoints.insert(endPoint);
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
                    
//                    currentIntersectionTypes.append(edgeEnd);
                }
                
                
//                currentIntersections.append(pHit);
                
                
                
//                if (edgeEnd == 1)
//                {
//                    // Close to end
//                    break;
//                }            
//                else
//                {
                    // Continue tracking. Shift by a small amount.
                    p0 = pHit.hitPoint() + tolVec;
//                }
                
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
        
        
//        intersections[edgeI].transfer(currentIntersections);
//        intersectionTypes[edgeI].transfer(currentIntersectionTypes);
    }
    
    Info << "ready" << nl;
    
    Info << edgeIntersections << nl;
    Info << pointIntersections << nl;
    Info << allCutEdges << nl;
    Info << allCutPoints.toc() << nl;
    
    
//    forAll(allCutEdges, i)
//    {
//        const label edgeI = allCutEdges[i];
//        const labelList intersectionTypesI = intersectionTypes[edgeI];
//        
//        forAll(intersectionTypesI, i)
//        {
//            const label intersectionType = intersectionTypesI[i];
//            
//            Info << "Type of " << intersections[edgeI][i] << " is " << intersectionType << nl;
//            
//            if(intersectionType == 0)
//            {
//                
//            }
//            else if(intersectionType == 1)
//            {
////                intersections[edgeI][i] = NULL;
//            }
//        } 
//    } 
//    
    
    
    
    
    
    
    
    
    
    
    List<label> nFacesCutEdges(mesh.nFaces(), 0);
    
    
    labelHashSet allCutFaces(mesh.nFaces());
    
    labelListList facesCutEdges(mesh.nFaces());
    
    forAll(allCutEdges, i)
    {
        const label cutEdgeI = allCutEdges[i];
        const labelList cutEdgeFaces = mesh.edgeFaces()[cutEdgeI];
        
        forAll(cutEdgeFaces, i)
        {
            const label cutEdgeFaceI = cutEdgeFaces[i];
            nFacesCutEdges[cutEdgeFaceI]++;
            if(nFacesCutEdges[cutEdgeFaceI] >= 2)
            {
                allCutFaces.insert(cutEdgeFaceI);
            }
            
            facesCutEdges[cutEdgeFaceI].append(cutEdgeI);
        }
    }
    
    const labelList owner = mesh.faceOwner();
    const labelList neighbour = mesh.faceNeighbour();
    const label nNeighbour = neighbour.size();
    
    List<label> nCellsCutFaces(mesh.nCells(), 0);
    labelHashSet allCutCells(mesh.nCells());
    labelListList cellsCutFaces(mesh.nCells());
    
//    labelListList alll = allCutFaces.toc();
    
//    labelList cutFaces = allCutFaces.toc();

    allCutFaces.shrink();
    
//    allCutFaces.clearStorage();
    
    
//    Info << "allCutFaces : " << allCutFaces << nl;
    
    forAll(allCutFaces, i)
    {
        const label cutFaceI = allCutFaces.toc()[i];
        const label cutFaceOwner = owner[cutFaceI];
        
        nCellsCutFaces[cutFaceOwner]++;
        cellsCutFaces[cutFaceOwner].append(cutFaceI);
        if(nCellsCutFaces[cutFaceOwner] >= 3)
        {
            allCutCells.insert(cutFaceOwner);
        }
        
        
        if(cutFaceI < nNeighbour)
        {
            const label cutFaceNeighbour = neighbour[cutFaceI];
            nCellsCutFaces[cutFaceNeighbour]++;
            cellsCutFaces[cutFaceNeighbour].append(cutFaceI);
            if(nCellsCutFaces[cutFaceNeighbour] >= 3)
            {
                allCutCells.insert(cutFaceNeighbour);
            }
        }
    }
    
    allCutCells.shrink();
    
//    Info  << cellsCutFaces << nl;
    
    
    forAll(allCutCells, i)
    {
        label cutCellI = allCutCells.toc()[i];
        labelList cutFaces = cellsCutFaces[cutCellI];
        
        forAll(cutFaces, i)
        {
            label firstFace = cutFaces[i];
            
            labelList cutEdges = facesCutEdges[firstFace];
            label firstEdge = cutEdges[0];
            
            List<pointIndexHit> cutPoints = intersections[firstEdge];
            
                Info <<  nl;
                
            forAll(cutPoints, i)
            {
                pointIndexHit firstPoint = cutPoints[i];
                Info << "this : " << firstPoint << nl;
            }
            
            
            
            if(cutEdges.size() >= 2)
            {
                for(int k = 1; k < cutEdges.size(); k++)
                {
                    label nextEdge = cutEdges[k];
                    
                    List<pointIndexHit> nextPoints = intersections[nextEdge];
                    
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
