/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistrianglebute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distrianglebuted in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cutSearcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Static Functions * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
        
//- Contruct from edge with edge weight
CutSearcher::CutSearcher
(
    const polyMesh& mesh,
    const triSurface& surf
) 
:
    mesh_(mesh),
    surf_(surf),
    cuts_(mesh.nPoints()*4),
    trianglesPerCell_(mesh.nCells()),
    alignedCos_(degToRad(89.0))
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

CutSearcher::~CutSearcher()
{
}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void CutSearcher::computeCuts()
{
    GeometryCut::setTriSurface(surf_);
    const vectorField& normals = surf_.faceNormals();
    

    labelList edgeLabels(mesh_.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
        
    triSurfaceSearch querySurf(surf_);
    const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

//  Find all cuts
    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];
        const edge e = mesh_.edges()[edgeI];
        const point pStart = mesh_.points()[e.start()] ;
        const point pEnd = mesh_.points()[e.end()] ;
        
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
                    cuts_.append(newCut);
                }
                else if (mag(pHit.hitPoint() - pEnd) < 0.01 * eMag)
                {
                    const label endPoint = e.end();
                    GeometryCut newCut(endPoint, pHit.index());
                    cuts_.append(newCut);
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
                    cuts_.append(newCut);
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
    cuts_.shrink();
    
//  Remove duplicate cuts
    labelList keepCuts;
    uniqueOrder(cuts_, keepCuts);
    PackedBoolList cutsToKeep(cuts_.size() , false);
    forAll(keepCuts, keepCutsI)
    {
        bool keep = true;
        label keepCut = keepCuts[keepCutsI];
        if(keepCutsI < keepCuts.size() - 1)
        {
            label nextKeepCut = keepCuts[keepCutsI + 1];
            
            if (cuts_[keepCut].isEqual(cuts_[nextKeepCut]))
                keep = false;
        }
        cutsToKeep[keepCut] = keep;
    }
    inplaceSubset(cutsToKeep, cuts_);
}



void CutSearcher::writeCuts()
{
    OFstream cutPointStream("cuts.obj");
    forAll(cuts_, cutI)
    {
        point addPoint;
        GeometryCut cut = cuts_[cutI];
        if(cut.isPoint())
        {
            label p = cut.geometry();
            addPoint = mesh_.points()[p];
        }
        else
        {
            label cutEdge = cut.geometry();
            scalar weight = cut.weight();
            const edge e = mesh_.edges()[cutEdge];
            const point pStart = mesh_.points()[e.start()];
            const point pEnd = mesh_.points()[e.end()];
            addPoint = pStart + (pEnd - pStart) * weight;
        }
        meshTools::writeOBJ(cutPointStream, addPoint);
    }
}


void CutSearcher::computeTrianglesPerCell()
{
    labelList edgeLabels(surf_.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
    Map<label> map(surf_.meshPointMap());
    
    treeBoundBox allBb(mesh_.points());
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
        treeDataFace(false, mesh_),
        allBb, // overall search domain
        8, // maxLevel
        10, // leafsize
        3.0 // duplicity
    );

    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];
        const edge e = surf_.edges()[edgeI];
        const point pStart = surf_.localPoints()[e.start()];
        const point pEnd = surf_.localPoints()[e.end()];
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
                
                labelList triangles = surf_.edgeFaces()[edgeI];
                
                forAll(triangles, i)
                {
                    label triangle = triangles[i];
                    trianglesPerCell_[mesh_.faceOwner()[faceI]].insert(triangle);
                    if (mesh_.isInternalFace(faceI))
                    {
                        trianglesPerCell_[mesh_.faceNeighbour()[faceI]].insert(triangle);
                    }
                } 
            
                const vector& area = mesh_.faceAreas()[pHit.index()];
                scalar typDim = Foam::sqrt(mag(area));
                if ((mag(pHit.hitPoint() - pEnd)/typDim) < SMALL)
                {
                    break;
                }
                p0 = pHit.hitPoint() + tolVec;
            }
            
        } while (pHit.hit());
    }
    
    
    labelList pointLabels(surf_.nPoints());
    forAll(pointLabels, i)
    {
        pointLabels[i] = i;
    }
    
    indexedOctree<treeDataCell> cellTree
    (
        treeDataCell(false, mesh_, polyMesh::FACEDIAGTETS),
        allBb, // overall search domain
        8, // maxLevel
        10, // leafsize
        3.0 // duplicity
    ); 
    
    forAll(pointLabels, pointI)
    {
        const point searchPoint = surf_.points()[pointI];
        const label cellI = cellTree.findInside(searchPoint);
        labelList triangles = surf_.pointFaces()[map[pointI]];
        if (cellI != -1)
        {
            forAll(triangles, i)
            {
                const label triangle = triangles[i];
                trianglesPerCell_[cellI].insert(triangle);
            }
        }
    }
}

//void agglomerateTriangles()
//{
//    List<labelHashSet> surfaces(triangles.size());
//    List<labelHashSet> points(triangles.size());
//    
//    forAll(triangles, triangleI)
//    {
//        labelList trianglePoints(surf[triangles[triangleI]]);
//        
//        surfaces[triangleI].insert(triangleI);
//        points[triangleI].insert(trianglePoints);
//    }
//    
//    bool foundConnected;
//    
//    do
//    {
//        foundConnected = false;
//        for (int i = 0; i < surfaces.size(); i++)
//        {
//            if (surfaces[i].size() > 0)
//            {
//                for (int j = i + 1; j < surfaces.size(); j++)
//                {
//                    if ((points[i] & points[j]).size() > 0)
//                    {
//                        points[i] += points[j];
//                        points[j].clear();
//                        surfaces[i] += surfaces[j];
//                        surfaces[j].clear();
//                        foundConnected = true;
//                    }
//                }
//            }
//        }
//    }
//    while (foundConnected);
//    
//    
//    label group = 0;
//    forAll(surfaces, surfaceI)
//    {
//        labelList triangles = surfaces[surfaceI].toc();
//        if (triangles.size() > 0)
//        {
//            forAll(triangles, triangleI)
//            {
//                label triangle = triangles[triangleI];
//                agglomeration[triangle] = group;
//            }
//            group++;
//        }
//    }
//}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

