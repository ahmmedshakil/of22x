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
    moveMeshPoints

Description
    Utility to move Points of a Mesh onto a near face.

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
