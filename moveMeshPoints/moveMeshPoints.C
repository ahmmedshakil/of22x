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
    argList::validArgs.append("tolerance");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");
    const fileName surfName = args[1];
    const scalar tol = args.argRead<scalar>(2);
    
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
    
    Info << "Find points near triSurface" << endl << endl;
    
    DynamicList<label> movePoints(mesh.nPoints());
    DynamicList<point> newLocations(mesh.nPoints());
    
    
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
        const scalar eTol = tol * eMag;
            
        bool foundStart = false;
        bool foundEnd = false;
        
        pointIndexHit nearHitStart;
        pointIndexHit nearHitEnd;
                
        while(true)
        {
            pointIndexHit pHit = tree.findLine(p0, p1);
            
            if(pHit.hit())
            {
                
                if (mag(pHit.hitPoint() - pStart) < eTol && !foundStart)
                {
                    nearHitStart = tree.findNearest(pStart, eTol);
                    foundStart = true;
                }
                else if (mag(pHit.hitPoint() - pEnd) < eTol)
                {
                    nearHitEnd = tree.findNearest(pEnd, eTol);
                    foundEnd = true;
                }
                
                p0 = pHit.hitPoint() + tolVec;
            }
            else
            {
                // No hit.
                break;
            }
        }
        
        if(foundStart)
        {
            movePoints.append(e.start());
            newLocations.append(nearHitStart.hitPoint());
        }
        
        if(foundEnd)
        {
            movePoints.append(e.end());
            newLocations.append(nearHitEnd.hitPoint());
        }
        
    }
    
    movePoints.shrink();
    newLocations.shrink();
    
    Info<< "Moving " << movePoints.size() << " points"  << endl << endl;
    
    polyTopoChange meshMod(mesh);
    forAll(movePoints, i)
    {
        label movePointI = movePoints[i];
        point newLocationI = newLocations[i];
        
        meshMod.modifyPoint(movePointI, newLocationI, -1, true);
    }
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (!overwrite)
    {
        runTime++;
    }
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing morphMesh to time " << runTime.timeName() << endl << endl;

    mesh.write(); 
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
