/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    refineWallLayer

Description
    Utility to refine cells next to patches.

    Takes a patchName and number of layers to refine. Works out cells within
    these layers and refines those in the wall-normal direction.

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
#include "cellSet.H"
#include "fvMeshSubset.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H" 
    argList::noParallel();
    argList::validArgs.append("input surfaceFile");
    argList::validArgs.append("feature angle");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
//#   include "createPolyMesh.H"
    #include "createNamedMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");
    const fileName surfName = args[1];
    const scalar featureAngle = args.argRead<scalar>(2);
    
    triSurface surf(surfName);
    surfaceFeatures set(surf, featureAngle);
    boolList borderEdge(surf.nEdges(), false);
    forAll(set.featureEdges(), i)
    {
        borderEdge[set.featureEdges()[i]] = true;
    }
    labelList faceRegion(surf.size());
    label nRegions = surf.markZones(borderEdge, faceRegion);
    forAll(surf, i)
    {
        surf[i].region() = faceRegion[i];
    }
    
    
    
    
    for( int regionI = 0; regionI < nRegions; regionI++)
    {
    
    
        pointField points = mesh.points();
        labelList edgeLabels(mesh.nEdges());
        forAll(edgeLabels, i)
        {
            edgeLabels[i] = i;
        }
        

        boolList includeMap(surf.size(), false);
        forAll(surf, faceI)
        {
            const labelledTri& f = surf[faceI];
            
            
            if (f.region() == regionI)
            {
                includeMap[faceI] = true;
            }
        }
        
        
        labelList pointMap;
        labelList faceMap;
        
        
        triSurface subSurf(
            surf.subsetMesh
            (
                includeMap,
                pointMap,
                faceMap
            )
        );
        
        
        
        
        
        triSurfaceSearch querySurf(subSurf);
        triSurfaceSearch searchSurf(subSurf);
        
        
        
        const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
        
        label nHits = 0;
        
        
        DynamicList<scalar> allCutEdgeWeights(mesh.nEdges());
        DynamicList<label> allCutEdges(mesh.nEdges()); 
        
        
        forAll(edgeLabels, i)
        {
            label edgeI = edgeLabels[i];
            const edge& e = mesh.edges()[edgeI];
            const point& pStart = points[e.start()] ;
            const point& pEnd = points[e.end()] ;
        
            pointIndexHit pHit = tree.findLine(pStart, pEnd);
            
            if(pHit.hit())
            {
                allCutEdges.append(edgeI);
                
                const vector eVec(pEnd - pStart);
                const vector pVec(pHit.hitPoint() - pStart);
                const scalar weight = mag(pVec)/mag(eVec);
                allCutEdgeWeights.append(weight);
                
                
                
            }
        }
          
          
        
        
        
        
        
        
        


        

        // Transfer DynamicLists to straight ones.
        scalarField cutEdgeWeights;
        cutEdgeWeights.transfer(allCutEdgeWeights);
        allCutEdgeWeights.clear();

        // Gets cuts across cells from cuts through edges.
        cellCuts cuts
        (
            mesh,
            labelList(0),       // cut vertices
            allCutEdges,        // cut edges
            cutEdgeWeights      // weight on cut edges
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
            
    }
    
    
    cellSet cells(mesh, "smallCells", mesh.nCells());
    

    List<scalar> volumes = mesh.cellVolumes();
    forAll(volumes, i)
    {
        scalar volumeI = volumes[i];
        if(volumeI > 1e-9)
        {
            cells.insert(i);
        }
    }
    
//    cells.instance() = oldInstance;
//    cells.write();
    
    
    if (!overwrite)
    {
        runTime++;
    }
    
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info<< "Writing refined morphMesh to time " << runTime.timeName() << endl;

    mesh.write(); 
    
    
    fvMeshSubset subsetter(mesh);
    
    subsetter.setLargeCellSubset(cells, -1, true);        
    subsetter.subMesh().setInstance(oldInstance);    
    subsetter.subMesh().write();
    
    
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
