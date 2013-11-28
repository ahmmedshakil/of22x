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

#include "geometryCut.H"
#include "triSurface.H"
#include "triSurfaceTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * Private Static Functions * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
        
//- Contruct from edge with edge weight
GeometryCut::GeometryCut
(
    const label geometry,
    const label triangle,
    const scalar weight
) 
:
    geometry_(geometry),
    triangle_(triangle),
    weight_(weight)
{}

//- Construct from point
GeometryCut::GeometryCut
(
    const label geometry,
    const label triangle
) 
:
    geometry_(geometry),
    triangle_(triangle),
    weight_(-1)
{}

//- Construct with default values
GeometryCut::GeometryCut() 
: 
    geometry_(-1),
    triangle_(-1),
    weight_(-1)
{}



// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

GeometryCut::~GeometryCut()
{
}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //
        

label GeometryCut::findNext(
    const triSurface surf,
    const List<GeometryCut> cuts,
    const labelList& otherCuts
)
{
    label nRuns = 100000;
    label nextCut = -1;
    labelHashSet visited(nRuns);
    label startTriangle = triangle_;
    
    labelList faceEdges = surf.faceEdges()[startTriangle];
    label thisEdge =  faceEdges[1];
    label lastEdge1 = faceEdges[0];
    label lastEdge2 = faceEdges[2];
    
    label lastTriangle = startTriangle;
    
    for(int i = 0; i < nRuns; i++)
    {
        label firstTriangle = 
            triSurfaceTools::otherFace(surf, lastTriangle, lastEdge1);
        label secondTriangle = 
            triSurfaceTools::otherFace(surf, lastTriangle, lastEdge2);
        
        if(visited[firstTriangle] && secondTriangle != -1)
        {
            thisEdge = lastEdge2;
            lastTriangle = secondTriangle;
        }
        else if (firstTriangle != -1)
        {
            thisEdge = lastEdge1;
            lastTriangle = firstTriangle;
        }
        
        triSurfaceTools::otherEdges(
            surf, 
            lastTriangle, 
            thisEdge, 
            lastEdge1, 
            lastEdge2
        );
        
//        if(nextTris[lastTri])
//        {
//            nextIntersection = lastTriangle;
//            break;
//        }
        
        forAll(otherCuts, otherCutI)
        {
            label cutI = otherCuts[otherCutI];
            GeometryCut cut = cuts[cutI];
            if(cut.triangle() == lastTriangle)
            {
                nextCut = cutI;
                i = nRuns;
            }
        }
        
        visited.insert(lastTriangle);
    }
    
    return nextCut;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

