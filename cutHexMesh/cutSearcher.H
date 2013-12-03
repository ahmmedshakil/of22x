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

Typedef
    Foam::CutSearcher

Description
    Searches for intersection between a triSurface and a mesh. 

SourceFiles
    cutSearcher.C

\*---------------------------------------------------------------------------*/

#ifndef cutSearcher_H
#define cutSearcher_H

#include "triSurface.H"
#include "polyMesh.H"
#include "meshTools.H"
#include "geometryCut.H"
#include "triSurfaceSearch.H"
#include "unitConversion.H"
#include "treeDataFace.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class CutSearcher
{    
    // Private data
    
        //- Cutting triSurface
        const triSurface& surf_;
        
        //- Mesh cut by triSurface
        const polyMesh& mesh_;
        
        //- List of GeometryCuts
        DynamicList<GeometryCut> cuts_;
        
        //- For every cell the triangles which are in the cell
        List<labelHashSet> trianglesPerCell_;
        
        //- Angle value of 89deg for parallel cuts check
        scalar alignedCos_;

public:

    //- Runtime type information
    ClassName("CutSearcher");
    

    // Constructors
        
        //- Contruct from mesh and triSurface
        CutSearcher(
            const polyMesh& mesh,
            const triSurface& surf
        );

    //- Destructor
    ~CutSearcher();
    

    // Member Functions

        // Access
        
        
        // Methods
        
            void computeCuts();
            
            void writeCuts();
            
            void computeTrianglesPerCell();
            
//            void agglomerateTriangles();
            
        // Operators
            
        
        // Write
        
        
};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

