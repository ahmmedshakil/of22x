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
    surfaceRefine

Description
    Refine edges, which are longer than specified length.

\*---------------------------------------------------------------------------*/


#include "triSurface.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("output surfaceFile");
    argList::validArgs.append("maximum length");
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const fileName outFileName = args[2];  
    const scalar maxLength = args.argRead<scalar>(3);

    Info<< "Reading surface from " << surfFileName << " ..." << endl;
    
//    Info<< maxLength << nl;

    triSurface surf1(surfFileName);
    
    pointField points = surf1.points();
    
    edgeList edges = surf1.edges();
    
    List<bool> refineEdges(edges.size(), false);
    
    labelListList edgeFaces = surf1.sortedEdgeFaces();
    
    
//    Info << edges << nl;
//    
//    Info << surf1.sortedEdgeFaces() << nl;
    
    List<bool> refineFaces(surf1.size(), false);
    
    forAll(edges, i)
    {   
        const edge edgeI = edges[i];
        
        if(edgeI.mag(points) > maxLength)
        {
            for(int k = 0; k < 2; k++)
            {
                label faceI = edgeFaces[i][k];
                refineFaces[faceI] = true;
            }
        }
    }
    
//    Info << refineEdges << nl;
    
    
//    List<bool> refineFaces(surf1.size(), false);
//    
//    forAll(surf1, i)
//    {
//        labelledTri faceI = surf1[i];
//             
//        Info << faceI << nl;   
////        if(getEdge
//        
//        
//        for(int k = 0; k < 3; k++)
//        {
//            label edgeI = faceI[k];
//            Info << edgeI << nl;
//            if(refineEdges[edgeI])
//            {
//                refineFaces[i] = true;
//            }
//        }
//    }
//    
////    Info << refineFaces << nl;
//    
    
    DynamicList<label> refineF(surf1.size());
    
    forAll(refineFaces, i)
    {
        if(refineFaces[i])
        {
            refineF.append(i);
        }
    }
    
    
    
    Info << refineF << nl;
    
//    triSurface refinedSurf = triSurfaceTools::greenRefine(surf1, refineEdges);        

    triSurface refinedSurf = triSurfaceTools::redGreenRefine(surf1, refineF);  

    Info<< "Writing refined surface to " << outFileName << " ..." << endl;

    refinedSurf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
