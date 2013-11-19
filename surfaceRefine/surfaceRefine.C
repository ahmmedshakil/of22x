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

    triSurface surf1(surfFileName);

    forAll(surf1, i)
    {
        labelledTri surfI = surf1[i];
        Info << surfI << nl;
    }
    
    triSurface surf2 = surf1;

    Info<< "Writing refined surface to " << outFileName << " ..." << endl;

    surf2.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
