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
    Foam::GeometryCut

Description
    Saves information about geometric cuts of a point or an edge with a
    triangle. Safes the point or edge label, as well as the triangle
    label. If the cut geometry is an edge, the position of the cut will be
    saved with an edge weight. 

SourceFiles
    geometryCut.C

\*---------------------------------------------------------------------------*/

#ifndef geometryCuts_H
#define geometryCuts_H

#include "DynamicList.H"
#include "typeInfo.H"
#include "triSurface.H"
#include "triSurfaceTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class GeometryCut
{    
    // Private data
    
        //- Index of cut geometry
        label geometry_;
        
        //- Index of cutting triangle
        label triangle_;
        
        //- Weight of cut on edge. -1 if point. 
        scalar weight_;
        
        static triSurface surf;

public:

    //- Runtime type information
    ClassName("GeometryCut");
    

    // Constructors
    
        //- Construct with default values
        GeometryCut();
        
        //- Contruct from edge with edge weight
        GeometryCut(
            const label geometry,
            const label triangle,
            const scalar weight
        );
        
        //- Construct from point
        GeometryCut(
            const label geometry,
            const label triangle
        );

    //- Destructor
    ~GeometryCut();
    

    // Member Functions

        // Access
        
            
            label geometry() const
            {
                return geometry_;
            }
            
            label triangle() const
            {
                return triangle_;
            }
            
            scalar weight() const
            {
                return weight_;
            }
            
            bool isEdge() const
            {
                if(weight_ == -1)
                    return false;
                else
                    return true;
            }
            
            bool isPoint() const
            {
                if(weight_ == -1)
                    return true;
                else
                    return false;
            }
            
        // Methods
            
            static void setTriSurface(const triSurface& triSurf);
            
            label findNext(
                const triSurface& surf,
                const List<GeometryCut>& cuts,
                const labelList& otherCuts
            );
            
            bool isEqual(
                const GeometryCut& otherCut
            );
            
        // Operators
            
            inline bool operator==(const GeometryCut &rhs) const
            {
                return
                    geometry_ == rhs.geometry()
                 && triangle_ == rhs.triangle()
                 && weight_ == rhs.weight();
            }

            inline bool operator!=(const GeometryCut &rhs) const
            {
                return !operator==(rhs);
            }

            inline bool operator<(const GeometryCut &rhs) const
            {
                bool lt = true;
                if (!(geometry_ < rhs.geometry()))
                {
                    if (geometry_ == rhs.geometry())
                    {
                        if (!(triangle_ < rhs.triangle()))
                        {
                            lt = false;
                        }
                    }
                    else
                    {
                        lt = false;
                    }
                }
                
                return lt;
            }
        
        // Write
        
            friend Ostream& operator<< (Ostream& os, const GeometryCut& cut)
            {
                if (os.format() == IOstream::ASCII)
                {
                    os << cut.geometry_ << token::SPACE 
                        << cut.triangle_ << token::SPACE << cut.weight_;
                }
                else
                {
                    os.write
                    (
                        reinterpret_cast<const char*>(&cut),
                        sizeof(GeometryCut)
                    );
                }

                // Check state of Ostream
                os.check("Ostream& operator<<(Ostream&, const PointIndexHit&)");

                return os;
            }
        
};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

