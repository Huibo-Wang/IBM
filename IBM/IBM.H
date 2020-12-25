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

Class
    Foam::IBM

Description
    Base class for immersed boundary models

SourceFiles
    IBM.C

\*---------------------------------------------------------------------------*/
#ifndef IBM_H
#define IBM_H

#include "fvMatrices.H"
#include "IOdictionary.H"
#include "scalarList.H"
#include "vectorList.H"
#include "fvMesh.H"
#include "Matrix.H"

namespace Foam
{

class IBM
{

private:

    //- Finite volume mesh reference
    const fvMesh& mesh_;

    //- Bool operator to clarify particle motion state
    bool motion_;

    //- Bool operator to clarity 2D(true) or 3D(false)
    bool dim_;

    //--------- circle(sphere) particle parameter ---------//
    //- Diameter [m]
    scalar diameter_;

    //- Center
    vector center_;

    //- Velocity
    vector velocity_;

    //- Velocity at old time
    vector velocityOld_;

    //- Angular velocity
    vector omega_;

    //- Angular velocity at old time
    vector omegaOld_;
    //--------- circle(sphere) particle parameter ---------//

public:

    //- Constructors
    IBM
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    ~IBM();

    //--------- Member functions ---------//

    //- The functions to get protected data

    inline const fvMesh& getMesh() const
    {
        return mesh_;
    }

    inline const bool& getMotion()
    {
        return motion_;
    }

    inline const bool& getDimention()
    {
        return dim_;
    }

    inline const scalar& getDiameter() const
    {
        return diameter_;
    }

    inline const vector& getCenter() const
    {
        return center_;
    }

    inline vector& getCenter()
    {
        return center_;
    }

    inline const vector& getVelocity() const
    {
        return velocity_;
    }

    inline vector& getVelocity()
    {
        return velocity_;
    }

    inline const vector& getVelocityOld() const
    {
        return velocityOld_;
    }

    inline vector& getVelocityOld()
    {
        return velocityOld_;
    }

    inline const vector& getOmega() const
    {
        return omega_;
    }

    inline vector& getOmega()
    {
        return omega_;
    }

    inline const vector& getOmegaOld() const
    {
        return omegaOld_;
    }

    inline vector& getOmegaOld()
    {
        return omegaOld_;
    }

    //- properties update
    inline void updateVelocityOld()
    {
        velocityOld_ = velocity_;
    }

    inline void updateOmegaOld()
    {
        omegaOld_ = omega_;
    }

    //- get particle volume
    scalar getVolume();

    //--------- Member functions ---------//
};

}
#endif // IMMERSEDBOUNDARYMETHOD_H

