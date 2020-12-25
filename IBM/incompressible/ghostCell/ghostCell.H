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
    Foam::IBM::ghostCell

Description
    Class for incompressible ghost cell IBM method

    cite: "A gohst-cell immersed boundary method for flow in complex geometry."

SourceFiles
    ghostCell.C

\*---------------------------------------------------------------------------*/
#ifndef ghostCell_H_
#define ghostCell_H_

#include "IBM.H"
#include "fvCFD.H"

namespace Foam
{

class ghostCell
:
    public IBM
{

private:

    // ------------ mesh part -------------//

    //- The width of mesh
    scalar h_;

    // ------------ mesh part -------------//

    // ------------ particle part -------------//

    //- The fluid density in specified temperature [kg/m3]
    scalar rhof_;

    //- particle density [kg/m3]
    scalar rhop_;

    //- Fluid cell indicator, marking only live fluid cells
    mutable volScalarField* fluidPtr_;

    //- Extend solid cell indicator, marking ghost and solid cells
    mutable volScalarField* solidExtPtr_;

    //- Ghost cell indicator, marking only ghost cells
    mutable labelList* ghostPtr_;

    //- Points at the immersed boundary (IB points)
    //  nearest to the IB cell centres
    mutable vectorList* ibPointsPtr_;

    //- Normals at IB points
    mutable vectorList* ibNormalsPtr_;

    //- Sampling points for IB cells
    mutable vectorList* ibImaginePointsPtr_;

    //- Interpolation cells for ghost cells
    mutable labelListList* ghostCellCellsPtr_;

    //- Interpolation matrices for Dirichlet BC at the IB
    mutable PtrList<scalarSquareMatrix>* dirichletMatricesPtr_;

    //- Interpolation matrices for Neumann BC at the IB
    mutable PtrList<scalarSquareMatrix>* neumannMatricesPtr_;

    //- The IBforce of particle
    vector force_;

    //- The torque of particle
    vector torque_;

    //- The Collisionforce of particle
    vector collisionForce;

    scalar gradValue_;

    // ------------ particle part -------------//

    // ------------ parameter calculation -------------//

    //- The parameter for calculation drag and lift force coefficient
    scalar rhoRef_;
    scalar magUInf_;
    scalar ARef_;
    vector dragDir_;
    vector liftDir_;

    //- Drag coefficient and Lift coefficient
    scalarList coeffs;
    // ------------ parameter calculation -------------//

public:

    // Constructors
    ghostCell
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    ~ghostCell();
        
    // Member functions

    //- The functions to get protected data

    inline const scalar& getH() const
    {
        return h_;
    }

    inline const scalar& getRhof() const
    {
        return rhof_;
    }

    inline const scalar& getRhop() const
    {
        return rhop_;
    }

    inline const vector& getIbForce() const
    {
        return force_;
    }

    inline const vector& getIbTorque() const
    {
        return torque_;
    }

    inline const scalarList& getCoeffs() const
    {
        return coeffs;
    }

    inline const labelList& getGhostCells() const
    {
        return *ghostPtr_;
    }

    //- The demand-driven function

    // ---------------- Properties calculation --------------- //

    //- Get the mass of fluent in sphere/circle
    scalar getFluMass();

    //- Get the mass of sphere/circle
    scalar getParMass();

    //- Get the mass difference between particle and fluid
    scalar getMass();

    //- Get the Moment of inertia of fluent in sphere/circle
    scalar getFluMomentInertia();

    //- Get the Moment of inertia of sphere/circle
    scalar getParMomentInertia();

    //- Get the moment inertia difference between particle and fluid
    scalar getMomentInertia();

    //- get the inertial moment for any particle
    tensor getJ();

    //- Get the eval Point Velocity
    inline vector evalPointVelocity(const vector& p) const
    {
        return getVelocity() + (getOmega() ^ (p - getCenter()));
    }

    //- Update data

    //- Update center and lagrange points coordinates
    void upCenter(const vector& displacement);

    // ---------------- Properties calculation --------------- //

    // ---------------- IBM ghost cell method implementation --------------- //

    //- Make fluid and solid extend cell indicator
    void makeFluidAndSolidExtCell() const;

    //- Make ghost cell indicator
    void makeGhostCell() const;

    //- Make immersed boundary points and normals
    void makeIbPointsAndNormals() const;

    //- Make interpolation cells to calculate the value in ghost cells
    void makeGhostCellCells() const;

    //- Find the cell with the nearest cell centre
    label findNearestCell(const point& location) const;

    //- Get the minmum width of grid
    scalar findMinGridWidth() const;

    //// ------Boundary evaluation matrices------ ////
    //- Make dirichlet interpolation matrices
    void makeDirichletMatrices() const;

    //- Impose dirichlet condition
    void imposeDirichletCondition(volVectorField& U);

    //- Make dirichlet interpolation matrices
    void makeNeumannMatrices() const;

    //- Impose dirichlet condition
    void imposeNeumannCondition(volScalarField& p);
    //// ------Boundary evaluation matrices------ ////

    //- solve the correction field
    void solve
    (
        volVectorField& U,
        volScalarField& p
    );

    //- Update velocity, angular velocity, position， 努塞尔数待定;
    void update
    (
        const scalar& deltaT,
        const vector& g
    );

    //- Level set function to calculate volume faction
    scalar Ix(const vector& x, const scalar& h);

    //- Clear all pointer data
    // initial function(velocity, enthalpy, force, torque etc)
    void clearOut();
    

    // ---------------- IBM ghost cell method implementation --------------- //


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End namespace Foam

#endif
