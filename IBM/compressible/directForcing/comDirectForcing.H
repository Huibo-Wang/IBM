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
    Foam::IBM::comDirectForcing

Description
    Class for compressible direct forcing IBM method

    cite: "A boundary condition-enforced immersed boundary method
    for compressible viscous flows."

SourceFiles
    comDirectForcing.C

\*---------------------------------------------------------------------------*/
#ifndef comDirectForcing_H_
#define comDirectForcing_H_

#include "IBM.H"
#include "fvCFD.H"

namespace Foam
{

class comDirectForcing
:
    public IBM
{

private:

    // ------------ mesh part -------------//

    //- The number of lagrangian points
    label number_;

    //- The width of mesh
    scalar h_;

    // ------------ mesh part -------------//

    // ------------ particle part -------------//

    //- The enthalpy in specified temperature [m2/s2]
    scalar enthalpy_;

    //- The fluid density in specified temperature [kg/m3]
    scalar rhof_;

    //- particle density [kg/m3]
    scalar rhop_;

    //- The left and right part for A
    //- A1
    scalarRectangularMatrix A1;
    //- A2
    scalarRectangularMatrix A2;

    //- List to store coordinates of lagrange points
    vectorList coordinates_;

    //- List of cells next to immersed boundary (IB cells)
    labelList ibCellList_;

    //- The IBvelocity of lagrange points
    vectorList velocityLagrange_;

    //- The IBenthalpy of lagrange points
    scalarList enthalpyLagrange_;

    //- The IBforce of particle
    vector force_;

    //- The torque of particle
    vector torque_;

    //- The Collisionforce of particle
    vector collisionForce;

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
    comDirectForcing
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const scalar& enthalpy,
        const scalar& rhof
    );

    ~comDirectForcing();
        
    // Member functions

    //- The functions to get protected data

    inline const label& getNumber()
    {
        return number_;
    }

    inline const scalar& getH() const
    {
        return h_;
    }

    inline const scalar& getEnthalpy() const
    {
        return enthalpy_;
    }

    inline const scalar& getRhof() const
    {
        return rhof_;
    }

    inline const scalar& getRhop() const
    {
        return rhop_;
    }

    inline const vectorList& getCoordinates() const
    {
        return coordinates_;
    }

    inline const labelList& getIbCellList() const
    {
        return ibCellList_;
    }

    inline const vectorList& getVelocityLagrange() const
    {
        return velocityLagrange_;
    }

    inline const scalarList& getEnthalpyLagrange() const
    {
        return enthalpyLagrange_;
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
    void upCoordinates(const vector& displacement);

    // ---------------- Properties calculation --------------- //

    // ---------------- IBM direct forcing implementation --------------- //

    //- Delta function
    // A.M. Roma, C.S. Peskin, M.J. Berger, An adaptive version of the immersed boundary method,
    // J. Comput. Phys. 153 (1999) 509–534.
    scalar deltaFunction(scalar r);

    //- find ibcell list
    void ibCellListSet();

    //- calculate deltaMatrices--the first part of eqn.10
    void calcA();

    //- Parameter (velocity or enthalpy correction)
    template<class fieldType>
    RectangularMatrix<fieldType> parameterCorrection
    (
        const GeometricField<fieldType, fvPatchField, volMesh>& eField,
        List<fieldType>& lField
    );

    //- Convert Euler field to lagrange field
    template<class EulerField, class LagrangeField>
    void EulerToLagrange
    (
        const EulerField& eField,
        LagrangeField& lField
    );

    //- Convert lagrange field to Euler field
    template<class EulerField, class LagrangeField>
    void LagrangeToEuler
    (
        const LagrangeField& lField,
        EulerField& eField
    );

    //- solve the correction field
    void solve
    (
        volVectorField& U,
        volScalarField& p,
        volScalarField& h,
        volScalarField& rho,
        volVectorField& drhoU,
        volScalarField& drho
    );

    //- Update velocity, angular velocity, position， 努塞尔数待定;
    void update
    (
        const scalar& deltaT,
        const vector& g
    );

    //- Level set function to calculate volume faction
    scalar Ix(const vector& x, const scalar& h);

    // initial function(velocity, enthalpy, force, torque etc)
    void clearOut();
    

    // ---------------- IBM direct forcing implementation --------------- //

    //// The list for old setting, maybe used in future
//        void velocityNextTimeStep(const scalar& deltaT, const vector& g, const volVectorField& U);

//        void displacementNextTimeStep(const scalar& deltaT);

//        //- The test to calculate u_p in next time
//        void testVelocityNextTimeStep
//        (
//            const volVectorField& U,
//            const volVectorField& Uold,
//            const scalar& deltaT,
//            const vector& g
//        );

//        void omegaNextTimeStep(const scalar& deltaT, const volVectorField& U);

//        //- The test to calculate omega_p in next time
//        void testOmegaNextTimeStep
//        (
//            const volVectorField& U,
//            const volVectorField& Uold,
//            const scalar& deltaT,
//            const vector& g
//        );

//        void initialization();

//        // calculate pressure force and viscous force
//        void calculateForce
//        (
//            const volVectorField& f,
//            const scalar& rhoRef_,
//            const scalar& magUInf_,
//            const scalar& Aref_
//        );

//        //- Calculate the collision force between particles
//        //void collision_twoParticles(const comDirectForcing& A);


//        void computeCyclicHits(const scalar L);

//        bool onMesh();

//        label patchHit(const label patchi, vector& R);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}   // End namespace Foam

#endif
