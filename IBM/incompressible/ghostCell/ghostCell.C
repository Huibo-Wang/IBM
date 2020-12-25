#include "ghostCell.H"

// Constructors

    Foam::ghostCell::ghostCell
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
    :
    IBM(mesh, dict),
    rhof_(readScalar(dict.lookup("rhof"))),
    rhop_(readScalar(dict.lookup("rhop"))),
    gradValue_(readScalar(dict.lookup("gradValue"))),
    rhoRef_(readScalar(dict.subDict("force").lookup("rhoRef"))),
    magUInf_(readScalar(dict.subDict("force").lookup("magUInf"))),
    ARef_((readScalar(dict.subDict("force").lookup("ARef")))),
    dragDir_(dict.subDict("force").lookupOrDefault("dragDir", vector(1,0,0))),
    liftDir_(dict.subDict("force").lookupOrDefault("liftDir", vector(0,1,0)))
    {
        const Vector<label>& directions = mesh.geometricD();
        if (mesh.nGeometricD() == 2)
        {
            if(directions[0] == -1 || directions[1] == -1)
            {
                FatalErrorIn("ghostCell::ghostCell()")
                    << "for 2-D simulations, we calculate using the x- and "
                    << "y- direction." << abort(FatalError);
            }

            h_ = Foam::sqrt(min(mesh.V().field())/mesh.bounds().span()[2]);
        }
        else
        {
            FatalErrorIn("ghostCell::ghostCell()")
                << "Current ghost cell method is just used in 2D"
                << abort(FatalError);
        }

        force_ = vector::zero;
        torque_ = vector::zero;
        collisionForce = vector::zero;

        // coeffs[0]:drag; coeffs[1]:lift
        coeffs = scalarList(2, 0.0);
    }

    Foam::ghostCell::~ghostCell()
    {
    }

// ----------------- Properties calculation ------------------ //
    //- Parameter calculation
    Foam::scalar Foam::ghostCell::getFluMass()
    {
        scalar Pi = constant::mathematical::pi;
        scalar mass;
        mass = this->rhof_ * 1.0/4.0*Pi*pow(this->getDiameter(),2);
        return mass;
    }

    Foam::scalar Foam::ghostCell::getParMass()
    {
        scalar Pi = constant::mathematical::pi;
        scalar mass;
        mass = this->rhop_ * 1.0/4.0*Pi*pow(this->getDiameter(),2);
        return mass;
    }

    Foam::scalar Foam::ghostCell::getMass()
    {
        return getParMass()-getFluMass();
    }

    Foam::scalar Foam::ghostCell::getFluMomentInertia()
    {
        scalar Pi = constant::mathematical::pi;
        scalar momentInertia;
        momentInertia = rhof_ * 1.0/32.0 * Pi * pow(this->getDiameter(), 4);
        return momentInertia;
    }

    Foam::scalar Foam::ghostCell::getParMomentInertia()
    {
        scalar Pi = constant::mathematical::pi;
        scalar momentInertia;
        momentInertia = rhop_ * 1.0/32.0 * Pi * pow(this->getDiameter(), 4);
        return momentInertia;
    }

    Foam::scalar Foam::ghostCell::getMomentInertia()
    {
        return getParMomentInertia()-getFluMomentInertia();
    }


    Foam::tensor Foam::ghostCell::getJ()
    {
        const fvMesh& mesh = getMesh();

        tensor J = tensor::zero;
        tensor I(1, 0, 0, 0, 1, 0, 0, 0, 1);
        const vectorField& C = mesh.cellCentres();
        forAll(C, celli)
        {
            vector r = C[celli] - getCenter();
            J = J + Ix(C[celli], h_) * getRhop() * mesh.V()[celli]
                  * ((r & r) * I - (r * r));
        }
        return J;
    }

    //- Update data

    void Foam::ghostCell::upCenter(const vector& displacement)
    {
        vector& center = getCenter();
        center += displacement;
    }

// ----------------- IBM ghost cell method implementation ------------------ //

    void Foam::ghostCell::makeFluidAndSolidExtCell() const
    {
        const fvMesh& mesh = getMesh();

//        if (debug)
//        {
//            InfoIn("void ghostCell::makeFluidAndSolidExtCell() const")
//                << "creating fluid and solid extend cell indicator "
//                << "for immersed boundary" << endl;
//        }

        if (fluidPtr_)
        {
            FatalErrorIn("void ghostCell::makeFluidAndSolidExtCell() const")
                << "fluid and solid extend cell indicator already exist"
                << "for immersed boundary" << abort(FatalError);
        }

        fluidPtr_ =
            new volScalarField
            (
                IOobject
                (
                    "fluid",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("1", dimless, 1)
            );

        scalarField& fluidI = fluidPtr_->internalField();

        solidExtPtr_ =
                new volScalarField
                (
                    IOobject
                    (
                        "solidExt",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("1", dimless, 1)
                );

        scalarField& solidExtI = solidExtPtr_->internalField();

        const vectorField& C = mesh.cellCentres();

        scalar R = getDiameter()/2.0;

        forAll(C, cellI)
        {
            scalar delta = mag(C[cellI] - getCenter());
            if(delta > R)
            {
                solidExtI[cellI] = 0.0;
            }
            else
            {
                fluidI[cellI] = 0.0;
            }
        }
    }

    void Foam::ghostCell::makeGhostCell() const
    {
//        if (debug)
//        {
//            InfoIn("void ghostCell::makeGhostCell() const")
//                << "creating ghost cell indicator "
//                << "for immersed boundary" << endl;
//        }

        if (ghostPtr_)
        {
            FatalErrorIn("void ghostCell::makeGhostCell() const")
                << "ghost cell indicator already exist"
                << "for immersed boundary" << abort(FatalError);
        }

        labelHashSet ghostCellSet;

        const fvMesh& mesh = getMesh();

        const unallocLabelList& owner = mesh.owner();
        const unallocLabelList& neighbour = mesh.neighbour();

        const scalarField& solidExtI = solidExtPtr_->internalField();

        forAll(neighbour, faceI)
        {
            if(mag(solidExtI[neighbour[faceI]] - solidExtI[owner[faceI]]) > SMALL)
            {
                if(solidExtI[owner[faceI]] > SMALL)
                {
                    if(!ghostCellSet.found(owner[faceI]))
                    {
                        ghostCellSet.insert(owner[faceI]);
                    }
                }
                else
                    if(solidExtI[neighbour[faceI]] > SMALL)
                    {
                        if(!ghostCellSet.found(neighbour[faceI]))
                        {
                            ghostCellSet.insert(neighbour[faceI]);
                        }
                    }
            }
        }

        ghostPtr_ = new labelList(ghostCellSet.toc());
        sort(*ghostPtr_);

        Pout << "Number of ghost cells: " << ghostPtr_->size() << endl;
    }

    void Foam::ghostCell::makeIbPointsAndNormals() const
    {
//        if (debug)
//        {
//            InfoIn
//            (
//                "void ghostCell::makeIbPointsAndNormals() const"
//            )   << "create immersed  boundary points and normals "
//                << "for immersed boundary " << endl;
//        }

        // It is an error to attempt to recalculate
        // if the pointer is already set
        if (ibPointsPtr_ || ibNormalsPtr_ || ibImaginePointsPtr_)
        {
            FatalErrorIn
            (
                "ghost::makeIbPointsAndNormals() const"
            )
                << "immersed boundary points and normals already exist"
                << "for immersed boundary " << abort(FatalError);
        }

        const labelList& ghostCell = *ghostPtr_;

        ibPointsPtr_ = new vectorList(ghostCell.size(), vector::zero);
        ibNormalsPtr_ = new vectorList(ghostCell.size(), vector::zero);
        ibImaginePointsPtr_ = new vectorList(ghostCell.size(), vector::zero);

        vectorList& ibPoints = *ibPointsPtr_;
        vectorList& ibNormals = *ibNormalsPtr_;
        vectorList& ibImaginePoints = *ibImaginePointsPtr_;

        const fvMesh& mesh = getMesh();
        const vectorField& C = mesh.cellCentres();

        const vector& center = getCenter();

        const scalar& diameter = getDiameter();
        scalar radii = diameter/2.0;

        forAll(ghostCell, i)
        {
            label cellI = ghostCell[i];
            vector delta = C[cellI] - center;

            ibNormals[i] = delta/mag(delta);
            ibPoints[i] = center + radii * ibNormals[i];
            ibImaginePoints[i] = 2 * ibPoints[i] - C[cellI];
        }

        Info << "ibNormals is : " << endl;
        forAll(ghostCell, cellI)
        {
            Info << ibNormals[cellI] << " " ;
        }
        Info << endl;

        Info << "ghostCell is : " << endl;
        forAll(ghostCell, cellI)
        {
            Info << C[ghostCell[cellI]] << " " ;
        }
        Info << endl;

        Info << "ibPoints is : " << endl;
        forAll(ghostCell, cellI)
        {
            Info << ibPoints[cellI] << " " ;
        }
        Info << endl;

        Info << "ibImaginePoints is : " << endl;
        forAll(ghostCell, cellI)
        {
            Info << ibImaginePoints[cellI] << " " ;
        }
        Info << endl;
    }

    void Foam::ghostCell::makeGhostCellCells() const
    {
//        if (debug)
//        {
//            InfoIn("void ghostCell::makeGhostCellCells() const")
//                << "Make interpolation cells to calculate the value in ghost cells"
//                << endl;
//        }

        // It is an error to attempt to recalculate
        // if the pointer is already set
        if
        (
            ghostCellCellsPtr_
        )
        {
            FatalErrorIn
            (
                "ghostCell::makeGhostCellCells() const"
            )   << "cell-cell addressing already exists"
                << abort(FatalError);
        }

        const labelList& ghostCell = *ghostPtr_;
        ghostCellCellsPtr_ = new labelListList(ghostCell.size());
        labelListList& cellCells = *ghostCellCellsPtr_;

        const fvMesh& mesh = getMesh();
        const vectorField& C = mesh.cellCentres();
        const vectorField& point = mesh.points();

        const labelListList& cellPoints = mesh.cellPoints();
        const labelListList& pointCells = mesh.pointCells();

        const volScalarField& fluid = *fluidPtr_;
        const vectorList& ibImaginePoints = *ibImaginePointsPtr_;

        Info << "ibImaginePoints: " << endl;
        forAll(ibImaginePoints, pointI)
        {
            Info << ibImaginePoints[pointI] << " " ;
        }
        Info << endl;

        Info << "cellcells' size is: " << cellCells.size() << endl;

        forAll(cellCells, lagCellI)
        {
            labelList selectedCells(2, -1);  // 2 cells to use for interpolation

            // find the cell where the imagine point is
            label imagineCellId = mesh.findCell(ibImaginePoints[lagCellI]);

            if(imagineCellId == -1)
            {
                FatalErrorIn
                (
                    "ghostCell::makeGhostCellCells() const"
                )   << "the function of findNearestCell is wrong"
                    << abort(FatalError);
            }

            const labelList& curCellPoints = cellPoints[imagineCellId];

            label pointId;
            scalar distance = GREAT;

            forAll(curCellPoints, pointI)
            {
                // The distance between imagine point and vortex
                scalar a = mag(point[curCellPoints[pointI]] - ibImaginePoints[lagCellI]);
                if(a < distance)
                {
                    pointId = pointI;
                    distance = a;
                }
            }
            labelList curCells(4, -1); // I will choose two cells in this list
            curCells = pointCells[pointId];

            // find two nearest fluid nodes in curCells
            scalar dis1 = GREAT;
            scalar dis2 = GREAT;

            forAll(curCells, i)
            {
                scalar cellI = curCells[i];
                scalar distance = mag(C[cellI] - ibImaginePoints[lagCellI]);
                if((distance < dis1) && (fluid[cellI] > SMALL))
                {
                    selectedCells[0] = cellI;
                    dis1 = distance;
                }
            }

            forAll(curCells, i)
            {
                scalar cellI = curCells[i];
                scalar distance = mag(C[cellI] - ibImaginePoints[lagCellI]);
                if((distance < dis2) && (fluid[cellI] > SMALL) && (cellI != selectedCells[0]))
                {
                    selectedCells[1] = cellI;
                    dis2 = distance;
                }
            }

            cellCells[lagCellI] = selectedCells;
        }
    }

    Foam::label Foam::ghostCell::findNearestCell(const point& location) const
    {
        const fvMesh& mesh = getMesh();
        const vectorField& C = mesh.cellCentres();

        scalar width = findMinGridWidth()/2.0;
        label nearestCellI = -1;

        forAll(C, cellI)
        {
            scalar dx = mag((location - C[cellI]).x());
            scalar dy = ((location - C[cellI]).y());
            if((dx <= width) && (dy <= width))
            {
                nearestCellI = cellI;
                return nearestCellI;
            }
        }
        return -1;
    }

    Foam::scalar Foam::ghostCell::findMinGridWidth() const
    {
        const fvMesh& mesh = getMesh();
        scalarField delta(mesh.V().field());

        const Vector<label>& directions = mesh.geometricD();

        if (directions[0] == -1 || directions[1] == -1)
        {
            FatalErrorIn("ghost::findMinGridWidth()")
                << "For 2-D simulations with the immersed boundary method "
                << "the geometry needs to be aligned with the z-direction.  "
                << "Having the x- or y-direction as empty is not allowed "
                << "because of the way the polynomials are expanded."
                << abort(FatalError);
        }

        scalar thickness = 0.0;

        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh.bounds().span()[dir];
                break;
            }
        }

        Info << "the thickNess is: " << thickness << endl;

        // Field created with mapping for IB cells only
        //delta = sqrt(scalarField(mesh.V().field(), ghostCells())/thickness);
        delta /= thickness;

        return min(sqrt(delta));
    }

    //// Boundary evaluation matrices
    void Foam::ghostCell::makeDirichletMatrices() const
    {
//        if (debug)
//        {
//            Info<< "ghostCell::makeDirichletMatrices() : "
//                << "making immersed boundary matrices"
//                << endl;
//        }

        // It is an error to attempt to recalculate
        // if the pointer is already set
        if (dirichletMatricesPtr_)
        {
            FatalErrorIn
            (
                "void ghostCell::makeDirichletMatrices()"
            )   << "immersed boundary least squares matrices already exist"
                << abort(FatalError);
        }

        // Get addressing
        const vectorList& imaginePoints = *ibImaginePointsPtr_;
        const vectorList& ibPoints = *ibPointsPtr_;
        const labelListList& ghostCellCells = *ghostCellCellsPtr_;

        dirichletMatricesPtr_ =
            new PtrList<scalarSquareMatrix>(imaginePoints.size());
        PtrList<scalarSquareMatrix>& dirichletMatrices = *dirichletMatricesPtr_;

        const fvMesh& mesh = getMesh();
        const vectorField& C = mesh.C().internalField();

        forAll(dirichletMatrices, lagCellI)
        {
            dirichletMatrices.set
            (
                lagCellI,
                new scalarSquareMatrix
                (
                    3,
                    3,
                    0.0
                )
            );
            scalarSquareMatrix& curMatrix = dirichletMatrices[lagCellI];

            for(label i = 0; i < 2; i++)
            {
                scalar X = C[ghostCellCells[lagCellI][i]].x();
                scalar Y = C[ghostCellCells[lagCellI][i]].y();

                label coeff = 0;
                curMatrix[i][coeff++] = 1;
                curMatrix[i][coeff++] = X;
                curMatrix[i][coeff++] = Y;
            }

            // boundary ib point
            curMatrix[2][0] = 1;
            curMatrix[2][1] = ibPoints[lagCellI].x();
            curMatrix[2][2] = ibPoints[lagCellI].y();
        }
    }

    void Foam::ghostCell::imposeDirichletCondition
    (
        volVectorField& U
    )
    {
        // Get addressing
        const labelList& ghostCells = *ghostPtr_;
        const labelListList& ghostCellCells = *ghostCellCellsPtr_;
        const vectorList& imaginePoints = *ibImaginePointsPtr_;
        const vectorList& ibPoints = *ibPointsPtr_;

        PtrList<scalarSquareMatrix>& dirichletMatrices = *dirichletMatricesPtr_;

        vectorField& UI = U.internalField();

        // Calculate Dirichlet values on ib points
        vectorList ibValues(ibPoints.size(), vector::zero);

        forAll(ghostCells, lagCellI)
        {
            ibValues[lagCellI] = evalPointVelocity(ibPoints[lagCellI]);

            vectorList source(3, vector::zero);
            vectorList a(3, vector::zero);  // a0, a1, a2

            for(label i = 0; i < 2; i++)
            {
                label cellI = ghostCellCells[lagCellI][i];
                source[i] = UI[cellI];
            }

            source[2] = ibValues[lagCellI];

            LUsolve(dirichletMatrices[lagCellI], source);

            a = source;

            vector imaginePointsVelocity
                = a[0] + imaginePoints[lagCellI].x()*a[1] + imaginePoints[lagCellI].y()*a[2];

            // the velocity on ghost cell center
            UI[ghostCells[lagCellI]] = 2*ibValues[lagCellI] - imaginePointsVelocity;
        }
    }

    void Foam::ghostCell::makeNeumannMatrices() const
    {
//        if (debug)
//        {
//            Info<< "ghostCell::makeNeumannMatrices() : "
//                << "making immersed boundary matrices"
//                << endl;
//        }

        // It is an error to attempt to recalculate
        // if the pointer is already set
        if (neumannMatricesPtr_)
        {
            FatalErrorIn
            (
                "void ghostCell::makeNeumannMatrices()"
            )   << "immersed boundary least squares matrices already exist"
                << abort(FatalError);
        }

        // Get addressing
        const vectorList& imaginePoints = *ibImaginePointsPtr_;
        const vectorList& ibNormals = *ibNormalsPtr_;
        const labelListList& ghostCellCells = *ghostCellCellsPtr_;

        neumannMatricesPtr_ =
            new PtrList<scalarSquareMatrix>(imaginePoints.size());
        PtrList<scalarSquareMatrix>& neumannMatrices = *neumannMatricesPtr_;

        const fvMesh& mesh = getMesh();
        const vectorField& C = mesh.C().internalField();

        forAll(neumannMatrices, lagCellI)
        {
            neumannMatrices.set
            (
                lagCellI,
                new scalarSquareMatrix
                (
                    3,
                    3,
                    0.0
                )
            );
            scalarSquareMatrix& curMatrix = neumannMatrices[lagCellI];

            for(label i = 0; i < 2; i++)
            {
                scalar X = C[ghostCellCells[lagCellI][i]].x();
                scalar Y = C[ghostCellCells[lagCellI][i]].y();

                label coeff = 0;
                curMatrix[i][coeff++] = 1;
                curMatrix[i][coeff++] = X;
                curMatrix[i][coeff++] = Y;
            }

            // boundary ib point
            curMatrix[2][0] = 0;
            curMatrix[2][1] = -ibNormals[lagCellI].y();
            curMatrix[2][2] = ibNormals[lagCellI].x();
        }
    }

    void Foam::ghostCell::imposeNeumannCondition
    (
        volScalarField& p
    )
    {
        // Get addressing
        const labelList& ghostCells = *ghostPtr_;
        const labelListList& ghostCellCells = *ghostCellCellsPtr_;
        const vectorList& imaginePoints = *ibImaginePointsPtr_;

        PtrList<scalarSquareMatrix>& neumannMatrices = *neumannMatricesPtr_;

        scalarField& pI = p.internalField();

        // neumann fixed grad
        const scalar& gradValue = gradValue_;

        forAll(ghostCells, lagCellI)
        {
            scalarList source(3, 0.0);
            scalarList a(3, 0.0);  // a0, a1, a2

            for(label i = 0; i < 2; i++)
            {
                label cellI = ghostCellCells[lagCellI][i];
                source[i] = pI[cellI];
            }

            source[2] = gradValue;

            LUsolve(neumannMatrices[lagCellI], source);

            a = source;

            scalar imaginePointsPressure
                = a[0] + imaginePoints[lagCellI].x()*a[1] + imaginePoints[lagCellI].y()*a[2];

            // the velocity on ghost cell center
            pI[ghostCells[lagCellI]] = imaginePointsPressure;
        }
    }

    void Foam::ghostCell::solve
    (
        volVectorField& U,
        volScalarField& p
    )
    {
        makeFluidAndSolidExtCell();
        makeGhostCell();
        makeIbPointsAndNormals();
        makeGhostCellCells();

        makeDirichletMatrices();
        makeNeumannMatrices();

        //- Correct velocity field
        imposeDirichletCondition(U);
        imposeNeumannCondition(p);

//        // ----------------  1. Velocity correction  ---------------- //
//        forAll(coordinates_, i)
//        {
//            velocityLagrange_[i] = evalPointVelocity(coordinates_[i]);
//        }

//        //- Get dU on boundary points
//        //RectangularMatrix<vector> Xu_B = parameterCorrection(U, velocityLagrange_);
//        List<vector> Xu_B(number_);

//        Xu_B = parameterCorrection(U, velocityLagrange_);

//        if(getDimention())
//        {
//            forAll(coordinates_, i)
//            {
//                Xu_B[i].z() = 0.0;
//                //Xu_B[i][0].z() = 0.0;
//            }
//        }

//        //- Get dU on Euler points
//        List<vector> Xu(ibCellList_.size(), vector::zero);
//        //RectangularMatrix<vector> Xu(ibCellList_.size(), 1, vector::zero);

//        forAll(ibCellList_, rowi)
//        {
//            forAll(coordinates_, columnj)
//            {
//                Xu[rowi] += A2[rowi][columnj] * Xu_B[columnj];
//                //Xu[rowi][0] += A2[rowi][columnj] * Xu_B[columnj][0];
//            }
//        }

//        if(getDimention())
//        {
//            forAll(coordinates_, i)
//            {
//                Xu[i].z() = 0.0;
//                //Xu[i][0].z() = 0.0;
//            }
//        }

//        forAll(ibCellList_, i)
//        {
//            label cellI = ibCellList_[i];
//            U[cellI] += Xu[i];
//            //U[cellI] += Xu[i][0];
//        }

//        // ----------------  2. Calculate force  ---------------- //
//        const scalar& deltaT = mesh.time().deltaT().value();
//        forAll(coordinates_, i)
//        {
//            vector r = vector(coordinates_[i] - getCenter());
//            force_ -= rhof_ * Xu_B[i] / deltaT ;
//            torque_ -= rhof_ * (r ^ Xu_B[i]) / deltaT;
//            //force_ -= rhof_ * Xu_B[i][0] / deltaT ;
//            //torque_ -= rhof_ * (r ^ Xu_B[i][0]) / deltaT;
//        }

//        scalar lagrangeVolume;
//        scalar Pi = constant::mathematical::pi;

//        if(getDimention())
//        {
//            lagrangeVolume = Pi * getDiameter() * h_/number_;
//        }
//        else
//        {
//            lagrangeVolume = (Pi * h_) * (3 * pow(getDiameter(), 2) + pow(h_, 2))/(3 * number_);
//        }
//        force_ *= lagrangeVolume;
//        torque_ *= lagrangeVolume;

//        Info << "the force is: " << force_ << endl;
//        Info << "the torque is: " << torque_ << endl;

//        if(getMotion())
//        {
//            if(getVelocity() == vector::zero)
//            {
//                magUInf_ = 1; // to prevent divide 0
//            }
//            else
//            {
//                magUInf_ = mag(getVelocity());
//            }
//        }

//        scalar pDyn = 0.5*rhoRef_*pow(magUInf_,2);

//        coeffs[0] = (force_ & dragDir_)/(ARef_*pDyn);   //drag coefficient
//        coeffs[1] = (force_ & liftDir_)/(ARef_*pDyn);   //lift coefficient
    }

    void Foam::ghostCell::update
    (
        const Foam::scalar& deltaT,
        const vector& g
    )
    {
        // Feng

        //- 1. velocity update
        const vector vt = getVelocity();
        const vector vtold = getVelocityOld();

        // The velocity to update
        vector& newVelocity = getVelocity();

        updateVelocityOld(); // update old velocity
        vector VpartA = force_;
        vector VpartB = getRhof()*getVolume()*(vt-vtold)/deltaT;
        vector VpartC = getVolume()*(getRhop()-getRhof())*g;

        vector vSum = (VpartA+VpartB+VpartC)*deltaT/getParMass();

        newVelocity.x() += vSum.x();
        newVelocity.y() += vSum.y();

        //- 2. omega update
        const vector& ot = getOmega();
        const vector& otold = getOmegaOld();

        // The omega to update
        vector& newOmega = getOmega();

        updateOmegaOld(); // update old omega
        vector OpartA = torque_;
        vector OpartB = getFluMomentInertia()*(ot-otold)/deltaT;

        vector oSum = (OpartA+OpartB)*deltaT/getParMomentInertia();

        newOmega.z() += oSum.z();

        if(!getDimention())
        {
            newVelocity.z() += vSum.z();

            newOmega.x() += vSum.x();
            newOmega.y() += vSum.y();
        }

        // update coordinates and center
        vector displacement = (getVelocity() + getVelocityOld())/2 * deltaT;

        vector & newcenter = getCenter();
        newcenter += displacement;
    }

    Foam::scalar Foam::ghostCell::Ix(const vector& x, const scalar& h)
    {
        scalar Ix, psix, Hs;

        // get psi(x)
        psix = mag(x - getCenter()) - getDiameter()/2;
        // get H(s)
        Hs = 0.5*(1+tanh(5*psix/h));
        //get I(x)
        Ix = 1 - Hs;
        return Ix;
    }

    void Foam::ghostCell::clearOut()
    {
        deleteDemandDrivenData(fluidPtr_);
        deleteDemandDrivenData(solidExtPtr_);
        deleteDemandDrivenData(ghostPtr_);
        deleteDemandDrivenData(ibPointsPtr_);
        deleteDemandDrivenData(ibNormalsPtr_);
        deleteDemandDrivenData(ibImaginePointsPtr_);
        deleteDemandDrivenData(ghostCellCellsPtr_);
        deleteDemandDrivenData(dirichletMatricesPtr_);
        deleteDemandDrivenData(neumannMatricesPtr_);

//        delete fluidPtr_;
//        delete solidExtPtr_;
//        delete ghostPtr_;
//        delete ibPointsPtr_;
//        delete ibNormalsPtr_;
//        delete ibImaginePointsPtr_;
//        delete ghostCellCellsPtr_;
        force_ = vector::zero;
        torque_ = vector::zero;
    }

