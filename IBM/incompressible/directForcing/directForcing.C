#include "directForcing.H"

// Constructors

    Foam::directForcing::directForcing
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
    :
    IBM(mesh, dict),
    rhof_(readScalar(dict.lookup("rhof"))),
    rhop_(readScalar(dict.lookup("rhop"))),
    rhoRef_(readScalar(dict.subDict("force").lookup("rhoRef"))),
    magUInf_(readScalar(dict.subDict("force").lookup("magUInf"))),
    ARef_((readScalar(dict.subDict("force").lookup("ARef")))),
    dragDir_(dict.subDict("force").lookupOrDefault("dragDir", vector(1,0,0))),
    liftDir_(dict.subDict("force").lookupOrDefault("liftDir", vector(0,1,0)))
    {
        const scalar& d = getDiameter();
        vector& center = getCenter();

        scalar radii = d/2.0;
        scalar Pi = constant::mathematical::pi;

        if(getDimention())
        {
            h_ = Foam::sqrt(min(mesh.V().field())/mesh.bounds().span()[2]);
            number_ = 2*Pi*radii/h_;

            coordinates_ = vectorList(number_, vector::zero);

            for(label i = 0; i < number_; i++)
            {
                scalar x = radii*cos(i*2*Pi/number_)+center.x();
                scalar y = radii*sin(i*2*Pi/number_)+center.y();
                scalar z = center.z();
                coordinates_[i] = vector(x,y,z);
            }
        }
        else
        {
            h_ = Foam::pow(min(mesh.V().field()), 1.0/3.0);    // uniform mesh
            number_ = Pi/3*(12*pow(radii,2)/pow(h_,2)+1);

            coordinates_ = vectorList(number_, vector::zero);
            scalar Phi = (1+sqrt(5.0))/2;

            for(label i = 0; i < number_; i++)
            {
                scalar theta = acos(1-2.0*i/number_);
                scalar phi = 2.0*i*Pi/Phi;

                scalar x = radii*cos(phi)*sin(theta)+center.x();
                scalar y = radii*sin(phi)*sin(theta)+center.y();
                scalar z = radii*cos(theta)+center.z();

                coordinates_[i] = vector(x, y, z);
            }
        }

        Pout << "num_points" << number_ << endl;

        velocityLagrange_ = vectorList(number_, vector::zero);

        force_ = vector::zero;
        torque_ = vector::zero;
        collisionForce = vector::zero;

        // coeffs[0]:drag; coeffs[1]:lift
        coeffs = scalarList(2, 0.0);
    }

    Foam::directForcing::~directForcing()
    {
    }

// ----------------- Properties calculation ------------------ //
    //- Parameter calculation
    Foam::scalar Foam::directForcing::getFluMass()
    {
        scalar Pi = constant::mathematical::pi;
        scalar mass;
        mass = this->rhof_ * 1.0/4.0*Pi*pow(this->getDiameter(),2);
        return mass;
    }

    Foam::scalar Foam::directForcing::getParMass()
    {
        scalar Pi = constant::mathematical::pi;
        scalar mass;
        mass = this->rhop_ * 1.0/4.0*Pi*pow(this->getDiameter(),2);
        return mass;
    }

    Foam::scalar Foam::directForcing::getMass()
    {
        return getParMass()-getFluMass();
    }

    Foam::scalar Foam::directForcing::getFluMomentInertia()
    {
        scalar Pi = constant::mathematical::pi;
        scalar momentInertia;
        momentInertia = rhof_ * 1.0/32.0 * Pi * pow(this->getDiameter(), 4);
        return momentInertia;
    }

    Foam::scalar Foam::directForcing::getParMomentInertia()
    {
        scalar Pi = constant::mathematical::pi;
        scalar momentInertia;
        momentInertia = rhop_ * 1.0/32.0 * Pi * pow(this->getDiameter(), 4);
        return momentInertia;
    }

    Foam::scalar Foam::directForcing::getMomentInertia()
    {
        return getParMomentInertia()-getFluMomentInertia();
    }


    Foam::tensor Foam::directForcing::getJ()
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

    void Foam::directForcing::upCoordinates(const vector& displacement)
    {
        vector& center = getCenter();
        center += displacement;
        forAll(coordinates_, i)
        {
            coordinates_[i] += displacement;
        }
    }

// ----------------- IBM direct forcing implementation ------------------ //

    //- Delta function
    // A.M. Roma, C.S. Peskin, M.J. Berger, An adaptive version of the immersed boundary method,
    // J. Comput. Phys. 153 (1999) 509–534.
    Foam::scalar Foam::directForcing::deltaFunction(Foam::scalar r)
    {
        const scalar deltaSupportRadius = 1.5;

        if (r > deltaSupportRadius)
        {
            return 0;
        }
        else
        {
            if (r <= 0.5)
                return (1 + sqrt(1 - 3*r*r))/3.0;
            else
                return (5 - 3*r - sqrt(1 - 3*(1-r)*(1-r)))/6.0;
        }
    }

    //- find ibCell list
    void Foam::directForcing::ibCellListSet()
    {
        const fvMesh& mesh = getMesh();

        ibCellList_ = labelList(mesh.nCells(), -1);

        scalar R = getDiameter()/2.0;
        scalar delta = 0.4 * R;
        scalar innerR = R - delta;
        scalar outerR = R + delta;

        const vectorField& C = mesh.cellCentres();

        label num = 0;

        forAll(C,cellI)
        {
            if (mag(C[cellI]-getCenter()) <= outerR && mag(C[cellI]-getCenter()) >= innerR)
            {
                ibCellList_[num++] = cellI;
            }
        }

        ibCellList_.resize(num);

        Pout << "Number of IB cells: " << ibCellList_.size() << endl;
    }

    void Foam::directForcing::calcA()
    {
        const fvMesh& mesh = getMesh();
        const vectorField& C = mesh.cellCentres();

        A1 = scalarRectangularMatrix(number_, ibCellList_.size());
        A2 = scalarRectangularMatrix(ibCellList_.size(), number_);

        scalar dx, dy, dz;
        scalar ds;

        forAll(coordinates_, rowi)
        {
            forAll(ibCellList_, columnj)
            {
                dx = mag(C[ibCellList_[columnj]].x() - coordinates_[rowi].x())/h_;
                dy = mag(C[ibCellList_[columnj]].y() - coordinates_[rowi].y())/h_;
                dz = mag(C[ibCellList_[columnj]].z() - coordinates_[rowi].z())/h_;

                if(getDimention())
                {
                    ds = constant::mathematical::pi * getDiameter() * h_/ number_;
                    A1[rowi][columnj] = deltaFunction(dx)*deltaFunction(dy);
                    A2[columnj][rowi] = deltaFunction(dx)*deltaFunction(dy)/pow(h_,2)*ds;
                }
                else
                {
                    ds = (constant::mathematical::pi * h_) * (3*pow(getDiameter(), 2)+pow(h_, 2))/(3*number_) ;
                    A1[rowi][columnj] = deltaFunction(dx)*deltaFunction(dy)*deltaFunction(dz);
                    A2[columnj][rowi] = deltaFunction(dx)*deltaFunction(dy)*deltaFunction(dz)/pow(h_,3)*ds;
                }
            }
        }
    }

    template<class fieldType>
    List<fieldType> Foam::directForcing::parameterCorrection
    //RectangularMatrix<fieldType> Foam::directForcing::parameterCorrection
    (
        const GeometricField<fieldType, fvPatchField, volMesh>& eField,
        List<fieldType>& lField
    )
    {
        //// invA
//        // calculate coefficient matrix A and B
//        scalarRectangularMatrix A(number_, number_);
//        RectangularMatrix<fieldType> B(number_, 1);
//        RectangularMatrix<fieldType> B1(number_, 1);
//        // B22: the right part of B2
//        RectangularMatrix<fieldType> B22(ibCellList_.size(), 1);

//        A = A1*A2;

//        forAll(coordinates_, i)
//        {
//            B1[i][0] = lField[i];
//        }
//        forAll(ibCellList_, j)
//        {
//            B22[j][0] = eField[ibCellList_[j]];
//        }

//        // get B
//        forAll(coordinates_, rowi)
//        {
//            forAll(ibCellList_, columnj)
//            {
//                B[rowi][0] -= A1[rowi][columnj] * B22[columnj][0];
//            }
//            B[rowi][0] += B1[rowi][0];
//        }

//        scalarRectangularMatrix invA = SVDinv(A);

//        RectangularMatrix<fieldType> Xu(number_, 1);

//        forAll(coordinates_, rowi)
//        {
//            forAll(coordinates_, columnj)
//            {
//                Xu[rowi][0] += invA[rowi][columnj] * B[columnj][0];
//            }
//        }

//        RectangularMatrix<fieldType> test(number_, 1);

//        forAll(coordinates_, rowi)
//        {
//            forAll(coordinates_, columnj)
//            {
//                test[rowi][0] += A[rowi][columnj] * Xu[columnj][0];
//            }
//            test[rowi][0] -= B[rowi][0];
//        }

//        Info << "AX-B= " << endl;
//        forAll(coordinates_, i)
//        {
//            Info << test[i][0] << " ";
//        }
//        Info << endl;

//        return Xu;

        //// Lusolve
//        // calculate coefficient matrix A and B
//        scalarSquareMatrix A(number_, number_);

//        List<fieldType> B(number_);
//        List<fieldType> B1(number_);
//        // B22: the right part of B2
//        List<fieldType> B22(ibCellList_.size());

//        forAll(coordinates_, rowi)
//        {
//            forAll(ibCellList_, columnj)
//            {
//                forAll(coordinates_, k)
//                {
//                    A[rowi][k] += A1[rowi][columnj] * A2[columnj][k];
//                }
//            }
//        }

//        scalarSquareMatrix testA = A;

//        forAll(coordinates_, i)
//        {
//            B1[i] = lField[i];
//        }
//        forAll(ibCellList_, j)
//        {
//            B22[j] = eField[ibCellList_[j]];
//        }

//        // get B
//        forAll(coordinates_, rowi)
//        {
//            forAll(ibCellList_, columnj)
//            {
//                B[rowi] -= A1[rowi][columnj] * B22[columnj];
//            }
//            B[rowi] += B1[rowi];
//        }

//        List<fieldType> Xu(number_);

//        List<fieldType> testB(B);

//        LUsolve(A, B);

//        List<fieldType> D(number_);

//        forAll(coordinates_, rowi)
//        {
//            forAll(coordinates_, columnj)
//            {
//                D[rowi] += testA[rowi][columnj] * B[columnj];
//            }
//            D[rowi] -= testB[rowi];
//        }

//        Info << "AX-B= " << endl;
//        forAll(coordinates_, i)
//        {
//            Info << D[i] << " ";
//        }

//        Info << endl;

//        Xu = B;

//        return Xu;

        //// jacobi iteration
        // calculate coefficient matrix A and B
        scalarSquareMatrix A(number_, number_);

        List<fieldType> B(number_);
        List<fieldType> B1(number_);
        // B22: the right part of B2
        List<fieldType> B22(ibCellList_.size());

        forAll(coordinates_, rowi)
        {
            forAll(ibCellList_, columnj)
            {
                forAll(coordinates_, k)
                {
                    A[rowi][k] += A1[rowi][columnj] * A2[columnj][k];
                }
            }
        }

        scalarSquareMatrix testA = A;

        forAll(coordinates_, i)
        {
            B1[i] = lField[i];
        }
        forAll(ibCellList_, j)
        {
            B22[j] = eField[ibCellList_[j]];
        }

        // get B
        forAll(coordinates_, rowi)
        {
            forAll(ibCellList_, columnj)
            {
                B[rowi] -= A1[rowi][columnj] * B22[columnj];
            }
            B[rowi] += B1[rowi];
        }

        List<fieldType> Xu(number_);

        List<fieldType> testB(B);

        scalar e = 10^-4;
        label N = 10000;

        getResult(e,N,A,B,Xu);

        List<fieldType> D(number_);

        forAll(coordinates_, rowi)
        {
            forAll(coordinates_, columnj)
            {
                D[rowi] += testA[rowi][columnj] * Xu[columnj];
            }
            D[rowi] -= testB[rowi];
        }

        Info << "AX-B= " << endl;
        forAll(coordinates_, i)
        {
            Info << D[i] << " ";
        }

        Info << endl;

        Info << "Xu: " << endl;
        forAll(Xu, i)
        {
            Info << Xu << " ";
        }
        Info << endl;

        Xu = B;

        return Xu;


    }

    template<class EulerField, class LagrangeField>
    void Foam::directForcing::EulerToLagrange
    (
        const EulerField& eField,
        LagrangeField& lField
    )
    {
        scalar dx, dy, dz;
        const vectorField& C = mesh_.cellCentres();

        for(label i = 0; i < coordinates_.size(); i++)
        {
            for(label j = 0; j < ibCellList_.size(); j++)
            {
                dx = mag(C[ibCellList_[j]].x()-coordinates_[i].x()) / h_;
                dy = mag(C[ibCellList_[j]].y()-coordinates_[i].y()) / h_;
                dz = mag(C[ibCellList_[j]].z()-coordinates_[i].z()) / h_;

                if(dim_)
                    lField[i] += deltaFunction(dx) * deltaFunction(dy) * eField[ibCellList_[j]];
                else
                    lField[i] += deltaFunction(dx) * deltaFunction(dy) * deltaFunction(dz) * eField[ibCellList_[j]];
            }
        }
    }

    template<class EulerField, class LagrangeField>
    void Foam::directForcing::LagrangeToEuler
    (
        const LagrangeField& lField,
        EulerField& eField
    )
    {
        scalar dx, dy, dz;
        // The volume for each lagrange point
        scalar lagrangeVolume;
        scalar Pi = constant::mathematical::pi;

        const vectorField& C = mesh_.cellCentres();

        for(label i = 0; i < ibCellList_.size(); i++)
        {
            for(label j = 0; j < coordinates_.size(); j++)
            {
                dx = mag(C[ibCellList_[i]].x()-coordinates_[j].x()) / h_;
                dy = mag(C[ibCellList_[i]].y()-coordinates_[j].y()) / h_;
                dz = mag(C[ibCellList_[i]].z()-coordinates_[j].z()) / h_;

                if(dim_)
                {
                    lagrangeVolume = Pi * getDiameter() * h_/number_;
                    eField[ibCellList_[i]] += 1/pow(h_,2) * deltaFunction(dx) * deltaFunction(dy) * lField[j];
                }
                else
                {
                    eField[ibCellList_[i]] += 1/pow(h_,2) * deltaFunction(dx) * deltaFunction(dy) * deltaFunction(dz) * lField[j];
                }
            }
        }
    }

    void Foam::directForcing::solve
    (
        volVectorField& U
    )
    {
        ibCellListSet();
        calcA();

        const fvMesh& mesh = getMesh();
        // ----------------  1. Velocity correction  ---------------- //
        forAll(coordinates_, i)
        {
            velocityLagrange_[i] = evalPointVelocity(coordinates_[i]);
        }

        //// Lusolve and jacobi iteration part
        //- Get dU on boundary points

        List<vector> Xu_B(number_);

        Xu_B = parameterCorrection(U, velocityLagrange_);

        if(getDimention())
        {
            forAll(coordinates_, i)
            {
                Xu_B[i].z() = 0.0;
                //Xu_B[i][0].z() = 0.0;
            }
        }

        //- Get dU on Euler points
        List<vector> Xu(ibCellList_.size(), vector::zero);
        //RectangularMatrix<vector> Xu(ibCellList_.size(), 1, vector::zero);

        forAll(ibCellList_, rowi)
        {
            forAll(coordinates_, columnj)
            {
                Xu[rowi] += A2[rowi][columnj] * Xu_B[columnj];
                //Xu[rowi][0] += A2[rowi][columnj] * Xu_B[columnj][0];
            }
        }

        if(getDimention())
        {
            forAll(coordinates_, i)
            {
                Xu[i].z() = 0.0;
                //Xu[i][0].z() = 0.0;
            }
        }

        forAll(ibCellList_, i)
        {
            label cellI = ibCellList_[i];
            U[cellI] += Xu[i];
            //U[cellI] += Xu[i][0];
        }

        // ----------------  2. Calculate force  ---------------- //
        const scalar& deltaT = mesh.time().deltaT().value();
        forAll(coordinates_, i)
        {
            vector r = vector(coordinates_[i] - getCenter());
            force_ -= rhof_ * Xu_B[i] / deltaT ;
            torque_ -= rhof_ * (r ^ Xu_B[i]) / deltaT;
            //force_ -= rhof_ * Xu_B[i][0] / deltaT ;
            //torque_ -= rhof_ * (r ^ Xu_B[i][0]) / deltaT;
        }



        //// invA
//        //- Get dU on boundary points
//        RectangularMatrix<vector> Xu_B = parameterCorrection(U, velocityLagrange_);

//        if(getDimention())
//        {
//            forAll(coordinates_, i)
//            {
//                //Xu_B[i].z() = 0.0;
//                Xu_B[i][0].z() = 0.0;
//            }
//        }

//        //- Get dU on Euler points
//        RectangularMatrix<vector> Xu(ibCellList_.size(), 1, vector::zero);

//        forAll(ibCellList_, rowi)
//        {
//            forAll(coordinates_, columnj)
//            {
//                //Xu[rowi] += A2[rowi][columnj] * Xu_B[columnj];
//                Xu[rowi][0] += A2[rowi][columnj] * Xu_B[columnj][0];
//            }
//        }

//        if(getDimention())
//        {
//            forAll(coordinates_, i)
//            {
//                //Xu[i].z() = 0.0;
//                Xu[i][0].z() = 0.0;
//            }
//        }

//        forAll(ibCellList_, i)
//        {
//            label cellI = ibCellList_[i];
//            //U[cellI] += Xu[i];
//            U[cellI] += Xu[i][0];
//        }

//        // ----------------  2. Calculate force  ---------------- //
//        const scalar& deltaT = mesh.time().deltaT().value();
//        forAll(coordinates_, i)
//        {
//            vector r = vector(coordinates_[i] - getCenter());
//            //force_ -= rhof_ * Xu_B[i] / deltaT ;
//            //torque_ -= rhof_ * (r ^ Xu_B[i]) / deltaT;
//            force_ -= rhof_ * Xu_B[i][0] / deltaT ;
//            torque_ -= rhof_ * (r ^ Xu_B[i][0]) / deltaT;
//        }

        //// invA

        scalar lagrangeVolume;
        scalar Pi = constant::mathematical::pi;

        if(getDimention())
        {
            lagrangeVolume = Pi * getDiameter() * h_/number_;
        }
        else
        {
            lagrangeVolume = (Pi * h_) * (3 * pow(getDiameter(), 2) + pow(h_, 2))/(3 * number_);
        }
        force_ *= lagrangeVolume;
        torque_ *= lagrangeVolume;

        Info << "the force is: " << force_ << endl;
        Info << "the torque is: " << torque_ << endl;

        if(getMotion())
        {
            if(getVelocity() == vector::zero)
            {
                magUInf_ = 1; // to prevent divide 0
            }
            else
            {
                magUInf_ = mag(getVelocity());
            }
        }

        scalar pDyn = 0.5*rhoRef_*pow(magUInf_,2);

        coeffs[0] = (force_ & dragDir_)/(ARef_*pDyn);   //drag coefficient
        coeffs[1] = (force_ & liftDir_)/(ARef_*pDyn);   //lift coefficient
    }

    void Foam::directForcing::update
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
        for(label i = 0; i < coordinates_.size(); i++)
        {
            coordinates_[i] += displacement;
        }
    }

    Foam::scalar Foam::directForcing::Ix(const vector& x, const scalar& h)
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

    void Foam::directForcing::clearOut()
    {
        force_ = vector::zero;
        torque_ = vector::zero;

        velocityLagrange_ = vectorList(number_, vector::zero);
    }

// ---------------- Iteration algorithm accomplishment --------------- //
    template<class fieldType>
    void Foam::directForcing::swapA
    (
        label i,
        label j,
        scalarSquareMatrix& A,
        List<fieldType>& B
    )
    {
        for(label column = 0; column < A.n(); column++)
        {
            scalar temp = A[i][column];
            A[i][column] = A[j][column];
            A[j][column] = temp;
        }

        fieldType temp = B[i];
        B[i] = B[j];
        B[j] = temp;
    }

    template<class fieldType>
    void Foam::directForcing::getResult
    (
        scalar e,    // precision
        label N,    // max iteration number
        scalarSquareMatrix& A,
        List<fieldType>& B,
        List<fieldType>& re
    )
    {
        label i,j,k;
        forAll(re, i)
        {
            re[i] = vector::zero;
        }

        for(i = 0; i < A.n(); i++)
        {
            if(fabs(A[i][i]) <= 1e-2)
            {
                for(j = 0; j < A.n(); j++)
                {
                    if(fabs(A[j][i]) > 1e-2)
                    {
                        swapA(i,j,A,B);
                        break;
                    }
                }
                if(j == A.n())
                {
                    Info << "" << endl;
                    FatalErrorIn("directForcing::getResult")
                        << "Invalid coefficient matrix."
                        << abort(FatalError);
                }
                i = 0;
            }
        }
        Info << "A: " << endl;
        for(i = 0; i < A.n(); i++)
        {
            for(j = 0; j < A.m(); j++)
            {
                Info << A[i][j] << " ";
            }
            Info << endl;
        }

        Info << "B" << endl;
        for(i = 0; i < B.size(); i++)
        {
            Info << B[i] << " ";
        }
        Info << endl;

        // Iteration
        k = 0;

        List<fieldType> x(A.n());

        for(i = 0; i < A.n(); i++)
        {
            x[i] = vector::zero;
        }

        while(k <= N)
        {
            k++;
            if(k > N)
            {
                Info << "迭代失败" << endl;
                break;
            }
            for(label i = 0; i < A.n(); i++)
            {
                re[i] = B[i];
                for(j = 0; j < A.n(); j++)
                {
                    if(j!=i)
                    {
                        re[i] = re[i] -A[i][j]*x[j];
                    }
                }
                re[i] = re[i] / A[i][i];
            }

            double maxXerror = 0.0;
            for(i = 0; i < A.n(); i++)
            {
                if(mag(x[i] - re[i]) > maxXerror)
                {
                    maxXerror = mag(x[i] - re[i]);
                }
            }
            if(maxXerror > e)
            {
                return;
            }
            for(i = 0; i < A.n(); i++)
            {
                x[i] = re[i];
            }
        }
    }

    //// The list for old setting, maybe used in future
//    void Foam::directForcing::velocityNextTimeStep(const Foam::scalar& deltaT, const Foam::vector& g, const Foam::volVectorField& U)
//    {
//        //old part
////        const scalar pi = 3.1415926;
////        velocityOld_ = velocity_;
////        vector a = getF_lag_sum();
////        vector b(vector::zero);
////        scalar h = Foam::sqrt(min(mesh_.V().field()));
////        const vectorField& C = mesh_.cellCentres();
////        for(label i = 0; i < ibCellList_.size(); i++)
////        {
////            b += Ix(C[ibCellList_[i]], h) * U[ibCellList_[i]] * getRhof() * pow(h,2)/deltaT;
////        }
////        vector c = getVolume() * (getRhop()-getRhof())*g;
////        vector sum = (a+b+c)*deltaT/getParMass();
////        velocity_.x() += sum.x();
////        velocity_.y() += sum.y();
////        Info << "the velocity is: " << velocity_ << endl;

//        //new part: feng
//        vector vt = velocity_;
//        vector vtold = velocityOld_;
//        velocityOld_ = velocity_;
//        vector a = getF_lag_sum();
//        vector b = getRhof()*getVolume()*(vt-vtold)/deltaT;
//        vector c = getVolume() * (getRhop()-getRhof())*g;
//        vector sum = (a+b+c)*deltaT/getParMass();
//        velocity_.x() += sum.x();
//        velocity_.y() += sum.y();
//        Info << "the velocity is: " << velocity_ << endl;
//    }

//    void Foam::directForcing::testVelocityNextTimeStep
//    (
//        const volVectorField& U,
//        const volVectorField& Uold,
//        const scalar& deltaT,
//        const vector& g
//    )
//    {
//        scalar Pi = constant::mathematical::pi;
//        scalar h = Foam::sqrt(min(mesh_.V().field()));
//        // Calculate A
//        vector A(vector::zero);
//        const vectorField& C = mesh_.cellCentres();

////        for(label i = 0; i < ibCellList_.size(); i++)
////        {
////            A += getRhof() * Ix(C[ibCellList_[i]], h) * U[ibCellList_[i]] * pow(h,2);
////            A -= getRhof() * Ix(C[ibCellList_[i]], h) * Uold[ibCellList_[i]] * pow(h,2);
////        }
//        A += deltaT * getVolume() * (getRhop() - getRhof()) * g;

//        Info << "A part is: " << A << endl;
//        // get lag part(the second term)
//        vector lagpart(vector::zero);
//        scalar  volumeLagrange = Pi*d_*h/number_;   // for 2D
//        for(label i = 0; i < coordinates_.size(); i++)
//        {
//            lagpart += getRhof() * velocityLagrange_[i] * volumeLagrange;
//        }

//        Info << "lagpart is: " << lagpart << endl;
//        // old velocity part;
//        velocityOld_ = velocity_;
//        // first part
//        vector velocitypart = velocityOld_ * getParMass();

//        Info << "velocitypart is: " << velocitypart << endl;
//        velocity_ = (velocitypart + lagpart + A) / (getParMass() + getFluMass());
//    }

//    void Foam::directForcing::displacementNextTimeStep(const Foam::scalar& deltaT)
//    {
//        vector displacement = (velocity_ + velocityOld_)/2* deltaT;
//        //coordinates_ += displacement;
//        center_ += displacement;
//        Info << "the center is: " << center_ << endl;
//        scalar maxa = 0;
//        for(label i = 0; i < coordinates_.size(); i++)
//        {
//            coordinates_[i] += displacement;
//            scalar a = mag(coordinates_[i]-center_);
//            if(a>maxa)
//                maxa = a;
//        }
//        Info << "the max is: " << maxa << endl;
//    }

//    void Foam::directForcing::omegaNextTimeStep(const Foam::scalar& deltaT, const Foam::volVectorField& U)
//    {
//        //old part
////        omega_ += torqueTotal / getMomentInertia() * deltaT;

//        //new part: feng
//        vector ot = omega_;
//        vector otold = omegaOld_;
//        omegaOld_ = omega_;
//        vector a = getTorque_lag_sum();
//        vector b = getFluMomentInertia()*(ot-otold)/deltaT;
//        vector sum = (a+b)*deltaT/getParMomentInertia();
//        Info << "a" << a << endl;
//        Info << "b" << b << endl;
//        omega_.z() += sum.z();

////        //newer part
////        const tensor J = getJ();
////        vector ot = omega_;
////        vector otold = omegaOld_;
////        omegaOld_ = omega_;
////        vector a = getTorque_lag_sum();
////        vector b = getRhof()/getRhop() *( (J & ((ot - otold)/deltaT)) + ot ^ (J & ot));
////        vector sum = (a + b - ot^(J&ot))/J*deltaT;
//////        omega_ = (1+getRhof()/getRhop()) * ot - getRhof()/getRhop() * otold
//////                + a / J;
////        omega_ = omega_ + sum;
//    }

//    void Foam::directForcing::testOmegaNextTimeStep
//    (
//        const volVectorField& U,
//        const volVectorField& Uold,
//        const scalar& deltaT,
//        const vector& g
//    )
//    {
//        scalar Pi = constant::mathematical::pi;
//        scalar h = Foam::sqrt(min(mesh_.V().field()));  //2D
//        // Calculate B
//        vector B(vector::zero);
//        const vectorField& C = mesh_.cellCentres();

//        for(label i = 0; i < ibCellList_.size(); i++)
//        {
//            B += getRhof() * Ix(C[ibCellList_[i]], h) * ((C[ibCellList_[i]] - center_) ^ U[ibCellList_[i]]) * pow(h,2);
//            B -= getRhof() * Ix(C[ibCellList_[i]], h) * ((C[ibCellList_[i]] - center_) ^ Uold[ibCellList_[i]]) * pow(h,2);
//        }

//        // get lag part(the second term)
//        vector lagpart(vector::zero);
//        scalar  volumeLagrange = Pi*d_*h/number_;   // for 2D
//        for(label i = 0; i < coordinates_.size(); i++)
//        {
//            lagpart += getRhof() * ((C[ibCellList_[i]] - center_) ^ velocityLagrange_[i]) * volumeLagrange;
//        }

//        // old velocity part;
//        omegaOld_ = omega_;
//        // first part
//        vector omegapart = omegaOld_ * getParMomentInertia();

//        omega_ = (omegapart + lagpart + B) / (getParMomentInertia() + getFluMomentInertia());

//    }

//    void Foam::directForcing::initialization()
//    {
//        if(motion_)
//        {
//            clearOut();
//        }
//        for(label i = 0; i < forceLagrange_.size(); i++)
//        {
//            forceLagrange_[i] = vector::zero;
//            velocityLagrange_[i] = vector::zero;
//            enthalpyLagrange_[i] = 0.0;
//        }
//    }

//    // calculate pressure force and viscous force
//    void Foam::directForcing::calculateForce
//    (
//        const volVectorField& f1,
//        const scalar& rhoRef_,
//        const scalar& magUInf_,
//        const scalar& Aref_
//    )
//    {
////        const dictionary forceDict = mesh_.foundObject<dictionary>("forceDict");

////        Foam::IFstream ifstream = Foam::IFstream("forceDict");
////        Foam::dictionary root(ifstream);

//        const scalar pi = 3.1415926;
//        // Read particle parameter
////        scalar rhoRef_ = Foam::readScalar(forceDict.lookup("rhoRef"));
////        scalar magUInf_ = Foam::readScalar(forceDict.lookup("magUInf"));
////        scalar Aref_ = Foam::readScalar(forceDict.lookup("Aref"));

//        vector dragDir_(1,0,0);
//        vector liftDir_(0,1,0);

//        const vectorField& C = mesh_.cellCentres();
//        F_lag_sum = vector::zero;
//        torque_lag_sum = vector::zero;
//        scalar volumeLagrange;

//        if(dim_)
//        {
//            volumeLagrange = pi*d_*h_/number_;
//        }
//        else
//        {
//            volumeLagrange = (pi * h_) * (3 * pow(d_, 2) + pow(h_, 2)) / (3 * number_) ;
//        }

//        vector F = vector::zero;
//        vector T = vector::zero;
//        for(label i = 0; i < forceLagrange_/*ibCellList_*/.size(); i++)
//        {
//            vector R =
//                coordinates_[i] - center_;

//            if(dim_)
//            {
//                F += forceMulti_[i] * pow(h_,2);
//                T += (R^forceMulti_[i]) * pow(h_,2);
//            }
//            else
//            {
//                F += /*f1[ibCellList_[i]]**/forceMulti_[i] * pow(h_,3);
//            }
//        }

////        //test part
////        for(label i = 0; i < ibCellList_.size(); i++)
////        {
////            vector R =
////                C[i] - center_;

////            if(dim_)
////            {
////                F += f1[ibCellList_[i]] * pow(h_,2);
////                T += (R^f1[ibCellList_[i]] * pow(h_,2));
////            }
////            else
////            {
////                F += f1[ibCellList_[i]] * pow(h_,3);
////            }
////        }
//        Info << "Flagsum" << F << endl;
//        F *= -1;
//        T *= -1;

//        if(Foam::Pstream::parRun())
//        {
//            Foam::vector forcePerCPU(F);

//            Foam::reduce(forcePerCPU, Foam::sumOp<Foam::vector>());

//            F = forcePerCPU;
//            Foam::Pstream::scatter(F);

//        }

//        F_lag_sum = F;
//        torque_lag_sum = T;
//        Pout << "Flagsum" << F_lag_sum << endl;
//        Pout << "Torquesum" << torque_lag_sum << endl;

//        scalar magu;
//        if(motion_)
//        {
//            if(velocity_ == vector::zero)
//            {
//                magu = 1;
//            }
//            else
//            {
//                magu = mag(velocity_);
//            }
//        }
//        else
//        {
//            magu = magUInf_;
//        }
//        scalar pDyn = 0.5*rhoRef_*magu*magu;

//        coeffs[0] = (F_lag_sum & dragDir_)/(Aref_*pDyn);
//        coeffs[1] = (F_lag_sum & liftDir_)/(Aref_*pDyn);
//        Info << "The drag force is " << coeffs[0] << endl;
//        Info << "The lift force is " << coeffs[1] << endl;
//    }

//    void Foam::directForcing::computeCyclicHits(const scalar L)
//    {
//        labelList patches;
//        labelList walls;

//        forAll(mesh_.boundary(), patchi)
//        {
//            if
//            (
//                (
//                    isA<cyclicFvPatch>(mesh_.boundary()[patchi])
//                 || isA<processorFvPatch>(mesh_.boundary()[patchi])
//                 || isA<processorCyclicFvPatch>(mesh_.boundary()[patchi])
//                )
//             && mesh_.boundary()[patchi].size() != 0
//            )
//            {
//                patches.append(patchi);
//            }
//            else if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
//            {
//                walls.append(patchi);
//            }
//        }

//        forAll(coordinates_, pti)
//        {
//            label celli = mesh_.findCell(coordinates_[pti]);

//            if(celli == -1)
//            {
//                vector delta(0, -L, 0);
//                coordinates_[pti] -= delta;
//            }
//        }

//        label celli = mesh_.findCell(center_);

//        if(celli == -1)
//        {
//            vector delta(0, -L, 0);
//            center_ -= delta;
//        }
//    }

//    bool Foam::directForcing::onMesh()
//    {
//        label nPtsOnMesh = 0;
//        forAll(coordinates_, pti)
//        {
//            vector pt = coordinates_[pti];
//            label celli = mesh_.findCell(pt);

//            if(celli != -1)
//            {
//                nPtsOnMesh++;
//            }
//        }

//        if (nPtsOnMesh == coordinates_.size())
//        {
//            return 1;
//        }
//        else
//        {
//            return 0;
//        }
//    }

//    Foam::label Foam::directForcing::patchHit
//    (
//        const label patchi,
//        vector &R
//    )
//    {
//        forAll(mesh_.boundaryMesh()[patchi], facei)
//        {
//            const vector& faceCentre =
//                mesh_.Cf().boundaryField()[patchi][facei];
//            vector norm =
//                mesh_.Sf().boundaryField()[patchi][facei]
//                /mesh_.magSf().boundaryField()[patchi][facei];

//            R = getCenter() - faceCentre;
//            if (mag(R) <= getDiameter()/2.0 && (velocity_ & norm) > 0)
//            {
//                return facei;
//            }
//        }
//        return -1;
//    }


