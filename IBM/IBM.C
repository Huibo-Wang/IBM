#include "IBM.H"

//- Constructors
    Foam::IBM::IBM
    (
        const fvMesh &mesh,
        const dictionary &dict
    )
    :
    mesh_(mesh),
    motion_(readBool(dict.lookup("motion"))),
    dim_(readBool(dict.lookup("dimension"))),
    diameter_(readScalar(dict.lookup("diameter"))),
    center_(dict.lookupOrDefault<vector>("center", vector::zero)),
    velocity_(dict.lookupOrDefault<vector>("velocity", vector::zero)),
    omega_(dict.lookupOrDefault<vector>("omega", vector::zero))
    {
        velocityOld_ = velocity_;
        omegaOld_ = omega_;
    }

    Foam::IBM::~IBM()
    {
    }

//- Member functions

    //- Parameter calculation
    Foam::scalar Foam::IBM::getVolume()
    {
        scalar Pi = constant::mathematical::pi;
        scalar volume;
        if(dim_)
        {
            volume = 1.0/4.0 * Pi * pow(this->getDiameter(), 2);
        }
        else
        {
            volume = 1.0/6.0*Pi*pow3(this->getDiameter());
        }
        return volume;
    }
