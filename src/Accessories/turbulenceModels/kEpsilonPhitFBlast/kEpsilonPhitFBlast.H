/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::RASModels::kEpsilonPhitFBlast

Group
    grpRASTurbulence

Description
    A version of the k-epsilon-phit-f model which is suitable for redistribution.
    Results should be equal to standard version.

SourceFiles
    kEpsilonPhitFBlast.C

SeeAlso
    kEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonPhitFBlast_H
#define kEpsilonPhitFBlast_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class kEpsilonPhitFBlast Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonPhitFBlast
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kEpsilonPhitFBlast(const kEpsilonPhitFBlast&) = delete;

        //- No copy assignment
        void operator=(const kEpsilonPhitFBlast&) = delete;


protected:

    // Protected Data

            Switch includeNu_;

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1a_;
            dimensionedScalar Ceps1b_;
            dimensionedScalar Ceps1c_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Cf1_;
            dimensionedScalar Cf2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar CT_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaPhit_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;

            //- Normalised wall-normal fluctuating velocity scale [-]
            volScalarField phit_;

            //- Elliptic relaxation factor [1/s]
            volScalarField f_;

//             //- Turbulent time scale [s]
//             volScalarField T_;


        // Bounding values

            dimensionedScalar phitMin_;
            dimensionedScalar fMin_;
            dimensionedScalar TMin_;
            dimensionedScalar L2Min_;


    // Protected Member Functions

        //- Update nut with the latest available k, phit, and T
        virtual void correctNut();

        //- Return the turbulent time scale, T
        tmp<volScalarField> Ts() const;

        //- Return the turbulent length scale, L
        tmp<volScalarField> Ls() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonPhitFBlast");


    // Constructors

        //- Construct from components
        kEpsilonPhitFBlast
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kEpsilonPhitFBlast() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k (LUU:Eq. 3)
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for epsilon (LUU:Eq. 4)
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    this->nut_/sigmaEps_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for phit (LUU:Eq. 17)
        tmp<volScalarField> DphitEff() const
        {
            auto tfld =
                tmp<volScalarField>::New("DphitEff", this->nut_/sigmaPhit_);

            if (includeNu_)
            {
                tfld.ref() += this->nu();
            }

            return tfld;
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the normalised wall-normal fluctuating velocity scale field
        virtual tmp<volScalarField> phit() const
        {
            return phit_;
        }

        //- Return the elliptic relaxation factor field
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonPhitFBlast.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
