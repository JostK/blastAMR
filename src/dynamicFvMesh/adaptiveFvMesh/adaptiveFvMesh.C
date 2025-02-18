/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |    Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "adaptiveFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "dimensionSets.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "sigFpe.H"
#include "fvMeshPolyRefiner.H"
#include "fvMeshHexRefiner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adaptiveFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        adaptiveFvMesh,
        IOobject
    );
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::adaptiveFvMesh::dynamicMeshDictIOobject
(
    const IOobject& io
)
{
    IOobject dictHeader
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    // defaultRegion (region0) gets loaded from constant, other ones get loaded
    // from constant/<regionname>. Normally we'd use polyMesh::dbDir() but we
    // haven't got a polyMesh yet ...
    return IOobject
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        (
            dictHeader.typeHeaderOk<IOdictionary>(true)
          ? IOobject::MUST_READ_IF_MODIFIED
          : IOobject::NO_READ
        ),
        IOobject::NO_WRITE
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::adaptiveFvMesh::readDict()
{
    const dictionary refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );
    refiner_->readDict(refineDict);
    error_->read(refineDict);

    if (!refineDict.found("correctFluxes"))
    {
        return;
    }

    List<Pair<word>> fluxVelocities = List<Pair<word>>
    (
        refineDict.lookup("correctFluxes")
    );
    // Rework into hashtable.
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }
}


void Foam::adaptiveFvMesh::updateMesh(const mapPolyMesh& map)
{
    // Update fluxes
    {
        const labelList& faceMap = map.faceMap();
        const labelList& reverseFaceMap = map.reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        // Estimate number of faces created

        bitSet masterFaces(nFaces());

        forAll(faceMap, facei)
        {
            const label oldFacei = faceMap[facei];

            if (oldFacei >= 0)
            {
                const label masterFacei = reverseFaceMap[oldFacei];

                if (masterFacei < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << facei << endl
                        << abort(FatalError);
                }
                else if (masterFacei != facei)
                {
                    masterFaces.set(masterFacei);
                }
            }
        }

        if (debug)
        {
            Pout<< "Found " << masterFaces.count() << " split faces " << endl;
        }

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );

        // Remove surfaceInterpolation to allow re-calculation on demand
        // This could be done in fvMesh::updateMesh but some dynamicFvMesh
        // might need the old interpolation fields (weights, etc).
        surfaceInterpolation::clearOut();

        forAllIters(fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            surfaceScalarField& phi = *iter();

            if (UName == "NaN")
            {
                Pout<< "Setting surfaceScalarField " << iter.key()
                    << " to NaN" << endl;

                sigFpe::fillNan(phi.primitiveFieldRef());

                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label facei = 0; facei < nInternalFaces(); ++facei)
            {
                const label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    phi[facei] = phiU[facei];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // face-from-masterface
                    phi[facei] = phiU[facei];
                }
            }

            // Recalculate new boundary faces.
            surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

            forAll(phiBf, patchi)
            {
                fvsPatchScalarField& patchPhi = phiBf[patchi];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchi];

                label facei = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    const label oldFacei = faceMap[facei];

                    if (oldFacei == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFacei] != facei)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    ++facei;
                }
            }

            // Update master faces
            for (const label facei : masterFaces)
            {
                if (isInternalFace(facei))
                {
                    phi[facei] = phiU[facei];
                }
                else
                {
                    const label patchi = boundaryMesh().whichPatch(facei);
                    const label i = facei - boundaryMesh()[patchi].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchi];

                    fvsPatchScalarField& patchPhi = phiBf[patchi];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }

    fvMesh::updateMesh(map);

    // Bugfix: update refiner object manually.
    if (refiner_.valid())
    {
        //fvMeshPolyRefiner* polyRefiner = dynamic_cast<fvMeshPolyRefiner*>(refiner_.get());
        //fvMeshHexRefiner* hexRefiner = dynamic_cast<fvMeshHexRefiner*>(refiner_.get());
        //if(polyRefiner != nullptr){
        //    polyRefiner->fvMeshPolyRefiner::updateMesh(map);
        //} else if (hexRefiner != nullptr){
        //    hexRefiner->fvMeshHexRefiner::updateMesh(map);
        //} else {
        //    FatalErrorInFunction
        //        << "fvMeshRefiner type " << refiner_->type() << " is not supported." << endl;
        //}
        refiner_->updateMesh(map);
    }
}


void Foam::adaptiveFvMesh::distribute
(
    const mapDistributePolyMesh& map
)
{
    refiner_->distribute(map);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveFvMesh::adaptiveFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshDict_(dynamicMeshDictIOobject(io)),
    //dynamicBlastFvMesh(io),
    error_(errorEstimator::New(*this, dynamicMeshDict())),
    refiner_(fvMeshRefiner::New(*this, dynamicMeshDict())),
    currentTimeIndex_(-1)
{
    // Read static part of dictionary
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adaptiveFvMesh::~adaptiveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adaptiveFvMesh::mapFields(const mapPolyMesh& mpm)
{
    dynamicFvMesh::mapFields(mpm);

    // Correct old-time volumes for refined/unrefined cells. We know at this
    // point that the points have not moved and the cells have only been split
    // or merged. We hope that dynamicMotionSolverListFvMesh::mapFields
    // does not use old-time volumes ...
    {
        const labelList& cellMap = mpm.cellMap();
        const labelList& reverseCellMap = mpm.reverseCellMap();

        // Each split cell becomes original + 7 additional
        labelList nSubCells(mpm.nOldCells(), 0);

        forAll(cellMap, celli)
        {
            const label oldCelli = cellMap[celli];
            if (oldCelli >= 0 && reverseCellMap[oldCelli] >= 0)
            {
                // Found master cell.
                nSubCells[oldCelli]++;
            }
        }

        // Start off from mapped old volumes
        scalarField correctedV0(this->V0());

        // Correct any split cells
        const auto& V = this->V();
        forAll(cellMap, celli)
        {
            const label oldCelli = cellMap[celli];
            if (oldCelli >= 0 && nSubCells[oldCelli] == 8)
            {
                // Found split cell. Use current volume instead of mapped
                // old one
                correctedV0[celli] = V[celli];
            }
        }

        const auto& cellsFromCells = mpm.cellsFromCellsMap();
        for (const auto& s : cellsFromCells)
        {
            // Reset unsplit cell
            const label celli = s.index();
            correctedV0[celli] = V[celli];
            //? Or sum up old volumes?
            //const labelList& oldCells = s.masterObjects();
            //for (const label oldCelli : oldCells)
            //{
            //}
        }

        setV0().field() = correctedV0;
    }
    
    // JK Also correct old-old-time volumes for refined/unrefined cells. We know at this
    // point that the points have not moved and the cells have only been split
    // or merged. We hope that dynamicMotionSolverListFvMesh::mapFields
    // does not use old-time volumes ...
    {
        const labelList& cellMap = mpm.cellMap();
        const labelList& reverseCellMap = mpm.reverseCellMap();

        // Each split cell becomes original + 7 additional
        labelList nSubCells(mpm.nOldCells(), 0);

        forAll(cellMap, celli)
        {
            const label oldCelli = cellMap[celli];
            if (oldCelli >= 0 && reverseCellMap[oldCelli] >= 0)
            {
                // Found master cell.
                nSubCells[oldCelli]++;
            }
        }

        // Start off from mapped old volumes
        scalarField correctedV00(this->V00());

        // Correct any split cells
        const auto& V = this->V();
        forAll(cellMap, celli)
        {
            const label oldCelli = cellMap[celli];
            if (oldCelli >= 0 && nSubCells[oldCelli] == 8)
            {
                // Found split cell. Use current volume instead of mapped
                // old one
                correctedV00[celli] = V[celli];
            }
        }

        const auto& cellsFromCells = mpm.cellsFromCellsMap();
        for (const auto& s : cellsFromCells)
        {
            // Reset unsplit cell
            const label celli = s.index();
            correctedV00[celli] = V[celli];
            //? Or sum up old volumes?
            //const labelList& oldCells = s.masterObjects();
            //for (const label oldCelli : oldCells)
            //{
            //}
        }
        
        // TODO setV00() does not exist, so you left me no choice
        const_cast<Foam::volScalarField::Internal&>(V00()).field() = correctedV00;
    }

    // Correct the flux for injected faces - these are the faces which have
    // no correspondence to the old mesh (i.e. added without a masterFace, edge
    // or point). An example is the internal faces from hexRef8.
    {
        const labelList& faceMap = mpm.faceMap();

        // JK TODO 
        mapNewInternalFaces<scalar>(this->Sf(), this->magSf(), faceMap);
        mapNewInternalFaces<vector>(this->Sf(), this->magSf(), faceMap);
        
//         mapNewInternalFaces<scalar>(faceMap);
//         mapNewInternalFaces<vector>(faceMap);

        // No oriented fields of more complex type
        mapNewInternalFaces<sphericalTensor>(faceMap);
        mapNewInternalFaces<symmTensor>(faceMap);
        mapNewInternalFaces<tensor>(faceMap);
    }
}

bool Foam::adaptiveFvMesh::firstUpdate()
{
    if (currentTimeIndex_ >= time().timeIndex())
    {
        return false;
    }
    currentTimeIndex_ = time().timeIndex();
    return true;
}

bool Foam::adaptiveFvMesh::update()
{
    if (!firstUpdate()) return false;
    bool changed = (refiner_->canRefine(true) || refiner_->canUnrefine(true)) && refine();
    reduce(changed, orOp<bool>());

    return changed;
}


bool Foam::adaptiveFvMesh::refine()
{
    // Re-read dictionary. Chosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    readDict();

    //- Update error
    error_->update();
    error_->error().correctBoundaryConditions();
    label nProtected = error_->protectPatches();
    Info << "Protecting " << returnReduce(nProtected, sumOp<label>())
        << " cells next to requested boundary patches." << endl;
    return refiner_->refine(error_->error(), error_->maxRefinement());
}


bool Foam::adaptiveFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    return
        dynamicFvMesh::writeObject(streamOpt, valid);
     //&& refiner_->write();
}


// ************************************************************************* //
