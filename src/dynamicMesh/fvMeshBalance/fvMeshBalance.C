/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "fvMeshBalance.H"
#include "decompositionMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "RefineBalanceMeshObject.H"
//#include "parcelCloud.H"
#include "preserveFaceZonesConstraint.H"
#include "singleProcessorFaceSetsConstraint.H"
#include "preservePatchesConstraint.H"
#include "preserveBafflesConstraint.H"

#include "faMeshesRegistry.H"
#include "parFvFieldDistributor.H"
#include "parPointFieldDistributor.H"
#include "fieldsDistributor.H"

using namespace Foam::decompositionConstraints;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshBalance, 0);

bool fvMeshBalance::balancing = false;
}

bool Foam::fvMeshBalance::isBalancing()
{
    return balancing;
}

bool Foam::fvMeshBalance::isBalancing(const polyMesh& mesh)
{
    if (balancing)
    {
        if (!mesh.db().foundObject<polyMesh>(mesh.name()))
        {
            return false;
        }
        if (&mesh == &mesh.db().lookupObject<polyMesh>(mesh.name()))
        {
            return !Pstream::parRun();
        }
    }
    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshBalance::fvMeshBalance(fvMesh& mesh)
:
    mesh_(mesh),
    decompositionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    modified_(false),
    constraintsDict_(decompositionDict_.subDictPtr("constraints")),
    //preserveFaceZonesDict_(nullptr),
    //singleProcessorFaceSetsDict_(nullptr),
    //preservePatchesDict_(nullptr),
    //preserveBafflesDict_(nullptr),
    distributor_(mesh_),
    balance_(true),
    allowableImbalance_(0.2)
{
    if (!constraintsDict_)
    {
        decompositionDict_.add("constraints", dictionary());
        constraintsDict_ = decompositionDict_.subDictPtr("constraints");
    }
    else
    {
        wordList toc(decompositionDict_.toc());
        forAll(toc, i)
        {
            if (decompositionDict_.isDict(toc[i]))
            {
                dictionary& dict(decompositionDict_.subDict(toc[i]));
                word type(dict.lookupOrDefault<word>("type", "none"));

                //if (type == preserveFaceZonesConstraint::typeName)
                //{
                //    preserveFaceZonesDict_ = &dict;
                //}
                //else if
                //(
                //    type == singleProcessorFaceSetsConstraint::typeName
                //)
                //{
                //    singleProcessorFaceSetsDict_ = &dict;
                //}
                //else if (type == preservePatchesConstraint::typeName)
                //{
                //    preservePatchesDict_ = &dict;
                //}
                //else if (type == preserveBafflesConstraint::typeName)
                //{
                //    preserveBafflesDict_ = &dict;
                //}
            }
        }
    }
}


Foam::fvMeshBalance::fvMeshBalance
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    decompositionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    modified_(false),
    constraintsDict_(decompositionDict_.subDictPtr("constraints")),
    //preserveFaceZonesDict_(nullptr),
    //singleProcessorFaceSetsDict_(nullptr),
    //preservePatchesDict_(nullptr),
    //preserveBafflesDict_(nullptr),
    distributor_(mesh_),
    balance_(false),
    allowableImbalance_(0.2)
{
    if (!constraintsDict_)
    {
        decompositionDict_.add("constraints", dictionary());
        constraintsDict_ = decompositionDict_.subDictPtr("constraints");
    }
    else
    {
        wordList toc(decompositionDict_.toc());
        forAll(toc, i)
        {
            if (decompositionDict_.isDict(toc[i]))
            {
                dictionary& dict(decompositionDict_.subDict(toc[i]));
                word type(dict.lookupOrDefault<word>("type", "none"));

                //if (type == preserveFaceZonesConstraint::typeName)
                //{
                //    preserveFaceZonesDict_ = &dict;
                //}
                //else if
                //(
                //    type == singleProcessorFaceSetsConstraint::typeName
                //)
                //{
                //    singleProcessorFaceSetsDict_ = &dict;
                //}
                //else if (type == preservePatchesConstraint::typeName)
                //{
                //    preservePatchesDict_ = &dict;
                //}
                //else if (type == preserveBafflesConstraint::typeName)
                //{
                //    preserveBafflesDict_ = &dict;
                //}
            }
        }
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshBalance::~fvMeshBalance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshBalance::read(const dictionary& balanceDict)
{
    if (!Pstream::parRun())
    {
        balance_ = false;
        return;
    }

    balance_ = balanceDict.lookupOrDefault("balance", true);

    if (!balance_)
    {
        return;
    }

    // Change decomposition method if entry is present
    if (balanceDict.found("method"))// || balanceDict.found("decomposer"))
    {
        word method(balanceDict.lookup("method"));
    //        (
    //            {"method"}
    //             {"decomposer", "method"},
    //        );
        decompositionDict_.set("method", method);

        // Remove optional coeffs dictionary since it would override entries and
        // not necessarily be overridden
        if (decompositionDict_.isDict(method + "Coeffs"))
        {
            decompositionDict_.remove(method + "Coeffs");
        }
    }
    decompositionDict_ <<= balanceDict;

    balanceDict.readIfPresent("allowableImbalance", allowableImbalance_);
}


void Foam::fvMeshBalance::addConstraint(const word& dictName, const dictionary& dict)
{
    // Add constraints dictionary
    if (!constraintsDict_->found(dictName))
    {
        modified_ = true;

        constraintsDict_->set(dictName, dict);

        // We need to apply the updated constraints, and this can only
        // be done when the decomposer is created
        if (decomposer_.valid())
        {
            decomposer_.clear();
        }
    }
}


//void Foam::fvMeshBalance::preserveFaceZone(const wordRe& zoneName)
//{
//    if (!preserveFaceZonesDict_)
//    {
//        constraintsDict_->add("faceZones", dictionary());
//        preserveFaceZonesDict_ = constraintsDict_->subDictPtr("faceZones");
//        preserveFaceZonesDict_->set
//        (
//            "type",
//            preserveFaceZonesConstraint::typeName
//        );
//    }
//    wordReList zones
//    (
//        preserveFaceZonesDict_->lookupOrDefault("zones", wordReList())
//    );
//
//    bool exists = false;
//    forAll(zones, i)
//    {
//        if (zones[i].match(zoneName))
//        {
//            exists = true;
//            break;
//        }
//    }
//    if (!exists)
//    {
//        modified_ = true;
//        zones.append(zoneName);
//        preserveFaceZonesDict_->set("zones", zones);
//
//        if (decomposer_.valid())
//        {
//            decomposer_.clear();
//        }
//    }
//}


//void Foam::fvMeshBalance::singleProcessorFaceSet
//(
//    const word& setName,
//    const label proc
//)
//{
//    //if (!preserveFaceZonesDict_)
//    //{
//    //    constraintsDict_->add("singleProcessorFaceSets", dictionary());
//    //    singleProcessorFaceSetsDict_ =
//    //        constraintsDict_->subDictPtr("singleProcessorFaceSets");
//    //    singleProcessorFaceSetsDict_->set
//    //    (
//    //        "type",
//    //        singleProcessorFaceSetsConstraint::typeName
//    //    );
//    //}
//    List<Tuple2<word, label>> setNameAndProcs
//    (
//        singleProcessorFaceSetsDict_->lookupOrDefault
//        (
//            "singleProcessorFaceSets",
//            List<Tuple2<word, label>>()
//        )
//    );
//
//    bool exists = false;
//    forAll(setNameAndProcs, i)
//    {
//        if (setName == setNameAndProcs[i].first())
//        {
//            exists = true;
//            break;
//        }
//    }
//    if (!exists)
//    {
//        modified_ = true;
//        setNameAndProcs.append(Tuple2<word, label>(setName, proc));
//        singleProcessorFaceSetsDict_->set
//        (
//            "setNameAndProcs",
//            setNameAndProcs
//        );
//
//        if (decomposer_.valid())
//        {
//            decomposer_.clear();
//        }
//    }
//}


//void Foam::fvMeshBalance::preservePatch(const wordRe& patchName)
//{
//    if (!preservePatchesDict_)
//    {
//        constraintsDict_->add("patches", dictionary());
//        preservePatchesDict_ = constraintsDict_->subDictPtr("patches");
//        preservePatchesDict_->set
//        (
//            "type",
//            preservePatchesConstraint::typeName
//        );
//    }
//    wordReList patches
//    (
//        preservePatchesDict_->lookupOrDefault("patches", wordReList())
//    );
//
//    bool exists = false;
//    forAll(patches, i)
//    {
//        if (patches[i].match(patchName))
//        {
//            exists = true;
//            break;
//        }
//    }
//    if (!exists)
//    {
//        modified_ = true;
//        patches.append(patchName);
//        preservePatchesDict_->set("patches", patches);
//
//        if (decomposer_.valid())
//        {
//            decomposer_.clear();
//        }
//    }
//}


//void Foam::fvMeshBalance::preserveBaffles()
//{
//    if (!preserveBafflesDict_)
//    {
//        modified_ = true;
//        constraintsDict_->add("baffles", dictionary());
//        preserveBafflesDict_ = constraintsDict_->subDictPtr("baffles");
//        preserveBafflesDict_->set
//        (
//            "type",
//            preserveBafflesConstraint::typeName
//        );
//    }
//}


void Foam::fvMeshBalance::makeDecomposer() const
{
    if (decomposer_.valid())
    {
        FatalErrorInFunction
            << "Decomposer already set" << endl
            << abort(FatalError);
    }
    decomposer_ = decompositionMethod::New(decompositionDict_);

    returnReduce(1, maxOp<label>());
    if (!decomposer_->parallelAware())
    {
        WarningInFunction
            << "You have selected decomposition method "
            << decomposer_->typeName
            << " which is not parallel aware." << endl;
    }
}


Foam::decompositionMethod& Foam::fvMeshBalance::decomposer() const
{
    if (!decomposer_.valid())
    {
        makeDecomposer();
    }
    return decomposer_();
}


bool Foam::fvMeshBalance::canBalance() const
{
    if (!balance_)
    {
        return false;
    }

    //First determine current level of imbalance - do this for all
    // parallel runs with a changing mesh, even if balancing is disabled
    label nGlobalCells = returnReduce(mesh_.nCells(), sumOp<label>());
    scalar idealNCells =
        scalar(nGlobalCells)/scalar(Pstream::nProcs());
    scalar localImbalance = mag(scalar(mesh_.nCells()) - idealNCells);
    scalar maxImbalance = returnReduce(localImbalance, maxOp<scalar>());
    scalar maxImbalanceRatio = maxImbalance/idealNCells;

    Info<<"Maximum imbalance = " << 100*maxImbalanceRatio << " %" << endl;

    if (debug)
    {
        Pout<< " localImbalance = "
            << 100.0*localImbalance/idealNCells << "%, "
            << "nCells = " << mesh_.nCells()
            << endl;
    }

    //If imbalanced, construct weighted coarse graph (level 0) with node
    // weights equal to their number of subcells. This partitioning works
    // as long as the number of level 0 cells is several times greater than
    // the number of processors.
    if (maxImbalanceRatio < allowableImbalance_)
    {
        return false;
    }

    // Decompose the mesh with uniform weights
    // The refinementHistory constraint is applied internally
    distribution_ = decomposer().decompose
    (
        mesh_,
        scalarField(mesh_.nCells(), 1.0)
    );

    // Check if distribution will improve anything
    labelList procLoadNew(Pstream::nProcs(), 0);
    forAll(distribution_, celli)
    {
        procLoadNew[distribution_[celli]]++;
    }
    reduce(procLoadNew, sumOp<labelList>());
    if (min(procLoadNew) == 0)
    {
        DebugInfo
            << "New distribtion results in a load of 0. Skipping" << endl;
        return false;
    }
    scalar averageLoadNew
    (
        scalar(sum(procLoadNew))/scalar(Pstream::nProcs())
    );
    scalar maxDevNew(max(mag(procLoadNew - averageLoadNew))/averageLoadNew);

    if (maxDevNew > maxImbalanceRatio*0.99)
    {
        Info
            << "    Not balancing because the new distribution does" << nl
            << "    not improve the load. Skipping" << nl
            << "    old imbalance: " << maxImbalanceRatio << nl
            << "    new imbalance: " << maxDevNew << nl
            << endl;
        return false;
    }

    return true;
}


Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::fvMeshBalance::distribute()
{
    PtrList<pointScalarField> pointScalarFields;
    PtrList<pointVectorField> pointVectorFields;
    PtrList<pointTensorField> pointTensorFields;
    PtrList<pointSphericalTensorField> pointSphTensorFields;
    PtrList<pointSymmTensorField> pointSymmTensorFields;

    {
    const Foam::HashTable<Foam::pointScalarField*>& pointFieldsTable(mesh_.lookupClass<Foam::pointScalarField>());
    pointScalarFields.setSize(pointFieldsTable.size());
    Foam::label i = 0;
    for (const auto* fieldPtr : pointFieldsTable)  
    {
        pointScalarFields.set(i++,  const_cast<Foam::pointScalarField*>(fieldPtr));  
    }
    }
    {
    const Foam::HashTable<Foam::pointVectorField*>& pointFieldsTable(mesh_.lookupClass<Foam::pointVectorField>());
    pointVectorFields.setSize(pointFieldsTable.size());
    Foam::label i = 0;
    for (const auto* fieldPtr : pointFieldsTable)  
    {
        pointVectorFields.set(i++,  const_cast<Foam::pointVectorField*>(fieldPtr));  
    }
    }
    {
    const Foam::HashTable<Foam::pointTensorField*>& pointFieldsTable(mesh_.lookupClass<Foam::pointTensorField>());
    pointTensorFields.setSize(pointFieldsTable.size());
    Foam::label i = 0;
    for (const auto* fieldPtr : pointFieldsTable)  
    {
        pointTensorFields.set(i++,  const_cast<Foam::pointTensorField*>(fieldPtr));  
    }
    }
    {
    const Foam::HashTable<Foam::pointSphericalTensorField*>& pointFieldsTable(mesh_.lookupClass<Foam::pointSphericalTensorField>());
    pointSphTensorFields.setSize(pointFieldsTable.size());
    Foam::label i = 0;
    for (const auto* fieldPtr : pointFieldsTable)  
    {
        pointSphTensorFields.set(i++,  const_cast<Foam::pointSphericalTensorField*>(fieldPtr));  
    }
    }
    {
    const Foam::HashTable<Foam::pointSymmTensorField*>& pointFieldsTable(mesh_.lookupClass<Foam::pointSymmTensorField>());
    pointSymmTensorFields.setSize(pointFieldsTable.size());
    Foam::label i = 0;
    for (const auto* fieldPtr : pointFieldsTable)  
    {
        pointSymmTensorFields.set(i++,  const_cast<Foam::pointSymmTensorField*>(fieldPtr));  
    }
    }
    
    // Check processors have meshes
    // - check for 'faces' file (polyMesh)
    // - check for 'faceLabels' file (faMesh)
//     boolList volMeshOnProc;
//     volMeshOnProc.setSize(UPstream::nProcs(), true);
//     boolList areaMeshOnProc;

    // Read handler on processors with a volMesh
//     refPtr<fileOperation> volMeshReadHandler = fileOperation::New(fileHandler(), volMeshOnProc, true);
//     refPtr<fileOperation> noReadHandler;
    
//     // All check if can read 'faces' file
//     volMeshOnProc = haveMeshFile
//     (
//         runTime,
//         volMeshMasterInstance/volMeshSubDir,
//         "faces"
//     );
//     
//     // Create 0 sized mesh to do all the generation of zero sized
//     // fields on processors that have zero sized meshes. Note that this is
//     // only necessary on master but since polyMesh construction with
//     // Pstream::parRun does parallel comms we have to do it on all
//     // processors
//     autoPtr<fvMeshSubset> subsetterPtr;
// 
//     // Missing a volume mesh somewhere?
//     if (volMeshOnProc.found(false))
//     {
//         // A zero-sized mesh with boundaries.
//         // This is used to create zero-sized fields.
//         subsetterPtr.reset(new fvMeshSubset(mesh, zero{}));
//         subsetterPtr().subMesh().init(true);
//         subsetterPtr().subMesh().globalData();
//         subsetterPtr().subMesh().tetBasePtIs();
//         subsetterPtr().subMesh().geometricD();
//     }
//     
Info << "HALLO 1" << endl;
    
    // Self-contained pointMesh for reading pointFields
    const pointMesh oldPointMesh(mesh_);
    refPtr<fileOperation> noWriteHandler; //TODO this is the problem
    parPointFieldDistributor pointDistributor
    (
        oldPointMesh,   // source mesh
        false,          // savePoints=false (ie, delay until later)
//         false           // Do not write
        noWriteHandler    // Do not write
    );
Info << "HALLO 2" << endl;
    IOobjectList objects = IOobjectList(mesh_, mesh_.time().timeName());
    
Info << "HALLO 3 objects"<< objects << endl;
    // pointFields
//     label nPointFields = 0;
/*
        #define doFieldReading(Storage)                                       \
        {                                                                     \
            fieldsDistributor::readFields                                     \
            (                                                                 \
                volMeshOnProc, noReadHandler, oldPointMesh,              \
                subsetterPtr, objects, Storage,                               \
                true                                \
            );                                                                \
            nPointFields += Storage.size();                                   \
        } 
*/
/*        #define doFieldReading(Storage)                                       \
        {                                                                     \
            fieldsDistributor::readFields                                     \
            (                                                                 \
                oldPointMesh,              \
                objects, Storage                               \
            );                                                                \
            nPointFields += Storage.size();                                   \
        }
*/
Info << "HALLO 4" << endl;
//     doFieldReading(pointScalarFields);
//     doFieldReading(pointVectorFields);
//     doFieldReading(pointSphTensorFields);
//     doFieldReading(pointSymmTensorFields);
//     doFieldReading(pointTensorFields);
//     #undef doFieldReading
// Info << "HALLO 5 nPointFields " <<  nPointFields << endl;
    // TODO only needed if pointFields are present
    pointDistributor.saveMeshPoints();
Info << "HALLO 6" << endl;
    
    
    
    
    
    
    
    //Correct values on all coupled patches
    correctBoundaries<volScalarField>();
    correctBoundaries<volVectorField>();
    correctBoundaries<volSphericalTensorField>();
    correctBoundaries<volSymmTensorField>();
    correctBoundaries<volTensorField>();

//     correctBoundaries<pointScalarField>();
//     correctBoundaries<pointVectorField>();
//     correctBoundaries<pointSphericalTensorField>();
//     correctBoundaries<pointSymmTensorField>();
//     correctBoundaries<pointTensorField>();

    blastMeshObject::preDistribute<fvMesh>(mesh_);
    
    // If faMeshesRegistry exists, it is also owned by the polyMesh and will
    // be destroyed by clearGeom() in fvMeshDistribute::distribute()
    //
    // Rescue faMeshesRegistry from destruction by temporarily moving
    // it to be locally owned.
    std::unique_ptr<faMeshesRegistry> faMeshesRegistry_saved
    (
        faMeshesRegistry::Release(mesh_)
    );

    Info<< "Distributing the mesh ..." << endl;
    balancing = true;
    autoPtr<mapDistributePolyMesh> map =
        distributor_.distribute(distribution_);
    balancing = false;
    
    // Restore ownership onto the polyMesh
    faMeshesRegistry::Store(std::move(faMeshesRegistry_saved));
    
    Info << "Successfully distributed mesh" << endl;
    label procLoadNew(mesh_.nCells());
    label overallLoadNew(returnReduce(procLoadNew, sumOp<label>()));
    scalar averageLoadNew(overallLoadNew/scalar(Pstream::nProcs()));

    scalar maxDevNew
    (
        returnReduce(mag(procLoadNew - averageLoadNew), maxOp<scalar>())
    );

    Info << "New max imbalance: " << maxDevNew/averageLoadNew*100.0 << "%"
        << endl;

    if (debug)
    {
        Pout<< " localImbalance = "
            << mag(procLoadNew - averageLoadNew)*100.0/averageLoadNew << "%, "
            << "Cells = " << procLoadNew
             << endl;
    }

    blastMeshObject::distribute<fvMesh>(mesh_, map());

    // TODO only needed if pointFields are present
    
    // Construct new pointMesh from distributed mesh
    const pointMesh& newPointMesh = pointMesh::New(mesh_);
    pointDistributor.resetTarget(newPointMesh, map());
    pointDistributor.distributeAndStore(pointScalarFields);
    pointDistributor.distributeAndStore(pointVectorFields);
    pointDistributor.distributeAndStore(pointSphTensorFields);
    pointDistributor.distributeAndStore(pointSymmTensorFields);
    pointDistributor.distributeAndStore(pointTensorFields);
    
    //Correct values on all coupled patches
    correctBoundaries<volScalarField>();
    correctBoundaries<volVectorField>();
    correctBoundaries<volSphericalTensorField>();
    correctBoundaries<volSymmTensorField>();
    correctBoundaries<volTensorField>();

//     correctBoundaries<pointScalarField>();
//     correctBoundaries<pointVectorField>();
//     correctBoundaries<pointSphericalTensorField>();
//     correctBoundaries<pointSymmTensorField>();
//     correctBoundaries<pointTensorField>();

    return map;
}

bool Foam::fvMeshBalance::write(const bool write) const
{
    if
    (
        balance_ && modified_ && write &&
        decompositionDict_.lookupOrDefault("writeDecomposeDict", false)
    )
    {
        modified_ = false;
        IOdictionary decomposeParDict
        (
            IOobject
            (
                "decomposeParDict",
                mesh_.time().caseSystem(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            decompositionDict_
        );
        return decomposeParDict.regIOobject::write();
    }
    return true;
}


// ************************************************************************* //
