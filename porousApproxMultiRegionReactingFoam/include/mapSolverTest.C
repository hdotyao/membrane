/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
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

Application
    chtMultiRegionFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "regionProperties.H"

#include "GeometricField.H"
#include "meshToMesh0.H"
#include "processorFvPatch.H"
// #include "MapMeshes.H"


template<class Type, class CombineOp>
void MapVolFields
(
    volScalarField&   fieldTarget,
    const volScalarField& fieldSource,
    const meshToMesh0& meshToMesh0Interp,
    const meshToMesh0::order& mapOrder,
    const CombineOp& cop
)
{


    // Interpolate field
    meshToMesh0Interp.interpolate
    (
        fieldTarget,
        fieldSource,
        mapOrder,
        cop
    );

}


template<template<class> class CombineOp>
void MapSubMesh
(
    volScalarField&   fieldTarget,
    const volScalarField& fieldSource,
    const meshToMesh0& meshToMesh0Interp,
    const meshToMesh0::order& mapOrder
)
{


    // Map volFields
    // ~~~~~~~~~~~~~
    MapVolFields<scalar>
    (
        fieldTarget,
        fieldSource,
        meshToMesh0Interp,
        mapOrder,
        CombineOp<scalar>()
    );

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{


    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"


    meshToMesh0::order mapOrder = meshToMesh0::INTERPOLATE;


/*     if (args.found("mapMethod"))
    {
        const word mapMethod(args["mapMethod"]);
        if (mapMethod == "mapNearest")
        {
            mapOrder = meshToMesh0::MAP;
        }
        else if (mapMethod == "interpolate")
        {
            mapOrder = meshToMesh0::INTERPOLATE;
        }
        else if (mapMethod == "cellPointInterpolate")
        {
            mapOrder = meshToMesh0::CELL_POINT_INTERPOLATE;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown mapMethod " << mapMethod << ". Valid options are: "
                << "mapNearest, interpolate and cellPointInterpolate"
                << exit(FatalError);
        }

        Info<< "Mapping method: " << mapMethod << endl;
    } */



    HashTable<word> patchMap;
    wordList cuttingPatches;

    IOdictionary mapFieldsDict
    (
        IOobject
        (
            "mapFieldsDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    mapFieldsDict.readEntry("patchMap", patchMap);
    mapFieldsDict.readEntry("cuttingPatches", cuttingPatches);
    

        // Create the interpolation scheme
    meshToMesh0 meshToMesh0Interp
    (
        fluidRegions[0],  //meshSource,
        fluidRegions[1],  //meshTarget,
        patchMap,
        cuttingPatches
    );

    MapSubMesh<eqOp>
    (
        Tmap[1],   //fieldTarget,
        T[0],   //fieldSource,
        meshToMesh0Interp,
        mapOrder
    );

    // Write field
    Tmap[1].write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
