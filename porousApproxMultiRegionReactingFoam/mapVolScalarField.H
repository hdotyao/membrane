#include "GeometricField.H"
#include "meshToMesh0.H"
#include "processorFvPatch.H"



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

