tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);
#include "IOdictionary.H"

IOdictionary membraneProperties//Dict Variable Name
(
    IOobject
    (
        "membraneProperties",
        runTime.constant(),   //Dict under Constant
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

if (Y.size())
{
    mvConvection = tmp<fv::convectionScheme<scalar>>
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );
}

{
    dimensionedScalar AbyV
    (
        membraneProperties.lookup("AbyV")
    );

    dimensionedScalar Qi
    (
        membraneProperties.lookup("Qi")
    );

    reaction.correct();
    Qdot = reaction.Qdot();
    volScalarField Yt
    (
        IOobject("Yt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Yt", dimless, 0)
    );

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            //source calculation ptrlist
            PtrList<volScalarField>& Sik = AbyV*W(i)*Qi*( p*thermo.W()*Yi/W(i) - p_Map* thermo.W()*Y_Map/W(i) );

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence.muEff(), Yi)
             ==
                reaction.R(Yi)
              + fvOptions(rho, Yi)
              + Sik
            );

//mol weight look up: /usr/OpenFOAM-v2206/tutorials/combustion/reactingFoam/laminar/counterFlowFlame2D

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    if (Y.size())
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
}
