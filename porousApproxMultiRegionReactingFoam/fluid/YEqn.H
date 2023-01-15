tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);


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

    reaction.correct();
    Qdot = reaction.Qdot();
    volScalarField Yt
    (
        IOobject("Yt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Yt", dimless, 0)
    );

    label rI = i ;
    const  meshToMesh0::order mapOrder = meshToMesh0::INTERPOLATE;
    const Foam::meshToMesh0& map = meshToMesh0Interp[rI];

    // get the source term for interregion mass transfer
    MapSubMesh<eqOp>
    (
        pMap[rI],   //fieldTarget,
        reactionFluid[1-rI].thermo().p(),    //fieldSource,
        map,
        mapOrder
    );
            
    forAll(Y, i)
    {

        if (i != inertIndex && composition.active(i))
        {

            volScalarField& Yi = Y[i];

            PtrList<volScalarField>& Ysrc = reactionFluid[1-rI].thermo().composition().Y();//Composition Y

            MapSubMesh<eqOp>
            (
                YMap[rI],   //fieldTarget,
                Ysrc[i],    //fieldSource,
                map,
                mapOrder
            );

            dimensionedScalar Qi("Qi",dimensionSet(0,-2,-1,0,1,0,0),Q[i]);
            dimensionedScalar sc(AbyV*Qi/thermo.composition().W(i));
            volScalarField    YSu(sc*pMap[rI]*MWavg[1-rI]*YMap[rI]);
            volScalarField    YSp(sc*p*MWavg[rI]);
    
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence.muEff(), Yi)
              + fvm::Sp(YSp,Yi)
             ==
                reaction.R(Yi)
              + fvOptions(rho, Yi)
              + YSu 
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;

            Smk += YSu - YSp*Yi;

        }
    }

    if (Y.size())
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
}
