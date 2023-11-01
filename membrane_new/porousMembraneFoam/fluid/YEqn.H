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

    pMap[rI] = p;//added

    // get the source term for interregion mass transfer
    Info<< "\nMapping pressure field " << endl;
    MapSubMesh<eqOp>
    (
        pMap[rI],   //fieldTarget,
        reactionFluid[1-rI].thermo().p(),    //fieldSource,
        map,
        mapOrder
    );
/*
X_.setSize(Y.size());

        forAll(Y, i)
        {
            X_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "X_" + Y[i].name(),
                        runTime.time().timeName(),
                        fluidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar(dimless, Zero)
                )
            );
        }
*/
    MWavgMap[rI] = MWavg[rI];
    Info<< "\nMapping average molecular weight field " << endl;
    MapSubMesh<eqOp>
    (
        MWavgMap[rI],   //fieldTarget,
        reactionFluid[1-rI].thermo().W(),    //fieldSource,
        map,
        mapOrder
    );

    Smk[rI] == dimensionedScalar("Smk", dimensionSet(1, -3, -1 , 0, 0), Zero) ;
    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
        Info<< "\nSolving for region number " <<  rI  << " (" <<fluidRegions[rI].name()  << ")" << endl;//rI??
        Info<< "& Species number " <<  i  << endl;
     
            volScalarField& Yi = Y[i];

            PtrList<volScalarField>& Ysrc = reactionFluid[1-rI].thermo().composition().Y();//Composition Y

//Info<< "\nmax MWAvg[" << rI << "] = " << max(MWavg[rI].internalField()) << endl;
//Info<< "\nmin MWAvg[" << rI << "] = " << min(MWavg[rI].internalField()) << endl;
//Info<< "\nMWAvg[" << rI << "] = " << MWavg[rI].internalField() << endl;
//Info<< "\nMWAvg[" << 1-rI << "](MWavgMap" << rI << ") = " << MWavgMap[rI].internalField() << endl;

/*            dimensionedScalar smMw
           (
                dimensionedScalar("smMw", dimensionSet(1, 0, 0 , 0, -1), 1.0e-5)
           );
           Info<< "\nsmMw = " << smMw << endl;*/
            YMap[rI] = Yi*MWavg[rI]/MWavgMap[rI];

            Info<< "\nMapping Y field " << endl;
            MapSubMesh<eqOp>
            (
                YMap[rI],   //fieldTarget,
                Ysrc[i],    //fieldSource,
                map,
                mapOrder
            );

            dimensionedScalar Qi("Qi",dimensionSet(-1, -1, 1, 0, 1, 0),Q[i]);
            //Info<< "\nQi = " << Qi << endl;
            dimensionedScalar sc(AbyV*Qi);// /thermo.composition().W(i)
            volScalarField    YSu(sc*pMap[rI]*MWavgMap[rI]*YMap[rI]);
            //Info<< "\nMWAvg[" << 1-rI << "] = " << MWavg[1-rI].internalField() << endl;
            //Info<< "\nYSu[" << i << "] of region " << rI << " = " << YSu.internalField() << endl;
            volScalarField    YSp(sc*p*MWavg[rI]);
            //Info<< "\nMWAvg[" << rI << "] = " << MWavg[rI].internalField() << endl;
            //Info<< "\nYSp[" << i << "] of region " << rI << " = " << YSp.internalField() << endl;

            Info<< "\nSolving YEqn " << endl;
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

            Smk[rI] += (YSu - YSp*Yi);

            //Info<< "\nSmk[" << rI << "] =  " << Smk.internalField() << endl;
/*
        const dimensionedScalar Wi
        (
            dimMass/dimMoles,
            thermo.composition().W(i)
        );

        X_[i] = MWavg[rI]*Y[i]/Wi;
*/
//calcMoleFractions();
        Info<< "\nYEqn solved " << endl;
        }
    }

    if (Y.size())
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
}

