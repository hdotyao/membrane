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

    dimensionedScalar WMi;
    if (i==0) {WMi = WM1;}
    else if (i==1){WMi = WM2;}
    else if (i==2){WMi = WM3;}

    volScalarField WMavg;

    forAll(Y[0],i)//cell
    {
        WMavg[i] = 0;//species
        forAll(Y,j) 
        {
            WMavg[i] += WMi(j)*Y[j][i];//todo
        }
    }

    forAll(Y, i)
    {
        dimensionedScalar Qi;
        if (i==0) {Qi = Q1;}
        else if (i==1) {Qi = Q2;}
        else if (i==2) {Qi = Q3;}

        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            volScalarField Sik = AbyV*WMi*Qi*(p*WMavg*Yi/WMi - p_Mapped* WMavg*Y_Mapped/WMi);

            forAll(Smk,i) //No need to forAll Smk
            {
                Smk[i] += Sik[i];
            }

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence.muEff(), Yi)
             ==
                reaction.R(Yi)
              + fvOptions(rho, Yi)
              + fvm::Su(Sik,Yi) //fvm::Sp, fvm::Susp ?
            );

//mol weight look: /usr/OpenFOAM-v2206/tutorials/combustion/reactingFoam/laminar/counterFlowFlame2D

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
