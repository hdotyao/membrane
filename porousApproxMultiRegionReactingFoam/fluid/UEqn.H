    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    dimensionedScalar AbyV
    (
        membraneProperties.lookup("AbyV")
    );

    dimensionedScalar Qi
    (
        membraneProperties.lookup("Qi")
    );
    
    PtrList<volVectorField>& Suk = -U * ( AbyV * W(i) * Qi*( p*thermo.W()*Y[i]/W(i) - p_Map* thermo.W()*Y_Map/W(i) ) );

    UEqn =
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
      + Suk
      // + Source Term due to porous area & Suk
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            ),
            mesh.solver(U.select(finalIter))
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
