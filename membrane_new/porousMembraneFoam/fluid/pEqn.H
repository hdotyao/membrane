//this pEqn is modified from the original chtMultiRegionFoam according to rhoPimpleFoam.
//bool closedVolume = p_rgh.needReference();
dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

rho = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

//surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::interpolate(rho)*fvc::flux(HbyA)
//      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
//  + phig
);

//MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf);//, MRF);
//Info<< "\npBefore p_DDtEqn = \n " << p_rghFluid[i] << endl;
{
    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))//comment out and add source in the RHS
       + fvc::div(phiHbyA)
       ==
        fvOptions(psi, p_rgh, rho.name()) 
        + Smk[i]
    );
//Info<< "\npAfter p_DDtEqn = \n " << p_rghFluid[i] << endl;
//Info<< "\nfvOptions term in pEqn = \n " << fvOptions(psi, p_rgh, rho.name()) << endl;        
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            p_rghDDtEqn
          - fvm::laplacian(rhorAUf, p_rgh)
        );

/*        p_rghEqn.setReference
        (
            pRefCell,
            compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
        );
*/
        p_rghEqn.solve
        (
            mesh.solver
            (
                p_rgh.select
                (
                    (
                        oCorr == nOuterCorr-1
                     && corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                )
            )
        );

/*Info<< "\nrho = \n " << rhoFluid[i] << endl;
Info<< "\nU = \n " << UFluid[i] << endl;
Info<< "\npsi = \n " << psi.internalField() << endl;
Info<< "\nHbyA = \n " << HbyA.internalField() << endl;
Info<< "\nrAU = \n " << rAU.internalField() << endl;
Info<< "\nphiHbyA = \n " << phiHbyA.internalField() << endl;
Info<< "\nrhorAUf = \n " << rhorAUf.internalField() << endl;
Info<< "\npAfter p_Eqn = \n " << p_rghFluid[i] << endl;*/
        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA //- rAU*fvc::grad(p);//Based on the code provided in "rhoMassSourceFoam/pEqn.H", the velocity is calculated using the relative velocity and the pressure gradient. Specifically, the code calculates the relative velocity using the product of the inverse of the A matrix and the gradient of the pressure, and then subtracts the product of the relative velocity and the inverse of the A matrix from the constrained HbyA vector.
            + rAU*fvc::reconstruct((/*phig +*/ p_rghEqn.flux())/rhorAUf); //the code calculates the velocity using the constrained HbyA vector and the pressure gradient. Here, HbyA is the constrained HbyA vector, rAU is the inverse of the A matrix, phig is the gravitational potential energy, p_rghEqn.flux() is the flux of the pressure equation, and rhorAUf is the product of the density, the inverse of the A matrix, and the face area.
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

}


if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}
//pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
/*if (closedVolume)
{
    if (!compressible)
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
    else
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /compressibility;
        thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
        rho = thermo.rho();
        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }
}
else
{*/
    thermo.correctRho(psi*p - psip0,  rhoMin, rhoMax);
//}

#include "rhoEqn.H"
#include "compressibleContinuityErrors.H"

rho = thermo.rho();

// Update pressure time derivative if needed
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
