forAll(fluidRegions, i)
{
    rhoThermo& thermo = thermoFluid[i];

    fv::options& fvOptions = fluidFvOptions[i];
    volScalarField& he = thermo.he();
    fvOptions.correct(he);
    thermo.correct();

}
