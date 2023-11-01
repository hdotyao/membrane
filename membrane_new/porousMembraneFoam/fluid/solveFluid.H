
if (finalIter)
{
    mesh.data::add("finalIteration", true);
}


#include "YEqn.H" //Solve the spices transport equation -> Here besides of the concentrations of the spices at the current time step the heat source due to chemical reaction is computed. The heat source therm is required in the energy equation. 

if (frozenFlow)
{
    #include "EEqn.H"
}
else
{

    if (oCorr == 0)
    {
//Info<< "\nrhoBefore[" << i << "] =\n " << rhoFluid[i] << endl;
        #include "rhoEqn.H" //Update the density with the help of the continuity equation
//Info<< "\nrhoAfter[" << i << "] =\n " << rhoFluid[i] << endl;   
    }


//Info<< "\nUBefore[" << i << "] =\n " << UFluid[i] << endl;
    #include "UEqn.H" //Solve the momentum equation; In general this U field does not satisfy the continuity equation
//Add mapping here (both Y and p) 
//(shell[0] to tubes[1])
//Info<< "\nUAfter[" << i << "] =\n " << UFluid[i] << endl;

    #include "EEqn.H" //Solve the energy equation -> Here the temperature at the new time step is computed

    if (!coupled)
    {
        Info<< "\nSolving for fluid region " << fluidRegions[i].name() << endl;

        // --- PISO loop
//Info<< "\npBefore pEqn [" << i << "] =\n " << p_rghFluid[i] << endl;
        for (int corr=0; corr<nCorr; corr++)
        {
            #include "pEqn.H"
        }
//Info<< "\npAfter pEqn[" << i << "] =\n " << p_rghFluid[i] << endl;

        turbulence.correct();

        rho = thermo.rho();
//Info<< "\npAfter rhoCorrection[" << i << "] =\n " << p_rghFluid[i] << endl;
    }
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
