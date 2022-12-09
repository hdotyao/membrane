if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)
{
    #include "EEqn.H"
}
else
{
    if (oCorr == 0)
    {
        #include "rhoEqn.H" //Update the density with the help of the continuity equation
    }

    #include "UEqn.H" //Solve the momentum equation; In general this U field does not satisfy the continuity equation
//Add mapping here (both Y and p) 
//(shell[0] to tubes[1])
    MapSubMesh<eqOp>
    (
        p_Map[1],   //fieldTarget,
        p[0],   //fieldSource,
        meshToMesh0Interp,
        mapOrder
    );

    MapSubMesh<eqOp>
    (
        Y_Map[1],   //fieldTarget,
        Y[0],   //fieldSource,
        meshToMesh0Interp,
        mapOrder
    );

//(tubes[1] to shell[0])
    MapSubMesh<eqOp>
    (
        p_Map[0],   //fieldTarget,
        p[1],   //fieldSource,
        meshToMesh0Interp,
        mapOrder
    );

    MapSubMesh<eqOp>
    (
        Y_Map[0],   //fieldTarget,
        Y[1],   //fieldSource,
        meshToMesh0Interp,
        mapOrder
    );

    #include "YEqn.H" //Solve the spices transport equation -> Here besides of the concentrations of the spices at the current time step the heat source due to chemical reaction is computed. The heat source therm is required in the energy equation. 
    #include "EEqn.H" //Solve the energy equation -> Here the temperature at the new time step is computed

    if (!coupled)
    {
        Info<< "\nSolving for fluid region " << fluidRegions[i].name() << endl;

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {
            #include "pEqn.H"
        }

        turbulence.correct();

        rho = thermo.rho();
    }
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
