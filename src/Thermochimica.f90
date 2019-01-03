
!-------------------------------------------------------------------------------
!
! Thermochimica, v 09.26.2012
! 
! Copyright (c) 2012 UT-Battelle, LLC All rights reserved.
! Redistribution and use, with or without modification, are permitted
! provided that the following conditions are met:
! 
! - Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! - Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! - Collection of administrative costs for redistribution of the source
! code or binary form is allowed. However, collection of a royalty or
! other fee in excess of good faith amount for cost recovery for such
! redistribution is prohibited.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER, THE DOE, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
! OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
! TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
! DAMAGE.
!
!
!
!> \file    Thermochimica.f90
!> \brief   The main thermochemical solver
!> \author  M.H.A. Piro
!> \date    July 17, 2012
!> \sa      CheckThermoInput.f90
!> \sa      InitThermo.f90
!> \sa      CheckSystem.f90
!> \sa      CompThermoData.f90
!> \sa      CheckThermoData.f90
!> \sa      LevelingSolver.f90
!> \sa      PostLevelingSolver.f90
!> \sa      GEMSolver.f90
!
!
!> \mainpage    Introduction to Thermochimica
!!
!! \section     Introduction 
!!  
!!          The purpose of Thermochimica is to compute the quantities of species
!!          and phases at thermodynamic equilibrium and various other 
!!          thermodynamic quantities.  For a detailed discussion on numerical 
!!          methods employed in this work and the principles of computational 
!!          thermodynamics, refer to the following literature:  
!!
!!          - M.H.A. Piro, S. Simunovic, T.M. Besmann, B.J. Lewis, W.T. Thompson,
!!            "The Thermochemistry Library Thermochimica," Computational Materials
!!            Science, 67 (2013) 266-272.
!!          - M.H.A. Piro and S. Simunovic, "Performance Enhancing Algorithms
!!            for Computing Thermodynamic Equilibria," CALPHAD, 39 (2012) 104-110.
!!          - M.H.A. Piro, "Computation of Thermodynamic Equilibria Pertinent to
!!            Nuclear Materials in Multi-Physics Codes," PhD Thesis, Royal 
!!            Military College of Canada (2011).
!!          - M.H.A. Piro, T.M. Besmann, S. Simunovic, B.J. Lewis, W.T. Thompson,
!!            "Numerical Verification of Equilibrium Computations in Nuclear 
!!            Fuel Performance Codes," Journal of Nuclear Materials, 414 (2011) 
!!            399-407.
!!
!!          The required input variables include:
!!          - dTemperature: double real scalar, represents temperature;
!!          - dPressure: double real scalar, represents absolute hydrostatic 
!!             pressure; 
!!          - dElementMass: double real vector [0:118], represents mass of each 
!!             chemical element;
!!          - cThermoInputUnits: character vector [1:3], represents the units of
!!             the input variables;
!!          - cThermoFileName: character scalar, represents data-file path and 
!!             name.
!!
!!          Useful information that Thermochimica provides as output include:
!!          - the identification of phases predicted to be stable at equilibrium
!!             (i.e., iAssemblage);
!!          - the number of moles (i.e., dMolesPhase) and the mole fraction 
!!             (i.e., dMolFraction) of each solution phase constituent;
!!          - the number of moles of each pure condensed phase 
!!             (i.e., dMolesPhase);
!!          - the chemical potentials of all species and phases 
!!             (i.e., dChemicalPotential); and  
!!          - the chemical potentials of the component elements 
!!             (i.e., dElementPotential). 
!!
!!          Thermochimica returns the variable INFOThermo (scalar integer) 
!!          indictating a successful computation or an error that has been 
!!          identified.  A description of each possible value of INFOThermo is 
!!          given in the ThermoDebug.f90 subroutine.
!!
!!          Any suggestions to improve documentation and/or programming are 
!!          always welcome and appreciated. Please send all suggestions to 
!!          Markus Piro at piromh@ornl.gov or markuspiro@gmail.com. 
!!
!!
!!
!! \section     Overview
!!
!!          Thermochimica is written in Fortran using the F90/F03 standard.  
!!          There are two calls required to make use of this library: 1) a call 
!!          to the data-file parser (i.e., ParseCSDataFile.f90), and 2) a call 
!!          to Thermochimica (i.e., Thermochimica.f90).  See thermo.f90 for an 
!!          example wrapper executable that calls the parser and Thermochimica 
!!          with input variables described above.  Finally, both the parser and 
!!          Thermochimica are reset by calling ResetThermoAll.f90.  
!!
!!          The main subroutines called by ParseCSDataFile.f90 are as follows:
!!
!!
!! <table border="1" width="800">
!! <tr>
!!    <td> <b> File name </td> <td> Description </b> </td>
!! </tr>
!! <tr>
!!    <td> ParseCSDataFile.f90 </td> 
!!    <td> Open the specified ChemSage data-file and return an error if 
!!         necessary.  </td>
!! </tr>
!! <tr>
!!    <td> ParseCSHeader.f90 </td> 
!!    <td> Parse the header section of a ChemSage data-file.  </td>
!! </tr>
!! <tr>
!!    <td> ParseCSDataBlock.f90 </td> 
!!    <td> Parse the data-block section of a ChemSage data-file.  </td>
!! </tr>
!! <tr>
!!    <td> ParseCSDataBlockGibbs.f90 </td> 
!!    <td> Parse the coefficients of the Gibbs energy equations in the 
!!          datablock section of a ChemSage data-file.  </td>
!! </tr>
!! <tr>
!!    <td> ModuleParseCS.f90 </td> 
!!    <td> Module containing variables used by the parser.  </td>
!! </tr>
!! </table>
!!
!!
!! The main subroutines called by Thermochimica.f90 are as follows:
!!
!! <table border="1" width="800">
!! <tr>
!!    <td> <b> File name </td> <td> Description </b> </td>
!! </tr>
!! <tr>
!!    <td> CheckThermoInput.f90 </td> 
!!    <td> Check the input variables to Thermochimica and return an error if 
!!          inappropriate.  </td>
!! </tr>
!! <tr>
!!    <td> InitThermo.f90 </td> 
!!    <td> Initialize Thermochimica (e.g., physical constants and numerical 
!!          tolerances).  </td>
!! </tr>
!! <tr>
!!    <td> CheckSystem.f90 </td> 
!!    <td> Check for consistency between the data-file and the current system. 
!!         </td>
!! </tr>
!! <tr>
!!    <td> CompThermoData.f90 </td> 
!!    <td> Compute thermodynamic data (e.g., standard Gibbs energies, etc.).    
!!          </td>
!! </tr>
!! <tr>
!!    <td> CheckThermoData.f90 </td> 
!!    <td> Check the thermodynamic data to ensure it is appropriate.  </td>
!! </tr>
!! <tr>
!!    <td> LevelingSolver.f90 </td> 
!!    <td> The Leveling Solver estimates the equilibrium phase assemblage 
!!          assuming that all solution phase constituents and pure condensed 
!!          phase may be initially treated as pure condensed phasess.  </td>
!! </tr>
!! <tr>
!!    <td> PostLevelingSolver.f90 </td> 
!!    <td> The Post-Leveling Solver further improves upon the estimates from 
!!          Leveling by including compositional dependent terms to solution 
!!          phase constituents.  </td>
!! </tr>
!! <tr>
!!    <td> GEMSolver.f90 </td> 
!!    <td> The Gibbs Energy Minimization (GEM) Solver computes thermodynamic 
!!          equilibrium including non-ideal mixing behaviour.  </td>
!! </tr>
!! </table>
!! 
!! 
!> \section     Style
!! 
!! All of the associated subroutines employ the following variable naming 
!! convention:
!!
!! <table border="1" width="800" align="center|left">
!! <tr>
!!    <td> <b> Prefix </td> <td> Description </b> </td>
!! </tr>
!! <tr>
!!    <td> i </td>
!!    <td> Integer variable </td>
!! </tr>   
!! <tr>
!!    <td> n </td>
!!    <td> Number of something (e.g., nElements refers to the number of 
!!          elements in the system) </td>
!! </tr>   
!! <tr>
!!    <td> d </td>
!!    <td> Double real variable (i.e., real(8)) </td>
!! </tr>   
!! <tr>
!!    <td> c </td>
!!    <td> Character variable </td>
!! </tr>   
!! <tr>
!!    <td> l </td>
!!    <td> Logical variable </td>
!! </tr>   
!! </table>
!!              
!!
!!
!> \section     Glossary
!!
!! The following gives a brief description of thermodynamics nomenclature used 
!!  in the soure code.      
!! <table border="1" width="800">
!! <tr>
!!    <td> <b> Term </td> <td> Description </b> </td>
!! </tr>
!! <tr>
!!    <td> Activity </td> 
!!    <td> A dimensionless quantity related to the chemical potential of a 
!!          substance and is represented by \f$ a_i \f$.  The activity is 
!!          equivalent to mole fraction for an ideal solution phase.  </td>
!! </tr>
!! <tr>
!!    <td> Activity Coefficient </td>
!!    <td> The activity coefficient accounts for the departure of a substance 
!!          from ideal behaviour and is represented by \f$ \gamma_i \f$. This 
!!          is related to the Partial molar excess Gibbs energy of mixing. </td>
!! </tr>
!! <tr>
!!    <td> Aqueous phase </td>
!!    <td> A particular solution phase where the solvent is water and many of 
!!          the solutes are electrically charged ions. </td>
!! </tr>   
!! <tr>
!!    <td> Chemical potential </td>
!!    <td> A measure of the effect on the Gibbs energy of the system by the 
!!          introduction of a substance. The chemical potential is defined as 
!!          \f$ \mu_i = \frac{\partial G_{sys}}{\partial n_i} 
!!          | _{T,P,n_{k \neq i}} \f$, which yields 
!!          \f$ \mu_i = g_i^{\circ} + RTln(a_i) \f$.   </td>
!! </tr>   
!! <tr>
!!    <td> Closed system </td>
!!    <td> A system that permits the exchange of heat and work with its 
!!          surroundings at constant mass. </td>
!! </tr>   
!! <tr>
!!    <td> Constituent </td>
!!    <td> A constituent of a solution phase refers to a particular species in 
!!          a particular phase. </td>
!! </tr>
!! <tr>
!!    <td> Element </td>
!!    <td> A chemical element, which is not to be confused with a nuclear fuel 
!!          element or an element of a vector/matrix. </td>
!! </tr>   
!! <tr>
!!    <td> Gibbs energy </td>
!!    <td> A thermodynamic potential measuring the maximum amount of useful work 
!!          obtainable from an isothermal-isobaric closed system.  The Gibbs 
!!          energy, represented as G, is defined as the difference between 
!!          enthalpy and the product of temperature and entropy 
!!          \f$ G = H - TS \f$. </td>
!! </tr>   
!! <tr>
!!    <td> Ion </td>
!!    <td> An atom or molecule where the number of electrons does not equal the 
!!          number of protons. </td>
!! </tr>   
!! <tr>
!!    <td> Isobaric </td>
!!    <td> A system at constant pressure. </td>
!! </tr>   
!! <tr>
!!    <td> Isothermal </td>
!!    <td> A system at constant temperature. </td>
!! </tr>   
!! <tr>
!!    <td> Molality </td>
!!    <td> Molality denotes the number of moles of a solute i per kilogram of 
!!          solvent (not solution). </td>
!! </tr> 
!! <tr>
!!    <td> Mole </td>
!!    <td> A quantity of mass measured as 6.02214179E23 atoms.  Equivalent to 
!!          gram-atom for a pure element. </td>
!! </tr>   
!! <tr>
!!    <td> Mole fraction </td>
!!    <td> The fraction of moles of a particular species in a particular 
!!          solution phase. </td>
!! </tr>   
!!
!! </table>
!!
!!
!! <table border="1" width="800">
!! <tr>
!!    <td> <b> Term </td> <td> Description </b> </td>
!! </tr>
!! <tr>
!!    <td> Partial molar excess Gibbs energy of mixing </td>
!!    <td> The partial molar excess Gibbs energy of mixing, represented as 
!!          \f$ g_i^{ex} \f$, represents the contribution to the chemical 
!!          potential term due to non-ideal behaviour.  This is the difference 
!!          between the real chemical potential of a substance and that if 
!!          assuming ideal mixing behaviour.  </td>
!! </tr>  
!! <tr>
!!    <td> Phase </td>
!!    <td> A body of matter that is uniform in chemical composition and physical 
!!          state. Phases are separated from one another by a physical 
!!          discontinuity.  A phase is not to be confused with a state of 
!!          matter.  For example, there are three different phases of uranium 
!!          in a solid state (orthogonal, tetragonal and body centred cubic). 
!!          </td>
!! </tr>   
!! <tr>
!!    <td> Phase assemblage </td>
!!    <td> A unique combination of phases predicted to be stable at 
!!          equilibrium. </td>
!! </tr> 
!! <tr>
!!    <td> Pure condensed phase </td>
!!    <td> A condensed phase with invariant stoichiometry and may be interpreted
!!          mathematically as containing a single species with unit 
!!          concentration. </td>
!! </tr> 
!! <tr>
!!    <td> Solution phase </td>
!!    <td> A phase containing a mixture of multiple species. A solution phase 
!!          can be in a gaseous, liquid or solid state.</td>
!! </tr> 
!! <tr>
!!    <td> Species </td>
!!    <td> A chemically distinct molecular entity.  For example, H2O has a 
!!          distinct chemical composition, but can be in gaseous, liquid or 
!!          solid phases.  This differs from the term constituent, which refers 
!!          to a particular species in a particular phase.  </td>
!! </tr> 
!! <tr>
!!    <td> Standard molar Gibbs energy </td>
!!    <td> The standard molar Gibbs energy of a pure species, represented as 
!!          \f$ g_i^{\circ} \f$, is the Gibbs energy of that species with unit 
!!          activity.  This quantity is computed using values from the ChemSage 
!!          data-file.</td>
!! </tr>
!!
!! <tr> 
!!    <td> State </td>
!!    <td> A state of matter distinguishes the distinct form that matter may 
!!          take on. This includes solid, liquid, gas and plasma...and for you 
!!          physicists, the Bose-Einstein condensate.  This is not to be 
!!          confused with the term "phase".  </td>
!! </tr>
!! </tr> 
!! <tr> 
!!    <td> Stoichiometry </td>
!!    <td> This refers to the relative amounts of atoms per formula mass of a 
!!          substance. </td>
!! </tr>
!! </tr> 
!! <tr> 
!!    <td> System </td>
!!    <td> A portion of the Universe with a perimeter defined by real or 
!!          imaginary boundaries.  </td>
!! </tr>
!!
!! </table>
!!
!!
!
!
! Revisions:
! ==========
! 
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   31/03/2011      M.H.A. Piro         Original code
!   21/10/2011      M.H.A. Piro         Clean up code: update modules, add 
!                                        error checking capability
!   04/26/2012      M.H.A. Piro         Add dOxygen documentation.
!
!
! Pertinent variables:
! ====================
!
! INFOThermo        An integer scalar identifying whether the program exits 
!                    successfully or if it encounters an error.  A description 
!                    for each error is given in ThermoDebug.f90.
!
!-------------------------------------------------------------------------------


subroutine Thermochimica
    
    USE ModuleThermoIO
    USE ModuleThermo

    implicit none


    ! Check the input variables:
    if (INFOThermo == 0) call CheckThermoInput

    ! Initialize Thermochimica:
    if (INFOThermo == 0) call InitThermo 

    ! Check that the system in the data-file is consistent with the input data 
    ! variables:
    if (INFOThermo == 0) call CheckSystem

    ! Compute thermodynamic data using the parameters from the specified 
    ! ChemSage data-file:
    if (INFOThermo == 0) call CompThermoData

    ! Check the thermodynamic database to ensure that it is appropriate:
    if (INFOThermo == 0) call CheckThermoData

    ! Estimate the equilibrium phase assemblage and other important properties 
    ! using the Leveling algorithm:
    if (INFOThermo == 0) call LevelingSolver  

    ! Improve estimates from the Leveling subroutine using the Post-Leveling 
    ! algorithm:
    !if (INFOThermo == 0) call PostLevelingSolver

    ! Compute the quantities of species and phases at thermodynamic equilibrium 
    ! using the GEM method:
    if (INFOThermo == 0) call GEMSolver

    ! Perform post-processing calculations of results:
    if (INFOThermo == 0) call PostProcess
    
    return
            
end subroutine Thermochimica
