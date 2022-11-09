subroutine PrintState() &
   bind(C, name="TCAPI_printState")

   USE ModuleThermo
   USE ModuleThermoIO
   USE ModuleGEMSolver

   implicit none

   integer :: i

   print *, "Temperature: ", dTemperature, " Pressure: ", dPressure

   print *, "Element Masses:"
   do i = 1, nElements
      print *, iElementSystem(i), " ", dMolesElement(i)
   end do

end subroutine PrintState
