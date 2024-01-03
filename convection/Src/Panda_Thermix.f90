!> @brief The Thermal-hydraulics Calculation of Helium in Pebble bed Reactor.
!! 
!> This Software was developed for the Thermal-hydraulics Calculation of helium in Pebble bed Reactor.
!> This module was developed for the main program.
!> @author XUDONGYU
!> @version 1.0
!> @date 
PROGRAM MAIN
    USE DriverTHInput        !< A module that processes input information.
    USE DriverTHOutput       !< A module that processes output information.
    USE DriverSteadyTHSolver !< A module for steady-state solution.

    IMPLICIT NONE
    
    CALL read_TH_Input()           !< A subroutine that read input information.
    CALL driver_Steady_TH_Solver() !< A subroutine that doing  steady state thermal hydraulic solution.
    CALL print_TH_Output()         !< A subroutine that print output information.

END PROGRAM MAIN