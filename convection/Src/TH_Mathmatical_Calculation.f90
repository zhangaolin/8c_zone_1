!> @brief A module that defines common mathematical functions
!!
!> This module is used to define the constants used to define variables and constants and functions required for mathematical calculations.
!> To avoid conflicts with other constant definitions and mathematical function definitions, the module's constant definitions should include the TH identifier whenever possible.
!!
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE THMathematicalCalculation
    USE GlobalTHConstants
    IMPLICIT NONE
    PUBLIC

CONTAINS
    !> @brief A subroutine that implements linear interpolation
    !! 
    !> @date 
    SUBROUTINE linear_Interpolation(y,x,x0,x1,y0,y1)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: y     !< The Calculation result of linear interpolation
        REAL(TH_KDUBLE), INTENT(IN)  :: x !< The substitution variable of linear interpolation
        REAL(TH_KDUBLE), INTENT(IN)  :: x0,y0 !< The first variable of linear interpolation
        REAL(TH_KDUBLE), INTENT(IN)  :: x1,y1 !< The second variable of linear interpolation
        y = y0+(x-x0)/(x1-x0)*(y1-y0)
    END SUBROUTINE linear_Interpolation         
    !> @brief A subroutine that Sums the types around the grid
    !! 
    !> @date 
    FUNCTION Summation_Surrounding_Mesh(i,j,k,data)  RESULT(result_sum)
        INTEGER,INTENT(IN) :: i !< The axial index of this mesh
        INTEGER,INTENT(IN) :: j !< The radial index of this mesh
        INTEGER,INTENT(IN) :: k !< The azimuthal index of this mesh
        INTEGER,INTENT(IN) :: data(:,:,:) !< The type of this mesh
        REAL(TH_KDUBLE) :: result_sum !< The Calculation result of Summation_Surrounding_Mesh
        result_sum = 0
        result_sum = result_sum + data(k,j,i) + data(k+1,j,i) + data(k,j+1,i) + data(k,j,i+1)
        result_sum = result_sum + data(k+1,j+1,i) + data(k+1,j,i+1) + data(k,j+1,i+1)
        result_sum = result_sum + data(k+1,j+1,i+1)
    END FUNCTION Summation_Surrounding_Mesh
END MODULE THMathematicalCalculation
