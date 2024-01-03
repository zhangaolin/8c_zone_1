!> @brief A module for steady-state solution
!! 
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE DriverSteadyTHSolver
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation
    USE GlobalTHVariables !< Global variables definition for thermal hydraulic calculation
    USE HeliumFlowSolver  !< A module for solving helium flow field
    USE HeliumTemperatureSolver !< A module for solving helium temperature field
    USE THMathematicalCalculation !< A module that defines common mathematical functions

    IMPLICIT NONE

CONTAINS
    !> @brief A subroutine that read input information
    !! 
    SUBROUTINE driver_Steady_TH_Solver()
        IMPLICIT NONE
        INTEGER :: i_iter !< Iteration number
        INTEGER :: i,j,k,temp_id !< Loop index
        INTEGER :: theta_num,radial_num,axial_num !< Number of mesh in theta, radial and axial direction
        REAL(TH_KDUBLE) :: max_error !< Maximum error
        LOGICAL :: is_converged,is_gas_tem_converged,is_cal_xkk
        
        i_iter = 0
        max_error = 1000.0
        is_gas_tem_converged = .FALSE.
        is_converged = .FALSE. 
        is_cal_xkk = .TRUE.
        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size)-1 !< 2
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size)-1 !< 29
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size)-1 !< 50

        T_solid(theta_num,radial_num+1,:) = T_solid(theta_num,radial_num,:)
        T_solid(theta_num,:,axial_num+1) = T_solid(theta_num,:,axial_num)
        DO i = 1, axial_num
            DO j = 1, radial_num
                DO k = 1, theta_num
                    T_fluid(k,j,i) = T_solid(k,j,i)
                    IF (Summation_Surrounding_Mesh(i,j,k,TH_mesh%th_mesh_type) .EQ. 0) T_fluid(k,j,i) = 0
                END DO
            END DO
        END DO

        DO i = 1,axial_num+1
            T_fluid(theta_num,radial_num+1,i) = T_solid(theta_num,radial_num,i)
        END DO
        DO j = 1,radial_num+1
            T_fluid(theta_num,j,axial_num+1) = T_solid(theta_num,j,axial_num)
        END DO
        T_fluid(1,:,:) = T_fluid(2,:,:)
        T_fluid(theta_num+1,:,:) = T_fluid(theta_num,:,:)
        
        
        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        WRITE(OUTUNIT,*) "xudongyu:T_solid"
        DO i = 1, axial_num+1
            DO j = 1, radial_num+1
                WRITE(OUTUNIT,*) i,j,T_solid(2,j,i)
            END DO
        END DO
        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        WRITE(OUTUNIT,*) "xudongyu:T_fluid"
        DO i = 1, axial_num+1
            DO j = 1, radial_num+1
                WRITE(OUTUNIT,*) i,j,T_fluid(2,j,i)
            END DO
        END DO
        CLOSE(OUTUNIT)
        
        CALL init_Flow_Variable()
        CALL init_Temperature_Variable()
        ! 网格之间的插值 call SETZT
        ! 与网格相关的净质量流量的确定 call BEZUG
        DO WHILE(NOT(is_converged))
            i_iter = i_iter + 1
            IF ((i_iter .EQ. Th_max_avg_gas_tem) .OR. is_gas_tem_converged) THEN
                is_converged = .TRUE.
            END IF
            CALL set_Flow_Resistance()
            CALL cal_Flow_Field(is_cal_xkk)
            is_cal_xkk = .FALSE.
            ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
            ! WRITE(OUTUNIT,*) "XUDONGYU: iter ",i_iter
            ! WRITE(OUTUNIT,*) "xudongyu:the finall pressure"
            ! DO i = 1, axial_num+1
            !     DO j = 1, radial_num+1
            !         WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
            !     END DO
            ! END DO
            ! WRITE(OUTUNIT,*) "xudongyu:the finall ifb"
            ! DO i = 1, axial_num+1
            !     DO j = 1, radial_num+1
            !         WRITE(OUTUNIT,*) i,j,th_mesh%th_mesh_type(2,j,i)
            !     END DO
            ! END DO
            ! CLOSE(OUTUNIT)

            CALL cal_Temperature_Field(is_gas_tem_converged)

            ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
            ! WRITE(OUTUNIT,*) "XUDONGYU: iter ",i_iter
            ! WRITE(OUTUNIT,*) "xudongyu:the finall t_fluid"
            ! DO i = 1, axial_num+1
            !     DO j = 1, radial_num+1
            !         WRITE(OUTUNIT,*) i,j,T_fluid(2,j,i)
            !     END DO
            ! END DO
            ! CLOSE(OUTUNIT)
            
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        WRITE(OUTUNIT,*) "xudongyu:the steady pressure"
        DO i = 1, axial_num+1
            DO j = 1, radial_num+1
                WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
            END DO
        END DO
        WRITE(OUTUNIT,*) "xudongyu:the steady T_fluid"
        DO i = 1, axial_num+1
            DO j = 1, radial_num+1
                WRITE(OUTUNIT,*) i,j,T_fluid(2,j,i)
            END DO
        END DO
        CLOSE(OUTUNIT)
        
    END SUBROUTINE driver_Steady_TH_Solver

END MODULE DriverSteadyTHSolver