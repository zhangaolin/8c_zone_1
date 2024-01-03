!> @brief A module for solving helium temperature field
!! 
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE HeliumTemperatureSolver
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation
    USE GlobalTHVariables !< Global variables definition for thermal hydraulic calculation
    USE HeliumPhysicalPropertiesCalculation !< The module for calculating helium physical properties
    USE HeliumCharacteristicNumberCalculation !< The module for calculating helium characteristic number
    IMPLICIT NONE

    REAL(TH_KDUBLE),ALLOCATABLE :: temp_t_fluid(:,:,:) !< The fluid temperature of the current iteration
    REAL(TH_KDUBLE),ALLOCATABLE :: lambda_turbulent(:,:,:) !< The transverse turbulent thermal conductivity
    REAL(TH_KDUBLE),ALLOCATABLE :: alga(:,:,:,:) !< An exponential term used to calculate temperature
    REAL(TH_KDUBLE),ALLOCATABLE :: temp_real_temperature_cal(:,:,:) !< The temp real variable in solving helium temperature field
    REAL(TH_KDUBLE) :: diameter_core !< The diameter of core
    REAL(TH_KDUBLE),ALLOCATABLE :: heat_flux_circumferential(:,:,:) !< The heat flux in circumferential direction
    REAL(TH_KDUBLE),ALLOCATABLE :: heat_flux_radial(:,:,:) !< The heat flux in radial direction
    REAL(TH_KDUBLE),ALLOCATABLE :: heat_flux_axial(:,:,:) !< The heat flux in axial direction
    REAL(TH_KDUBLE),ALLOCATABLE :: field_source(:, :, :) !< The source term of the temperature field

CONTAINS
    !> @brief A subroutine that initialize the flow variable
    !! 
    !> @date 
    SUBROUTINE init_Temperature_Variable()
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        INTEGER :: i,j,k,n !< The index of mesh
        INTEGER :: temp_id !< The material id of mesh
        REAL(TH_KDUBLE) :: temp_real

        temp_real = 0.0
        ! ALLOCATE(diameter_core(Th_mesh%radial_mesh_num))
        diameter_core = 0.0
        DO j = 1, Th_mesh%radial_mesh_num
            temp_real = temp_real + Th_mesh%radial_mesh_size(j)
            DO i = 1, Th_mesh%axial_mesh_num
                DO k = 1, Th_mesh%circumferential_mesh_num
                    temp_id = Th_mesh%th_mesh_material_id(k,j,i)
                    temp_id = Channels(temp_id)%channel_type
                    IF (temp_id .EQ. 1) THEN
                        diameter_core = 2*temp_real
                    ELSE
                        CYCLE
                    END IF
                END DO
            END DO
        END DO

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        ALLOCATE(temp_t_fluid(theta_num,radial_num,axial_num))
        temp_t_fluid = 0.0
        ALLOCATE(temp_real_temperature_cal(theta_num,radial_num,axial_num))
        temp_real_temperature_cal = 0.0
        ALLOCATE(lambda_turbulent(theta_num,radial_num,axial_num))
        lambda_turbulent = 0.0
        ALLOCATE(alga(theta_num,radial_num,axial_num,2)) 
        alga = 0.0

        ALLOCATE(heat_flux_circumferential(theta_num,radial_num,axial_num))
        ALLOCATE(heat_flux_radial(theta_num,radial_num,axial_num))
        ALLOCATE(heat_flux_axial(theta_num,radial_num,axial_num))
        heat_flux_circumferential = 0.0
        heat_flux_radial = 0.0
        heat_flux_axial = 0.0

        ALLOCATE(field_source(theta_num,radial_num,axial_num))
        field_source = 0.0

    END SUBROUTINE init_Temperature_Variable
    !> @brief A subroutine that calculate temperature fields
    !! 
    !> @note This subroutine has the same functionality as the <b>GASTEM</b> 
    !> @date 
    SUBROUTINE cal_Temperature_Field(is_gas_tem_converged)
        !开始迭代
        !计算温度计算常数 
        LOGICAL, INTENT(INOUT) :: is_gas_tem_converged !< The flag of the convergence of the temperature and pressure field
        INTEGER :: i_iter
        INTEGER :: i,j,k,n,temp_id !< The index of mesh
        INTEGER :: ifque
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        REAL(TH_KDUBLE) :: temp_real,fluid_temperature_ratio,max_fluid_temperature_ratio
        LOGICAL :: is_converged

        i_iter = 0
        is_converged = .FALSE.
        fluid_temperature_ratio = 0.0
        max_fluid_temperature_ratio = 0.0
        ifque = 0
        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51
        
        CALL cal_Temperature_Constants()

        !> Store the fluid temperature before calculation
        !> Determine the fine-mesh location of the fluid mass source and add the corresponding th_mesh_type to +20 from the original value
        DO i = 1, axial_num-1
            DO j= 1, radial_num-1
                DO k = 1, theta_num-1
                    temp_t_fluid(k,j,i) = T_fluid(k,j,i)
                    temp_id = Th_mesh%th_mesh_type(k,j,i)
                    IF (temp_id .EQ. 0.0) THEN
                        CYCLE
                    ELSE
                        temp_id = th_mesh%th_fine_mesh_material_id(k,j,i)
                        temp_real = Channels(temp_id)%source_mass_flow
                        IF (temp_real .GT. 0.0) THEN
                            th_mesh%th_mesh_type(k,j,i) = Th_mesh%th_mesh_type(k,j,i) + 20
                        ELSE
                            CYCLE
                        END IF
                    END IF
                END DO
            END DO
        END DO
        field_source = 0.0
        !> Start the iteration of the temperature field calculation
        DO WHILE(NOT(is_converged))
            i_iter = i_iter + 1
            fluid_temperature_ratio = 0.0
            max_fluid_temperature_ratio = 0.0
            DO i = 1, axial_num-1
                DO j = 1, radial_num-1
                    DO k = 2, theta_num-1
                        temp_id = Summation_Surrounding_Mesh(i,j,k,TH_mesh%th_mesh_type)
                        IF (temp_id .EQ. 0.0) THEN
                            CYCLE
                        ELSE
                            IF (i .EQ. 40 .AND. j .EQ. 13) THEN
                                WRITE(*,*)
                            END IF
                            temp_real = T_fluid(k,j,i)
                            CALL cal_Temperature_Mesh(i,j,k,ifque,temp_real)
                            IF (temp_real .NE. 0.0) THEN
                                fluid_temperature_ratio = ABS(T_fluid(k,j,i)/temp_real-1.0)
                            END IF
                            T_fluid(k,j,i) = temp_real
                            IF (fluid_temperature_ratio .GE. max_fluid_temperature_ratio) THEN
                                max_fluid_temperature_ratio = fluid_temperature_ratio
                            END IF
                        END IF
                    END DO
                END DO
            END DO
            CALL is_Temperature_Field_Coverges(is_converged,i_iter,max_fluid_temperature_ratio,ifque)
            ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
            ! WRITE(OUTUNIT,*) "xudongyu: i_iter ",i_iter
            ! WRITE(OUTUNIT,*) "xudongyu: PREL ",max_fluid_temperature_ratio
            ! WRITE(OUTUNIT,*) "xudongyu: T_fluid "
            ! DO i = 1, axial_num
            !     DO j = 1, radial_num
            !         WRITE(OUTUNIT,*) i,j,T_fluid(2,j,i)
            !     END DO
            ! END DO
            ! CLOSE(OUTUNIT)
        END DO
        T_fluid(1,:,:) = T_fluid(2,:,:)
        T_fluid(theta_num,:,:) = T_fluid(theta_num-1,:,:)

        n = 0
        fluid_temperature_ratio = 0.0
        max_fluid_temperature_ratio = 0.0
        DO i = 1,axial_num-1
            DO j = 1,radial_num-1
                DO k = 2, theta_num-1
                    temp_id = Summation_Surrounding_Mesh(i,j,k,TH_mesh%th_mesh_type)
                    IF (temp_id .EQ. 0) THEN
                        CYCLE
                    ELSE
                        temp_id = th_mesh%th_mesh_type(k,j,i)
                        IF (temp_id .GT. 10) THEN
                            Th_mesh%th_mesh_type(k,j,i) = Th_mesh%th_mesh_type(k,j,i) - 20
                        END IF
                        IF (T_fluid(k,j,i) .NE. 0.0) THEN
                            n = n + 1
                            fluid_temperature_ratio = fluid_temperature_ratio+ABS(temp_t_fluid(k,j,i)/T_fluid(k,j,i)-1.0)
                        END IF
                    END IF
                END DO
            END DO
        END DO
        fluid_temperature_ratio = fluid_temperature_ratio/n
        IF (fluid_temperature_ratio .LT. Th_epsi_avg_gas_tem) THEN
            is_gas_tem_converged = .TRUE.
        ELSE
            is_gas_tem_converged = .FALSE.
        END IF
        ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: final temp_t_fluid "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,temp_t_fluid(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: final T_fluid "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,T_fluid(2,j,i)
        !     END DO
        ! END DO
        ! CLOSE(OUTUNIT)
    END SUBROUTINE cal_Temperature_Field
    !> @brief A subroutine that calculate the constants variable used in the temperature fields calcualtion
    !! 
    !> @note This subroutine has the same functionality as the <b>GASKON</b> 
    !> @date 
    SUBROUTINE cal_Temperature_Constants()
        ! 计算温度计算所需要的常数
        INTEGER :: i,j,k,n
        INTEGER :: temp_id
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        REAL(TH_KDUBLE) :: temp_real

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51
        WRITE(*,*) "ALPHAK",Radial_mass_flow(2,27,41)
        CALL cal_Alphak()
        WRITE(*,*) "XLTURB",Radial_mass_flow(2,27,41)
        CALL cal_Transverse_Turbulent_Conductivity()
        WRITE(*,*) "ALGAS",Radial_mass_flow(2,27,41)
        CALL cal_Algas()
        WRITE(*,*) "GASKON",Radial_mass_flow(2,27,41)
        DO i = 1, axial_num
            Do j = 1, radial_num
                DO k = 1, theta_num
                    temp_id = Th_mesh%th_mesh_type(k,j,i)
                    IF (temp_id .EQ. 0) THEN
                        heat_flux_circumferential(k,j,i) = 0.0
                        heat_flux_radial(k,j,i) = 0.0
                        heat_flux_axial(k,j,i) = 0.0
                    ELSE 
                        heat_flux_circumferential(k,j,i) = HALF*Circumferential_mass_flow(k,j,i)*Gas_heat_capacity
                        heat_flux_radial(k,j,i) = HALF*Radial_mass_flow(k,j,i)*Gas_heat_capacity
                        heat_flux_axial(k,j,i) = HALF*Axial_mass_flow(k,j,i)*Gas_heat_capacity
                    END IF
                END DO
            END DO
        END DO

        DO i = 1, axial_num - 1
            DO j = 1, radial_num - 1
                DO k = 1, theta_num - 1
                    temp_id = Summation_Surrounding_Mesh(i,j,k,TH_mesh%th_mesh_type)
                    IF (temp_id .EQ. 0) THEN
                        Circumferential_mass_flow(k,j,i) = 0.0
                        Radial_mass_flow(k,j,i) = 0.0
                        Axial_mass_flow(k,j,i) = 0.0
                        temp_real_temperature_cal(k,j,i) = 0.0
                    ELSE 
                        Circumferential_mass_flow(k,j,i) = 0.0
                        Radial_mass_flow(k,j,i) = 0.0
                        Axial_mass_flow(k,j,i) = 0.0
                        IF (lambda_turbulent(k,j,i) .GT. 0.0 .OR. lambda_turbulent(k,j,i+1) .GT. 0.0) THEN
                            temp_real = 2.0*th_mesh%radial_fine_mesh_size(j)/(Th_mesh%th_mesh_area(k,j,i,5)*(lambda_turbulent(k,j,i)+lambda_turbulent(k,j,i+1)))
                            Radial_mass_flow(k,j,i) = 1.0/temp_real
                        END IF
                        IF (lambda_turbulent(k,j,i) .GT. 0.0 .OR. lambda_turbulent(k,j+1,i) .GT. 0.0) THEN
                            temp_real = 2.0*th_mesh%axial_fine_mesh_size(i)/(Th_mesh%th_mesh_area(k,j,i,6)*(lambda_turbulent(k,j,i)+lambda_turbulent(k,j+1,i)))
                            Axial_mass_flow(k,j,i) = 1.0/temp_real
                        END IF
                        IF (i .EQ. 1 .OR. i .EQ. (axial_num-1)) THEN
                            Radial_mass_flow(k,j,i) = 2.0*Radial_mass_flow(k,j,i)
                        END IF
                        IF (j .EQ. 1 .OR. j .EQ. (radial_num-1)) THEN
                            Axial_mass_flow(k,j,i) = 2.0*Axial_mass_flow(k,j,i)
                        END IF
                    END IF
                END DO
            END DO
        END DO

        ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: heat_flux_radial "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,heat_flux_radial(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: heat_flux_axial "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,heat_flux_axial(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: rho "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,temp_real_temperature_cal(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: Radial_mass_flow "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,Radial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: Axial_mass_flow "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,Axial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        ! CLOSE(OUTUNIT)

    END SUBROUTINE cal_Temperature_Constants
    !> @brief A subroutine that calculate alpha
    !! 
    !> @note This subroutine has the same functionality as the <b>ALPHAK</b> 
    !> @date 2023-11-08
    SUBROUTINE cal_Alphak()
        ! 遍历所有节块，判断流道类型计算努塞尔数，然后计算指数尽速所需要的临时参数
        !DROHR,ALPHA-Coefficient_heat_transition;QTV-Qtv;LTV-Ltv;XGEO-Xgeo;T-T_solid;TFL-T_fluid;IFB,FZQ,FRQ,DZ,EPSIL-THMesh;STROM-temp_real_temperature_cal;DHYD-Hydraulic_diameter;RHO-Helium_density;MZ-Axial_mass_flow;MR-Radial_mass_flow;XVER-是周向长度的累计，需要和流道以及尺寸编号对应-sum_axial_mesh;ALFISO-Coefficient_heat_transition;
        INTEGER :: i,j,k,n
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        INTEGER :: temp_id 
        REAL(TH_KDUBLE) :: temp_real !< The temp real variable in solving helium flow field
        REAL(TH_KDUBLE) :: thermal_conductivity,relo,nusselt !< The thermal conductivity, reynolds number and nusselt number
        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        
        temp_real_temperature_cal = 0.0

        DO i = 2, axial_num-1
            DO j = 2, radial_num-1
                DO k = 2, theta_num-1
                    temp_id = Th_mesh%th_mesh_type(k,j,i)
                    IF (temp_id .EQ. 0) THEN
                    ELSE IF (temp_id .EQ. 1) THEN
                        temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                        temp_real = T_fluid(k,j,i) + T_fluid(k-1,j,i) + T_fluid(k,j-1,i) + T_fluid(k,j,i-1)
                        temp_real = temp_real + T_fluid(k-1,j-1,i) + T_fluid(k-1,j,i-1) + T_fluid(k,j-1,i-1)
                        temp_real = temp_real + T_fluid(k-1,j-1,i-1)
                        temp_real = temp_real + T_solid(k,j,i)
                        temp_real = temp_real + T_solid(k-1,j,i) + T_solid(k,j-1,i) + T_solid(k,j,i-1)
                        temp_real = temp_real + T_solid(k-1,j-1,i) + T_solid(k-1,j,i-1) + T_solid(k,j-1,i-1)
                        temp_real = temp_real + T_solid(k-1,j-1,i-1)
                        temp_real = temp_real / 16
                        ! WRITE(*,*) k,j,i
                        ! WRITE(*,*) temp_real
                        CALL cal_Helium_Conductivity(thermal_conductivity,Gas_pressure,temp_real)
                        ! WRITE(*,*) k,j,i,temp_real,thermal_conductivity
                        CALL cal_Helium_Re(relo,temp_real,Helium_density(k,j,i),Channels(temp_id)%channel_type,Th_mesh%th_mesh_porosity(k,j,i), &
                                       &  Compositions(temp_id)%diameter_fuel_element,Channels(temp_id)%hydraulic_diameter,Th_mesh%th_mesh_area(k,j,i,3), &
                                       &  Th_mesh%th_mesh_area(k,j,i,2),Th_mesh%th_mesh_area(k,j,i,1),Axial_mass_flow(k,j,i),Radial_mass_flow(k,j,i),Circumferential_mass_flow(k,j,i))
                        ! WRITE(*,*) thermal_conductivity
                        relo = relo*(1-th_mesh%th_mesh_porosity(k,j,i))*1.5
                        CALL cal_Helium_Nu_Pebble_Bed(nusselt,Prandtl_number,relo,th_mesh%th_mesh_porosity(k,j,i))
                        ! WRITE(*,*) nusselt
                        nusselt = nusselt*thermal_conductivity/Compositions(temp_id)%diameter_fuel_element
                        nusselt = nusselt*Th_mesh%th_mesh_area(k,j,i,3)*Th_mesh%axial_fine_mesh_size(i)*(1-th_mesh%th_mesh_porosity(k,j,i))*6.0/Compositions(temp_id)%diameter_fuel_element
                        ! WRITE(*,*) Compositions(temp_id)%diameter_fuel_element,Th_mesh%th_mesh_area(k,j,i,3),Th_mesh%axial_fine_mesh_size(i),th_mesh%th_mesh_porosity(k,j,i)
                        temp_real_temperature_cal(k,j,i) = nusselt
                    ELSE IF (temp_id .EQ. 2) THEN
                        DO n = 1, SIZE(material_channel2)
                            IF (material_channel2(n) .EQ. Th_mesh%th_fine_mesh_material_id(k,j,i)) THEN
                                temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                                temp_real = T_fluid(k,j,i) + T_fluid(k-1,j,i) + T_fluid(k,j-1,i) + T_fluid(k,j,i-1)
                                temp_real = temp_real + T_fluid(k-1,j-1,i) + T_fluid(k-1,j,i-1) + T_fluid(k,j-1,i-1)
                                temp_real = temp_real + T_fluid(k-1,j-1,i-1)
                                temp_real = temp_real + T_solid(k,j,i)
                                temp_real = temp_real + T_solid(k-1,j,i) + T_solid(k,j-1,i) + T_solid(k,j,i-1)
                                temp_real = temp_real + T_solid(k-1,j-1,i) + T_solid(k-1,j,i-1) + T_solid(k,j-1,i-1)
                                temp_real = temp_real + T_solid(k-1,j-1,i-1)
                                temp_real = temp_real / 16
                                CALL cal_Helium_Conductivity(thermal_conductivity,Gas_pressure,temp_real)
                                CALL cal_Helium_Re(relo,temp_real,Helium_density(k,j,i),Channels(temp_id)%channel_type,Th_mesh%th_mesh_porosity(k,j,i), &
                                               &  Compositions(temp_id)%diameter_fuel_element,Channels(temp_id)%hydraulic_diameter,Th_mesh%th_mesh_area(k,j,i,3), &
                                               &  Th_mesh%th_mesh_area(k,j,i,2),Th_mesh%th_mesh_area(k,j,i,1),Axial_mass_flow(k,j,i),Radial_mass_flow(k,j,i),Circumferential_mass_flow(k,j,i))
                                temp_real = Channels(temp_id)%hydraulic_diameter/th_mesh%sum_axial_mesh(k,j,n)
                                CALL cal_Helium_Nu_Vertical_Pipes(nusselt,Prandtl_number,relo,temp_real)
                                nusselt = nusselt*thermal_conductivity/Channels(temp_id)%hydraulic_diameter
                                nusselt = 4*nusselt*Th_mesh%th_mesh_area(k,j,i,3)*Th_mesh%axial_fine_mesh_size(i)*(th_mesh%th_mesh_porosity(k,j,i))/Channels(temp_id)%hydraulic_diameter
                                temp_real_temperature_cal(k,j,i) = nusselt
                            ELSE
                                CYCLE
                            END IF
                        END DO
                    ELSE IF (temp_id .EQ. 5) THEN
                        CYCLE
                    END IF
                END DO
            END DO
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: strom "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,temp_real_temperature_cal(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)

    END SUBROUTINE cal_Alphak
    !> @brief A subroutine that calculate alpha
    !! 
    !> @note This subroutine has the same functionality as the <b>ALGAS</b> 
    !> @date 
    SUBROUTINE cal_Algas()
        ! 遍历所有节块，判断流道类型计算努塞尔数，然后计算指数尽速所需要的临时参数
        !alga-alga;RHO-Helium_density;MZ-Axial_mass_flow;MR-Radial_mass_flow;
        INTEGER :: i,j,k,n
        INTEGER :: temp_id
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        REAL(TH_KDUBLE) :: temp_real !< The temp real variable in solving helium flow field
        REAL(TH_KDUBLE) :: sum_mass_flow
        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51


        DO i = 1, axial_num
            DO j = 1, radial_num
                DO k = 1, theta_num
                    IF (k .EQ. 2 .AND. i .EQ. 21 .AND. j .EQ. 7) THEN
                        ! WRITE(*,*) ""
                    END IF
                    sum_mass_flow = ABS(Axial_mass_flow(k,j,i)) + ABS(Radial_mass_flow(k,j,i)) + ABS(Circumferential_mass_flow(k,j,i))
                    ! WRITE(*,*) k,j,i
                    ! WRITE(*,*) Axial_mass_flow(k,j,i),Radial_mass_flow(k,j,i),Circumferential_mass_flow(k,j,i)
                    temp_real = temp_real_temperature_cal(k,j,i)

                    IF (sum_mass_flow .EQ. 0.0) THEN
                        temp_real = 0.0
                    ELSE 
                        temp_real = temp_real / sum_mass_flow / Gas_heat_capacity
                    END IF

                    IF (temp_real .GT. 30) THEN
                        temp_real_temperature_cal(k,j,i) = temp_real
                    ELSE
                        temp_real_temperature_cal(k,j,i) = temp_real
                        alga(k,j,i,1) = EXP(-temp_real)
                        temp_real = -HALF*temp_real
                        alga(k,j,i,2) = EXP(temp_real)
                    END IF
                END DO
            END DO
        END DO

        ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: alga 1 "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,alga(2,j,i,1)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: alga 2 "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,alga(2,j,i,2)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: rho "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,temp_real_temperature_cal(2,j,i)
        !     END DO
        ! END DO
        ! CLOSE(OUTUNIT)
    END SUBROUTINE cal_Algas
    !> @brief A subroutine that Calculate transverse turbulent thermal conductivity
    !! 
    !> @note This subroutine has the same functionality as the <b>ALGAS</b> 
    !> @date 
    SUBROUTINE cal_Transverse_Turbulent_Conductivity()
    ! 遍历所有节块，计算横向湍流热导率
    !LAMTUR-lambda_turbulent;TFL-T_fluid;RHO-Helium_density;IFB,KOM,FZQ,FRQ-THMESH;MZ-Axial_mass_flow;MR-Radial_mass_flow;XKON-Additional_pressure_drop;
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        INTEGER :: i,j,k,n !< The index of mesh
        INTEGER :: temp_id !< The material id of mesh
        REAL(TH_KDUBLE) :: temp_real !< The temp real variable in solving helium flow field
        REAL(TH_KDUBLE) :: velocity,thermal_conductivity
        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        lambda_turbulent = 0.0
        DO i = 1, axial_num
            DO j = 1, radial_num
                DO k = 1, theta_num
                    temp_id = Th_mesh%th_mesh_type(k,j,i)
                    IF (k .EQ. 2) THEN
                        ! WRITE(*,*) ""
                    END IF
                    IF (temp_id .EQ. 0) THEN
                        CYCLE
                    ELSE IF (temp_id .EQ. 1) THEN
                        temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                        temp_real = diameter_core
                        temp_real = 8.0*(2.0-(1.0-2.0*Compositions(temp_id)%diameter_fuel_element/temp_real)**2)
                        CALL cal_Scalar_Velocity(velocity,Helium_density(k,j,i),th_mesh%th_mesh_area(k,j,i,3),th_mesh%th_mesh_area(k,j,i,2), &
                                                & th_mesh%th_mesh_area(k,j,i,1),Axial_mass_flow(k,j,i),Radial_mass_flow(k,j,i),Circumferential_mass_flow(k,j,i))
                        lambda_turbulent(k,j,i) = Helium_density(k,j,i)*Gas_heat_capacity*velocity*Compositions(temp_id)%diameter_fuel_element/temp_real
                    ELSE IF (temp_id .EQ. 2) THEN
                        CYCLE
                    ELSE IF (temp_id .EQ. 5) THEN
                        temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                        temp_real = T_fluid(k,j,i) + T_fluid(k-1,j,i) + T_fluid(k,j-1,i) + T_fluid(k,j,i-1)
                        temp_real = temp_real + T_fluid(k-1,j-1,i) + T_fluid(k-1,j,i-1) + T_fluid(k,j-1,i-1)
                        temp_real = temp_real + T_fluid(k-1,j-1,i-1)
                        temp_real = temp_real / 8
                        CALL cal_Helium_Conductivity(thermal_conductivity,Gas_pressure,temp_real)
                        lambda_turbulent(k,j,i) = thermal_conductivity*(1.0+Channels(temp_id)%additional_pressure_drop)
                    END IF
                END DO
            END DO
        END DO

        ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: lambda_turbulent "
        ! DO i = 1, axial_num
        !     DO j = 1, radial_num
        !         WRITE(OUTUNIT,*) i,j,lambda_turbulent(2,j,i)
        !     END DO
        ! END DO
        ! CLOSE(OUTUNIT)
    END SUBROUTINE cal_Transverse_Turbulent_Conductivity
    !> @brief A subroutine that calculate temperature fields
    !! 
    !> @note This subroutine has the same functionality as the <b>GASTEM</b> 
    !> @date 
    SUBROUTINE cal_Temperature_Mesh(i,j,k,ifque,mesh_t_fluid)
        ! REAL(TH_KDUBLE), INTENT(INOUT) :: temperature !< The temperature of this element 
        INTEGER, INTENT(IN) :: i,j,k !< The index of this element
        ! REAL(TH_KDUBLE), INTENT(IN) :: source_mass(6)!< The Mass volume source
        INTEGER, INTENT(IN) :: ifque !< The flag of the mass source
        REAL(TH_KDUBLE),INTENT(INOUT) :: mesh_t_fluid !< The fluid temperature of the current mesh
        INTEGER :: temp_id,n,m
        REal(TH_KDUBLE) :: temp_real
        INTEGER :: num_face,num_direction
        REAL(TH_KDUBLE) :: mesh_t_solid !< The solid temperature of the current mesh
        REAL(TH_KDUBLE),ALLOCATABLE :: surrounding_t_fluid(:) !< The surrounding fluid temperature
        REAL(TH_KDUBLE),ALLOCATABLE :: surrounding_t_solid(:) !< The surrounding solid temperature
        REAL(TH_KDUBLE),ALLOCATABLE :: thermal_conductance(:) !< The thermal conductance of the mesh
        REAL(TH_KDUBLE),ALLOCATABLE :: heat_flux(:,:) !< The heat flux of the mesh
        REAL(TH_KDUBLE),ALLOCATABLE :: centre_exponent(:),boundary_exponent(:) !< The exponent of the mesh
        REAL(TH_KDUBLE),ALLOCATABLE :: heat_transfer_efficiency(:) !< The heat transfer efficiency of the mesh
        REAL(TH_KDUBLE),ALLOCATABLE :: tef1(:,:),tef2(:,:),tef3(:,:)
        REAL(TH_KDUBLE),ALLOCATABLE :: tfh,tfe,tfa
        INTEGER,ALLOCATABLE :: if1(:)
        REAL(TH_KDUBLE) :: numerator,denominator

        num_face = 4
        num_direction = 2

        IF (NOT(ALLOCATED(surrounding_t_fluid))) THEN
            ALLOCATE(surrounding_t_fluid(num_face)) !< pb
            ALLOCATE(surrounding_t_solid(num_face)) !> pc
            ALLOCATE(if1(num_face))
        END IF

        IF(NOT(ALLOCATED(thermal_conductance))) THEN
            temp_id = 2*num_direction
            ALLOCATE(thermal_conductance(temp_id))
            temp_id = 2**num_direction
            ALLOCATE(heat_flux(temp_id,num_direction))
            ALLOCATE(centre_exponent(temp_id))
            ALLOCATE(boundary_exponent(temp_id))
            ALLOCATE(heat_transfer_efficiency(temp_id))
            ALLOCATE(tef1(temp_id,num_direction))
            ALLOCATE(tef2(temp_id,num_direction))
            ALLOCATE(tef3(temp_id,num_direction))
        END IF

        !> Left side
        IF (j .EQ. 1) THEN
            surrounding_t_fluid(1) = T_fluid(k,j,i)
            surrounding_t_solid(1) = surrounding_t_fluid(1)
        ELSE
            surrounding_t_fluid(1) = T_fluid(k,j-1,i)
            surrounding_t_solid(1) = T_solid(k,j-1,i)
        END IF
        !> Top side
        IF (i .EQ. 1) THEN
            surrounding_t_fluid(2) = T_fluid(k,j,i)
            surrounding_t_solid(2) = surrounding_t_fluid(2)
        ELSE
            surrounding_t_fluid(2) = T_fluid(k,j,i-1)
            surrounding_t_solid(2) = T_solid(k,j,i-1)
        END IF
        !> Right side
        surrounding_t_fluid(3) = T_fluid(k,j+1,i)
        surrounding_t_solid(3) = T_solid(k,j+1,i)
        !> Bottom side
        surrounding_t_fluid(4) = T_fluid(k,j,i+1)
        surrounding_t_solid(4) = T_solid(k,j,i+1)
        !> Solid temperature of the current mesh
        mesh_t_solid = T_solid(k,j,i)
        !> The thermal conductance of the mesh
        thermal_conductance(1) = Radial_mass_flow(k,j,i)
        thermal_conductance(2) = Axial_mass_flow(k,j,i)
        thermal_conductance(3) = Radial_mass_flow(k,j+1,i)
        thermal_conductance(4) = Axial_mass_flow(k,j,i+1)
        !> The heat flux of the mesh
        heat_flux(1,1) = heat_flux_radial(k,j,i)
        heat_flux(1,2) = heat_flux_axial(k,j,i)
        heat_flux(2,1) = -heat_flux_radial(k,j+1,i)
        heat_flux(2,2) = heat_flux_axial(k,j+1,i)
        heat_flux(3,1) = heat_flux_radial(k,j,i+1)
        heat_flux(3,2) = -heat_flux_axial(k,j,i+1)
        heat_flux(4,1) = -heat_flux_radial(k,j+1,i+1)
        heat_flux(4,2) = -heat_flux_axial(k,j+1,i+1)
        !> The exponent of the mesh
        centre_exponent(1) = alga(k,j,i,1)
        centre_exponent(2) = alga(k,j+1,i,1)
        centre_exponent(3) = alga(k,j,i+1,1)
        centre_exponent(4) = alga(k,j+1,i+1,1)
        boundary_exponent(1) = alga(k,j,i,2)
        boundary_exponent(2) = alga(k,j+1,i,2)
        boundary_exponent(3) = alga(k,j,i+1,2)
        boundary_exponent(4) = alga(k,j+1,i+1,2)
        !> The heat transfer efficiency of the mesh
        heat_transfer_efficiency(1) = temp_real_temperature_cal(k,j,i)
        heat_transfer_efficiency(2) = temp_real_temperature_cal(k,j+1,i)
        heat_transfer_efficiency(3) = temp_real_temperature_cal(k,j,i+1)
        heat_transfer_efficiency(4) = temp_real_temperature_cal(k,j+1,i+1)

        DO n = 1, num_face
            if1(n) = 0
            IF (centre_exponent(n) .EQ. 1.) THEN
                centre_exponent(n) = 1.0
                boundary_exponent(n) = 1.0
                if1(n) = 2
                heat_transfer_efficiency(n) = 1.0
            ELSE
                CYCLE
            END IF
        END DO

        tef1(1,1) = surrounding_t_fluid(1)-surrounding_t_solid(1)
        tef1(1,2) = surrounding_t_fluid(2)-surrounding_t_solid(2)
        tef1(2,1) = surrounding_t_fluid(3)-surrounding_t_solid(3)
        tef1(2,2) = surrounding_t_fluid(2)-surrounding_t_solid(2)
        tef1(3,1) = surrounding_t_fluid(1)-surrounding_t_solid(1)
        tef1(3,2) = surrounding_t_fluid(4)-surrounding_t_solid(4)
        tef1(4,1) = surrounding_t_fluid(3)-surrounding_t_solid(3)
        tef1(4,2) = surrounding_t_fluid(4)-surrounding_t_solid(4)
        tef2(1,1) = surrounding_t_solid(1)
        tef2(1,2) = surrounding_t_solid(2)
        tef2(2,1) = surrounding_t_solid(3)
        tef2(2,2) = surrounding_t_solid(2)
        tef2(3,1) = surrounding_t_solid(1)
        tef2(3,2) = surrounding_t_solid(4)
        tef2(4,1) = surrounding_t_solid(3)
        tef2(4,2) = surrounding_t_solid(4)
        tef3(1,1) = (mesh_t_solid-surrounding_t_solid(1))*thermal_conductance(1)
        tef3(1,2) = (mesh_t_solid-surrounding_t_solid(2))*thermal_conductance(1)
        tef3(2,1) = (mesh_t_solid-surrounding_t_solid(3))*thermal_conductance(2)
        tef3(2,2) = (mesh_t_solid-surrounding_t_solid(2))*thermal_conductance(2)
        tef3(3,1) = (mesh_t_solid-surrounding_t_solid(1))*thermal_conductance(3)
        tef3(3,2) = (mesh_t_solid-surrounding_t_solid(4))*thermal_conductance(3)
        tef3(4,1) = (mesh_t_solid-surrounding_t_solid(3))*thermal_conductance(4)
        tef3(4,2) = (mesh_t_solid-surrounding_t_solid(4))*thermal_conductance(4)

        DO n = 1,2**num_direction
            DO m = 1,num_direction
                temp_real = heat_flux(n,m)
                IF (temp_real .LE. 0.0) THEN
                    CYCLE
                ELSE
                    tfe = tef1(n,m)+tef2(n,m)
                    IF (if1(n) .EQ. 2) THEN
                        tfh = tfe
                        tfa = tfe
                        tef2(n,m) = tfe - mesh_t_solid
                    ELSE
                        IF (ifque .NE. 1) THEN
                            tef2(n,m) = tef1(n,m)*centre_exponent(n)+tef3(n,m)*(centre_exponent(n)-1.0)
                            tfa = tef2(n,m) + mesh_t_solid
                        ELSE
                            tfh = tef1(n,m)*boundary_exponent(n)-tef3(n,m)*(boundary_exponent(n)-1.0)
                            tfh = tfh + HALF*(tef2(n,m)+mesh_t_solid)
                            tef2(n,m) = tef1(n,m)*centre_exponent(n)+tef3(n,m)*(centre_exponent(n)-1.0)
                            tfa = tef2(n,m) + mesh_t_solid
                        END IF
                    END IF
                    IF (ifque .EQ. 1) THEN
                        CALL cal_Field_Source(i,j,k,n,m,tfe,tfh,tfa,heat_flux(n,m))
                    END IF
                END IF
            END DO
        END DO

        numerator = 0.0
        denominator = 0.0
        DO n = 1,2**num_direction
            tfe = surrounding_t_fluid(n)-mesh_t_solid
            numerator = numerator + tfe*thermal_conductance(n)
            denominator = denominator + thermal_conductance(n)
        END DO
        DO n = 1,2**num_direction
            DO m = 1, num_direction
                IF(heat_flux(n,m) .LE. 0.0) THEN
                    CYCLE
                ELSE
                    numerator = numerator + tef2(n,m)*heat_flux(n,m)
                    denominator = denominator + heat_flux(n,m)
                END IF
            END DO
        END DO

        temp_id = Summation_Surrounding_Mesh(i,j,k,th_mesh%th_mesh_type)
        IF (temp_id .GE. 20) THEN
            DO n = 0,1
                DO m = 0,1
                    temp_id = Th_mesh%th_fine_mesh_material_id(k,j+n,i+m)
                    IF (temp_id .NE. 0) THEN 
                        temp_real = Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j,i)*Gas_heat_capacity
                        numerator = numerator + HALF*HALF*(Channels(temp_id)%source_temperature_in-mesh_t_solid)*temp_real
                        denominator = denominator + HALF*HALF*temp_real
                    END IF
                END DO
            END DO
        END IF
        
        ! temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
        ! IF (temp_id .NE. 0) THEN 
        !     temp_real = Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j,i)*Gas_heat_capacity
        !     numerator = numerator + HALF*HALF*(Channels(temp_id)%source_temperature_in-mesh_t_solid)*temp_real
        !     denominator = denominator + HALF*HALF*temp_real
        ! END IF
        ! temp_id = Th_mesh%th_fine_mesh_material_id(k,j+1,i)
        ! IF (temp_id .NE. 0) THEN 
        !     temp_real = Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j,i)*Gas_heat_capacity
        !     numerator = numerator + HALF*HALF*(Channels(temp_id)%source_temperature_in-mesh_t_solid)*temp_real
        !     denominator = denominator + HALF*HALF*temp_real
        ! END IF
        ! temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i+1)
        ! IF (temp_id .NE. 0) THEN 
        !     temp_real = Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j,i)*Gas_heat_capacity
        !     numerator = numerator + HALF*HALF*(Channels(temp_id)%source_temperature_in-mesh_t_solid)*temp_real
        !     denominator = denominator + HALF*HALF*temp_real
        ! END IF
        ! temp_id = Th_mesh%th_fine_mesh_material_id(k,j+1,i+1)
        ! IF (temp_id .NE. 0) THEN 
        !     temp_real = Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j,i)*Gas_heat_capacity
        !     numerator = numerator + HALF*HALF*(Channels(temp_id)%source_temperature_in-mesh_t_solid)*temp_real
        !     denominator = denominator + HALF*HALF*temp_real
        ! END IF
        IF (denominator .GT. 0.0) THEN
            mesh_t_fluid = mesh_t_solid + numerator/denominator
        ELSE
            mesh_t_fluid = mesh_t_solid
        END IF

    END SUBROUTINE cal_Temperature_Mesh
    !> @brief A subroutine that calculate temperature fields
    !! 
    !> @note This subroutine has the same functionality as the <b>QUELL</b> 
    !> @date 
    SUBROUTINE cal_Field_Source(i,j,k,n,m,tfe,tfh,tfa,heat_flux)
        INTEGER,INTENT(IN) :: i,j,k,n,m
        REAL(TH_KDUBLE),INTENT(IN) :: tfe,tfh,tfa
        REAL(TH_KDUBLE),INTENT(IN) :: heat_flux
        INTEGER :: field_source_i,field_source_j,field_source_k
        
        IF (n .EQ. 1) THEN
            IF (m .EQ. 1) THEN
                field_source_i = i
                field_source_j = j-1
            ELSE 
                field_source_i = i-1
                field_source_j = j
            END IF
        ELSE IF (n .EQ. 2) THEN
            IF (m .EQ. 1) THEN
                field_source_i = i
                field_source_j = j+1
            ELSE 
                field_source_i = i-1
                field_source_j = j
            END IF
        ELSE IF (n .EQ. 3) THEN
            IF (m .EQ. 1) THEN
                field_source_i = i
                field_source_j = j-1
            ELSE 
                field_source_i = i+1
                field_source_j = j
            END IF
        ELSE IF (n .EQ. 4) THEN
            IF (m .EQ. 1) THEN
                field_source_i = i
                field_source_j = j+1
            ELSE 
                field_source_i = i+1
                field_source_j = j
            END IF
        END IF

        field_source(k,j,i) = field_source(k,j,i) + heat_flux*(tfh-tfa)
        field_source(k,field_source_j,field_source_i) = field_source(k,field_source_j,field_source_i) + heat_flux*(tfe-tfh)

    END SUBROUTINE
    !> @brief A subroutine that calculate temperature fields
    !! 
    !> @note This subroutine has the same functionality as the <b>GASTEM</b> 
    !> @date 
    SUBROUTINE is_Temperature_Field_Coverges(is_converged,i_iter,max_fluid_temperature_ratio,ifque)
        LOGICAL,INTENT(INOUT) :: is_converged !< The flag of the convergence of the temperature field
        INTEGER,INTENT(IN) :: i_iter !< The iteration number
        REAL(TH_KDUBLE),INTENT(IN) :: max_fluid_temperature_ratio !< The maximum ratio of the fluid temperature
        INTEGER,INTENT(INOUT) :: ifque !< The flag of the mass source
        

        is_converged = ifque /= 0

        ! IF (is_converged .AND. (max_fluid_temperature_ratio > Th_epsi_gas_tem) .AND. (i_iter < Th_max_iter_gas_tem)) THEN
        !     is_converged = .FALSE.
        !     ifqu = 0
        ! ELSE
        !     ifqu = 1
        ! END IF

        ! IF (ifqu .EQ. 0) THEN
        !     is_converged = .FALSE.
        ! ELSE
        !     is_converged = .TRUE.
        ! END IF

        IF ((max_fluid_temperature_ratio .GT. Th_epsi_gas_tem) .AND. (i_iter .LT. Th_max_iter_gas_tem)) THEN
            is_converged = .FALSE.
            ifque = 0
        ELSE
            ifque = 1
        END IF
    END SUBROUTINE is_Temperature_Field_Coverges
    !> @brief A subroutine that calculate temperature fields
    !! 
    !> @note This subroutine has the same functionality as the <b>GASTEM</b> 
    !> @date 
    ! SUBROUTINE cal_Temperature_Field_Elem4(temperature,i,j,k,source_mass)
    !     REAL(TH_KDUBLE), INTENT(INOUT) :: temperature !< The temperature of this element 
    !     INTEGER, INTENT(IN) :: i,j,k !< The index of this element
    !     REAL(TH_KDUBLE), INTENT(IN) :: source_mass(6)!< The Mass volume source
    !     !进行温度场的求解 尽量把多维数组的索引变成多个一维数组的索引
    !     !通过定义numface来确定一维数组的大小
    !     !数据计算还是从轴向到径向再到周向，这样先写完二维再升级成三维的时候会好一点
    !     !通过mesh计算单个网格的温度，通过field计算所有区域网格的温度，并将温度计算之后进行更新
    !     ! I,N-局部索引编号;PB,MB-内部函数临时参数无需传递;PP-求解得到的温度;QUEL-source_mass-Source_mass_flow质量源;TQUE-Source_temperature_in;FELD-feld;IFQUE-判断参数;ALGA-alga;T-T_solid;TFL-T_fluid;RHO-Helium_density;MZ-Axial_mass_flow;MR-Radial_mass_flow;XKR-Radial_w;XKZ-Axial_mass_flow；
    !     !> cal_Flow_Mesh_Channel
    !     INTEGER :: temp_id,n,m
    !     INTEGER :: num_face,num_direction
    !     REAL(TH_KDUBLE),ALLOCATABLE :: mesh_centre_t_fluid(:),mesh_boundary_t_fluid(:) !< The fluid temperature of the current mesh
    !     REAL(TH_KDUBLE) :: mesh_t_solid !< The solid temperature of the current mesh
    !     REAL(TH_KDUBLE),ALLOCATABLE :: surrounding_t_fluid(:) !< The surrounding fluid temperature
    !     REAL(TH_KDUBLE),ALLOCATABLE :: surrounding_t_solid(:) !< The surrounding solid temperature
    !     REAL(TH_KDUBLE),ALLOCATABLE :: mesh_mass_flow(:) !< The mass flow of the current mesh
    !     REAL(TH_KDUBLE),ALLOCATABLE :: forward_mass_flow(:) !< The forward mass flow
    !     REAL(TH_KDUBLE),ALLOCATABLE :: mesh_heat_flux(:) !< The heat flux of the current mesh
    !     REAL(TH_KDUBLE),ALLOCATABLE :: forward_heat_flux_circumferential(:) !< The forward circumferential heat flux
    !     REAL(TH_KDUBLE),ALLOCATABLE :: forward_radial_heat_flux_radial(:) !< The forward radial heat flux
    !     REAL(TH_KDUBLE),ALLOCATABLE :: forward_axial_heat_flux_axial(:) !< The forward axial heat flux
    !     REAL(TH_KDUBLE) :: mesh_centre_exponent,mesh_boundary_exponent !< The exponent of the current mesh
    !     REAL(TH_KDUBLE),ALLOCATABLE :: forward_centre_exponent(:) !< THe forward centre exponent
    !     REAL(TH_KDUBLE),ALLOCATABLE :: forward_boundary_exponent(:) !< The forward boundary exponent
    !     REAL(TH_KDUBLE):: mesh_heat_transfer_efficiency !< The heat transfer efficiency of the current mesh
    !     REAL(TH_KDUBLE),ALLOCATABLE :: forward_heat_transfer_efficiency(:) !< The forward heat transfer efficiency

    !     mesh_centre_t_fluid = 0.0 !> The fluid temperature of the current mesh
    !     mesh_boundary_t_fluid = 0.0 !> The fluid temperature of the current mesh
    !     num_face = 6 !> 1-left side; 2-top side; 3-right side; 4-bottom side; 5-backward; 6-forward side
    !     num_direction = 3 !> 1-radial direction; 2-axial direction; 3-circumferential direction

        


    !     IF (NOT(ALLOCATED(surrounding_t_fluid))) THEN
    !         ALLOCATE(surrounding_t_fluid(num_face)) !< pb
    !         ALLOCATE(surrounding_t_solid(num_face)) !> pc
    !         temp_id = 2*num_face
    !         ALLOCATE(mesh_centre_t_fluid(temp_id))
    !         ALLOCATE(mesh_boundary_t_fluid(temp_id))
    !     END IF

    !     !> Left side
    !     IF (j .EQ. 1) THEN
    !         surrounding_t_fluid(1) = T_fluid(k,j,i)
    !         surrounding_t_solid(1) = surrounding_t_fluid(1)
    !     ELSE
    !         surrounding_t_fluid(1) = T_fluid(k,j-1,i)
    !         surrounding_t_solid(1) = T_solid(k,j-1,i)
    !     END IF
    !     !> Top side
    !     IF (i .EQ. 1) THEN
    !         surrounding_t_fluid(2) = T_fluid(k,j,i)
    !         surrounding_t_solid(2) = surrounding_t_fluid(2)
    !     ELSE
    !         surrounding_t_fluid(2) = T_fluid(k,j,i-1)
    !         surrounding_t_solid(2) = T_solid(k,j,i-1)
    !     END IF
    !     !> Right side
    !     surrounding_t_fluid(3) = T_fluid(k,j+1,i)
    !     surrounding_t_solid(3) = T_solid(k,j+1,i)
    !     !> Bottom side
    !     surrounding_t_fluid(4) = T_fluid(k,j,i+1)
    !     surrounding_t_solid(4) = T_solid(k,j,i+1)
    !     !> Solid temperature of the current mesh
    !     mesh_t_solid = T_solid(k,j,i)

    !     IF (NOT(ALLOCATED(mesh_mass_flow))) THEN
    !         ALLOCATE(mesh_mass_flow(num_direction))
    !         ALLOCATE(forward_mass_flow(2*num_direction))
    !         ALLOCATE(mesh_heat_flux(num_face,num_direction))
    !     END IF

    !     !> The circumferential mass flow of the current mesh
    !     forward_mass_flow(1) = Radial_mass_flow(k,j,i)
    !     !> The radial mass flow of the current mesh
    !     forward_mass_flow(2) = Axial_mass_flow(k,j,i)
    !     !> The axial mass flow of the current mesh
    !     forward_mass_flow(3) = Circumferential_mass_flow(k,j,i)
    !     !> The forward circumferential mass flow
    !     forward_mass_flow(4) = Radial_mass_flow(k,j+1,i)
    !     !> The forward radial mass flow
    !     forward_mass_flow(5) = Axial_mass_flow(k,j,i+1)
    !     !> The forward axial mass flow
    !     forward_mass_flow(6) = Circumferential_mass_flow(k+1,j,i)

    !     mesh_heat_flux(1,1) = heat_flux_radial(k,j,i)
    !     mesh_heat_flux(1,2) = heat_flux_axial(k,j,i)
    !     mesh_heat_flux(1,3) = heat_flux_circumferential(k,j,i)
    !     mesh_heat_flux(2,1) = -heat_flux_radial(k,j+1,i)
    !     mesh_heat_flux(2,2) = heat_flux_axial(k,j+1,i)
    !     mesh_heat_flux(2,3) = heat_flux_circumferential(k,j+1,i)
    !     mesh_heat_flux(3,1) = heat_flux_radial(k,j,i+1)
    !     mesh_heat_flux(3,2) = -heat_flux_axial(k,j,i+1)
    !     mesh_heat_flux(3,3) = heat_flux_circumferential(k,j,i+1)
    !     mesh_heat_flux(4,1) = -heat_flux_radial(k,j+1,i+1)
    !     mesh_heat_flux(4,2) = -heat_flux_axial(k,j+1,i+1)
    !     mesh_heat_flux(4,3) = heat_flux_circumferential(k,j+1,i+1)
    !     mesh_heat_flux(5,1) = heat_flux_radial(k+1,j,i)
    !     mesh_heat_flux(5,2) = heat_flux_axial(k+1,j,i)
    !     mesh_heat_flux(5,3) = -heat_flux_circumferential(k+1,j,i)
    !     mesh_heat_flux(6,1) = -heat_flux_radial(k+1,j+1,i)
    !     mesh_heat_flux(6,2) = heat_flux_axial(k+1,j+1,i)
    !     mesh_heat_flux(6,3) = -heat_flux_circumferential(k+1,j+1,i)



    !     IF (NOT(ALLOCATED(forward_heat_flux_circumferential))) THEN
    !         temp_id = 2**num_direction !> C(num_dircetion,1)->C(num_dircetion,2)->...->C(num_dircetion,num_dircetion)
    !         ALLOCATE(forward_heat_flux_circumferential(temp_id))
    !         ALLOCATE(forward_radial_heat_flux_radial(temp_id))
    !         ALLOCATE(forward_axial_heat_flux_axial(temp_id))
    !         ALLOCATE(forward_centre_exponent(temp_id))
    !         ALLOCATE(forward_boundary_exponent(temp_id))
    !         ALLOCATE(forward_heat_transfer_efficiency(temp_id))
    !     END IF


    !     !> The circumferential mass flow of the current mesh
    !     mesh_mass_flow(1) = Circumferential_mass_flow(k,j,i)
    !     !> The radial mass flow of the current mesh
    !     mesh_mass_flow(2) = Radial_mass_flow(k,j,i)
    !     !> The axial mass flow of the current mesh
    !     mesh_mass_flow(3) = Axial_mass_flow(k,j,i)
    !     !> The forward circumferential mass flow
    !     forward_mass_flow(1) = Circumferential_mass_flow(k+1,j,i)
    !     !> The forward radial mass flow
    !     forward_mass_flow(2) = Radial_mass_flow(k,j+1,i)
    !     !> The forward axial mass flow
    !     forward_mass_flow(3) = Axial_mass_flow(k,j,i+1)
    !     !> The circumferential heat flux of the current mesh
    !     mesh_heat_flux(1) = heat_flux_circumferential(k,j,i)
    !     !> The radial heat flux of the current mesh
    !     mesh_heat_flux(2) = heat_flux_radial(k,j,i)
    !     !> The axial heat flux of the current mesh
    !     mesh_heat_flux(3) = heat_flux_axial(k,j,i)
    !     !> The forward circumferential heat flux
    !     forward_heat_flux_circumferential(1) = -heat_flux_circumferential(k,j,i)
    !     forward_heat_flux_circumferential(2) = -heat_flux_circumferential(k+1,j,i)
    !     forward_axial_heat_cicumferential(3) = heat_flux_circumferential(k,j+1,i)
    !     forward_axial_heat_cicumferential(4) = heat_flux_circumferential(k,j,i+1)
    !     forward_heat_flux_circumferential(5) = -heat_flux_circumferential(k+1,j+1,i)
    !     forward_heat_flux_circumferential(6) = -heat_flux_circumferential(k+1,j,i+1)
    !     forward_heat_flux_circumferential(7) = heat_flux_circumferential(k,j+1,i+1)
    !     forward_heat_flux_circumferential(8) = -heat_flux_circumferential(k+1,j+1,i+1)
    !     !> The forward radial heat flux
    !     forward_radial_heat_flux_radial(1) = heat_flux_radial(k+1,j,i)
    !     forward_radial_heat_flux_radial(2) = -heat_flux_radial(k,j+1,i)
    !     forward_radial_heat_flux_radial(3) = heat_flux_radial(k,j,i+1)
    !     forward_radial_heat_flux_radial(4) = -heat_flux_radial(k+1,j+1,i)
    !     forward_radial_heat_flux_radial(5) = heat_flux_radial(k+1,j,i+1)
    !     forward_radial_heat_flux_radial(6) = -heat_flux_radial(k,j+1,i+1)
    !     forward_radial_heat_flux_radial(7) = -heat_flux_radial(k+1,j+1,i+1)
    !     !> The forward axial heat flux
    !     forward_axial_heat_flux_axial(1) = heat_flux_axial(k+1,j,i)
    !     forward_axial_heat_flux_axial(2) = heat_flux_axial(k,j+1,i)
    !     forward_axial_heat_flux_axial(3) = -heat_flux_axial(k,j,i+1)
    !     forward_axial_heat_flux_axial(4) = heat_flux_axial(k+1,j+1,i)
    !     forward_axial_heat_flux_axial(5) = -heat_flux_axial(k+1,j,i+1)
    !     forward_axial_heat_flux_axial(6) = -heat_flux_axial(k,j+1,i+1)
    !     forward_axial_heat_flux_axial(7) = -heat_flux_axial(k+1,j+1,i+1)
    !     !> The centre exponent of the current mesh
    !     mesh_centre_exponent = alga(k,j,i,1)
    !     !> The boundary exponent of the current mesh
    !     mesh_boundary_exponent = alga(k,j,i,2)
    !     !> The forward centre exponent
    !     forward_centre_exponent(1) = alga(k+1,j,i,1)
    !     forward_centre_exponent(2) = alga(k,j+1,i,1)
    !     forward_centre_exponent(3) = alga(k,j,i+1,1)
    !     forward_centre_exponent(4) = alga(k+1,j+1,i,1)
    !     forward_centre_exponent(5) = alga(k+1,j,i+1,1)
    !     forward_centre_exponent(6) = alga(k,j+1,i+1,1)
    !     forward_centre_exponent(7) = alga(k+1,j+1,i+1,1)
    !     !> The forward boundary exponent
    !     forward_boundary_exponent(1) = alga(k+1,j,i,2)
    !     forward_boundary_exponent(2) = alga(k,j+1,i,2)
    !     forward_boundary_exponent(3) = alga(k,j,i+1,2)
    !     forward_boundary_exponent(4) = alga(k+1,j+1,i,2)
    !     forward_boundary_exponent(5) = alga(k+1,j,i+1,2)
    !     forward_boundary_exponent(6) = alga(k,j+1,i+1,2)
    !     forward_boundary_exponent(7) = alga(k+1,j+1,i+1,2)
    !     !> The heat transfer efficiency of the current mesh
    !     mesh_heat_transfer_efficiency = temp_real_temperature_cal(k,j,i)
    !     !> The forward heat transfer efficiency
    !     forward_heat_transfer_efficiency(1) = temp_real_temperature_cal(k+1,j,i)
    !     forward_heat_transfer_efficiency(2) = temp_real_temperature_cal(k,j+1,i)
    !     forward_heat_transfer_efficiency(3) = temp_real_temperature_cal(k,j,i+1)
    !     forward_heat_transfer_efficiency(4) = temp_real_temperature_cal(k+1,j+1,i)
    !     forward_heat_transfer_efficiency(5) = temp_real_temperature_cal(k+1,j,i+1)
    !     forward_heat_transfer_efficiency(6) = temp_real_temperature_cal(k,j+1,i+1)
    !     forward_heat_transfer_efficiency(7) = temp_real_temperature_cal(k+1,j+1,i+1)

    !     !> Start calculating the mesh temperature
    !     !> The calculation of the current mesh
    !     mesh_centre_t_fluid(1) = 0.0
    !     DO n = 1,2**num_direction
    !         DO m = 1, num_direction 
    !             IF (mesh_heat_flux(n,m) .LE. 0.0) THEN
                    
    !             ELSE
                    
    !             END IF
    !         END DO

    !     END DO


    ! END SUBROUTINE cal_Temperature_Field_Elem4
END MODULE HeliumTemperatureSolver