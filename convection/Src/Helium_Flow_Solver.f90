!> @brief A module for solving helium flow field
!! 
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE HeliumFlowSolver
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation
    USE GlobalTHVariables !< Global variables definition for thermal hydraulic calculation
    USE HeliumPhysicalPropertiesCalculation !< Physical property parameter calculation of helium for thermal hydraulic calculation
    USE HeliumCharacteristicNumberCalculation !< Characteristic number calculation of helium for thermal hydraulic calculation
    IMPLICIT NONE

    REAL(TH_KDUBLE),ALLOCATABLE :: Axial_w(:,:,:) !< The axial interphase friction factor 
    REAL(TH_KDUBLE),ALLOCATABLE :: Radial_w(:,:,:) !< The radial interphase friction factor 
    REAL(TH_KDUBLE),ALLOCATABLE :: Circumferential_w(:,:,:) !< The circumferential interphase friction factor 
    REAL(TH_KDUBLE),ALLOCATABLE :: rogg(:,:,:) !< The gravitational pressure drop
    REAL(TH_KDUBLE),ALLOCATABLE :: sum_rogg(:,:,:) !< The sum of gravitational pressure drop
    REAL(TH_KDUBLE),ALLOCATABLE :: sum_w(:,:,:) !< The sum of interphase friction factor
    REAL(TH_KDUBLE),ALLOCATABLE :: temp_p_fluid(:,:,:) !< The temp real variable in solving helium flow field
    REAL(TH_KDUBLE),ALLOCATABLE :: temp_real_flow_cal(:,:,:) !< The temp real variable in solving helium flow field

    INTEGER,ALLOCATABLE :: material_channel1(:) !< The material id of channel type 1
    REAL(TH_KDUBLE),ALLOCATABLE :: volume_channel1(:) !< The volume of channel type 1

    INTEGER,ALLOCATABLE :: num_channel5 !< The number of channel type 5
    INTEGER,ALLOCATABLE :: material_channel5(:) !< The material id of channel type 5
    INTEGER,ALLOCATABLE :: theta_channel5(:) !< The theta index of channel type 5
    INTEGER,ALLOCATABLE :: radial_channel5(:) !< The radial index of channel type 5
    INTEGER,ALLOCATABLE :: axial_channel5(:) !< The axial index of channel type 5
    INTEGER,ALLOCATABLE :: num_mesh_channel5(:) !< The number of mesh of channel type 5
    REAL(TH_KDUBLE),ALLOCATABLE :: volume_channel5(:) !< The volume of channel type 5

CONTAINS
    !> @brief A subroutine that initialize the flow variable
    !! 
    !> @date 
    SUBROUTINE init_Flow_Variable()
        !遍历所有流道，判断流道类型
        !计算密度 CALL cal_Helium_Density()
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        INTEGER :: i,j,k,n !< The index of mesh
        INTEGER :: temp_id !< The material id of mesh

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        ALLOCATE(rogg(theta_num,radial_num,axial_num))
        rogg = TH_REAL_ZERO
        ALLOCATE(temp_p_fluid(theta_num,radial_num,axial_num))
        temp_p_fluid = TH_REAL_ZERO
        ALLOCATE(Axial_w(theta_num,radial_num,axial_num))
        Axial_w = TH_REAL_ZERO
        ALLOCATE(Radial_w(theta_num,radial_num,axial_num))
        Radial_w = TH_REAL_ZERO
        ALLOCATE(Circumferential_w(theta_num,radial_num,axial_num))
        Circumferential_w = TH_REAL_ZERO

        IF (NOT(ALLOCATED(temp_real_flow_cal))) THEN
            ALLOCATE(temp_real_flow_cal(theta_num,radial_num,axial_num))
        END IF
        temp_real_flow_cal = TH_REAL_ZERO
        
        !> @todo Determine the number of each channel type
        n = 0
        j = 0
        k = 0
        DO i = 1, SIZE(Channels)
            IF (Channels(i)%channel_type .EQ. 1) THEN
                n = n + 1
            END IF
            IF (Channels(i)%channel_type .EQ. 2) THEN
                j = j + 1
            END IF
            IF (Channels(i)%channel_type .EQ. 5) THEN
                k = k + 1
            END IF
        END DO

        !> @todo Allocate the memory for channel type 1
        ALLOCATE(material_channel1(n))
        ALLOCATE(volume_channel1(n))
        material_channel1 = 0
        volume_channel1 = 0
        !> @todo Allocate the memory for th_mesh%sum_axial_mesh
        ALLOCATE(th_mesh%sum_axial_mesh(theta_num,radial_num,j))
        th_mesh%sum_axial_mesh = 0.0
        !> @todo Allocate the memory for channel type 2
        ALLOCATE(material_channel2(j))
        ALLOCATE(volume_channel2(j))
        ALLOCATE(upper_channel2(theta_num,radial_num,j))
        ALLOCATE(lower_channel2(theta_num,radial_num,j))
        material_channel2 = 0
        volume_channel2 = 0
        upper_channel2 = 0
        lower_channel2 = 0
        !> @todo Allocate the memory for channel type 5
        ALLOCATE(material_channel5(k))
        ALLOCATE(theta_channel5(k))
        ALLOCATE(radial_channel5(k))
        ALLOCATE(axial_channel5(k))
        ALLOCATE(num_mesh_channel5(k))
        ALLOCATE(volume_channel5(k))
        num_channel5 = 0
        material_channel5 = 0
        theta_channel5 = 0
        radial_channel5 = 0
        axial_channel5 = 0
        num_mesh_channel5 = 0
        volume_channel5 = 0.0
        !> @todo Allocate the memory for sum_rogg and sum_w
        ALLOCATE(sum_rogg(theta_num,radial_num,j))
        ALLOCATE(sum_w(theta_num,radial_num,j))
        sum_rogg = 0
        sum_w = 0
        !> @todo Determine the material id of each channel type
        j = 0
        n = 0
        k = 0
        DO i = 1, SIZE(Channels)
            IF (Channels(i)%channel_type .EQ. 1) THEN
                n = n + 1
                material_channel1(n) = Compositions(i)%composition_id
            END IF
            IF (Channels(i)%channel_type .EQ. 2) THEN
                j = j + 1
                material_channel2(j) = Compositions(i)%composition_id
            END IF
            IF (Channels(i)%channel_type .EQ. 5) THEN
                k = k + 1
                material_channel5(k) = Compositions(i)%composition_id
            END IF
        END DO
        !> @todo Determine the volume of channel type 1
        DO j = 2, radial_num-1
            DO k = 2, theta_num-1
                DO n = 1, SIZE(material_channel1)
                    DO i = 2, axial_num-1
                        temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                        IF (temp_id .EQ. material_channel1(n)) THEN
                            volume_channel1(n) = volume_channel1(n) + Th_mesh%th_mesh_volume(k,j,i)
                        END IF
                    END DO
                END DO
            END DO
        END DO
        !> @todo Change the unit of source mass flow of channel type 1
        DO i = 1, SIZE(material_channel1)
            Channels(material_channel1(i))%source_mass_flow = Channels(material_channel1(i))%source_mass_flow / volume_channel1(i)
        END DO
        !> @todo Determine the volume of channel type 2
        !> @todo Determine the upper and lower boundary of channel type 2
        DO j = 2, radial_num-1
            DO k = 2, theta_num-1
                DO n = 1, SIZE(material_channel2)
                    DO i = 2, axial_num-1
                        temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                        IF (temp_id .EQ. material_channel2(n)) THEN
                            volume_channel2(n) = volume_channel2(n) + Th_mesh%th_mesh_volume(k,j,i)
                            IF (Th_mesh%th_fine_mesh_material_id(k,j,i-1) .NE. material_channel2(n)) THEN
                                upper_channel2(k,j,n) = i
                            END IF
                            IF (Th_mesh%th_fine_mesh_material_id(k,j,i+1) .NE. material_channel2(n)) THEN
                                lower_channel2(k,j,n) = i
                            END IF
                        END IF
                    END DO
                END DO
            END DO
        END DO
        !> @todo Change the unit of source mass flow of channel type 2
        DO i = 1, SIZE(material_channel2)
            Channels(material_channel2(i))%source_mass_flow = Channels(material_channel2(i))%source_mass_flow / volume_channel2(i)
        END DO
        DO j = 2, radial_num-1
            DO k = 2, theta_num-1
                DO n = 1, SIZE(material_channel2)
                    DO i = upper_channel2(k,j,n),lower_channel2(k,j,n)
                        IF (i .EQ. 0) THEN
                            CYCLE
                        ELSE
                            Th_mesh%sum_axial_mesh(k,j,n) = Th_mesh%sum_axial_mesh(k,j,n) + Th_mesh%axial_fine_mesh_size(i)
                        END IF
                    END DO
                END DO
            END DO
        END DO
        !> @todo Determine the volume of channel type 5
        !> @todo Determine the number of mesh of channel type 5
        !> @todo Determine the theta, radial and axial index of channel type 5
        DO n = 1, SIZE(Compositions)
            IF (Channels(n)%channel_type .EQ. 5) THEN
                DO i = 2, axial_num-1
                    DO j = 2, radial_num-1
                        DO k = 2, theta_num-1
                            temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                            IF (temp_id .EQ. n) THEN
                                temp_id = Th_mesh%th_fine_mesh_material_id(k,j-1,i)
                                IF (temp_id .EQ. 0) THEN
                                    num_channel5 = num_channel5 + 1
                                    theta_channel5(num_channel5) = k
                                    radial_channel5(num_channel5) = j
                                    axial_channel5(num_channel5) = i
                                    num_mesh_channel5(num_channel5) = num_mesh_channel5(num_channel5)+1
                                    volume_channel5(num_channel5) = volume_channel5(num_channel5) + Th_mesh%th_mesh_volume(k,j,i)
                                ELSE
                                    IF (Channels(temp_id)%channel_type .NE. 5) THEN
                                        num_channel5 = num_channel5 + 1
                                        theta_channel5(num_channel5) = k
                                        radial_channel5(num_channel5) = j
                                        axial_channel5(num_channel5) = i
                                        num_mesh_channel5(num_channel5) = num_mesh_channel5(num_channel5)+1
                                        volume_channel5(num_channel5) = volume_channel5(num_channel5) + Th_mesh%th_mesh_volume(k,j,i)
                                    ELSE
                                        num_mesh_channel5(num_channel5) = num_mesh_channel5(num_channel5)+1
                                        volume_channel5(num_channel5) = volume_channel5(num_channel5) + Th_mesh%th_mesh_volume(k,j,i)
                                    END IF
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            ELSE    
                CYCLE
            END IF
        END DO
        !> @todo Change the unit of source mass flow of channel type 5
        DO n = 1, num_channel5
            i = axial_channel5(n)
            j = radial_channel5(n)
            k = theta_channel5(n)
            temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
            Channels(temp_id)%source_mass_flow = Channels(temp_id)%source_mass_flow / volume_channel5(n)
        END DO


    END SUBROUTINE init_Flow_Variable
    !> @brief A subroutine that set the flow resistance of all flow channels
    !! 
    !> @note This subroutine has the same functionality as the <b>vorp</b> subroutine
    !> @date 
    SUBROUTINE set_Flow_Resistance()
        !遍历所有流道，判断流道类型
        ! rogg-rogg;tfl-T_fluid;RHO-Helium_density;IFB,DZ,KOM-Th_mesh;MZ-Axial_mass_flow;MR-Radial_mass_flow;XKR-Radial_w;XKZ-Axial_mass_flow;STROM-temp_real_flow_cal;SUMRG-sum_rogg
        !计算密度 CALL cal_Helium_Density()
        INTEGER :: i,j,k,n,temp_id !< The index of mesh
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction
        REAL(TH_KDUBLE) :: temp_real !< The temp real variable in solving helium flow field
        
        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        !> Variable initialization
        Helium_density = TH_REAL_ZERO
        Axial_w = TH_REAL_ZERO
        Radial_w = TH_REAL_ZERO
        Circumferential_w = TH_REAL_ZERO
        rogg = TH_REAL_ZERO
        sum_rogg = TH_REAL_ZERO
        sum_w = TH_REAL_ZERO
        Axial_mass_flow = TH_REAL_ZERO
        Radial_mass_flow = TH_REAL_ZERO
        Circumferential_mass_flow = TH_REAL_ZERO
        temp_p_fluid = TH_REAL_ZERO

        DO i = 1, axial_num
            DO j = 1, radial_num
                DO k = 1, theta_num
                    temp_id = Th_mesh%th_mesh_type(k,j,i)
                    IF (temp_id .EQ. 0) THEN
                        IF (k .NE. 1) Circumferential_w(k-1,j,i) = TH_REAL_INFINITY
                        IF (j .NE. 1) Radial_w(k,j-1,i) = TH_REAL_INFINITY
                        IF (i .NE. 1) Axial_w(k,j,i-1) = TH_REAL_INFINITY
                        Circumferential_w(k,j,i) = TH_REAL_INFINITY
                        Radial_w(k,j,i) = TH_REAL_INFINITY
                        Axial_w(k,j,i) = TH_REAL_INFINITY
                    ELSE IF (temp_id .EQ. 1) THEN
                        CYCLE
                    ELSE IF (temp_id .EQ. 2) THEN
                        IF (k .NE. 1) Circumferential_w(k-1,j,i) = TH_REAL_INFINITY
                        IF (j .NE. 1) Radial_w(k,j-1,i) = TH_REAL_INFINITY
                        Circumferential_w(k,j,i) = TH_REAL_INFINITY
                        Radial_w(k,j,i) = TH_REAL_INFINITY
                    ELSE IF (temp_id .EQ. 5) THEN
                        CYCLE
                    END IF
                END DO
            END DO
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*)"xudongyu:axial_w"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,Axial_w(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:radial_w"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,Radial_w(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:circumferential_w"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,Circumferential_w(1,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)

        theta_num = theta_num - 1 !< 2
        radial_num = radial_num - 1 !< 29
        axial_num = axial_num - 1 !< 50

        DO i = 2, axial_num
            DO j = 2, radial_num
                DO k = 2, theta_num
                    temp_real = T_fluid(k,j,i) + T_fluid(k-1,j,i) + T_fluid(k,j-1,i) + T_fluid(k,j,i-1)
                    temp_real = temp_real + T_fluid(k-1,j-1,i) + T_fluid(k-1,j,i-1) + T_fluid(k,j-1,i-1)
                    temp_real = temp_real + T_fluid(k-1,j-1,i-1)
                    temp_real = temp_real / 8
                    IF (i .EQ. 13 .AND. j .EQ. 2) THEN
                        WRITE(*,*) T_fluid(k,j,i), T_fluid(k-1,j,i), T_fluid(k,j-1,i), T_fluid(k,j,i-1)
                        WRITE(*,*) T_fluid(k-1,j-1,i), T_fluid(k-1,j,i-1), T_fluid(k,j-1,i-1)
                        WRITE(*,*) T_fluid(k-1,j-1,i-1)
                        WRITE(*,*) temp_real
                    END IF
                    CALL cal_Helium_Density(Helium_density(k,j,i),gas_pressure,temp_real)
                END DO
            END DO
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*)"xudongyu:Hslium density"
        ! DO i = 1,axial_num+1
        !     DO j  = 1,radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Helium_density(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)

        axial_num = axial_num - 1 !< 49
        DO i = 2, axial_num
            DO j = 2, radial_num
                DO k = 1, theta_num
                    rogg(k,j,i) = (Helium_density(k,j,i)*Th_mesh%axial_fine_mesh_size(i)+Helium_density(k,j,i+1)*Th_mesh%axial_fine_mesh_size(i+1))
                    rogg(k,j,i) = HALF * TH_GRAVITY * rogg(k,j,i)
                    temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                    IF (temp_id .EQ. 0) THEN
                        CYCLE
                    END IF
                    IF (Channels(temp_id)%channel_type .EQ. 2) THEN
                        DO n = 1, SIZE(material_channel2)
                            IF (temp_id .EQ. material_channel2(n)) THEN
                                sum_rogg(k,j,n) = sum_rogg(k,j,n) + rogg(k,j,i)
                            ELSE
                                CYCLE
                            END IF
                        END DO
                    END IF
                END DO
            END DO
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*)"xudongyu:rogg"
        ! DO i = 1,axial_num+2
        !     DO j  = 1,radial_num+1
        !         WRITE(OUTUNIT,*) i,j,rogg(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:summg"
        ! DO i = 1,SIZE(material_channel2)
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,sum_rogg(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)

    END SUBROUTINE set_Flow_Resistance
    !> @brief A subroutine that calculate pressure and flow fields
    !! 
    !> @note This subroutine has the same functionality as the <b>STROEM</b> 
    !> @date 
    SUBROUTINE cal_Flow_Field(is_cal_xkk)
        LOGICAL,INTENT(INOUT) :: is_cal_xkk !< A logical variable that indicates whether to calculate the additional pressure drop
        REAL(TH_KDUBLE) :: min_pressure !< The minimum pressure value in the field
        INTEGER :: i_iter,i_extrapolation
        LOGICAL :: is_converged = .FALSE.
        INTEGER :: i,j,k,n,temp_id !< The index of mesh
        INTEGER :: theta_num,radial_num,axial_num !< The number of mesh in each direction

        
        i_iter = 0
        i_extrapolation = 0
        is_converged = .FALSE.

        !开始迭代
        !计算流体计算常数 CALL cal_Flow_Friction_Factor()
        DO WHILE(NOT(is_converged))
            i_iter = i_iter + 1
            ! IF (i_iter .EQ. 50) THEN
            !     WRITE(*,*) ""
            ! END IF
            min_pressure = TH_REAL_INFINITY
            ! WRITE(*,*) "iter:",i_iter
            CALL cal_Flow_Friction_Factor(is_cal_xkk)
            is_cal_xkk = .TRUE.
            CALL cal_Flow_Field_Channel5(min_pressure)
            CALL cal_Flow_Field_Channel1(min_pressure)
            CALL cal_Flow_Field_Channel2(min_pressure)
            CALL update_Flow_Field(min_pressure,i_iter,i_extrapolation)
            CALL is_Flow_Field_Coverges(is_converged,i_iter)
        END DO


        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51
        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*)"xudongyu:final pressure"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:final Axial_mass_flow"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,Axial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:final Radial_mass_flow"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,Radial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)
    END SUBROUTINE cal_Flow_Field
    !> @brief A subroutine that calculate friction factor required for pressure field calculation
    !! 
    !> @note This subroutine has the same functionality as the <b>STRLAM</b> 
    !> @date 
    SUBROUTINE cal_Flow_Friction_Factor(is_cal_xkk)
        LOGICAL,INTENT(IN) :: is_cal_xkk !< A logical variable 
        INTEGER :: i,j,k,n,temp_id
        INTEGER :: theta_num,radial_num,axial_num
        REAL(TH_KDUBLE) :: temp_real
        REAL(TH_KDUBLE) :: viscosity,velocity,relo
        

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        ! ALLOCATE(temp_real_flow_cal(theta_num,radial_num,axial_num))
        ! temp_real_flow_cal = TH_REAL_ZERO
        !遍历所有流道，判断流道类型
        !计算动力粘度 CALL cal_Helium_Viscosity()
        !计算雷洛数
        !计算标量速度，然后计算每个节块的阻力
        ! IFXKK用于控制跳过计算，想把它优化掉；TFL-T_fluid;RHO-Helium_density;IFB,FZQ,FRQ,DR,DZ,KOM,EPSIL-TH_mesh;DHYD-Hydraulic_diameter;T-T_solid;MZ-Axial_mass_flow;MR-Radial_mass_flow;XKON-Additional_pressure_drop;SUMXK-sum_w;QTV-Qtv;LTV-Ltv;XGEO-Xgeo;
        sum_w = 0.0
        ! theta_num = theta_num - 1
        ! radial_num = radial_num - 1
        ! DO n = 1, SIZE(material_channel2)
        !     DO j = 2, radial_num
        !         DO k  = 2, theta_num
        !             sum_w(k,j,n) = 0
        !         END DO
        !     END DO
        ! END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*)"xudongyu:Radial_mass_flow"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,Radial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:Axial_mass_flow"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,Axial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)

        theta_num = theta_num - 1 !< 2
        radial_num = radial_num - 1 !< 29
        axial_num = axial_num - 1 !< 50
        IF(is_cal_xkk) THEN
            DO i = 2, axial_num
                DO j = 2, radial_num
                    DO k = 2, theta_num
                        temp_id = th_mesh%th_fine_mesh_material_id(k,j,i)
                        IF (Channels(temp_id)%channel_type .EQ. 0) THEN
                            CYCLE
                        ELSE IF (Channels(temp_id)%channel_type .EQ. 1) THEN
                            temp_real = T_fluid(k,j,i) + T_fluid(k-1,j,i) + T_fluid(k,j-1,i) + T_fluid(k,j,i-1)
                            temp_real = temp_real + T_fluid(k-1,j-1,i) + T_fluid(k-1,j,i-1) + T_fluid(k,j-1,i-1)
                            temp_real = temp_real + T_fluid(k-1,j-1,i-1)
                            temp_real = temp_real / 8

                            CALL cal_Helium_Viscosity(viscosity,temp_real)
                            temp_real = temp_real * 8
                            temp_real = temp_real + T_solid(k,j,i)
                            temp_real = temp_real + T_solid(k-1,j,i) + T_solid(k,j-1,i) + T_solid(k,j,i-1)
                            temp_real = temp_real + T_solid(k-1,j-1,i) + T_solid(k-1,j,i-1) + T_solid(k,j-1,i-1)
                            temp_real = temp_real + T_solid(k-1,j-1,i-1)
                            temp_real = temp_real / 16
                            CALL cal_Helium_Re(relo,temp_real,Helium_density(k,j,i),Channels(temp_id)%channel_type,Th_mesh%th_mesh_porosity(k,j,i), &
                                           &  Compositions(temp_id)%diameter_fuel_element,Channels(temp_id)%hydraulic_diameter,Th_mesh%th_mesh_area(k,j,i,3), &
                                           &  Th_mesh%th_mesh_area(k,j,i,2),Th_mesh%th_mesh_area(k,j,i,1),Axial_mass_flow(k,j,i),Radial_mass_flow(k,j,i),Circumferential_mass_flow  (k,j,i))
                            ! relo = 1.5*relo*(1-Th_mesh%th_mesh_porosity(k,j,i))

                            relo = relo**0.9
                            relo = 142.22+3.663*relo
                            temp_real_flow_cal(k,j,i) = 0.75/Compositions(temp_id)%diameter_fuel_element*relo*viscosity
                            temp_real_flow_cal(k,j,i) = temp_real_flow_cal(k,j,i)*(1-Th_mesh%th_mesh_porosity(k,j,i))**2
                            temp_real_flow_cal(k,j,i) = 1.5*temp_real_flow_cal(k,j,i)/Th_mesh%th_mesh_porosity(k,j,i)**3/Compositions(temp_id)%diameter_fuel_element
                            ! IF (relo .LT. 2.0) THEN
                            !     nus = 4.03
                            ! ELSE

                            ! END IF
                        ELSE IF (Channels(temp_id)%channel_type .EQ. 2) THEN
                            CALL cal_Scalar_Velocity(velocity,Helium_density(k,j,i),th_mesh%th_mesh_area(k,j,i,3),th_mesh%th_mesh_area(k,j,i,2), &
                                                    & th_mesh%th_mesh_area(k,j,i,1),Axial_mass_flow(k,j,i),Radial_mass_flow(k,j,i),Circumferential_mass_flow(k,j,i))
                            IF (velocity .EQ. 0) THEN
                                velocity = 0.0002
                            END IF
                            temp_real = T_fluid(k,j,i) + T_fluid(k-1,j,i) + T_fluid(k,j-1,i) + T_fluid(k,j,i-1)
                            temp_real = temp_real + T_fluid(k-1,j-1,i) + T_fluid(k-1,j,i-1) + T_fluid(k,j-1,i-1)
                            temp_real = temp_real + T_fluid(k-1,j-1,i-1)
                            temp_real = temp_real + T_solid(k,j,i)
                            temp_real = temp_real + T_solid(k-1,j,i) + T_solid(k,j-1,i) + T_solid(k,j,i-1)
                            temp_real = temp_real + T_solid(k-1,j-1,i) + T_solid(k-1,j,i-1) + T_solid(k,j-1,i-1)
                            temp_real = temp_real + T_solid(k-1,j-1,i-1)
                            temp_real = temp_real / 16
                            CALL cal_Helium_Re(relo,temp_real,Helium_density(k,j,i),Channels(temp_id)%channel_type,Th_mesh%th_mesh_porosity(k,j,i), &
                                           &  Compositions(temp_id)%diameter_fuel_element,Channels(temp_id)%hydraulic_diameter,Th_mesh%th_mesh_area(k,j,i,3), &
                                           &  Th_mesh%th_mesh_area(k,j,i,2),Th_mesh%th_mesh_area(k,j,i,1),Axial_mass_flow(k,j,i),Radial_mass_flow(k,j,i),Circumferential_mass_flow  (k,j,i))
                            CALL cal_Helium_AxialRe(relo)
                            temp_real_flow_cal(k,j,i) = HALF*velocity*Helium_density(k,j,i)/(th_mesh%th_mesh_porosity(k,j,i))**2*(relo/Channels(temp_id)%hydraulic_diameter+Channels    (temp_id)%additional_pressure_drop)
                        ELSE IF (Channels(temp_id)%channel_type .EQ. 5) THEN
                            CYCLE
                        END IF
                    END DO 
                END DO
            END DO
        END IF

        ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*)"xudongyu:th_mesh_porosity"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,th_mesh%th_mesh_porosity(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:xkk"
        ! DO i = 1,axial_num
        !     DO j  = 1,radial_num
        !         WRITE(OUTUNIT,*) i,j,temp_real_flow_cal(2,j,i)
        !     END DO
        ! END DO
        ! CLOSE(OUTUNIT)

        DO i = 2, axial_num
            DO j = 2, radial_num
                DO k = 2, theta_num
                    temp_id = th_mesh%th_fine_mesh_material_id(k,j,i)
                    IF (Channels(temp_id)%channel_type .EQ. 0) THEN
                        CYCLE
                    ELSE
                        IF (Radial_w(k,j,i) .NE. TH_REAL_INFINITY) THEN 
                            Radial_w(k,j,i) = temp_real_flow_cal(k,j,i)/Helium_density(k,j,i)/th_mesh%th_mesh_area(k,j,i,2)*th_mesh%radial_fine_mesh_size(j)
                            Radial_w(k,j,i) = Radial_w(k,j,i) + temp_real_flow_cal(k,j+1,i)/Helium_density(k,j+1,i)/Th_mesh%th_mesh_area(k,j+1,i,2)*Th_mesh%radial_fine_mesh_size(j+1)
                            Radial_w(k,j,i) = HALF*Radial_w(k,j,i)
                        END IF
                        IF (Axial_w(k,j,i) .NE. TH_REAL_INFINITY) THEN
                            Axial_w(k,j,i) = temp_real_flow_cal(k,j,i)/Helium_density(k,j,i)/Th_mesh%th_mesh_area(k,j,i,3)*Th_mesh%axial_fine_mesh_size(i)
                            Axial_w(k,j,i) = Axial_w(k,j,i) + temp_real_flow_cal(k,j,i+1)/Helium_density(k,j,i+1)/Th_mesh%th_mesh_area(k,j,i+1,3)*Th_mesh%axial_fine_mesh_size(i+1)
                            Axial_w(k,j,i) = HALF*Axial_w(k,j,i)
                            ! IF(i .EQ. 34) WRITE(*,*) i,j,Axial_w(k,j,i)
                            ! IF(i .EQ. 34) WRITE(*,*) temp_real_flow_cal(k,j,i),Helium_density(k,j,i),Th_mesh%th_mesh_area(k,j,i,3),Th_mesh%axial_fine_mesh_size(i)
                            ! IF(i .EQ. 34) WRITE(*,*) temp_real_flow_cal(k,j,i+1),Helium_density(k,j,i+1),Th_mesh%th_mesh_area(k,j,i+1,3),Th_mesh%axial_fine_mesh_size(i+1)
                        END IF
                        IF (Channels(temp_id)%channel_type .EQ. 2) THEN
                            DO n = 1, SIZE(material_channel2)
                                IF (temp_id .EQ. material_channel2(n)) THEN
                                    sum_w(k,j,n) = sum_w(k,j,n) + Axial_w(k,j,i)
                                END IF
                            END DO
                        END IF
                    END IF
                END DO
            END DO
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*)"xudongyu:Radial_w"
        ! DO i = 1,axial_num+1
        !     DO j  = 1,radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Radial_w(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:Axial_w"
        ! DO i = 1,axial_num+1
        !     DO j  = 1,radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Axial_w(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:frq"
        ! DO i = 1,axial_num+1
        !     DO j  = 1,radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Th_mesh%th_mesh_area(2,j,i,2)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:dr"
        ! DO j = 1, radial_num+1
        !     WRITE(OUTUNIT,*) 1,j,Th_mesh%radial_fine_mesh_size(j)
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:fzq"
        ! DO i = 1,axial_num+1
        !     DO j  = 1,radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Th_mesh%th_mesh_area(2,j,i,3)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:dz"
        ! DO j = 1, axial_num+1
        !     WRITE(OUTUNIT,*) 1,j,Th_mesh%axial_fine_mesh_size(j)
        ! END DO
        ! WRITE(OUTUNIT,*)"xudongyu:sumw"
        ! DO i = 1,SIZE(material_channel2)
        !     DO j  = 1,radial_num+1
        !         WRITE(OUTUNIT,*) i,j,sum_w(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)
        

    END SUBROUTINE cal_Flow_Friction_Factor
    !> @brief A subroutine that calculate pressure and flow fields for channel type 5
    !! 
    !> @note 
    !> @date 2023-10-26
    SUBROUTINE cal_Flow_Mesh_Channel(numerator,denominator,rogp,xkp,i,j,k,direction)
        REAL(TH_KDUBLE), INTENT(INOUT) :: numerator,denominator,rogp,xkp
        INTEGER, INTENT(IN) :: i,j,k,direction
        INTEGER :: temp_id,temp_i,n
        !> the left mesh of this mesh 
        IF (direction .EQ. 1) THEN
            numerator = numerator + P_fluid(k,j-1,i)/Radial_w(k,j-1,i)
            denominator = denominator + 1.0/Radial_w(k,j-1,i)
            xkp = Radial_w(k,j-1,i)
            IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j-1,i)
            ! WRITE(*,*) Radial_w(k,j-1,i),P_fluid(k,j-1,i),numerator,denominator,direction
        !> the upper mesh of this mesh
        ELSE IF (direction .EQ. 2) THEN
            temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i-1)
            IF (temp_id .EQ. 0) THEN
                numerator = numerator + (P_fluid(k,j,i-1)+rogg(k,j,i-1))/Axial_w(k,j,i-1)
                denominator = denominator + 1.0/Axial_w(k,j,i-1)
                xkp = Axial_w(k,j,i-1)
                IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j,i-1) + rogg(k,j,i-1)
                ! WRITE(*,*) Axial_w(k,j,i-1),P_fluid(k,j,i-1),rogg(k,j,i-1),numerator,denominator,direction
            ELSE
                IF (Channels(temp_id)%channel_type .EQ. 2) THEN
                    DO n = 1, SIZE(material_channel2)
                        IF (temp_id .EQ. material_channel2(n)) THEN
                            temp_i = upper_channel2(k,j,n)-1
                            numerator = numerator + (P_fluid(k,j,temp_i)+sum_rogg(k,j,n)+rogg(k,j,temp_i))/(sum_w(k,j,n)+Axial_w(k,j,temp_i))
                            denominator = denominator + 1.0/(sum_w(k,j,n)+Axial_w(k,j,temp_i))
                            xkp = sum_w(k,j,n)+Axial_w(k,j,temp_i)
                            IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j,temp_i)+sum_rogg(k,j,n)+rogg(k,j,temp_i)
                            ! WRITE(*,*) denominator,sum_w(k,j,n),Axial_w(k,j,temp_i),P_fluid(k,j,temp_i),direction
                            ! WRITE(*,*) numerator,sum_rogg(k,j,n),rogg(k,j,temp_i),temp_i
                        ELSE
                            CYCLE
                        END IF
                    END DO
                ELSE
                    numerator = numerator + (P_fluid(k,j,i-1)+rogg(k,j,i-1))/Axial_w(k,j,i-1)
                    denominator = denominator + 1.0/Axial_w(k,j,i-1)
                    xkp = Axial_w(k,j,i-1)
                    IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j,i-1) + rogg(k,j,i-1)
                    ! WRITE(*,*) Axial_w(k,j,i-1),P_fluid(k,j,i-1),numerator,denominator,direction
                END IF
            END IF
        !> the right mesh of this mesh
        ELSE IF (direction .EQ. 3) THEN
            numerator = numerator + P_fluid(k,j+1,i)/Radial_w(k,j,i)
            denominator = denominator + 1.0/Radial_w(k,j,i)
            xkp = Radial_w(k,j,i)
            IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j+1,i)
            ! WRITE(*,*) Radial_w(k,j,i),P_fluid(k,j+1,i),numerator,denominator,direction
        !> the lower mesh of this mesh
        ELSE IF (direction .EQ. 4) THEN
            temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i+1)
            IF (temp_id .EQ. 0) THEN
                numerator = numerator + (P_fluid(k,j,i+1)-rogg(k,j,i))/Axial_w(k,j,i)
                denominator = denominator + 1.0/Axial_w(k,j,i)
                xkp = Axial_w(k,j,i)
                IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j,i+1)-rogg(k,j,i)
                ! WRITE(*,*) Axial_w(k,j,i),P_fluid(k,j,i+1),rogg(k,j,i),numerator,denominator,direction
            ELSE 
                IF (Channels(temp_id)%channel_type .EQ. 2) THEN
                    DO n = 1, SIZE(material_channel2)
                        IF (temp_id .EQ. material_channel2(n)) THEN
                            temp_i = lower_channel2(k,j,n)+1
                            numerator = numerator + (P_fluid(k,j,temp_i)-sum_rogg(k,j,n)-rogg(k,j,i))/(sum_w(k,j,n)+Axial_w(k,j,i))
                            denominator = denominator + 1.0/(sum_w(k,j,n)+Axial_w(k,j,i))
                            xkp = sum_w(k,j,n)+Axial_w(k,j,i)
                            IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j,temp_i)-sum_rogg(k,j,n)-rogg(k,j,i)
                            ! WRITE(*,*) denominator,sum_w(k,j,n),Axial_w(k,j,i),P_fluid(k,j,temp_i),direction
                            ! WRITE(*,*) numerator,sum_rogg(k,j,n),rogg(k,j,i),temp_i
                        ELSE
                            CYCLE
                        END IF
                    END DO 
                ELSE
                    numerator = numerator + (P_fluid(k,j,i+1)-rogg(k,j,i))/Axial_w(k,j,i)
                    denominator = denominator + 1.0/Axial_w(k,j,i)
                    xkp = Axial_w(k,j,i)
                    IF (xkp .LT. TH_REAL_INFINITY) rogp = P_fluid(k,j,i+1)-rogg(k,j,i)
                    ! WRITE(*,*) Axial_w(k,j,i),P_fluid(k,j,i+1),rogg(k,j,i),numerator,denominator,direction
                END IF
            END IF
        END IF

    END SUBROUTINE cal_Flow_Mesh_Channel
    !> @brief A subroutine that calculate pressure and flow fields for channel type 5
    !! 
    !> @note This subroutine has the same functionality as the <b>ELEM3A</b> 
    !> @date 2023-10-26
    SUBROUTINE cal_Flow_Field_Channel5(min_pressure)
        !计算流体压力分布和质量流量
        ! J用判断网格位置；PNN是计算得到的网格压力，在程序外部看来和P没有区别因此优化为内部的变量；IFMK用于判断边界，能不能用I和N来确定；OVERL-Th_max_sor;IFMXY用于控制跳过计算想优化掉；XZSUM没有用上优化掉；XNSUM1没有用上优化掉；ROGG-rogg;STZUK-Source_mass_flow;MZ-Axial_mass_flow;MR-Radial_mass_flow;P-P_fluid;XKR-Radial_w;XKZ-Axial_w;FZQ,DZ,IFB,KOM-TH_MEsh;IOVER-位置有关的变量；IUVER-位置有关的变量；SUMXK-sum_w;SUMRG-sum_rogg;IPAR,JPAR,NPAR-位置有关的变量；XKP-局部变量直接在程序内部定义；ROGP-局部变量直接在程序内部定义;
        REAL(TH_KDUBLE),INTENT(INOUT) :: min_pressure !< The minimum pressure value in the field
        INTEGER :: i,j,k,n,m,num,temp_id,temp_i
        INTEGER :: num_mesh
        INTEGER :: theta_num,radial_num,axial_num
        REAL(TH_KDUBLE) :: pressure,numerator,denominator,mass_flow(4)
        REAL(TH_KDUBLE) :: temp_real
        REAL(TH_KDUBLE),ALLOCATABLE :: rogp(:),xkp(:)

        
        DO num = 1, num_channel5
            i = axial_channel5(num)
            j = radial_channel5(num)
            k = theta_channel5(num)
            num_mesh = num_mesh_channel5(num)
            ! WRITE(*,*) "XUDONGYU"
            ! WRITE(*,*) i,j,k,num_mesh
            IF (ALLOCATED(rogp)) DEALLOCATE(rogp)
            IF (ALLOCATED(xkp)) DEALLOCATE(xkp)
            ALLOCATE(rogp(2*(num_mesh+1)))
            ALLOCATE(xkp(2*(num_mesh+1)))

            pressure = 0.0
            numerator = 0.0
            denominator = 0.0
            mass_flow = 0.0
            rogp = 0.0
            xkp = 0.0
            !> calculate the pressure field
            IF (num_mesh .EQ. 1) THEN
                !> the left mesh of this mesh
                CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(1),xkp(1),i,j,k,1)
                !> the upper mesh of this mesh
                CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2),xkp(2),i,j,k,2)
                !> the right mesh of this mesh
                CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(4),xkp(4),i,j,k,3)
                !> the lower mesh of this mesh
                CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(3),xkp(3),i,j,k,4)
            ELSE
                DO m = 1, num_mesh
                    !> the first mesh of this channel
                    IF (m .EQ. 1) THEN
                        !> the left mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(1),xkp(1),i,j,k,1)
                        !> the upper mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2),xkp(2),i,j,k,2)
                        !> the lower mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(3),xkp(3),i,j,k,4)
                    !> the last mesh of this channel
                    ELSE IF (m .EQ. num_mesh) THEN
                        !> the upper mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2*m),xkp(2*m),i,j+m-1,k,2)
                        !> the right mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2*m+2),xkp(2*m+2),i,j+m-1,k,3)
                        !> the lower mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2*m+1),xkp(2*m+1),i,j+m-1,k,4)
                    !> the middle mesh of this channel
                    ELSE
                        !> the upper mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2*m),xkp(2*m),i,j+m-1,k,2)
                        !> the lower mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2*m+1),xkp(2*m+1),i,j+m-1,k,4)
                    END IF
                END DO
            END IF
            temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
            numerator = numerator + Channels(temp_id)%source_mass_flow * volume_channel5(num)
            pressure = numerator / denominator
            IF (pressure .LT. min_pressure) min_pressure = pressure
            ! WRITE(*,*) numerator,denominator,pressure
            ! WRITE(*,*) Channels(temp_id)%source_mass_flow,volume_channel5(num)
            DO m = 1,num_mesh
                P_fluid(k,j+m-1,i) = pressure
            END DO

            !> calculate the mass field
            temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
            !> the only one mesh of this channel
            IF (num_mesh .EQ. 1) THEN
                mass_flow = 0.0
                pressure = P_fluid(k,j,i)
                IF (xkp(1) .LT. TH_REAL_INFINITY) mass_flow(1) = (rogp(1)-pressure)/xkp(1)
                IF (xkp(2) .LT. TH_REAL_INFINITY) mass_flow(2) = (rogp(2)-pressure)/xkp(2)
                IF (xkp(3) .LT. TH_REAL_INFINITY) mass_flow(4) = (rogp(3)-pressure)/xkp(3)
                mass_flow(3) = -mass_flow(1)-mass_flow(2)-mass_flow(4)-Channels(temp_id)%source_mass_flow*volume_channel5(num)

                Radial_mass_flow(k,j,i) = HALF * (mass_flow(1)-mass_flow(3))
                Axial_mass_flow(k,j,i) = HALF * (mass_flow(2)-mass_flow(4))
            ELSE
                DO m = 1,num_mesh
                    mass_flow = 0.0
                    pressure = P_fluid(k,j+m-1,i)
                    !> the first mesh of this channel
                    IF (m .EQ. 1) THEN
                        IF (xkp(1) .LT. TH_REAL_INFINITY) mass_flow(1) = (rogp(1)-pressure)/xkp(1)
                        IF (xkp(2) .LT. TH_REAL_INFINITY) mass_flow(2) = (rogp(2)-pressure)/xkp(2)
                        IF (xkp(3) .LT. TH_REAL_INFINITY) mass_flow(4) = (rogp(3)-pressure)/xkp(3)
                        mass_flow(3) = -mass_flow(1)-mass_flow(2)-mass_flow(4)-Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j+m-1,i)
                        temp_real = -mass_flow(3)

                        Radial_mass_flow(k,j,i) = HALF * (mass_flow(1)-mass_flow(3))
                        Axial_mass_flow(k,j,i) = HALF * (mass_flow(2)-mass_flow(4))

                    !> the last mesh of this channel
                    ELSE IF (m .EQ. num_mesh) THEN
                        mass_flow(1) = temp_real
                        IF (xkp(2*m) .LT. TH_REAL_INFINITY) mass_flow(2) = (rogp(2*m)-pressure)/xkp(2*m)
                        IF (xkp(2*m+1) .LT. TH_REAL_INFINITY) mass_flow(4) = (rogp(2*m+1)-pressure)/xkp(2*m+1)
                        mass_flow(3) = -mass_flow(1)-mass_flow(2)-mass_flow(4)-Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j+m-1,i)
                        temp_real = -mass_flow(3)

                        Radial_mass_flow(k,j+m-1,i) = HALF * (mass_flow(1)-mass_flow(3))
                        Axial_mass_flow(k,j+m-1,i) = HALF * (mass_flow(2)-mass_flow(4))

                    !> the middle mesh of this channel
                    ELSE
                        mass_flow(1) = temp_real
                        IF (xkp(2*m) .LT. TH_REAL_INFINITY) mass_flow(2) = (rogp(2*m)-pressure)/xkp(2*m)
                        IF (xkp(2*m+1) .LT. TH_REAL_INFINITY) mass_flow(4) = (rogp(2*m+1)-pressure)/xkp(2*m+1)
                        mass_flow(3) = -mass_flow(1)-mass_flow(2)-mass_flow(4)-Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j+m-1,i)
                        temp_real = -mass_flow(3)

                        Radial_mass_flow(k,j+m-1,i) = HALF * (mass_flow(1)-mass_flow(3))
                        Axial_mass_flow(k,j+m-1,i) = HALF * (mass_flow(2)-mass_flow(4))
                    END IF
                END DO
            END IF
        END DO


        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size)-1 !< 2
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size)-1 !< 29
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size)-1 !< 50
        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: channel 5 presssure"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: channel 5 radial_mass_flow"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Radial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: channel 5 axial_mass_flow"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Axial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)
    END SUBROUTINE cal_Flow_Field_Channel5

    !> @brief A subroutine that calculate pressure and flow fields
    !! 
    !> @date 
    SUBROUTINE cal_Flow_Field_Channel1(min_pressure)
        !计算流体压力分布和质量流量
        ! IFXKK用于控制跳过计算，想把它优化掉；P-P_fluid;PP是计算得到的压力；KX-w;rog2/4-rogg;strom没有什么用；MX-Radial_mass_flow;MY-Axial_mass_flow;STZU-Source_mass_flow x V;IFMXY计算控制参数
        REAL(TH_KDUBLE),INTENT(INOUT) :: min_pressure !< The minimum pressure value in the field
        INTEGER :: i,j,k,n,m,num,temp_id,temp_i
        INTEGER :: num_mesh
        INTEGER :: theta_num,radial_num,axial_num
        REAL(TH_KDUBLE) :: pressure,numerator,denominator,mass_flow(4)
        REAL(TH_KDUBLE) :: temp_real
        REAL(TH_KDUBLE) :: rogp(4),xkp(4)

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size)-1 !< 2
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size)-1 !< 29
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size)-1 !< 50

        DO i = 2,axial_num
            DO j = 2,radial_num
                DO k = 2,theta_num
                    numerator = 0.0
                    denominator = 0.0
                    mass_flow = 0.0
                    rogp = 0.0
                    xkp = 0.0
                    temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                    IF (Channels(temp_id)%channel_type .EQ. 1) THEN
                        !> the left mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(1),xkp(1),i,j,k,1)
                        !> the upper mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(2),xkp(2),i,j,k,2)
                        !> the right mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(3),xkp(3),i,j,k,3)
                        !> the lower mesh of this mesh
                        CALL cal_Flow_Mesh_Channel(numerator,denominator,rogp(4),xkp(4),i,j,k,4)
                        !> calculate the pressure field
                        pressure = numerator / denominator
                        IF (pressure .LT. min_pressure) min_pressure = pressure
                        P_fluid(k,j,i) = pressure
                        !> calculate the mass field
                        DO n = 1, 4
                            IF (xkp(n) .LT. TH_REAL_INFINITY) mass_flow(n) = (rogp(n)-pressure)/xkp(n)
                        END DO
                        Radial_mass_flow(k,j,i) = HALF * (mass_flow(1)-mass_flow(3))
                        Axial_mass_flow(k,j,i) = HALF * (mass_flow(2)-mass_flow(4))
                    END IF
                END DO
            END DO
        END DO

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size)-1 !< 2
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size)-1 !< 29
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size)-1 !< 50
        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: channel 1 presssure"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: channel 1 radial_mass_flow"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Radial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: channel 1 axial_mass_flow"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Axial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)
    END SUBROUTINE cal_Flow_Field_Channel1
    !> @brief A subroutine that calculate pressure and flow fields
    !! 
    !> @date 
    SUBROUTINE cal_Flow_Field_Channel2(min_pressure)
        !计算流体压力分布和质量流量,
        ! DIFF是最小压力，我打算放在外面判断；OVREL-Th_max_sor;IFMXY计算控制参数;IOVER-位置有关的变量；IUVER-位置有关的变量;ROGG-rogg;MZ-Axial_mass_flow;mr-Radial_mass_flow;p-P_fluid;XKZ-Axial_w;
        REAL(TH_KDUBLE),INTENT(INOUT) :: min_pressure !< The minimum pressure value in the field
        INTEGER :: i,j,k,n 
        INTEGER :: upper_index,lower_index 
        REAL(TH_KDUBLE) :: upper_pressure,lower_pressure
        INTEGER :: theta_num,radial_num,axial_num
        REAL(TH_KDUBLE) :: pressure,mass_flow(1)
        REAL(TH_KDUBLE) :: rogp(1),xkp(1)

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size)-1 !< 2
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size)-1 !< 29
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size)-1 !< 50

        DO n = 1, SIZE(material_channel2)
            DO j = 2, radial_num
                DO k = 2, theta_num
                    upper_index = upper_channel2(k,j,n)
                    IF (upper_index .EQ. 0) CYCLE
                    lower_index = lower_channel2(k,j,n)
                    upper_pressure = P_fluid(k,j,upper_index-1)
                    lower_pressure = P_fluid(k,j,lower_index+1)
                    rogp(1) = rogg(k,j,upper_index-1)
                    xkp(1) = Axial_w(k,j,upper_index-1)
                    DO i = upper_index,lower_index
                        rogp(1) = rogp(1) + rogg(k,j,i)
                        xkp(1) = xkp(1) + Axial_w(k,j,i)
                    END DO
                    mass_flow(1) = (upper_pressure+rogp(1)-lower_pressure)/xkp(1)
                    DO i = upper_index,lower_index
                        Radial_mass_flow(k,j,i) = 0.0
                        Axial_mass_flow(k,j,i) = mass_flow(1)
                        pressure = upper_pressure + rogg(k,j,i-1) - Axial_w(k,j,i-1)*mass_flow(1)
                        upper_pressure = pressure
                        P_fluid(k,j,i) = pressure
                        IF (pressure .LT. min_pressure) min_pressure = pressure
                    END DO

                END DO
            END DO
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: channel 2 presssure"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu:channel 2 radial_mass_flow"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Radial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: channel 2 axial_mass_flow"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,Axial_mass_flow(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)
    END SUBROUTINE cal_Flow_Field_Channel2

      SUBROUTINE update_Flow_Field(min_pressure,i_iter,i_extrapolation)
        REAL(TH_KDUBLE),INTENT(INOUT) :: min_pressure !< The minimum pressure value in the field
        INTEGER, INTENT(IN) :: i_iter
        INTEGER, INTENT(INOUT) :: i_extrapolation
        INTEGER :: i,j,k,n,m,num,temp_id,temp_i
        INTEGER :: theta_num,radial_num,axial_num

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size)-1 !< 2
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size)-1 !< 29
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size)-1 !< 50

        !> update the pressure field
        DO i = 2, axial_num
            DO j = 2, radial_num
                DO k = 2, theta_num
                    P_fluid(k,j,i) = P_fluid(k,j,i) - min_pressure
                END DO
            END DO
        END DO
        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: i_iter",i_iter
        ! WRITE(OUTUNIT,*) "xudongyu: diff",min_pressure
        ! WRITE(OUTUNIT,*) "xudongyu: diff presssure"
        ! DO i = 1, axial_num+1
        !     DO j = 1, radial_num+1
        !         WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)
        !> define the extrapolation pressure field
        IF (i_extrapolation .EQ. 4) THEN
            DO i = 2, axial_num
                DO j = 2, radial_num
                    DO k = 2, theta_num
                        temp_p_fluid(k,j,i) = P_fluid(k,j,i)
                    END DO
                END DO
            END DO
        END IF
        i_extrapolation = i_extrapolation + 1
        !> the pressure field is updated by external pressure propulsion
        IF (i_iter .GE. 30) THEN
            IF (i_extrapolation .GE. 10) THEN
                i_extrapolation = 0
                DO i = 2, axial_num
                    DO j = 2, radial_num
                        DO k = 2, theta_num
                            temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                            IF (Channels(temp_id)%channel_type .EQ. 0) THEN
                                CYCLE
                            ELSE
                                P_fluid(k,j,i) = P_fluid(k,j,i) + (P_fluid(k,j,i)-temp_p_fluid(k,j,i))*The_extropolation_factor
                            END IF
                        END DO
                    END DO
                END DO
                OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
                ! WRITE(*,*) "xudongyu: vom1",The_extropolation_factor
                ! WRITE(OUTUNIT,*) "xudongyu: strom presssure"
                ! DO i = 1, axial_num+1
                !     DO j = 1, radial_num+1
                !         WRITE(OUTUNIT,*) i,j,P_fluid(2,j,i)
                !     END DO
                ! END DO
                ! WRITE(OUTUNIT,*) "xudongyu: strom"
                ! DO i = 1, axial_num+1
                !     DO j = 1, radial_num+1
                !         WRITE(OUTUNIT,*) i,j,temp_p_fluid(2,j,i)
                !     END DO
                ! END DO
                CLOSE(OUTUNIT)
            END IF
        END IF
    END SUBROUTINE update_Flow_Field
    !> @brief A subroutine that calculate pressure and flow fields
    !! 
    !> @date 
    SUBROUTINE is_Flow_Field_Coverges(is_coverged,i_iter)
        ! 用于版判断计算是否收敛
        LOGICAL, INTENT(INOUT) :: is_coverged
        INTEGER, INTENT(IN) :: i_iter
        INTEGER :: i,j,k,temp_id,temp_i
        INTEGER :: theta_num,radial_num,axial_num
        INTEGER, ALLOCATABLE :: ifqrow(:,:)
        REAL(TH_KDUBLE) :: smax,srowm,sdiff,tdiff,srow,ifqr,srwov,rowstr,asr,difrow,smaxa,ira,srowv
        REAL(TH_KDUBLE) :: temp_real
        REAL(TH_KDUBLE), ALLOCATABLE :: strbez(:,:) 

        smax = 0.0
        srowm = 0.0
        sdiff = 0.0
        tdiff = 0.0
        smaxa = 1.0

        theta_num = SIZE(Th_mesh%circumferential_fine_mesh_size)-1 !< 2
        radial_num = SIZE(Th_mesh%radial_fine_mesh_size)-1 !< 29
        axial_num = SIZE(Th_mesh%axial_fine_mesh_size)-1 !< 50

        ALLOCATE(ifqrow(theta_num,axial_num))
        ALLOCATE(strbez(theta_num,axial_num))
        ifqrow = 1
        strbez = 0.0

        DO k = 2,theta_num
            temp_real = 0.0
            DO i = 2,axial_num
                DO j = 2,radial_num
                    temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                    IF (temp_id .EQ. 0) THEN
                        CYCLE
                    ELSE
                        IF (Channels(temp_id)%source_mass_flow .GT. 0.0) THEN
                            ifqrow(k,i) = 2
                            temp_real = temp_real + Channels(temp_id)%source_mass_flow*Th_mesh%th_mesh_volume(k,j,i)
                        END IF
                    END IF
                END DO
                strbez(k,i) = temp_real
            END DO
        END DO

        DO k = 2, theta_num
            DO i = 2, axial_num
                srow = 0.0
                ifqr = ifqrow(k,i)
                IF (ifqr .EQ. 1) THEN
                    srowv = strbez(k,i)
                    rowstr = 0.0
                    DO j = 2, radial_num
                        srow = srow + Axial_mass_flow(k,j,i)
                        asr = ABS(srow)
                        IF (asr .GT. rowstr) THEN
                            rowstr = asr
                        END IF
                    END DO
                    IF (rowstr .GT. smax) THEN
                        smax = rowstr
                    END IF
                    difrow = ABS((srowv-srow)/smaxa)
                    IF (difrow .GT. srowm) THEN
                        srowm = difrow
                        ira = i
                    ELSE
                        CYCLE
                    END IF
                ELSE
                    CYCLE
                END IF
            END DO

            sdiff = ABS(smaxa-smax)
            tdiff = sdiff/smax
            sdiff = sdiff/smaxa
            smaxa = smax

            IF(i_iter .GE. Th_max_iter_gas_mass_flow) THEN
                is_coverged = .TRUE.
            ELSE IF (srowm .LT. Th_epsi_gas_mass_flow) THEN
                IF (i_iter .GT. 3) THEN
                    is_coverged = .TRUE.
                ELSE IF (sdiff .GT. 1 .OR. tdiff .GT. 1) THEN
                    is_coverged = .TRUE.
                END IF
            END IF
        END DO

    END SUBROUTINE is_Flow_Field_Coverges


END MODULE HeliumFlowSolver