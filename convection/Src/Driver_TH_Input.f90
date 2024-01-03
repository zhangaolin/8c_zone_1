!> @brief A module that processes input information
!! 
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE DriverTHInput
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation
    USE GlobalTHVariables !< Global variables definition for thermal hydraulic calculation
    USE, INTRINSIC  :: ISO_FORTRAN_ENV !< The module provides access to the ISO Fortran environment and various constant values defined by the Fortran standard

    IMPLICIT NONE
    PUBLIC
    INTEGER, PARAMETER, PRIVATE     :: n_section = 3 !< The number of sections of input file
    INTEGER, PARAMETER, PRIVATE     :: n_keyword = 20 !< The number of keywords of each section
    CHARACTER(LEN=TH_MAX_WORD_LEN)  :: inp_section(n_section) = CHAR_SENTINEL !< The name of each section
    CHARACTER(LEN=TH_MAX_WORD_LEN)  :: th_card(n_keyword) = CHAR_SENTINEL !< The name of each keyword
    CHARACTER(LEN=TH_MAX_WORD_LEN)  :: th_other_card(n_keyword) = CHAR_SENTINEL 
    INTEGER, PARAMETER, PRIVATE  :: max_int_parameter = 3000 !< The maximum number of int variable in each line
    INTEGER, PARAMETER, PRIVATE  :: max_real_parameter = 1000 !< The maximum number of real variable in each line
    INTEGER, PARAMETER, PRIVATE  :: max_char_parameter = 100 !< The maximum number of character variable in each line
    INTEGER, PARAMETER, PRIVATE  :: max_log_parameter = 100 !< The maximum number of logical variable in each line

    INTEGER :: dummy_int(max_int_parameter) !< The int variable in each line
    REAL(TH_KDUBLE) :: dummy_real(max_real_parameter) !< The real variable in each line
    LOGICAL    :: dummy_log(max_log_parameter) !< The logical variable in each line
    CHARACTER(len=TH_MAX_WORD_LEN) :: dummy_char(max_char_parameter) !< The character variable in each line

    CHARACTER(LEN=TH_MAX_WORD_LEN)  :: words(TH_MAX_WORDS) !< The words of each line
    CHARACTER(LEN=TH_MAX_WORD_LEN)  :: keyword !< The key word in the input file
    CHARACTER(LEN=TH_MAX_LINE_LEN)  :: aline !< The line of the input file
    INTEGER  :: nword !< The number of words in each line
CONTAINS
    !> @brief A subroutine that set the section keyword information
    !! 
    SUBROUTINE set_Section_Keyword()
        inp_section(1:n_section) =   ['CASENAME:        ',    & !< The input of case name
                                  &   'TH:              ',    & !< The input of thermal calculation variable
                                  &   'TH_OTHER:        '] !< The input of other variable
    END SUBROUTINE
    !> @brief A subroutine that set the card keyword information
    !! 
    SUBROUTINE set_Card_Keyword()
        th_card(1:n_keyword) =   ["Th_control        ",    &
                            &     "Th_sor            ",    &
                            &     "Gas_physical_properties",    &
                            &     "Mesh_num          ",    &
                            &     "Axial_mesh_size   ",    &
                            &     "Axial_fine_mesh_num",    &
                            &     "Radial_mesh_size  ",    &
                            &     "Radial_fine_mesh_num",    &
                            &     "Circumferential_mesh_size",    &
                            &     "Circumferential_fine_mesh_num",    &
                            &     "Composition_num   ",    &
                            &     "Composition       ",    &
                            &     "Channel           ",    &
                            &     "N1_confg          ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  "]
        th_other_card(1:n_keyword) =   [     "T_soild           ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  ",    &
                            &     "                  "]
    END SUBROUTINE set_Card_Keyword
    !> @brief A function that determine the identifier in the input card
    !! 
    FUNCTION is_Keyword(list,key)  result(is_true)
        CHARACTER(TH_MAX_WORD_LEN), INTENT(IN)  :: list(:) !< The list of strings to be compared
        CHARACTER(TH_MAX_WORD_LEN), INTENT(IN)  :: key !< The string to be compared
        LOGICAL :: is_true !< The logical variable used to record whether a match is made
        INTEGER :: i !< The loop variable

        is_true = .FALSE.
        DO i = 1, SIZE(list)
            IF (TRIM(ADJUSTL(list(i))) == TRIM(ADJUSTL(key))) THEN
                is_true = .TRUE.
                EXIT
            END IF
        END DO
    END FUNCTION is_Keyword
    !> @brief A subroutine that split the string for each line by space or other identifier
    !! 
    SUBROUTINE split_String(string,awords,n,sub)
        CHARACTER(LEN=*), INTENT(IN)  :: string !< The character string to be split
        CHARACTER(LEN=*), INTENT(INOUT)  :: awords(TH_MAX_WORDS) !< The split character array
        INTEGER, INTENT(OUT)  :: n !< The number of characters after splitting
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL  :: sub !< An identifier used to split a string. The default is a space
        CHARACTER(LEN=LEN(string))  :: new_string !< The string after removing the leading space
        CHARACTER(LEN=1)  :: str !< The identifier used to split a string
        INTEGER :: ibeg !< The beginning position of the string to be split
        INTEGER :: iend !< The end position of the string to be split
        INTEGER :: i !< The loop variable

        IF (PRESENT(sub)) THEN
            str = sub
        ELSE
            str = " "
        END IF
        new_string = ADJUSTL(string)
        ibeg = 0
        iend = 0
        n = 0
        DO i = 1, LEN_TRIM(new_string)
            IF ((ibeg .EQ. 0) .AND. (new_string(i:i) /= str)) THEN
                ibeg = i
            END IF
            IF (ibeg .GT. 0) THEN
                IF (new_string(i:i) .EQ. str) THEN
                    iend = i - 1
                END IF
                IF (i .EQ. LEN_TRIM(new_string) .AND. str .EQ. " ") THEN
                    iend = i
                END IF
                IF (iend .GT. 0) THEN 
                    n = n + 1
                    IF (iend - ibeg+1 .GT. len(awords(n))) THEN
                    END IF
                    awords(n) = ADJUSTL(new_string(ibeg:iend))
                    ibeg = 0
                    iend = 0
                END IF
            END IF
        END DO
    END SUBROUTINE split_String
    !> @brief A subroutine that read casename information
    !! 
    SUBROUTINE read_Casename(file_in,file_out)
        INTEGER, INTENT(IN)  :: file_in !< The unit of input file
        INTEGER, INTENT(IN)  :: file_out !< The unit of output file
        INTEGER  :: io_error !< The error code of reading the input file
        CHARACTER(LEN=max_char_parameter)  :: char_tmp !< The character variable in each line

        BACKSPACE(file_in, IOSTAT=io_error) !< The file pointer is moved back one line
        
        READ(unit=file_in, FMT="(A)", IOSTAT=io_error) aline
        CALL split_String(aline,words,nword)
        if (nword .GE. 2) READ(UNIT=aline, FMT=*, IOSTAT=io_error) keyword, char_tmp
        if (nword .GE. 2) CASENAME = char_tmp

    END SUBROUTINE read_Casename
    !> @brief A subroutine that read thermal information
    !! 
    SUBROUTINE read_TH(file_in,file_out)
        INTEGER, INTENT(IN)  :: file_in !< The unit of input file
        INTEGER, INTENT(IN)  :: file_out !< The unit of output file
        CHARACTER(TH_MAX_WORD_LEN) :: aword !< The word in each line
        INTEGER  :: io_error !< The error code of reading the input file
        INTEGER  :: i,j,k !< The loop variable
        CHARACTER(LEN=max_char_parameter)  :: char_tmp !< The character variable in each line
        INTEGER  :: theta_index,radial_index,axial_index !< The number of mesh in each direction
        INTEGER  :: temp_id !< The identifier of the material

        DO 
            READ(UNIT=file_in,FMT='(A)',IOSTAT=io_error) aline
            IF(ALINE(1:TH_MAX_LINE_LEN) .EQ. " ") THEN
                CYCLE
            END IF
            IF (io_error .EQ. IOSTAT_END) THEN
                EXIT
            END IF
            READ(UNIT=aline,FMT=*,IOSTAT=io_error) aword
            IF (is_Keyword(inp_section, aword)) THEN
                BACKSPACE(file_in, IOSTAT=io_error)
                EXIT
            END IF
            IF (is_Keyword(th_card, aword)) THEN
                SELECT CASE(TRIM(ADJUSTL(aword)))
                CASE("Th_control")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1:3),dummy_real(1:4)
                    IF(dummy_int(1) .GT. TH_INT_ZERO)Th_max_iter_gas_tem = dummy_int(1)
                    IF(dummy_int(2) .GT. TH_INT_ZERO)Th_max_iter_gas_mass_flow = dummy_int(2)
                    IF(dummy_int(3) .GT. TH_INT_ZERO)Th_max_avg_gas_tem = dummy_int(3)
                    IF(dummy_real(1) .GT. TH_REAL_ZERO)Th_epsi_gas_tem = dummy_real(1)
                    IF(dummy_real(2) .GT. TH_REAL_ZERO)Th_epsi_gas_mass_flow = dummy_real(2)
                    IF(dummy_real(3) .GT. TH_REAL_ZERO)Th_epsi_avg_gas_tem = dummy_real(3)
                    IF(dummy_real(4) .GT. TH_REAL_ZERO)The_extropolation_factor = dummy_real(4)

                CASE("Th_sor")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_real(1:2),dummy_int(1)
                    IF(dummy_real(1) .GT. TH_REAL_ZERO)Th_max_sor = dummy_real(1)
                    IF(dummy_real(2) .GT. TH_REAL_ZERO)Th_min_sor = dummy_real(2)
                    IF(dummy_int(1) .GT. TH_INT_ZERO) Th_max_sor_change_time = dummy_int(1)

                CASE("Gas_physical_properties")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_real(1:4)
                    Gas_heat_capacity = dummy_real(1)
                    Gas_prandtl_number = dummy_real(2)
                    Gas_pressure = dummy_real(3)
                    Prandtl_number = dummy_real(4)
                    IF (dummy_real(1) .EQ. -1) Gas_heat_capacity = 5195
                    IF (dummy_real(2) .EQ. -1) Gas_prandtl_number = 0.66
                    
                CASE("Mesh_num")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1:3)
                    CALL Th_mesh%alloc(dummy_int(1),dummy_int(2),dummy_int(3)) !< 49 20 1 -> 49 20 1
                    
                CASE("Axial_mesh_size")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_real(1:Th_mesh%axial_mesh_num)
                    Th_mesh%axial_mesh_size(1:Th_mesh%axial_mesh_num) = dummy_real(1:Th_mesh%axial_mesh_num)/TH_REAL_HUNDREAD

                CASE("Axial_fine_mesh_num")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1:Th_mesh%axial_mesh_num)
                    Th_mesh%axial_fine_mesh_num(1:Th_mesh%axial_mesh_num) = dummy_int(1:Th_mesh%axial_mesh_num)

                CASE("Radial_mesh_size")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_real(1:Th_mesh%radial_mesh_num)
                    Th_mesh%radial_mesh_size(1:Th_mesh%radial_mesh_num) = dummy_real(1:Th_mesh%radial_mesh_num)/TH_REAL_HUNDREAD

                CASE("Radial_fine_mesh_num")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1:Th_mesh%radial_mesh_num)
                    Th_mesh%radial_fine_mesh_num(1:Th_mesh%radial_mesh_num) = dummy_int(1:Th_mesh%radial_mesh_num)

                CASE("Circumferential_mesh_size")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_real(1:Th_mesh%circumferential_mesh_num)
                    Th_mesh%circumferential_mesh_size(1:Th_mesh%circumferential_mesh_num) = dummy_real(1:Th_mesh%circumferential_mesh_num)

                CASE("Circumferential_fine_mesh_num")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1:Th_mesh%circumferential_mesh_num)
                    Th_mesh%circumferential_fine_mesh_num(1:Th_mesh%circumferential_mesh_num) = dummy_int(1:Th_mesh%circumferential_mesh_num)

                CASE("Composition_num")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1)
                    ALLOCATE(Compositions(dummy_int(1)))
                    ALLOCATE(Channels(dummy_int(1)))

                CASE("Composition")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1:5),dummy_real(1:11),dummy_int(6),dummy_real(12),dummy_int(7)
                    Compositions(dummy_int(1))%composition_id = dummy_int(1)
                    Compositions(dummy_int(1))%is_heat_exchange = dummy_int(2)
                    Compositions(dummy_int(1))%method_cal_heat_capacity = dummy_int(3)
                    Compositions(dummy_int(1))%method_cal_thermal_conductivity = dummy_int(4)
                    Compositions(dummy_int(1))%direction_radiation = dummy_int(5)
                    Compositions(dummy_int(1))%ntvar = dummy_real(1)
                    Compositions(dummy_int(1))%vol_fraction_soild_material = dummy_real(2)
                    Compositions(dummy_int(1))%heat_capacity = dummy_real(3)
                    Compositions(dummy_int(1))%thermal_conductivity = dummy_real(4)
                    Compositions(dummy_int(1))%lam0 = dummy_real(5)
                    Compositions(dummy_int(1))%coefficient_heat_radiation = dummy_real(6)
                    Compositions(dummy_int(1))%coefficient_heat_radiation_wall = dummy_real(7)
                    Compositions(dummy_int(1))%horizontal_gap = dummy_real(8)
                    Compositions(dummy_int(1))%composition_input_temperature = dummy_real(9)
                    Compositions(dummy_int(1))%composition_power_density = dummy_real(10)
                    Compositions(dummy_int(1))%is_gas_streaming = dummy_int(6)
                    Compositions(dummy_int(1))%diameter_fuel_element = dummy_real(12)/TH_REAL_HUNDREAD
                    Compositions(dummy_int(1))%number_radial_mesh = dummy_int(7)

                CASE("Channel")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(2:3),dummy_real(1:6),dummy_int(4)
                    Channels(dummy_int(1))%is_convective_heat_source = dummy_int(2)
                    Channels(dummy_int(1))%channel_type = dummy_int(3)
                    Channels(dummy_int(1))%pressure_begin = dummy_real(1)
                    Channels(dummy_int(1))%additional_pressure_drop = dummy_real(2)*TH_REAL_HUNDREAD !< /cm -> /m
                    Channels(dummy_int(1))%heat_transition_coefficient = dummy_real(3)
                    Channels(dummy_int(1))%hydraulic_diameter = dummy_real(4)/TH_REAL_HUNDREAD
                    Channels(dummy_int(1))%source_mass_flow = dummy_real(5)
                    Channels(dummy_int(1))%source_temperature_in = dummy_real(6)
                    Channels(dummy_int(1))%is_forced_cooling = dummy_int(4)

                CASE("N1_confg")
                    READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1),dummy_int(2:Th_mesh%radial_mesh_num+1)
                    Th_mesh%th_mesh_material_id(dummy_int(1),1:Th_mesh%radial_mesh_num,1) = dummy_int(2:Th_mesh%radial_mesh_num+1)
                    DO i = 1, Th_mesh%axial_mesh_num-1
                        READ(UNIT=file_in,FMT=*,iostat=io_error) dummy_int(1),dummy_int(2:Th_mesh%radial_mesh_num+1)
                        Th_mesh%th_mesh_material_id(dummy_int(1),1:Th_mesh%radial_mesh_num,i+1) = dummy_int(2:Th_mesh%radial_mesh_num+1)
                    END DO

                END SELECT
            END IF
        END DO

        CALL Th_mesh%init() !< Initialize the mesh information

        theta_index = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_index = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_index = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        ALLOCATE(T_fluid(theta_index,radial_index,axial_index+1))
        T_fluid = TH_REAL_ZERO
        ALLOCATE(T_solid(theta_index,radial_index,axial_index))
        T_solid = TH_REAL_ZERO
        ALLOCATE(Helium_density(theta_index,radial_index,axial_index))
        Helium_density = TH_REAL_ZERO
        ALLOCATE(Axial_mass_flow(theta_index,radial_index,axial_index))
        Axial_mass_flow = TH_REAL_ZERO
        ALLOCATE(Radial_mass_flow(theta_index,radial_index,axial_index))
        Radial_mass_flow = TH_REAL_ZERO
        ALLOCATE(Circumferential_mass_flow(theta_index,radial_index,axial_index))
        Circumferential_mass_flow = TH_REAL_ZERO
        ALLOCATE(P_fluid(theta_index,radial_index,axial_index))
        P_fluid = TH_REAL_ZERO

        DO i = 2, axial_index-1
            DO j = 2, radial_index-1
                DO k = 2, theta_index-1
                    temp_id = Th_mesh%th_fine_mesh_material_id(k,j,i)
                    Th_mesh%th_mesh_type(k,j,i) = Channels(temp_id)%channel_type
                    Th_mesh%th_mesh_porosity(k,j,i) = 1-Compositions(temp_id)%vol_fraction_soild_material
                    IF (Channels(temp_id)%channel_type .EQ. 1 .OR. Channels(temp_id)%channel_type .EQ. 2 .OR. Channels(temp_id)%channel_type .EQ. 5) THEN
                        IF (Channels(temp_id)%pressure_begin .LT. TH_REAL_ZERO) THEN
                            P_fluid(k,j,i) = Gas_pressure
                        ELSE
                            P_fluid(k,j,i) = Channels(temp_id)%pressure_begin
                        END IF
                    END IF
                END DO
            END DO
        END DO

        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu:th_mesh_type"
        ! DO i = 1, axial_index
        !     DO j = 1, radial_index
        !         WRITE(OUTUNIT,*) i,j,Th_mesh%th_mesh_type(2,j,i)
        !     END DO
        ! END DO
        CLOSE(OUTUNIT)

    END SUBROUTINE read_TH
    !> @brief A subroutine that read other information
    !! 
    SUBROUTINE read_THOther(file_in,file_out)
        INTEGER, INTENT(IN)  :: file_in !< The unit of input file
        INTEGER, INTENT(IN)  :: file_out !< The unit of output file
        CHARACTER(TH_MAX_WORD_LEN) :: aword !< The word in each line
        INTEGER  :: io_error !< The error code of reading the input file
        INTEGER  :: i !< The loop variable
        CHARACTER(LEN=max_char_parameter)  :: char_tmp !< The character variable in each line   
        INTEGER  :: theta_index,radial_index,axial_index !< The number of mesh in each direction

        theta_index = SIZE(Th_mesh%circumferential_fine_mesh_size) !< 3
        radial_index = SIZE(Th_mesh%radial_fine_mesh_size) !< 30
        axial_index = SIZE(Th_mesh%axial_fine_mesh_size) !< 51

        DO 
            READ(UNIT=file_in,FMT='(A)',IOSTAT=io_error) aline
            IF(ALINE(1:TH_MAX_LINE_LEN) .EQ. " ") THEN
                CYCLE
            END IF
            IF (io_error .EQ. IOSTAT_END) THEN
                EXIT
            END IF
            READ(UNIT=aline,FMT=*,IOSTAT=io_error) aword
            IF (is_Keyword(inp_section, aword)) THEN
                BACKSPACE(file_in, IOSTAT=io_error)
                EXIT
            END IF
            IF (is_Keyword(th_other_card, aword)) THEN
                SELECT CASE(TRIM(ADJUSTL(aword)))
                    CASE("T_soild")
                        READ(UNIT=aline,FMT=*,iostat=io_error) keyword,dummy_int(1),dummy_real(1:radial_index-1)
                        T_solid(dummy_int(1)+1,1:radial_index-1,1) = dummy_real(1:radial_index-1)
                        DO i = 1, axial_index-2
                            READ(UNIT=file_in,FMT=*,iostat=io_error) dummy_int(1),dummy_real(1:radial_index-1)
                            T_solid(dummy_int(1)+1,1:radial_index-1,i+1) = dummy_real(1:radial_index-1)
                        END DO

                END SELECT
            END IF
        END DO
        DO i = 1, axial_index
            T_solid(1,1:radial_index,i) = T_solid(2,1:radial_index,i)
            T_solid(theta_index,1:radial_index,i) = T_solid(theta_index-1,1:radial_index,i)
        END DO
    END SUBROUTINE read_THOther
    !> @brief A subroutine that drives the reading of input card information
    !! 
    SUBROUTINE driver_Plain_Read()
        CHARACTER(LEN=TH_MAX_WORD_LEN)  :: section_name !< The name of each section
        INTEGER  :: io_error !< The error code of reading the input file

        OPEN(UNIT=INUNIT, FILE=INPUT_FILE)

        DO
            READ(UNIT=INUNIT, FMT='(A)', IOSTAT=io_error) aline
            IF (aline(1:TH_MAX_LINE_LEN) .EQ. " ") THEN 
                CYCLE
            END IF
            IF (io_error .EQ. IOSTAT_END) EXIT
            READ(UNIT=ALINE, FMT=*, IOSTAT=io_error) section_name
            IF (is_Keyword(inp_section, section_name)) THEN
                SELECT CASE(TRIM(ADJUSTL(section_name)))
                CASE("CASENAME:")
                    CALL read_Casename(INUNIT,OUTUNIT)
                CASE("TH:")
                    CALL read_TH(INUNIT,OUTUNIT)
                CASE("TH_OTHER:")
                    CALL read_THOther(INUNIT,OUTUNIT)
                CASE("END:")
                    EXIT
                END SELECT
            END IF
        END DO

        CLOSE(INUNIT)
    END SUBROUTINE driver_Plain_Read
    !> @brief A subroutine that read input information
    !! 
    SUBROUTINE read_TH_Input()
        OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,STATUS="REPLACE",POSITION="APPEND")
        WRITE(OUTUNIT,*) "xudongyu"
        CLOSE(OUTUNIT)
        CALL set_Section_Keyword()
        CALL set_Card_Keyword()
        CALL driver_Plain_Read()
    END SUBROUTINE read_TH_Input
END MODULE DriverTHInput