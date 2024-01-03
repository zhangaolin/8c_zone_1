 !> @brief A module for the mesh of thermal hydraulic calculations
!! 
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE THMeshHeader
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation

    IMPLICIT NONE

    TYPE AxialMesh
    END TYPE AxialMesh

    TYPE RadialMesh
    END TYPE RadialMesh

    TYPE CircumferentialMesh
    END TYPE CircumferentialMesh
    !> @class THMesh
    !> @brief The type of mesh of thermal hydraulic calculations
    !> @date 
    TYPE THMesh
        INTEGER          :: th_mesh_num              !< Integer variable: the numbers of thermial mesh 
        INTEGER          :: axial_mesh_num           !< Integer variable: the numbers of axial mesh 
        INTEGER          :: radial_mesh_num          !< Integer variable: the numbers of radial mesh 
        INTEGER          :: circumferential_mesh_num !< Integer variable: the numbers of circumferential mesh 
        INTEGER,ALLOCATABLE  :: axial_fine_mesh_num(:) !< One-dimensional integer variable array: the numbers of axial fine mesh
        INTEGER,ALLOCATABLE  :: radial_fine_mesh_num(:) !< One-dimensional integer variable array: the numbers of radial fine mesh
        INTEGER,ALLOCATABLE  :: circumferential_fine_mesh_num(:) !< One-dimensional integer variable array: the numbers of circumferential fine mesh
        INTEGER,ALLOCATABLE  :: th_mesh_type(:,:,:) !< Three-dimensional integer variable array: thermial mesh type
        INTEGER,ALLOCATABLE  :: th_mesh_material_id(:,:,:) !< Three-dimensional integer variable array: the material number in the thermial mesh -th_mesh_material_id(theta,radial,axial)
        INTEGER,ALLOCATABLE  :: th_fine_mesh_material_id(:,:,:) !< Three-dimensional integer variable array: the material number in the thermial mesh -th_mesh_material_id(theta,radial,axial)
        REAL(TH_KDUBLE),ALLOCATABLE  :: th_mesh_porosity(:,:,:) !< Three-dimensional integer variable array: the porosity in the thermial mesh 
        REAL(TH_KDUBLE),ALLOCATABLE  :: axial_mesh_size(:)           !< One-dimensional real variable array: axial mesh size
        REAL(TH_KDUBLE),ALLOCATABLE  :: radial_mesh_size(:)          !< One-dimensional real variable array: radial mesh size
        REAL(TH_KDUBLE),ALLOCATABLE  :: circumferential_mesh_size(:) !< One-dimensional real variable array: circumferential mesh size
        REAL(TH_KDUBLE),ALLOCATABLE  :: axial_fine_mesh_size(:) !< One-dimensional real variable array: axial fine mesh size
        REAL(TH_KDUBLE),ALLOCATABLE  :: radial_fine_mesh_size(:) !< One-dimensional real variable array: radial fine mesh size
        REAL(TH_KDUBLE),ALLOCATABLE  :: circumferential_fine_mesh_size(:) !< One-dimensional real variable array: circumferential fine mesh size
        REAL(TH_KDUBLE),ALLOCATABLE  :: sum_axial_mesh(:,:,:) !< Sum of axial mesh dimensions for a particular flow channel
        REAL(TH_KDUBLE),ALLOCATABLE  :: th_mesh_volume(:,:,:)        !< Three-dimensional real variable array: thermial mesh volume
        REAL(TH_KDUBLE),ALLOCATABLE  :: th_mesh_area(:,:,:,:)        !< Four-dimensional real variable array: thermial mesh area 拆开
        
    CONTAINS
        PROCEDURE, PUBLIC  :: alloc      => allocate_size    !< @copydoc THMesh::allocate_size       
        PROCEDURE, PUBLIC  :: init       => initialize       !< @copydoc THMesh::initialize
        PROCEDURE, PUBLIC  :: cal_volume => calculate_volume !< @copydoc THMesh::calculate_volume
        PROCEDURE, PUBLIC  :: cal_area   => calculate_area   !< @copydoc THMesh::calculate_area
        PROCEDURE, PRIVATE :: cre_fmesh  => create_fine_mesh !< @copydoc THMehs::create_fine_mesh
        PROCEDURE, PUBLIC  :: free       => free_memory      !< @copydoc THMesh::free_memory
        ! 计算面积好像应该在计算体积之前，面积完全没必要存这么多数据我感觉，至少轴向上过去面积都是一样的，周向上也是，我想用type这个变量来确定这个网格的类型，类似于ifbr，通过判断这个type来确定计算采用的方式
    END TYPE THMesh

CONTAINS
!> @brief Mesh size allocate for thermal hydraulic computation
    !> @ingroup THMesh
    !> @date 
    SUBROUTINE allocate_size(this, axial_mesh_num, radial_mesh_num, circumferential_mesh_num)
        CLASS(THMesh), INTENT(INOUT)  :: this                     !< THMesh variable: the mesh of thermal hydraulic calculations
        INTEGER,       INTENT(IN)     :: axial_mesh_num           !< Integer variable: the numbers of axial mesh 
        INTEGER,       INTENT(IN)     :: radial_mesh_num          !< Integer variable: the numbers of radial mesh 
        INTEGER,       INTENT(IN)     :: circumferential_mesh_num !< Integer variable: the numbers of circumferential 
        INTEGER                       :: num_face = 6             !< Integer variable: the number of face of a discrete element

        this%axial_mesh_num = axial_mesh_num !< 49
        this%radial_mesh_num = radial_mesh_num !< 20
        this%circumferential_mesh_num = circumferential_mesh_num !< 1

        this%th_mesh_num = axial_mesh_num * radial_mesh_num * circumferential_mesh_num

        ALLOCATE(this%th_mesh_material_id(circumferential_mesh_num,radial_mesh_num,axial_mesh_num))
        this%th_mesh_material_id = 0
        ALLOCATE(this%axial_mesh_size(axial_mesh_num))
        this%axial_mesh_size = 0
        ALLOCATE(this%axial_fine_mesh_num(axial_mesh_num))
        this%axial_fine_mesh_num = 0
        ALLOCATE(this%radial_mesh_size(radial_mesh_num))
        this%radial_mesh_size = 0
        ALLOCATE(this%radial_fine_mesh_num(radial_mesh_num))
        this%radial_fine_mesh_num = 0
        ALLOCATE(this%circumferential_mesh_size(circumferential_mesh_num))
        this%circumferential_mesh_size = 0
        ALLOCATE(this%circumferential_fine_mesh_num(circumferential_mesh_num))
        this%circumferential_fine_mesh_num = 0
        
        
    END SUBROUTINE allocate_size
    !> @brief Mesh initialization for thermal hydraulic computation
    !> @ingroup THMesh
    !> @date 
    SUBROUTINE initialize(this)
        CLASS(THMesh), INTENT(INOUT)  :: this                     !< THMesh variable: the mesh of thermal hydraulic calculations
        
        CALL this%cre_fmesh() !< Create the fine mesh data
        CALL this%cal_area() !< Calculate the area of each mesh
        CALL this%cal_volume() !< Calculate the volume of each mesh
    END SUBROUTINE initialize

    SUBROUTINE calcilate_sum_size(this,num_channel2)
        CLASS(THMesh), INTENT(IN OUT)  :: this !< THMesh variable: the mesh of thermal
        INTEGER, INTENT(IN) :: num_channel2 !< Integer variable: the number of flow channel
        INTEGER  :: theta_index,radial_index,axial_index !< Integer variable: the index of theta, radial and axial
    END SUBROUTINE calcilate_sum_size

    !> @brief Calculate the volume of each mesh
    !> @ingroup THMesh
    !> @date 
    SUBROUTINE calculate_volume(this)
        CLASS(THMesh), INTENT(IN OUT)  :: this !< THMesh variable: the mesh of thermal 
        INTEGER  :: theta_index,radial_index,axial_index !< Integer variable: the index of theta, radial and axial
        REAL(TH_KDUBLE),ALLOCATABLE  :: delta_r(:),middle_r(:) !< One-dimensional real variable array: the difference of radius and the middle radius
        INTEGER :: i,j,k,n,temp_id !< Integer variable: the index of i,j,k,n and temp_id

        theta_index = SIZE(this%circumferential_fine_mesh_size)
        radial_index = SIZE(this%radial_fine_mesh_size)
        axial_index = SIZE(this%axial_fine_mesh_size)

        ALLOCATE(this%th_mesh_volume(theta_index,radial_index,axial_index))
        this%th_mesh_volume  = TH_REAL_ZERO

        DO i = 2, axial_index-1
            DO j = 2, radial_index-1
                DO k = 2, theta_index-1
                    this%th_mesh_volume(k,j,i) = this%th_mesh_area(k,j,i,3)*this%axial_fine_mesh_size(i)
                END DO
            END DO
        END DO
        
    END SUBROUTINE calculate_volume
    !> @brief Calculate the area of each mesh
    !> @ingroup THMesh
    !> @date 
    SUBROUTINE calculate_area(this)
        CLASS(THMesh), INTENT(IN OUT)  :: this !< THMesh variable: the mesh of thermal 
        REAL(TH_KDUBLE),ALLOCATABLE  :: delta_r(:),middle_r(:) !< One-dimensional real variable array: the difference of radius and the middle radius
        INTEGER                       :: num_face = 6             !< Integer variable: the number of face of a discrete element
        INTEGER  :: theta_index,radial_index,axial_index !< Integer variable: the index of theta, radial and axial
        INTEGER :: i,j,k,n,temp_id !< Integer variable: the index of i,j,k,n and temp_id

        theta_index = SIZE(this%circumferential_fine_mesh_size)
        radial_index = SIZE(this%radial_fine_mesh_size)
        axial_index = SIZE(this%axial_fine_mesh_size)

        ALLOCATE(this%th_mesh_area(theta_index,radial_index,axial_index,num_face))
        this%th_mesh_area  = TH_REAL_ZERO
        ALLOCATE(middle_r(radial_index))
        middle_r = TH_REAL_ZERO
        DO j = 2, radial_index-1
            middle_r(j) = HALF * (this%radial_fine_mesh_size(j)+this%radial_fine_mesh_size(j-1)) + middle_r(j-1)

        END DO

        DO i = 2, axial_index-1
            DO j = 2, radial_index-1
                DO k = 2, theta_index-1
                    !< the area in the cricumferential
                    this%th_mesh_area(k,j,i,1) = this%radial_fine_mesh_size(j)*this%axial_fine_mesh_size(i)
                    !< the area in the radial 
                    this%th_mesh_area(k,j,i,2) = 2.0*TH_PI*middle_r(j)*this%axial_fine_mesh_size(i)*this%circumferential_fine_mesh_size(k)/TH_ANGLE_CIRCUMFERENCE
                    ! IF (i .EQ. 4) WRITE(*,*) i,j,this%th_mesh_area(k,j,i,2)
                    ! IF (i .EQ. 4) WRITE(*,*) middle_r(j),this%axial_fine_mesh_size(i)
                    ! IF (i .EQ. 4) WRITE(*,*) this%circumferential_fine_mesh_size(k)
                    !< the area in the axial
                    this%th_mesh_area(k,j,i,3) = 2.0*TH_PI*middle_r(j)*this%radial_fine_mesh_size(j)*this%circumferential_fine_mesh_size(k)/TH_ANGLE_CIRCUMFERENCE
                    !< the area in the cricumferential
                    this%th_mesh_area(k,j,i,4) = 0.0
                    !< the area in the radial
                    this%th_mesh_area(k,j,i,5) = TH_PI*middle_r(j)*(this%axial_fine_mesh_size(i)+this%axial_fine_mesh_size(i+1))*this%circumferential_fine_mesh_size(k)/TH_ANGLE_CIRCUMFERENCE
                    !< the area in the axial
                    ! this%th_mesh_area(k,j,i,6) = 2.0*TH_PI*middle_r(j+1)*this%radial_fine_mesh_size(j+1)*this%circumferential_fine_mesh_size(k)/TH_ANGLE_CIRCUMFERENCE
                    ! this%th_mesh_area(k,j,i,6) = HALF*(this%th_mesh_area(k,j,i,3)+this%th_mesh_area(k,j,i,6))
                END DO
            END DO
        END DO

        DO j = 1, radial_index
            DO k = 1, theta_index
                this%th_mesh_area(k,j,1,5) =TH_PI*middle_r(j)*this%axial_fine_mesh_size(2)
            END DO
        END DO

        DO i = 1, axial_index
            DO j = 1, radial_index - 1
                DO k = 1, theta_index
                    this%th_mesh_area(k,j,i,6) = HALF*(this%th_mesh_area(k,j,i,3)+this%th_mesh_area(k,j+1,i,3))
                END DO
            END DO
        END DO
        this%th_mesh_area(:,:,1,6) = this%th_mesh_area(:,:,2,6)
        this%th_mesh_area(:,:,axial_index,6) = this%th_mesh_area(:,:,axial_index-1,6)

        ! OPEN(UNIT=OUTUNIT,FILE=OUTPUT_FILE,POSITION="APPEND")
        ! WRITE(OUTUNIT,*) "xudongyu: frq1 "
        ! DO i = 1, axial_index
        !     DO j = 1, radial_index
        !         WRITE(OUTUNIT,*) i,j,this%th_mesh_area(2,j,i,5)
        !     END DO
        ! END DO
        ! WRITE(OUTUNIT,*) "xudongyu: fzq1 "
        ! DO i = 1, axial_index
        !     DO j = 1, radial_index
        !         WRITE(OUTUNIT,*) i,j,this%th_mesh_area(2,j,i,6)
        !     END DO
        ! END DO
        ! CLOSE(OUTUNIT)

    END SUBROUTINE calculate_area
    !> @brief Create the fine mesh data
    !> @ingroup THMesh
    !> @date 
    SUBROUTINE create_fine_mesh(this)
        CLASS(THMesh), INTENT(IN OUT)  :: this !< THMesh variable: the mesh of thermal 
        INTEGER  :: i,j,k,ia,ja,ka !< Integer variable: the loop variable
        INTEGER  :: theta_index,radial_index,axial_index !< Integer variable: the index of theta, radial and axial
        REAL(TH_KDUBLE)  :: diifer_raidial_area,radius,radius1,radius2 !< Real variable: the difference of raidial area and the radius

        i = 0
        j = 1
        k = 1
        DO i = 1, this%axial_mesh_num
            j = j + this%axial_fine_mesh_num(i)
        END DO
        ALLOCATE(this%axial_fine_mesh_size(j+1)) !< 49 -> 51 2-50 = 49
        this%axial_fine_mesh_size(1) = TH_REAL_ZERO
        this%axial_fine_mesh_size(j+1) = TH_REAL_ZERO
        DO i = 1, this%axial_mesh_num
            DO j = 1, this%axial_fine_mesh_num(i)
                k = k + 1
                this%axial_fine_mesh_size(k) = this%axial_mesh_size(i)/ABS(this%axial_fine_mesh_num(i))
            END DO
        END DO

        j = 1
        K = 1
        diifer_raidial_area = 0.0
        radius = 0.0
        radius1 = 0.0
        radius2 = 0.0
        DO i = 1, this%radial_mesh_num
            j = j + this%radial_fine_mesh_num(i)
        END DO
        ALLOCATE(this%radial_fine_mesh_size(j+1)) !< 28 -> 30 2-29 = 28
        this%radial_fine_mesh_size(1) = TH_REAL_ZERO
        this%radial_fine_mesh_size(j+1) = TH_REAL_ZERO 
        DO i = 1, this%radial_mesh_num
            IF (i .EQ. 1) THEN
                diifer_raidial_area = this%radial_mesh_size(i)**2/ABS(this%radial_fine_mesh_num(i))
                radius = this%radial_mesh_size(i)
            ELSE
                radius = radius+this%radial_mesh_size(i)
                diifer_raidial_area = 2*radius*this%radial_mesh_size(i)-this%radial_mesh_size(i)**2 !< r^2-(r-b)^2=2rb-b^2
                ! diifer_raidial_area = this%radial_mesh_size(i)**2-this%radial_mesh_size(i-1)**2
                diifer_raidial_area = diifer_raidial_area/ABS(this%radial_fine_mesh_num(i))
            END IF
            DO j = 1, this%radial_fine_mesh_num(i)
                k = k + 1
                ! this%radial_fine_mesh_size(k) = this%radial_mesh_size(i)/ABS(this%radial_fine_mesh_num(i))
                radius2 = SQRT(radius1**2+diifer_raidial_area)
                this%radial_fine_mesh_size(k) = radius2 - radius1
                radius1 = radius2
            END DO
        END DO

        j = 1
        k = 1
        DO i = 1, this%circumferential_mesh_num
            j = j + this%axial_fine_mesh_num(i)
        END DO
        ALLOCATE(this%circumferential_fine_mesh_size(j+1)) !< 1 -> 3 2-2 = 1
        this%circumferential_fine_mesh_size(1) = 0
        this%circumferential_fine_mesh_size(j+1) = 0
        DO i = 1, this%circumferential_mesh_num
            DO j = 1, this%circumferential_fine_mesh_num(i)
                k = k + 1
                this%circumferential_fine_mesh_size(k) = this%circumferential_mesh_size(i)/ABS(this%circumferential_fine_mesh_num(i))
            END DO
        END DO
        
        theta_index = SIZE(this%circumferential_fine_mesh_size) !< 3
        radial_index = SIZE(this%radial_fine_mesh_size) !< 30
        axial_index = SIZE(this%axial_fine_mesh_size) !< 51
        
        ALLOCATE(this%th_fine_mesh_material_id(theta_index,radial_index,axial_index)) !< 3 30 51
        this%th_fine_mesh_material_id(:,:,:) = TH_INT_ZERO
        ! this%th_fine_mesh_material_id(:,radial_index,:) = TH_INT_ZERO
        ! this%th_fine_mesh_material_id(:,:,axial_index) = TH_INT_ZERO

        ALLOCATE(this%th_mesh_type(theta_index,radial_index,axial_index))
        this%th_mesh_type(:,:,:) = TH_INT_ZERO
        ! this%th_mesh_type(:,radial_index,:) = TH_INT_ZERO
        ! this%th_mesh_type(:,:,axial_index) = TH_INT_ZERO

        ALLOCATE(this%th_mesh_porosity(theta_index,radial_index,axial_index))
        this%th_mesh_porosity = 0.39

        axial_index = 1
        DO i = 1, this%axial_mesh_num
            DO ia = 1, this%axial_fine_mesh_num(i)
                axial_index = axial_index + 1
                radial_index = 1
                DO j = 1, this%radial_mesh_num
                    DO ja = 1, this%radial_fine_mesh_num(j)
                        radial_index = radial_index + 1
                        theta_index = 1
                        DO k = 1, this%circumferential_mesh_num
                            DO ka = 1, this%circumferential_fine_mesh_num(k)
                                theta_index = theta_index + 1
                                this%th_fine_mesh_material_id(theta_index,radial_index,axial_index) = this%th_mesh_material_id(k,j,i)
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE create_fine_mesh
    !> @brief Free the memory of the thermal hydraulic computing mesh
    !> @ingroup THMesh
    !> @date 
    SUBROUTINE free_memory(this)
        CLASS(THMesh), INTENT(IN OUT)  :: this !< THMesh variable: the mesh of thermal 

        DEALLOCATE(this%axial_mesh_size)
        DEALLOCATE(this%radial_mesh_size)
        DEALLOCATE(this%circumferential_mesh_size)
        DEALLOCATE(this%th_mesh_volume)
        DEALLOCATE(this%th_mesh_area)

    END SUBROUTINE free_memory
END MODULE THMeshHeader