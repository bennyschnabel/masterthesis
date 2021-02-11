PROGRAM main
        !-------
        ! Calculate Mean Intercept Length
        !-------
        
        USE aux_routines
        USE ISO_FORTRAN_ENV
        USE calculate_mil

        IMPLICIT NONE

        ! Parameter
        REAL(KIND = real64), PARAMETER :: pi = 3.141593
        REAL(KIND = real64), PARAMETER :: r = 1

        ! Internal variables 
        CHARACTER(LEN = 3) :: version
        CHARACTER(LEN = 6) :: time
        CHARACTER(LEN = 8) :: date
        CHARACTER(LEN = 100) :: n2s
        CHARACTER(LEN = 100) :: fileName
        CHARACTER(LEN = 100) :: fileNameExport

        INTEGER(KIND = int32) :: k
        INTEGER(KIND = int64) :: i, j
        INTEGER(KIND = int64) :: noOrientations   
        INTEGER(KIND = int64) :: noRepetitions
        INTEGER(KIND = int64) :: fun1, fun2, fun3, fun4
        INTEGER(KIND = int64) :: status
        INTEGER(KIND = int64), DIMENSION(3) :: dims
        
        REAL(KIND = real64) :: theta
        REAL(KIND = real64) :: phi
        REAL(KIND = real64) :: mil
        REAL(KIND = real64) :: rad2deg
        REAL(KIND = real64) :: start, finish
        REAL(KIND = real64), DIMENSION(3) :: spcng
        REAL(KIND = real64), DIMENSION(3,1) :: n
        REAL(KIND = real64), DIMENSION(3,1)  :: p0, p1
        REAL(KIND = real64), DIMENSION (:,:,:), ALLOCATABLE :: array


        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: arrayData
        INTEGER(KIND = int64), DIMENSION(:,:,:), ALLOCATABLE :: array1

        !----------------------------------------
        ! User input
        !----------------------------------------

        ! File name
        fileName = 'Knochenprobe2_1.2mm_1.vtk'
        ! Number of randomly generated orientation (positiv integer, minimum 9)
        noOrientations = 20
        ! Number of repetitions (positiv integer)
        noRepetitions = 1
        ! VTK import ('mat', 'geb') 'mat'...Matlab, 'geb'... Gebert
        version = 'mat'

        !----------------------------------------
        !
        !----------------------------------------

        CALL CPU_TIME(start)

        fun1 = 5
        fun2 = 10
        fun3 = 15
        fun4 = 20

        WRITE(n2s,*) noOrientations
        

        SELECT CASE(version)
                CASE('mat')
                        CALL getSize(fun3, fileName, dims)
                        ALLOCATE(array1(dims(1), dims(2), dims(3)))
                        CALL importVTK(fun4, fileName, dims, array1)
                CASE('geb')
                        CALL read_vtk(fun1, fileName, array, dims, spcng)
                CASE DEfAULT
        END SELECT
        
        DO i = 1, noRepetitions

                CALL DATE_AND_TIME(date, time)
                fileNameExport = fileName(1:LEN_TRIM(fileName)-4) // '_' // TRIM(ADJUSTL(n2s)) // '_' // date // time // '.dat'

                OPEN(UNIT = fun2, FILE = fileNameExport, IOSTAT = status)
                WRITE(fun2,'(I7)') noOrientations

                DO j = 1, noOrientations
                        ! Generate random value for angle theta and phi
                        k = 12
                        CALL RANDOM_SEED(SIZE = k)
                        CALL RANDOM_NUMBER(theta)
                        theta = theta * pi
                        CALL RANDOM_NUMBER(phi)
                        phi = 2 * phi * pi

                        ! Origin
                        p0 = RESHAPE((/0, 0, 0 /), (/ 3, 1 /))

                        ! Spherical coordinates to Cartesian coordinates
                        p1(1,1) = r * SIN(theta) * COS(phi)
                        p1(2,1) = r * SIN(theta) * SIN(phi)
                        p1(3,1) = r * COS(theta)

                        ! Direction vector
                        n = (1 / SQRT(SUM((p1 - p0)**2))) * (p1 - p0)
                        
                        CALL test(n, dims, array1, mil)
                        WRITE(*,*) 'j: ',j , '/ ', noOrientations, 'theta: ', rad2deg(theta), &
                                'phi: ', rad2deg(phi), 'mil: ', mil
                        WRITE(fun2,'(1f6.4,A1,1f6.4,A1,3f8.4)') theta, ',', phi, ',', mil
                END DO
                CLOSE(fun2)
                
                CALL mil_fabric_tensor(fileNameExport)
        END DO

        CALL CPU_TIME(finish)
        PRINT '("Time = ",f8.3," Minutes")',(finish-start)/60
END PROGRAM main

FUNCTION rad2deg(angleInRadians)
        USE ISO_FORTRAN_ENV

        IMPLICIT NONE

        ! Parameter
        REAL(KIND = real64), PARAMETER :: pi = 3.141593

        ! External variables
        REAL(KIND = real64) :: angleInRadians

        ! Internal variables
        REAL(KIND = real64) :: rad2deg

        rad2deg = (180/pi) * angleInRadians

END FUNCTION rad2deg

SUBROUTINE getSize(fun, fileName, dims)
        USE ISO_FORTRAN_ENV

        IMPLICIT NONE

        ! External variables
        CHARACTER(LEN = 100), INTENT(IN) :: fileName
        INTEGER(KIND = int64), INTENT(IN) :: fun
        INTEGER(KIND = int64), DIMENSION(3), INTENT(OUT) :: dims

        ! Internal variables
        CHARACTER(len = 100) :: dummy
        INTEGER(KIND = int64) :: status

        OPEN(UNIT=fun, FILE=fileName, IOSTAT=status, STATUS='old')

        ! First two lines contain the version number and a title
        READ(fun, '(A)') dummy
        READ(fun, '(A)') dummy
        ! The third line contains ASCII or BINARY
        READ(fun, '(A)') dummy
        ! Fourth line contains the type of dataset
        READ(fun, '(A)') dummy
        ! Dimensions
        READ(fun, '(A10,1X,I10,1X,I10,1X,I10)') dummy, dims

        CLOSE(fun)
END SUBROUTINE getSize
