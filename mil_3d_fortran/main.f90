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
        CHARACTER(LEN = 100) :: fileNameExportTensor
        CHARACTER(LEN = 100) :: formatString

        INTEGER(KIND = int32) :: k
        INTEGER(KIND = int64) :: i, j
        INTEGER(KIND = int64) :: noOrientations   
        INTEGER(KIND = int64) :: noRepetitions
        INTEGER(KIND = int64) :: fun1, fun2, fun3, fun4, fun5, fun6, fun7
        INTEGER(KIND = int64) :: status
        INTEGER(KIND = int64) :: dpi
        INTEGER(KIND = int64) :: xD
        INTEGER(KIND = int64) :: nnDmax
        INTEGER(KIND = int64), DIMENSION(3) :: nnD
        INTEGER(KIND = int64), DIMENSION(3) :: dims
        
        REAL(KIND = real64) :: theta
        REAL(KIND = real64) :: phi
        REAL(KIND = real64) :: mil
        REAL(KIND = real64) :: rad2deg
        REAL(KIND = real64) :: start, finish
        REAL(KIND = real64) :: domainSize
        REAL(KIND = real64) :: delta 
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
        noOrientations = 200
        ! VTK import ('mat', 'geb') 'mat'...Matlab, 'geb'... Gebert
        version = 'mat'
        ! DPI of image stack
        dpi = 1814
        ! Domiain size per subvolume [mm]
        domainSize = 0.6

        ! REMOVE???
        ! Number of repetitions (positiv integer)
        noRepetitions = 1

        !----------------------------------------
        ! Acutal calculation
        !----------------------------------------
        !Test
        
        

        CALL CPU_TIME(start)

        ! Open file values
        fun1 = 5
        fun2 = 10
        fun3 = 15
        fun4 = 20
        fun5 = 25
        fun6 = 30
        fun7 = 35
                
        ! Selected version for VTK import
        SELECT CASE(version)
                CASE('mat')
                        CALL getSize(fun3, fileName, dims)
                        ALLOCATE(array1(dims(1), dims(2), dims(3)))
                        CALL importVTK(fun4, fileName, dims, array1)
                CASE('geb')
                        CALL read_vtk(fun1, fileName, array, dims, spcng)
                CASE DEfAULT
        END SELECT

        ! Decomposition of the volume into subvolumes
        delta = 25.4 / REAL(DPI)
        WRITE(*,*) delta
        xD = NINT(domainSize / delta)
        WRITE(*,*) xD
        DO i = 1, 3
                nnD(i) = NINT(REAL(dims(i)) / REAL(xD))
                WRITE(*,*) nnD(i)
        END DO
        nnDmax = nnd(1) * nnd(2) * nnd(3)
        WRITE(*,*) nnDmax

        ! Export fabric tensor H, csv file
        fileNameExportTensor = fileName(1:LEN_TRIM(fileName)-4) // '_H.csv'
        OPEN(UNIT = fun7, FILE = fileNameExportTensor, IOSTAT = status)
        WRITE(fun7,*) 'DomainNo; H11; H12; H13; H21; H22; H23; H31; H32; H33'
        CLOSE(fun7)

        !MIL calculation loop
        DO i = 1, nnDmax
                WRITE(*,*) '======================================'

                CALL DATE_AND_TIME(date, time)
                WRITE(n2s,*) noOrientations
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
                        formatString = '(I5, A1, I5, A1, I5, A1, I5, A9, F7.3, A7, F7.3, A7, F8.3)'
                        WRITE(*,formatString) i, '/', nnDmax, ':', j, '/', noOrientations, ', theta: ', rad2deg(theta), &
                                ', phi: ', rad2deg(phi), ', MIL: ', mil
                        !WRITE(*,*) 'j: ',j , '/ ', noOrientations, 'theta: ', rad2deg(theta), &
                        !        'phi: ', rad2deg(phi), 'mil: ', mil
                        WRITE(fun2,'(1f6.4,A1,1f6.4,A1,3f8.4)') theta, ',', phi, ',', mil
                END DO
                CLOSE(fun2)
                
                CALL mil_fabric_tensor(fileNameExport,fileNameExportTensor, i, fun5, fun6)
                WRITE(*,*) '======================================'
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
