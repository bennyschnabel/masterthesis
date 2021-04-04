PROGRAM main
	!-------
	! Calculate the Mean Intercept Length of a .vtk Volume
	! 
	! Author: Benjamin Schnabel
	!-------

	USE aux_routines
	USE calculate_mil
	USE ISO_FORTRAN_ENV

	IMPLICIT NONE

	! Parameter
	INTEGER(KIND = int64), PARAMETER :: fun1 = 5
	INTEGER(KIND = int64), PARAMETER :: fun2 = 10
	INTEGER(KIND = int64), PARAMETER :: fun3 = 15
	INTEGER(KIND = int64), PARAMETER :: fun4 = 20
	INTEGER(KIND = int64), PARAMETER :: fun5 = 25
	INTEGER(KIND = int64), PARAMETER :: fun6 = 30
	INTEGER(KIND = int64), PARAMETER :: fun7 = 35
	INTEGER(KIND = int64), PARAMETER :: fun8 = 40
	REAL(KIND = real64), PARAMETER :: pi = 3.141593
	REAL(KIND = real64), PARAMETER :: r = 1

	! Internal variables 
	CHARACTER(LEN = 100) :: n2s
	CHARACTER(LEN = 200) :: fileName
	CHARACTER(LEN = 200) :: fileNameExportDat
	CHARACTER(LEN = 200) :: fileNameExportVtk
	CHARACTER(LEN = 200) :: fileNameExportTensor
	CHARACTER(LEN = 200) :: fileNameExportBVTV
	CHARACTER(LEN = 200) :: formatString

	INTEGER(KIND = int32) :: p
	INTEGER(KIND = int64) :: i, j, k, l, m
	INTEGER(KIND = int64) :: noOrientations
	INTEGER(KIND = int64) :: threshold
	INTEGER(KIND = int64) :: status
	INTEGER(KIND = int64) :: dpi
	INTEGER(KIND = int64) :: xD
	INTEGER(KIND = int64) :: nnDmax
	INTEGER(KIND = int64) :: BV
	INTEGER(KIND = int64) :: TV
	INTEGER(KIND = int64), DIMENSION(3) :: nnD
	INTEGER(KIND = int64), DIMENSION(3) :: dims
	INTEGER(KIND = int64), DIMENSION(3) :: dimsSub

	REAL(KIND = real64) :: theta
	REAL(KIND = real64) :: phi
	REAL(KIND = real64) :: mil
	REAL(KIND = real64) :: rad2deg
	REAL(KIND = real64) :: start, finish
	REAL(KIND = real64) :: domainSize
	REAL(KIND = real64) :: delta
	REAL(KIND = real64) :: BVTV
	REAL(KIND = real64), DIMENSION(3) :: spcng
	REAL(KIND = real64), DIMENSION(3,1) :: n
	REAL(KIND = real64), DIMENSION(3,1)  :: p0, p1
	REAL(KIND = real64), DIMENSION (:,:,:), ALLOCATABLE :: array
	REAL(KIND = real64), DIMENSION (:,:,:), ALLOCATABLE :: arraySub

	!NEW 
	INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: arrayH
	INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: idx
	INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: idy
	INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: idz
	INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: indexX1
	INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: indexX2
	INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: indexX3

	!----------------------------------------
	! User input
	!----------------------------------------

	! File name
	fileName = 'Knochenprobe_2.vtk'
	! Number of randomly generated orientation (positiv integer, minimum 10)
	noOrientations = 600
	! DPI of image stack (positiv integer)
	dpi = 1814
	! Domiain size per subvolume [mm] (positiv integer)
	domainSize = 0.6
	! Threshold (positiv integer)
	threshold = 14250

	!----------------------------------------
	! Acutal calculation
	!----------------------------------------

	! Track calculation time
	CALL CPU_TIME(start)
	
	! Titlescreen
	WRITE(*,*) '======================================'
	WRITE(*,*) 'Mean Intercept Length - Calculation'
	WRITE(*,*) ''
	WRITE(*,*) 'Version: 423423'
	WRITE(*,*) 'Author: Benjamin Schnabel'
	WRITE(*,*) '======================================'

	! Selected version for VTK import
	WRITE(*,*) 'File import started'
	CALL read_vtk(fun1, fileName, array, dims, spcng)
	WRITE(*,*) 'File import done'

	WRITE(*,*) 'MAXVAL:', MAXVAL(array)
	WRITE(*,*) 'MINVAL:', MINVAL(array)
	
	! Decomposition of the volume into subvolumes
	WRITE(*,*) 'Dims:', dims
	delta = 25.4 / REAL(DPI)
	WRITE(*,*) 'Delta: ', delta
	xD = NINT(domainSize / delta)
	WRITE(*,*) 'xD: ', xD
	DO i = 1, 3
		nnD(i) = FLOOR(REAL(dims(i)) / REAL(xD))
		WRITE(*,*) 'nnD(i): ', nnD(i)
	END DO
	nnDmax = nnD(1) * nnD(2) * nnD(3)
	WRITE(*,*) 'nnDmax: ', nnDmax
	
	! New dims of subvolume
	dimsSub(1) = xD
	dimsSub(2) = xD
	dimsSub(3) = xD
	
	! Todo
	ALLOCATE(arrayH(xD))
	ALLOCATE(idx(xD))
	ALLOCATE(idy(xD))
	ALLOCATE(idz(xD))
	ALLOCATE(indexX1(2 * nnDmax))
	ALLOCATE(indexX2(2 * nnDmax))
	ALLOCATE(indexX3(2 * nnDmax))
	
	DO i = 1, xD
		arrayH(i) = i
	END DO
	
	l = 1
	DO i = 1, nnD(1)
		idx = arrayH + (i - 1) * xD
		DO j = 1, nnD(2)
			idy = arrayH + (j - 1) * xD
			DO k = 1, nnD(3)
				idz = arrayH + (k - 1) * xD
				indexX1(l) = idx(1)
				indexX1(l+1) = idx(SIZE(idx))
				indexX2(l) = idy(1)
				indexX2(l+1) = idy(SIZE(idy))
				indexX3(l) = idz(1)
				indexX3(l+1) = idz(SIZE(idz))
				l = l + 2
			END DO
		END DO
	END DO
	
	DEALLOCATE(arrayH)
	DEALLOCATE(idx)
	DEALLOCATE(idy)
	DEALLOCATE(idz)
	
	ALLOCATE(arraySub(xD, xD, xD))
	
	! Export fabric tensor H, csv file
	fileNameExportBVTV = fileName(1:LEN_TRIM(fileName)-4) // '_bvtv.csv'
	fileNameExportTensor = fileName(1:LEN_TRIM(fileName)-4) // '_M.dat'
        
        OPEN(UNIT = fun4, FILE = fileNameExportBVTV, IOSTAT = status)
	WRITE(fun4,*) 'DomainNo;BVTV'
	CLOSE(fun4)
	
	! MIL calculation loop
	l = 1
	DO i = 0, (nnDmax - 1)
		! Create file for export 
		WRITE(n2s,*) i
		fileNameExportDat = fileName(1:LEN_TRIM(fileName)-4) // '_' // TRIM(adjustl(n2s)) // '.dat'
		fileNameExportVtk = fileName(1:LEN_TRIM(fileName)-4) // '_' // TRIM(adjustl(n2s)) // '.vtk'
		
		WRITE(*,*) '------------------'
		WRITE(*,'(I3,A1,I3)') i, '/', nnDmax
		
		! Create subvolume of input
		arraySub = array(indexX1(l) : indexX1(l+1), indexX2(l) : indexX2(l+1), indexX3(l) : indexX3(l+1))
		threshold = INT((MAXVAL(arraySub) - MINVAL(arraySub)) * (6.0_real64 / 7.0_real64), KIND = int64)
		CALL thresholding(dimsSub, arraySub, threshold)
		l = l + 2

		!CALL write_vtk(fun3, fileNameExportVtk, arraySub, spcng, dimsSub)

		WRITE(n2s,*) noOrientations
		OPEN(UNIT = fun2, FILE = fileNameExportDat, IOSTAT = status)
		
		DO m = 1, noOrientations
			! Generate random value for angle theta and phi
			p = 12
			CALL RANDOM_SEED(SIZE = p)
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

			CALL test(n, dimsSub, INT(arraySub, KIND = int64), mil)
			formatString = '(I5, A1, I5, A1, I5, A1, I5, A10, F7.3, A7, F7.3, A7, F8.3)'
			WRITE(*,formatString) i, '/', (nnDmax - 1), ':', m, '/', noOrientations, ' - theta: ', rad2deg(theta), &
			', phi: ', rad2deg(phi), ', MIL: ', mil
			WRITE(fun2,'(1f6.4,A1,1f6.4,A1,3f14.4)') theta, ',', phi, ',', mil
		END DO
		
		CLOSE(fun2)
		
		! Bone volume fraction (BV/TV)
		TV = SIZE(arraySub)
		BV = TV - SUM(arraySub, MASK = arraySub.GT.0)
		BVTV = REAL(BV) / REAL(TV)
		
		OPEN(UNIT = fun4, FILE = fileNameExportBVTV, IOSTAT = status, ACCESS = 'append')
		formatString = '(I6, A1, F7.5)'
		WRITE(fun4,formatString) i, ';', BVTV
		CLOSE(fun4)
		
		CALL mil_fabric_tensor(fileNameExportDat,fileNameExportTensor, i, fun5, fun6, fun8, noOrientations)
		WRITE(*,*) '------------------'
	END DO
	
	DEALLOCATE(arraySub)
	DEALLOCATE(indexX1)
	DEALLOCATE(indexX2)
	DEALLOCATE(indexX3)
	
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

SUBROUTINE thresholding(dims, array, threshold)
	USE ISO_FORTRAN_ENV
	
	IMPLICIT NONE
	
	! External variables
	INTEGER(KIND = int64), INTENT(IN) :: threshold
	INTEGER(KIND = int64), DIMENSION(3), INTENT(IN) :: dims
	REAL(KIND = real64), DIMENSION(dims(1),dims(2),dims(3)), INTENT(INOUT) :: array
	
	! Internal variables
	INTEGER(KIND = int64) :: i, j, k
	
	DO i = 1, dims(1)
		DO j = 1, dims(2)
			DO k = 1, dims(3)
				IF (array(i,j,k) >= threshold) THEN
					array(i,j,k) = 1
				ELSE
					array(i,j,k) = 0
				END IF
			END DO
		END DO
	END DO
END SUBROUTINE thresholding
