PROGRAM fabric

	USE ISO_FORTRAN_ENV
	USE numeric_funcs
	USE mkl95_lapack
	USE mkl95_precision

	IMPLICIT NONE
	
	! Parameter
	INTEGER(KIND = int64), PARAMETER :: fun1 = 5
	INTEGER(KIND = int64), PARAMETER :: fun2 = 10
	REAL(KIND = real64), PARAMETER :: abserr = 1.0e-09
	
	! Internal variables
	CHARACTER(LEN = 100) :: dummy
	CHARACTER(LEN = 200) :: fileName
	CHARACTER(LEN = 200) :: fileNameExport
	CHARACTER(LEN = 200) :: formatString1
	CHARACTER(LEN = 200) :: formatString2
	CHARACTER(LEN = 500) :: dummyLong
	
	INTEGER(KIND = int64) :: i, j, k, l
	INTEGER(KIND = int64) :: noLines
	INTEGER(KIND = int64) :: domainNo
	INTEGER(KIND = int32) :: status
	INTEGER(KIND = int32) :: it_num, rot_num
	
	REAL(KIND = real64), DIMENSION (3) :: eigenvalues
	REAL(KIND = real64), DIMENSION (3,3) :: tensorM
	REAL(KIND = real64), DIMENSION (3,3) :: tensorMInv
	REAL(KIND = real64), DIMENSION (3,3) :: tensorMInvR
	REAL(KIND = real64), DIMENSION (3,3) :: eig
	REAL(KIND = real64), DIMENSION (3,3) :: evecs
	REAL(KIND = real64), DIMENSION (3,3) :: D05
	REAL(KIND = real64), DIMENSION (3,3) :: eigenvectors
	REAL(KIND = real64), DIMENSION (3,3) :: tensorH
	
	!----------------------------------------
	! User input
	!----------------------------------------
	
	fileName = 'Knochenprobe_2_M.dat'
	noLines = 312
	
	! dapow
	
	fileNameExport = fileName(1:LEN_TRIM(fileName)-6) // '_H.csv'
	OPEN(UNIT = fun2, FILE = fileNameExport, IOSTAT = status)
	WRITE(fun2,*) 'DomainNo;H11;H12;H13;H21;H22;H23;H31;H32;H33'
	CLOSE(fun2)
	
	OPEN(UNIT = fun1, FILE = fileName, IOSTAT = status, status = 'old')
	
	formatString1 = '(I6, 9(A1, E15.7))'
	formatString2 = '(I6, 9(A1, F11.5))'
	
	DO i = 1, noLines
		READ(fun1, formatString1, IOSTAT = status) domainNo, dummy, tensorM(1,1), dummy, tensorM(1,2), dummy, tensorM(1,3), &
		dummy, tensorM(2,1), dummy, tensorM(2,2), dummy, tensorM(2,3), dummy, tensorM(3,1), dummy, tensorM(3,2), dummy, tensorM(3,3)
		
		WRITE(*,*) 'MIL tensor M:'
		DO j = 1, SIZE(tensorM, 1)
			WRITE(*,'(20G12.4)') tensorM(j,:)
		END DO
		
		WRITE(*,*) 'Inverse'
		! Calculate fabric tensor F
                CALL inverse(tensorM, tensorMInv, 3_int64)
                
                ! Important to get a 'symmetric matric' for jacobi algorithm
		
		!eig =  nint(tensorMInv * 1000.0) / 1000.0
		
                !CALL jacobi(eig, evecs, abserr, 3_int64)
                
                tensorMInvR =  nint(tensorMInv * 1000.0) / 1000.0
                !WRITE(*,*) 'Jacobi1'
                !CALL Jacobi1(tensorMInvR,3,eigenvalues,eigenvectors, 25)
               ! CALL jacobi_eigenvalue ( 3, tensorMInvR, 100, eigenvectors, eigenvalues, it_num, rot_num )
                !WRITE(*,*) eigenvalues
		
		CALL find_eigens(eigval, eigvec, nx, nev, which)
		
		DO k = 1, 3
			DO l = 1, 3
				D05(k, l) = 0
				IF (k == l) THEN
					!D05(k, l) = sqrt(eig(k, l))
		!			D05(k, l) = sqrt(eigenvalues(l))
				END IF
			END DO
		END DO
		
		!tensorH = MATMUL(MATMUL(evecs, D05), TRANSPOSE(evecs))
		!tensorH = MATMUL(MATMUL(eigenvectors, D05), TRANSPOSE(eigenvectors))
		
		!WRITE(*,*) 'Fabric tensor H:'
		!DO j = 1, SIZE(tensorH, 1)
		!	WRITE(*,'(20G12.4)') tensorH(j,:)
		!END DO
		
		!OPEN(UNIT = fun2, FILE = fileNameExport, IOSTAT = status, ACCESS = 'append')
		!WRITE(fun2, formatString2) domainNo, ';', tensorH(1,1), ';', tensorH(1,2), ';', tensorH(1,3), ';', tensorH(2,1), &
		!';', tensorH(2,2), ';', tensorH(2,3), ';', tensorH(3,1), ';', tensorH(3,2), ';', tensorH(3,3)
		!CLOSE(fun2)
	END DO
	CLOSE(fun1)

END PROGRAM fabric
