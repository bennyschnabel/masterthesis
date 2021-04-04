MODULE numeric_funcs
	USE ISO_FORTRAN_ENV
	
	IMPLICIT NONE

	CONTAINS

                SUBROUTINE gauss_2(a,b,x,n)
                !===========================================================
                ! Solutions to a system of linear equations A*x=b
                ! Method: Gauss elimination (with scaling and pivoting)
                ! Alex G. (November 2009)
                !-----------------------------------------------------------
                ! input ...
                ! a(n,n) - array of coefficients for matrix A
                ! b(n) - array of the right hand coefficients b
                ! n - number of equations (size of matrix A)
                ! output ...
                ! x(n) - solutions
                ! coments ...
                ! the original arrays a(n,n) and b(n) will be destroyed
                ! during the calculation
                !===========================================================
                IMPLICIT NONE

                INTEGER(KIND = int64) :: n
                REAL(KIND = real64) :: a(n,n), b(n), x(n)
                REAL(KIND = real64) :: s(n)
                REAL(KIND = real64) :: c, pivot, store
                INTEGER(KIND = int64) :: i, j, k, l

                ! step 1: begin forward elimination
                DO k=1, n-1
                        ! step 2: "scaling"
                        ! s(i) will have the largest element from row i
                        DO i=k,n ! loop over rows
                                s(i) = 0.0
                                DO j=k,n ! loop over elements of row i
                                        s(i) = max(s(i),ABS(a(i,j)))
                                END DO
                        END DO
                        ! step 3: "pivoting 1"
                        ! find a row with the largest pivoting element
                        pivot = ABS(a(k,k)/s(k))
                        l=k
                        DO j=k+1,n
                                IF(abs(a(j,k)/s(j)) > pivot) THEN
                                        pivot = abs(a(j,k)/s(j))
                                        l=j
                                END IF
                        END DO
                        ! Check IF the system has a sigular matrix
                        IF(pivot == 0.0) THEN
                                write(*,*) 'The matrix is sigular'
                                RETURN
                        END IF
                        ! step 4: "pivoting 2" interchange rows k and l (IF needed)
                        IF (l /= k) THEN
                                DO j=k,n
                                        store = a(k,j)
                                        a(k,j) = a(l,j)
                                        a(l,j) = store
                                END DO
                                store = b(k)
                                b(k) = b(l)
                                b(l) = store
                        END IF
                        ! step 5: the elimination (after scaling and pivoting)
                        DO i=k+1,n
                                c=a(i,k)/a(k,k)
                                a(i,k) = 0.0
                                b(i)=b(i)- c*b(k)
                                DO j=k+1,n
                                        a(i,j) = a(i,j)-c*a(k,j)
                                END DO
                        END DO
                END DO
                ! step 6: back substiturion
                x(n) = b(n)/a(n,n)
                DO i=n-1,1,-1
                        c=0.0
                        DO j=i+1,n
                                c= c + a(i,j)*x(j)
                        END DO
                        x(i) = (b(i)- c)/a(i,i)
                END DO
                END SUBROUTINE gauss_2

                SUBROUTINE jacobi(a,x,abserr,n)
                !===========================================================
                ! Evaluate eigenvalues and eigenvectors
                ! of a real symmetric matrix a(n,n): a*x = lambda*x 
                ! method: Jacoby method for symmetric matrices 
                ! Alex G. (December 2009)
                !-----------------------------------------------------------
                ! input ...
                ! a(n,n) - array of coefficients for matrix A
                ! n      - number of equations
                ! abserr - abs tolerance [sum of (off-diagonal elements)^2]
                ! output ...
                ! a(i,i) - eigenvalues
                ! x(i,j) - eigenvectors
                ! comments ...
                !===========================================================
                IMPLICIT NONE
                INTEGER(KIND = int64) :: i, j, k, n, o
                REAL(KIND = real64) :: a(n,n),x(n,n)
                REAL(KIND = real64) :: abserr, b2, bar
                REAL(KIND = real64) :: beta, coeff, c, s, cs, sc

                ! initialize x(i,j)=0, x(i,i)=1
                ! *** the array operation x=0.0 is specIFic for Fortran 90/95
                x = 0.0
                DO i=1,n
                        x(i,i) = 1.0
                END DO

                ! find the sum of all off-diagonal elements (squared)
                b2 = 0.0
                DO i=1,n
                        DO j=1,n
                                IF (i.ne.j) b2 = b2 + a(i,j)**2
                        END DO
                END DO

                IF (b2 <= abserr) RETURN

                ! average for off-diagonal elements /2
                bar = 0.5*b2/FLOAT(n*n)
                o = 0

                DO WHILE (b2.gt.abserr)
                        DO i=1,n-1
                                DO j=i+1,n
                                        IF (a(j,i)**2 <= bar) cycle  ! DO not touch small elements
                                        b2 = b2 - 2.0*a(j,i)**2
                                        bar = 0.5*b2/float(n*n)
                                        ! calculate coefficient c and s for Givens matrix
                                        beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
                                        coeff = 0.5*beta/sqrt(1.0+beta**2)
                                        s = sqrt(max(0.5+coeff,0.0))
                                        c = sqrt(max(0.5-coeff,0.0))
                                        ! recalculate rows i and j
                                        DO k=1,n
                                                cs =  c*a(i,k)+s*a(j,k)
                                                sc = -s*a(i,k)+c*a(j,k)
                                                a(i,k) = cs
                                                a(j,k) = sc
                                        END DO
                                        ! new matrix a_{k+1} from a_{k}, and eigenvectors 
                                        DO k=1,n
                                                cs =  c*a(k,i)+s*a(k,j)
                                                sc = -s*a(k,i)+c*a(k,j)
                                                a(k,i) = cs
                                                a(k,j) = sc
                                                cs =  c*x(k,i)+s*x(k,j)
                                                sc = -s*x(k,i)+c*x(k,j)
                                                x(k,i) = cs
                                                x(k,j) = sc
                                        END DO
                                END DO
                        END DO
                        o = o + 1
                        
                        IF (o >= 10000) THEN
                        	b2 = abserr
                        END IF
                        
                        !IF (o >= 50) THEN
                        !	abserr = 1.0e-06
                        !	WRITE(*,*) 'abserr = 1.0e-06'
                	!ELSEIF (o >= 70) THEN
                !		abserr = 1.0e-03
                !		WRITE(*,*) 'abserr = 1.0e-03'
        	!	ELSEIF (o>= 100) THEN
        	!		abserr = 1.0
        	!		WRITE(*,*) 'abserr = 1.0'
                !	END IF
                END DO
                RETURN
                END SUBROUTINE jacobi
                
                SUBROUTINE inverse(a,c,n)
                !============================================================
                ! Inverse matrix
                ! Method: Based on DOolittle LU factorization for Ax=b
                ! Alex G. December 2009
                !-----------------------------------------------------------
                ! input ...
                ! a(n,n) - array of coefficients for matrix A
                ! n      - dimension
                ! output ...
                ! c(n,n) - inverse matrix of A
                ! comments ...
                ! the original matrix a(n,n) will be destroyed 
                ! during the calculation
                !===========================================================
                IMPLICIT NONE
                INTEGER(KIND = int64) :: n
                REAL(KIND = real64) :: a(n,n), c(n,n)
                REAL(KIND = real64) :: L(n,n), U(n,n), b(n), d(n), x(n)
                REAL(KIND = real64) :: coeff
                INTEGER(KIND = int64) :: i, j, k

                ! step 0: initialization for matrices L and U and b
                ! Fortran 90/95 aloows such operations on matrices
                L=0.0
                U=0.0
                b=0.0

                ! step 1: forward elimination
                DO k=1, n-1
                        DO i=k+1,n
                                coeff=a(i,k)/a(k,k)
                                L(i,k) = coeff
                                DO j=k+1,n
                                        a(i,j) = a(i,j)-coeff*a(k,j)
                                END DO
                        END DO
                END DO

                ! Step 2: prepare L and U matrices 
                ! L matrix is a matrix of the elimination coefficient
                ! + the diagonal elements are 1.0
                DO i=1,n
                        L(i,i) = 1.0
                END DO
                ! U matrix is the upper triangular part of A
                DO j=1,n
                        DO i=1,j
                                U(i,j) = a(i,j)
                        END DO
                END DO

                ! Step 3: compute columns of the inverse matrix C
                DO k=1,n
                        b(k)=1.0
                        d(1) = b(1)
                        ! Step 3a: Solve Ld=b using the forward substitution
                        DO i=2,n
                                d(i)=b(i)
                                DO j=1,i-1
                                        d(i) = d(i) - L(i,j)*d(j)
                                END DO
                        END DO
                        ! Step 3b: Solve Ux=d using the back substitution
                        x(n)=d(n)/U(n,n)
                        DO i = n-1,1,-1
                                x(i) = d(i)
                                DO j=n,i+1,-1
                                        x(i)=x(i)-U(i,j)*x(j)
                                END DO
                                x(i) = x(i)/u(i,i)
                        END DO
                        ! Step 3c: fill the solutions x(n) into column k of C
                        DO i=1,n
                                c(i,k) = x(i)
                        END DO
                        b(k)=0.0
                END DO
                END SUBROUTINE inverse
END MODULE numeric_funcs
