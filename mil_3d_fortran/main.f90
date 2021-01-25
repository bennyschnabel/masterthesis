PROGRAM main
	!===========================================================
	! Calculate Mean Intercept Length
	!-----------------------------------------------------------
	! input ...
	! fileName - File with data
	! numberOfOrientations - 	Number of randomly generated 
	!							orientation (positiv integer, 
	!							minimum 9)
	! numberOfRepetitions - Number of repetitions (positiv integer)
	! output ...
	! M - MIL tensor
	!===========================================================
	
	IMPLICIT NONE
	
	INTERFACE
		FUNCTION getSize(fileName) RESULT(dimensions)
			CHARACTER(len=100), INTENT(IN)			:: fileName
			INTEGER, DIMENSION(3)	 				:: dimensions
		END FUNCTION getSize
		
		FUNCTION importVTK(fileName, a, b, c) RESULT(x)
			CHARACTER(len=100), INTENT(IN)			:: fileName
			INTEGER, INTENT(IN)						:: a, b, c
			INTEGER, DIMENSION (:), ALLOCATABLE	:: x  
		END FUNCTION importVTK
	END INTERFACE
	
	CHARACTER(len=6) 						:: time
	CHARACTER(len=8) 						:: date
	CHARACTER(len=20) 						:: fmt
	CHARACTER(len=100)						:: fileName, fileNameExport, testFile
	INTEGER 								:: numberOfOrientations, numberOfRepetitions, increment, ii, jj, status
	INTEGER, DIMENSION(3)					:: dimensions
	INTEGER, DIMENSION (:), ALLOCATABLE		:: fileData
	INTEGER, DIMENSION (:,:,:), ALLOCATABLE	:: voxel
	REAL 									:: theta, phi, thetaR, phiR, calculatemil, mil, start, finish, rad2deg
	REAL, PARAMETER 						:: pi = 3.141593, r = 1
	REAL, DIMENSION(3,1)					:: p0, p1, n
	
	CALL CPU_TIME(start)
	
	! ----------------------
	! User input area
	!
	!-----------------------
	
	! File name
	fileName = 'Knochenprobe2_1mm_1.vtk'
	! Number of randomly generated orientation (positiv integer, minimum 9)
	numberOfOrientations = 1000
	! Number of repetitions (positiv integer)
	numberOfRepetitions = 1
	! Distance between two created lines (positiv integer) 
	! TODO
	increment = 1
	
	
	! ----------------------
	! TODO
	!
	!-----------------------
	
	WRITE(fmt,*) numberOfOrientations
	CALL DATE_AND_TIME(date, time)
	fileNameExport = fileName(1:LEN_TRIM(fileName)-4) // '_' // TRIM(ADJUSTL(fmt)) // '_' // date // time // '.dat'
	
	dimensions = getSize(fileName)
	
	ALLOCATE(fileData(dimensions(1) * dimensions(2) * dimensions(3)))
	fileData = importVTK(fileName, dimensions(1), dimensions(2), dimensions(3))
	voxel = RESHAPE(fileData, (/ dimensions(1), dimensions(2), dimensions(3) /))
	DEALLOCATE(fileData)
	
	DO ii = 1, numberOfRepetitions
	
		WRITE(*,*) 'ii: ', ii, '/', numberOfRepetitions
	
		!WRITE(fmt,*) numberOfOrientations
		CALL DATE_AND_TIME(date, time)
		fileNameExport = fileName(1:LEN_TRIM(fileName)-4) // '_' // TRIM(ADJUSTL(fmt)) // '_' // date // time // '.dat'
			
		OPEN(UNIT = 1, FILE = fileNameExport, IOSTAT = status)
		! Write number of orientations into dat file
		WRITE(1,'(I7)') numberOfOrientations
		
		DO jj = 1, numberOfOrientations
			CALL RANDOM_NUMBER(thetaR)
			theta = thetaR * pi
			CALL RANDOM_NUMBER(phiR)
			phi = phiR * 2 * pi
			! Origin
			p0 = RESHAPE((/0, 0, 0 /), (/ 3, 1 /))
			! Spherical coordinates to Cartesian coordinates
			p1(1,1) = r * SIN(theta) * COS(phi)
			p1(2,1) = r * SIN(theta) * SIN(phi)
			p1(3,1) = r * COS(theta)
			
			n = (1 / SQRT(SUM((p1 - p0)**2))) * (p1 - p0)
			mil = calculatemil(n, REAL(dimensions(1)), REAL(dimensions(2)), REAL(dimensions(3)), voxel)
			WRITE(*,*) 'jj: ',jj , '/ ', numberOfOrientations, 'theta: ', rad2deg(theta), 'phi: ', rad2deg(phi), 'mil: ', mil
			
			WRITE(1,'(1f6.4,A1,1f6.4,A1,3f8.4)') theta, ',', phi, ',', mil
		END DO
		
		! Close export file
		CLOSE(1)

		CALL mil_tensor(fileNameExport)
	
	END DO
	
	DEALLOCATE(voxel)
	
	CALL CPU_TIME(finish)
	PRINT '("Time = ",f8.3," Minutes")',(finish-start)/60
	
END PROGRAM	main

FUNCTION calculatemil(n, a, b, c, voxel)
	IMPLICIT NONE
	
	INTERFACE
		FUNCTION rot3Axis(a, b, phi) RESULT(rot)
			REAL, INTENT(IN)					:: phi
			REAL, DIMENSION(4,4)				:: rot
			REAL, DIMENSION(3,1), INTENT(IN) 	:: a, b
		END FUNCTION rot3Axis
		
		FUNCTION bresenhamAlgorithm(P1,P2,d) RESULT(XYZ)
			INTEGER, INTENT(IN)					:: d
			INTEGER, DIMENSION(3,1), INTENT(IN)	:: P1, P2
			INTEGER, DIMENSION(3*d)				:: XYZ
		END FUNCTION bresenhamAlgorithm
	END INTERFACE
	
	INTEGER									:: d, ii, jj, kk, ll, mm, nn, l
	INTEGER									:: counter, counterL, counterU
	INTEGER, DIMENSION(3,1)					:: ps, pQ1, pQ2, pQ3, pQ4, pO1, pO2, pO3, pO4
	INTEGER, DIMENSION(3,1)					:: pQ12, pQ34, pO12, pO34
	INTEGER, DIMENSION(3,1)					:: q, o, SP, EP
	INTEGER, DIMENSION(:), ALLOCATABLE		:: xQ12, x1Q12, x2Q12, x3Q12
	INTEGER, DIMENSION(:), ALLOCATABLE		:: xQ34, x1Q34, x2Q34, x3Q34
	INTEGER, DIMENSION(:), ALLOCATABLE		:: xO12, x1O12, x2O12, x3O12
	INTEGER, DIMENSION(:), ALLOCATABLE		:: xO34, x1O34, x2O34, x3O34
	INTEGER, DIMENSION(:), ALLOCATABLE		:: xQ, x1Q, x2Q, x3Q
	INTEGER, DIMENSION(:), ALLOCATABLE		:: xO, x1O, x2O, x3O
	INTEGER, DIMENSION(:), ALLOCATABLE		:: x, x1, x2, x3
	INTEGER, DIMENSION(:), ALLOCATABLE		:: xse, x1se, x2se, x3se
	INTEGER, DIMENSION(:), ALLOCATABLE		:: indizes
	
	REAL									:: calculatemil, h, cv, a, b, c, dr, ru
	REAL, PARAMETER							:: pi = 3.141593
	REAL, DIMENSION(3,1)					:: n, m, v0
	REAL, DIMENSION(4,1)					:: pQ1l4, pQ2l4, pQ3l4, pQ4l4
	
	INTEGER, DIMENSION(nint(a), nint(b), nint(c))	:: voxel
	
	h = 0.0
	cv = 0.0
	
	! TODO
	!IF (n(1,1) == 1 .AND. n(2,1) == 0 .AND. n(3,1) == 0) THEN
	!ELSE
	
	
	! Room diagonal
	!dr = SQRT(real(a)**2 + real(b)**2 + real(c)**2)
	dr = SQRT(a**2 + b**2 + c**2)
	!WRITE(*,*) 'dr: ', dr
	! Radius of the sphere (Half room diagonal)
	ru = 0.5 * dr
	!WRITE(*,*) 'ru: ', ru
	! Origin of the sphere (center of the read data)
	m = RESHAPE((/ a/2, b/2, c/2 /), (/ 3, 1 /))
	!WRITE(*,*) 'm: ', m
	! Point on the sphere in direction n
	ps = NINT(m + ru * n)
	!WRITE(*,*) 'ps: ', ps
	! First directional vector (direction 1, 0)
	v0 = RESHAPE((/ 1.0, 0.0, (-n(1,1) / n(3,1)) /), (/ 3, 1 /))
	v0 = (1 / SQRT(SUM(v0**2))) * v0
	!WRITE(*,*) 'v0: ', v0
	
	! Boundary point 1 on the plane Q
	pQ1 = NINT(ps + ru * v0)
	!WRITE(*,*) 'pQ1: ', pQ1
	pQ1l4 = RESHAPE((/ REAL(pQ1(1,1)), REAL(pQ1(2,1)), REAL(pQ1(3,1)), 1.0 /), (/ 4,1 /))
	!WRITE(*,*) 'pQ1l4: ', pQ1l4
	
	! Boundary point 2 on the plane Q
	pQ2l4 = MATMUL(rot3Axis(m, n, pi/2), pQ1l4)
	!WRITE(*,*) 'pQ2l4: ', pQ2l4
	pQ2 = RESHAPE(NINT((/ pQ2l4(1,1), pQ2l4(2,1), pQ2l4(3,1) /)), (/ 3,1 /))
	!WRITE(*,*) 'pQ2: ', pQ2
	
	! Boundary point 3 on the plane Q
	pQ3l4 = MATMUL(rot3Axis(m, n, pi), pQ1l4)
	!WRITE(*,*) 'pQ3l4: ', pQ3l4
	pQ3 = RESHAPE(NINT((/ pQ3l4(1,1), pQ3l4(2,1), pQ3l4(3,1) /)), (/ 3,1 /))
	!WRITE(*,*) 'pQ3: ', pQ3
	
	! Boundary point 4 on the plane Q
	pQ4l4 = MATMUL(rot3Axis(m, n, (3 * pi) / 2), pQ1l4)
	!WRITE(*,*) 'pQ4l4: ', pQ4l4
	pQ4 = RESHAPE(NINT((/ pQ4l4(1,1), pQ4l4(2,1), pQ4l4(3,1) /)), (/ 3,1 /))
	!WRITE(*,*) 'pQ4: ', pQ4
	
	! Boundary point 1 on the plane O
	pO1 = pQ1 - NINT(dr * n)
	!WRITE(*,*) 'pO1: ', pO1
	! Boundary point 2 on the plane O
	pO2 = pQ2 - NINT(dr * n)
	!WRITE(*,*) 'pO2: ', pO2
	! Boundary point 3 on the plane O
	pO3 = pQ3 - NINT(dr * n)
	!WRITE(*,*) 'pO3: ', pO3
	! Boundary point 4 on the plane O
	pO4 = pQ4 - NINT(dr * n)
	!WRITE(*,*) 'pO4: ', pO4
	
	! Creation of planes Q and O
	! 
	d = MAXVAL(ABS(pQ2 - pQ1) + 1)
	ALLOCATE(xQ12(3 * d))
	xQ12 = bresenhamAlgorithm(pQ1,pQ2,d)
	ALLOCATE(x1Q12(d))
	ALLOCATE(x2Q12(d))
	ALLOCATE(x3Q12(d))
	DO ii = 1, d
		x1Q12(ii) = xQ12(ii)
		x2Q12(ii) = xQ12(d + ii)
		x3Q12(ii) = xQ12(2 * d + ii)
	END DO
	DEALLOCATE(xQ12)
	
	! Second
	d = MAXVAL(ABS(pQ4 - pQ3) + 1)
	ALLOCATE(xQ34(3 * d))
	xQ34 = bresenhamAlgorithm(pQ3,pQ4,d)
	ALLOCATE(x1Q34(d))
	ALLOCATE(x2Q34(d))
	ALLOCATE(x3Q34(d))
	DO ii = 1, d
		x1Q34(ii) = xQ34(ii)
		x2Q34(ii) = xQ34(d + ii)
		x3Q34(ii) = xQ34(2 * d + ii)
	END DO
	DEALLOCATE(xQ34)
	
	! Third
	d = MAXVAL(ABS(pO1 - pO2) + 1)
	ALLOCATE(xO12(3 * d))
	xO12 = bresenhamAlgorithm(pO1,pO2,d)
	ALLOCATE(x1O12(d))
	ALLOCATE(x2O12(d))
	ALLOCATE(x3O12(d))
	DO ii = 1, d
		x1O12(ii) = xO12(ii)
		x2O12(ii) = xO12(d + ii)
		x3O12(ii) = xO12(2 * d + ii)
	END DO
	DEALLOCATE(xO12)
	
	! Fourth
	d = MAXVAL(ABS(pO3 - pO4) + 1)
	ALLOCATE(xO34(3 * d))
	xO34 = bresenhamAlgorithm(pO3,pO4,d)
	ALLOCATE(x1O34(d))
	ALLOCATE(x2O34(d))
	ALLOCATE(x3O34(d))
	DO ii = 1, d
		x1O34(ii) = xO34(ii)
		x2O34(ii) = xO34(d + ii)
		x3O34(ii) = xO34(2 * d + ii)
	END DO
	DEALLOCATE(xO34)
	
	l = MIN(SIZE(x1Q12), SIZE(x1Q34), SIZE(x1O12), SIZE(x1O34))
	
	DO jj = 1, (l - 1)
		pQ12 = RESHAPE((/x1Q12(jj), x2Q12(jj), x3Q12(jj) /), (/ 3, 1 /))
		pQ34 = RESHAPE((/x1Q34(l - jj), x2Q34(l - jj), x3Q34(l - jj) /), (/ 3, 1 /))
		pO12 = RESHAPE((/x1O12(jj), x2O12(jj), x3O12(jj) /), (/ 3, 1 /))
		pO34 = RESHAPE((/x1O34(l - jj), x2O34(l - jj), x3O34(l - jj) /), (/ 3, 1 /))
		
		! First
		d = MAXVAL(ABS(pQ12 - pQ34) + 1)
		ALLOCATE(xQ(3 * d))
		xQ = bresenhamAlgorithm(pQ12,pQ34,d)
		ALLOCATE(x1Q(d))
		ALLOCATE(x2Q(d))
		ALLOCATE(x3Q(d))
		DO ii = 1, d
			x1Q(ii) = xQ(ii)
			x2Q(ii) = xQ(d + ii)
			x3Q(ii) = xQ(2 * d + ii)
		END DO
		DEALLOCATE(xQ)
		
		! Second
		d = MAXVAL(ABS(pO12 - pO34) + 1)
		ALLOCATE(xO(3 * d))
		xO = bresenhamAlgorithm(pO12,pO34,d)
		ALLOCATE(x1O(d))
		ALLOCATE(x2O(d))
		ALLOCATE(x3O(d))
		DO ii = 1, d
			x1O(ii) = xO(ii)
			x2O(ii) = xO(d + ii)
			x3O(ii) = xO(2 * d + ii)
		END DO
		DEALLOCATE(xO)
		
		DO kk = 1, SIZE(x1Q)
			q = RESHAPE((/x1Q(kk), x2Q(kk), x3Q(kk) /), (/ 3, 1 /))
			o = RESHAPE((/x1O(kk), x2O(kk), x3O(kk) /), (/ 3, 1 /))
			
			! First
			d = MAXVAL(ABS(q - o) + 1)
			ALLOCATE(x(3 * d))
			x = bresenhamAlgorithm(q,o,d)
			ALLOCATE(x1(d))
			ALLOCATE(x2(d))
			ALLOCATE(x3(d))
			DO ii = 1, d
				x1(ii) = x(ii)
				x2(ii) = x(d + ii)
				x3(ii) = x(2 * d + ii)
			END DO
			DEALLOCATE(x)
			
			counter = 0
			DO ll = 1, SIZE(x1)
				IF ((x1(ll) >= 1 .AND. x1(ll) <= a) .AND. (x2(ll) >= 1 .AND. x2(ll) <= b) .AND. (x3(ll) >=1 .AND. x3(ll) <= c)) THEN
					counter = counter + 1
				END IF				
			END DO
			
			
			IF (counter > 2) THEN
				ALLOCATE(indizes(counter))
				mm = 1
				DO ll = 1, SIZE(x1)
					IF ((x1(ll) >= 1 .AND. x1(ll) <= a) .AND. (x2(ll) >= 1 .AND. x2(ll) <= b) .AND. (x3(ll) >=1 .AND. x3(ll) <= c)) THEN
						indizes(mm) = ll
						mm = mm + 1
					END IF
				END DO
				
				!WRITE(*,*) 'indizes', indizes
				
				counterL = indizes(1)
				counterU = indizes(counter)
				
				SP = RESHAPE((/ x1(counterL), x2(counterL), x3(counterL) /), (/ 3, 1 /))
				
				EP = RESHAPE((/ x1(counterU), x2(counterU), x3(counterU) /), (/ 3, 1 /))
				h = h + SQRT(SUM(REAL(EP - SP)**2))
				
				! Bresenham
				d = MAXVAL(ABS(SP - EP) + 1)
				ALLOCATE(xse(3 * d))
				xse = bresenhamAlgorithm(SP,EP,d)
				ALLOCATE(x1se(d))
				ALLOCATE(x2se(d))
				ALLOCATE(x3se(d))
				DO ii = 1, d
					x1se(ii) = xse(ii)
					x2se(ii) = xse(d + ii)
					x3se(ii) = xse(2 * d + ii)
				END DO
				DEALLOCATE(xse)
				
				! Calculation of the intersection points with the cancellous bone
				DO nn = 1, (SIZE(x1se) - 1)
					IF (voxel(x1se(nn), x2se(nn), x3se(nn)) == 0 .AND. voxel(x1se(nn + 1), x2se(nn + 1), x3se(nn + 1)) ==1) THEN
						cv = cv + 1.0
					END IF
				END DO
				
				DEALLOCATE(x1se)
				DEALLOCATE(x2se)
				DEALLOCATE(x3se)
				DEALLOCATE(indizes)
			END IF
			
			
			DEALLOCATE(x1)
			DEALLOCATE(x2)
			DEALLOCATE(x3)
		END DO
		
		DEALLOCATE(x1Q)
		DEALLOCATE(x2Q)
		DEALLOCATE(x3Q)
		
		DEALLOCATE(x1O)
		DEALLOCATE(x2O)
		DEALLOCATE(x3O)
	END DO
	
	DEALLOCATE(x1Q12)
	DEALLOCATE(x2Q12)
	DEALLOCATE(x3Q12)
	
	DEALLOCATE(x1Q34)
	DEALLOCATE(x2Q34)
	DEALLOCATE(x3Q34)
	
	DEALLOCATE(x1O12)
	DEALLOCATE(x2O12)
	DEALLOCATE(x3O12)
	
	DEALLOCATE(x1O34)
	DEALLOCATE(x2O34)
	DEALLOCATE(x3O34)

	!END IF
		
	calculatemil = h / cv
END FUNCTION calculatemil

FUNCTION rot3Axis(a, b, phi) RESULT(rot)
	IMPLICIT NONE
	
	INTERFACE
		FUNCTION translationMatrixPositiv(t) RESULT(M)
			REAL, DIMENSION(3,1), INTENT(IN)	:: t
			REAL, DIMENSION(4,4) 				:: M
		END FUNCTION translationMatrixPositiv
		
		FUNCTION translationMatrixNegativ(t) RESULT(M)
			REAL, DIMENSION(3,1), INTENT(IN)	:: t
			REAL, DIMENSION(4,4) 				:: M
		END FUNCTION translationMatrixNegativ
		
		FUNCTION rotationMatrixX2Positiv(b, d) RESULT(R)
			REAL, INTENT(IN) 					:: d
			REAL, DIMENSION(3,1), INTENT(IN) 	:: b
			REAL, DIMENSION(4,4) 				:: R
		END FUNCTION rotationMatrixX2Positiv
		
		FUNCTION rotationMatrixX2Negativ(b, d) RESULT(R)
			REAL, INTENT(IN) 					:: d
			REAL, DIMENSION(3,1), INTENT(IN) 	:: b
			REAL, DIMENSION(4,4) 				:: R
		END FUNCTION rotationMatrixX2Negativ
		
		FUNCTION rotationMatrixX3Positiv(b, d) RESULT(R)
			REAL, INTENT(IN) 					:: d
			REAL, DIMENSION(3,1), INTENT(IN) 	:: b
			REAL, DIMENSION(4,4) 				:: R
		END FUNCTION rotationMatrixX3Positiv
		
		FUNCTION rotationMatrixX3Negativ(b, d) RESULT(R)
			REAL, INTENT(IN) 					:: d
			REAL, DIMENSION(3,1), INTENT(IN) 	:: b
			REAL, DIMENSION(4,4) 				:: R
		END FUNCTION rotationMatrixX3Negativ
		
		FUNCTION rotationMatrixX3(alpha) RESULT(R)
			REAL, INTENT(IN) 					:: alpha
			REAL, DIMENSION(4,4) 				:: R
		END FUNCTION rotationMatrixX3
	END INTERFACE
	
	REAL								:: d
	REAL, INTENT(IN)					:: phi
	REAL, DIMENSION(4,4)				:: rot, Tp, Tn, Rx2p, Rx2n, Rx3p, Rx3n, Rx3
	REAL, DIMENSION(3,1), INTENT(IN) 	:: a, b
	
	d = sqrt(b(1,1)**2 + b(2,1)**2)
	
	Tp = translationMatrixPositiv(a)
	Tn = translationMatrixNegativ(a)
	Rx2p = rotationMatrixX2Positiv(b, d)
	Rx2n = rotationMatrixX2Negativ(b, d)
	Rx3p = rotationMatrixX3Positiv(b, d)
	Rx3n = rotationMatrixX3Negativ(b, d)
	Rx3 = rotationMatrixX3(phi)
	rot = MATMUL(MATMUL(MATMUL(Tp,Rx3p), MATMUL(Rx2p, Rx3)), MATMUL(MATMUL(Rx2n, Rx3n), Tn))
END FUNCTION rot3Axis

FUNCTION translationMatrixPositiv(t) RESULT(M)
	! The translation matrix
	IMPLICIT NONE
	
	INTEGER								:: ii, jj
	REAL, DIMENSION(3,1), INTENT(IN)	:: t
	REAL, DIMENSION(4,4)				:: M
	
	DO ii = 1,4
		DO jj = 1,4
			M(ii,jj) = 0.0
			IF (ii == jj) M (ii,jj) = 1.0
		END DO
	END DO
	
	M(1,4) = t(1,1)
	M(2,4) = t(2,1)
	M(3,4) = t(3,1)
END FUNCTION translationMatrixPositiv

FUNCTION translationMatrixNegativ(t) RESULT(M)
	! The translation matrix
	IMPLICIT NONE
	
	INTEGER								:: ii, jj
	REAL, DIMENSION(3,1), INTENT(IN)	:: t
	REAL, DIMENSION(4,4)				:: M
	
	DO ii = 1,4
		DO jj = 1,4
			M(ii,jj) = 0.0
			IF (ii == jj) M (ii,jj) = 1.0
		END DO
	END DO
	
	M(1,4) = -t(1,1)
	M(2,4) = -t(2,1)
	M(3,4) = -t(3,1)
END FUNCTION translationMatrixNegativ

FUNCTION rotationMatrixX2Positiv(b, d) RESULT(R)
	! The rotation matrix
	IMPLICIT NONE
	
	INTEGER								:: ii, jj
	REAL, INTENT(IN)					:: d
	REAL, DIMENSION(3,1), INTENT(IN)	:: b
	REAL, DIMENSION(4,4)				:: R
	
	DO ii = 1,4
		DO jj = 1,4
			R(ii,jj) = 0.0
			IF (ii == jj) R (ii,jj) = 1.0
		END DO
	END DO
	
	R(1,1) = b(3,1)
	R(1,3) = d
	R(3,1) = -d
	R(3,3) = b(3,1)
END FUNCTION rotationMatrixX2Positiv

FUNCTION rotationMatrixX2Negativ(b, d) RESULT(R)
	! The rotation matrix
	IMPLICIT NONE
	
	INTEGER								:: ii, jj
	REAL, INTENT(IN)					:: d
	REAL, DIMENSION(3,1), INTENT(IN)	:: b
	REAL, DIMENSION(4,4)				:: R
	
	DO ii = 1,4
		DO jj = 1,4
			R(ii,jj) = 0.0
			IF (ii == jj) R (ii,jj) = 1.0
		END DO
	END DO
	
	R(1,1) = b(3,1)
	R(1,3) = -d
	R(3,1) = d
	R(3,3) = b(3,1)
END FUNCTION rotationMatrixX2Negativ

FUNCTION rotationMatrixX3Positiv(b, d) RESULT(R)
	! The rotation matrix
	IMPLICIT NONE
	
	INTEGER								:: ii, jj
	REAL, INTENT(IN)					:: d
	REAL, DIMENSION(3,1), INTENT(IN)	:: b
	REAL, DIMENSION(4,4)				:: R
	
	DO ii = 1,4
		DO jj = 1,4
			R(ii,jj) = 0.0
			IF (ii == jj) R (ii,jj) = 1.0
		END DO
	END DO
	
	R(1,1) = (1/d) * b(1,1)
	R(1,2) = - (1/d) * b(2,1)
	R(2,1) = (1/d) * b(2,1)
	R(2,2) = (1/d) * b(1,1)
END FUNCTION rotationMatrixX3Positiv

FUNCTION rotationMatrixX3Negativ(b, d) RESULT(R)
	! The rotation matrix
	IMPLICIT NONE
	
	INTEGER								:: ii, jj
	REAL, INTENT(IN)					:: d
	REAL, DIMENSION(3,1), INTENT(IN)	:: b
	REAL, DIMENSION(4,4)				:: R
	
	DO ii = 1,4
		DO jj = 1,4
			R(ii,jj) = 0.0
			IF (ii == jj) R (ii,jj) = 1.0
		END DO
	END DO
	
	R(1,1) = (1/d) * b(1,1)
	R(1,2) = (1/d) * b(2,1)
	R(2,1) = - (1/d) * b(2,1)
	R(2,2) = (1/d) * b(1,1)
END FUNCTION rotationMatrixX3Negativ

FUNCTION rotationMatrixX3(alpha) RESULT(R)
	! The rotation matrix
	IMPLICIT NONE
	
	INTEGER								:: ii, jj
	REAL, INTENT(IN)					:: alpha
	REAL, DIMENSION(4,4)				:: R
	
	DO ii = 1,4
		DO jj = 1,4
			R(ii,jj) = 0.0
			IF (ii == jj) R (ii,jj) = 1.0
		END DO
	END DO
	
	R(1,1) = cos(alpha)
	R(1,2) = -sin(alpha)
	R(2,1) = sin(alpha)
	R(2,2) = cos(alpha)
END FUNCTION rotationMatrixX3

FUNCTION bresenhamAlgorithm(P1,P2,d) RESULT(XYZ)
	IMPLICIT NONE
	
	INTEGER								:: ii, x1, y1, z1, x2, y2, z2, dx, dy, dz
	INTEGER								:: ax, ay, az, sx, sy, sz, x, y, z, idx
	INTEGER, INTENT(IN)					:: d
	INTEGER, DIMENSION(3,1), INTENT(IN)	:: P1, P2
	INTEGER, DIMENSION(d)				:: XR, YR, ZR
	INTEGER, DIMENSION(3*d)				:: XYZ
	REAL								:: xd, yd, zd
		
	DO ii = 1, d
		XR(ii) = 0
	END DO
	YR = XR
	ZR = XR
	
	x1 = P1(1,1)
	y1 = P1(2,1)
	z1 = P1(3,1)
	x2 = P2(1,1)
	y2 = P2(2,1)
	z2 = P2(3,1)
	dx = x2 - x1
	dy = y2 - y1
	dz = z2 - z1
	ax = ABS(dx) * 2
	ay = ABS(dy) * 2
	az = ABS(dz) * 2
	sx = SIGN(1, dx)
	sy = SIGN(1, dy)
	sz = SIGN(1, dz)
	x = x1
	y = y1
	z = z1
	idx = 1
	
	IF (ax >= MAX(ay, az)) THEN
		! x1 DOminant
		yd = ay - (ax / 2)
		zd = az - (ax / 2)
		DO WHILE(.TRUE.)
			XR(idx) = x
			YR(idx) = y
			ZR(idx) = z
			idx = idx + 1
			IF (x == x2) EXIT
			IF (yd >= 0) THEN
				! Move along y
				y = y + sy
				yd = yd - ax
			END IF
			IF (zd >= 0) THEN
				! Move along z
				z = z + sz
				zd = zd - ax
			END IF
			! Move along x
			x = x + sx
			yd = yd + ay
			zd = zd + az
		END DO
	ELSE IF (ay >= MAX(ax, az)) THEN
		! x2 DOminant
		xd = ax - (ay / 2)
		zd = az - (ay / 2)
		DO WHILE(.TRUE.)
			XR(idx) = x
			YR(idx) = y
			ZR(idx) = z
			idx = idx + 1
			IF (y == y2) EXIT
			IF(xd >= 0) THEN
				! Move along x
				x = x + sx
				xd = xd - ay
			END IF
			IF (zd >= 0) THEN
				! Move along z
				z = z + sz
				zd = zd - ay
			END IF
			! Move along y
			y = y + sy
			xd = xd + ax
			zd = zd + az
		END DO
	ELSE IF (az >= MAX(ax, ay)) THEN
		! x3 DOminant
		xd = ax - (az / 2)
		yd = ay - (az / 2)
		DO WHILE(.TRUE.)
			XR(idx) = x
			YR(idx) = y
			ZR(idx) = z
			idx = idx + 1
			IF (z == z2) EXIT
			IF (xd >= 0) THEN
				x = x + sx
				xd = xd - az
			END IF
			IF (yd >= 0) THEN
				y = y + sy
				yd = yd - az
			END IF
			z = z + sz
			xd = xd + ax
			yd = yd + ay
		END DO
	END IF
	
	DO ii = 1, d
		XYZ(ii) = XR(ii)
		XYZ(d + ii) = YR(ii)
		XYZ(2 * d + ii) = ZR(ii)
	END DO
END FUNCTION bresenhamAlgorithm

FUNCTION deg2rad(angleInDegrees)
	IMPLICIT NONE
	
	REAL					:: angleInDegrees, deg2rad
	REAL, PARAMETER 		:: pi = 3.141593
	
	deg2rad = (pi/180) * angleInDegrees
	
END FUNCTION deg2rad

FUNCTION rad2deg(angleInRadians)
	IMPLICIT NONE
	
	REAL					:: angleInRadians, rad2deg
	REAL, PARAMETER 		:: pi = 3.141593
	
	rad2deg = (180/pi) * angleInRadians
	
END FUNCTION rad2deg

FUNCTION getSize(fileName) RESULT(dimensions)
	IMPLICIT NONE

	CHARACTER(len=100), INTENT(IN)			:: fileName
	CHARACTER(len=100) 						:: dummy
	INTEGER 								:: status
	INTEGER, DIMENSION(3)	 				:: dimensions
	
	OPEN(UNIT=8, FILE=fileName, IOSTAT=status, STATUS='old')
	
	! First two lines contain the version number and a title
	READ(8, '(A)') dummy
	READ(8, '(A)') dummy
	!  The third line contains ASCII or BINARY
	READ(8, '(A)') dummy
	! Fourth line contains the type of dataset
	READ(8, '(A)') dummy
	! Dimensions
	READ(8, '(A10,1X,I10,1X,I10,1X,I10)') dummy, dimensions
	
	CLOSE(8)
END FUNCTION getSize

FUNCTION importVTK(fileName, a, b, c) RESULT(x)
	!IMPLICIT NONE
	
	CHARACTER(len=100), INTENT(IN)			:: fileName
	
	CHARACTER(len=100) 						:: dummy
	CHARACTER(len=20) 						:: dataset, fmt
	INTEGER 								:: status, np, ln
	INTEGER, INTENT(IN)						:: a, b, c
	INTEGER, DIMENSION(:), ALLOCATABLE 		:: x
	INTEGER, DIMENSION(a * b,c)				:: datas
	
	OPEN(UNIT=9, FILE=fileName, IOSTAT=status, STATUS='old')
	
	! 1
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 2
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 3
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 4
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 5
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 6
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 7
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 8
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 9
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	! 10
	READ(9, '(A)') dummy
	!WRITE(*,*) dummy
	
	ALLOCATE(x(a * b * c))
		
	WRITE(fmt,*) a * b
	DO ln = 1, c
		READ(9, '(I1,' // ADJUSTL(fmt) // '(2X, I1))') datas(:,ln)
	END DO
	x = RESHAPE((datas), (/ a * b * c /))
	WRITE(*,*) SIZE(x)
	!WRITE(*,*) x
	CLOSE(9)
END FUNCTION importVTK

SUBROUTINE mil_tensor(fileName)
	!===========================================================
	! Calculate MIL Tensor M
	!-----------------------------------------------------------
	! input ...
	! a(n,n) - array of coefficients for matrix A
	! b(n) - array of the right hand coefficients b
	! n - number of equations (size of matrix A)
	! output ...
	! x(n) - solutions
	!===========================================================
	
	IMPLICIT NONE
	
	CHARACTER(len = 10)						:: dummy
	CHARACTER(len=100), INTENT(IN)			:: fileName
	INTEGER									:: status, n, ii
	INTEGER, PARAMETER						:: q = 13
	REAL, PARAMETER							:: abserr = 1.0e-09
	REAL, DIMENSION(3)						:: center, evals
	REAL, DIMENSION(9)						:: u
	REAL, DIMENSION(10)						:: v
	REAL, DIMENSION(3,3)					:: eig, evecs, M, H
	REAL, DIMENSION(4,4)					:: A, T, R
	REAL, DIMENSION(:), ALLOCATABLE			:: theta, phi, mil, x1, x2, x3, d2
	REAL, DIMENSION(:,:), ALLOCATABLE		:: D
	
	OPEN(UNIT=q, FILE=fileName, IOSTAT=status, STATUS='old')
	READ(q,*) n
	
	ALLOCATE(theta(n))
	ALLOCATE(phi(n))
	ALLOCATE(mil(n))
	ALLOCATE(x1(n))
	ALLOCATE(x2(n))
	ALLOCATE(x3(n))
	
	DO ii = 1, n
		READ(q, '(1f6.4,A1,1f6.4,A1,3f8.4)') theta(ii), dummy, phi(ii), dummy, mil(ii)
		!WRITE(*,*) 'theta: ', theta(ii), 'phi: ', phi(ii), 'mil: ', mil(ii)
		
		! Convert from spherical coordinates to cartesian coordinates
		x1(ii) = MIL(ii) * SIN(theta(ii)) * COS(phi(ii))
		x2(ii) = MIL(ii) * SIN(theta(ii)) * SIN(phi(ii))
		x3(ii) = MIL(ii) * COS(theta(ii))
		!WRITE(*,*) 'x1(ii) :', x1(ii), 'x2(ii): ', x2(ii), 'x3(ii): ', x3(ii)
	END DO
	
	DEALLOCATE(theta)
	DEALLOCATE(phi)
	DEALLOCATE(mil)
	
	! fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx +
	! Hy + 2Iz + J = 0 and A + B + C = 3 constraint removing one extra
	! parameter
	ALLOCATE(D(n,9))
	ALLOCATE(d2(n))
	
	DO ii = 1, n
		D(ii, 1) = x1(ii) * x1(ii) + x2(ii) * x2(ii) - 2 * x3(ii) * x3(ii)
		D(ii, 2) = x1(ii) * x1(ii) + x3(ii) * x3(ii) - 2 * x2(ii) * x2(ii)
		D(ii, 3) = 2 * x1(ii) * x2(ii)
		D(ii, 4) = 2 * x1(ii) * x3(ii)
		D(ii, 5) = 2 * x2(ii) * x3(ii)
		D(ii, 6) = 2 * x1(ii)
		D(ii, 7) = 2 * x2(ii)
		D(ii, 8) = 2 * x3(ii)
		D(ii, 9) = 1 + 0 * x1(ii)
		
		d2(ii) = x1(ii) * x1(ii) + x2(ii) * x2(ii) + x3(ii) * x3(ii)
	END DO
	
	CALL gauss_2(MATMUL(TRANSPOSE(D), D), MATMUL(TRANSPOSE(D), d2), u, 9)
	
	! Find the ellipsoid parameters and convert back to the conventional algebraic form
	v(1) = u(1) + u(2) - 1
	v(2) = u(1) - 2 * u(2) - 1
	v(3) = u(2) - 2 * u(1) - 1
	v(4:10) = u(3:9)
	
	! Form the algebraic form of the ellipsoid
	A = RESHAPE((/ v(1), v(4), v(5), v(7), v(4), v(2), v(6), v(8), v(5), v(6), v(3), v(9), v(7), v(8), v(9), v(10) /), (/ 4, 4/))
	
	! Find the center of the ellipsoid
	CALL gauss_2(-A(1:3, 1:3), v(7:9), center, 3)
	
	! Form the corresponding translation matrix
	T = 0
	FORALL(ii = 1 : 4) T(ii,ii) = 1
	T(4, 1:3) = center(1:3)
	!WRITE(*,*) T
	R = MATMUL(MATMUL(T, A), TRANSPOSE(T))
	!WRITE(*,*) R
	
	!WRITE(*,*) R(1:3, 1:3) / (-R(4,4))
	
	eig = R(1:3, 1:3) / (-R(4,4))
		
	CALL jacobi(eig, evecs, abserr, 3)
	
	DO ii = 1, 3
		evals(ii) = eig(ii, ii)
	END DO
	
	
	
	
	IF (ABS(v(10)) > 1E-6) THEN
		v = -v / v(10)
	END IF
		
	WRITE(*,*) '-----------------'
	
	M = RESHAPE((/ v(1), v(4), v(5), v(4), v(2), v(6), v(5), v(6), v(3) /), (/ 3, 3/))
	
	WRITE(*,*) 'MIL tensor M: '
	DO ii = 1, SIZE(M, 1)
		WRITE(*,'(20G12.4)') M(ii,:)
	END DO
	
	H = M**(-1/2)
	
	WRITE(*,*) 'Fabric tensor H: '
	DO ii = 1, SIZE(H, 1)
		WRITE(*,'(20G12.4)') H(ii,:)
	END DO 
	
	DEALLOCATE(D)
	DEALLOCATE(d2)
	DEALLOCATE(x1)
	DEALLOCATE(x2)
	DEALLOCATE(x3)
	
END SUBROUTINE mil_tensor

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

	INTEGER				:: n
	REAL				:: a(n,n), b(n), x(n)
	REAL				:: s(n)
	REAL				:: c, pivot, store
	INTEGER				:: i, j, k, l

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
	INTEGER					:: i, j, k, n
	REAL					:: a(n,n),x(n,n)
	REAL					:: abserr, b2, bar
	REAL					:: beta, coeff, c, s, cs, sc

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
	END DO
	RETURN
END SUBROUTINE jacobi