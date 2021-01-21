PROGRAM main
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
	
	CHARACTER(len=100)						:: fileName, fileNameExport
	INTEGER 								:: numberOfOrientations, increment, ii, status
	INTEGER, DIMENSION(3)					:: dimensions
	INTEGER, DIMENSION (:), ALLOCATABLE		:: fileData
	INTEGER, DIMENSION (:,:,:), ALLOCATABLE	:: voxel
	REAL 									:: theta, phi, thetaR, phiR, calculateMIL, MIL, start, finish
	REAL, PARAMETER 						:: pi = 3.141593, r = 1
	REAL, DIMENSION(3,1)					:: p0, p1, n
	
	CALL CPU_TIME(start)
	
	! Number of ranDOmly generated orientation (positiv integer, minimum 9)
	numberOfOrientations = 100
	! Distance between two created lines (positiv integer)
	increment = 1
	! File name
	fileName = 'Knochenprobe2_1mm_1.vtk'
	
	fileNameExport = fileName(1:LEN_TRIM(fileName)-4) // '.dat'
	
	dimensions = getSize(fileName)
	ALLOCATE(fileData(dimensions(1) * dimensions(2) * dimensions(3)))
	WRITE(*,*) dimensions
	fileData = importVTK(fileName, dimensions(1), dimensions(2), dimensions(3))
	voxel = RESHAPE(fileData, (/ dimensions(1), dimensions(2), dimensions(3) /))
	DEALLOCATE(fileData)
	
	OPEN(UNIT = 1, FILE = fileNameExport, IOSTAT = status)
	
	DO ii = 1, numberOfOrientations
		CALL RANDOM_NUMBER(thetaR)
		theta = thetaR * pi
		!theta = pi / 4
		CALL RANDOM_NUMBER(phiR)
		phi = phiR * 2 * pi
		!phi = pi / 2
		! Origin
		p0 = RESHAPE((/0, 0, 0 /), (/ 3, 1 /))
		! Spherical coordinates to Cartesian coordinates
		p1(1,1) = r * SIN(theta) * COS(phi)
		p1(2,1) = r * SIN(theta) * SIN(phi)
		p1(3,1) = r * COS(theta)
		
		n = (1 / SQRT(SUM((p1 - p0)**2))) * (p1 - p0)
		MIL = calculateMIL(n, REAL(dimensions(1)), REAL(dimensions(2)), REAL(dimensions(3)), voxel)
		WRITE(*,*) 'ii: ',ii , '/ ', numberOfOrientations, 'theta: ', theta, 'phi: ', phi, 'MIL: ', MIL
		
		WRITE(1,'(1f6.4,A1,1f6.4,A1,3f8.4)') theta, ',', phi, ',', MIL
	END DO
	
	DEALLOCATE(voxel)
	CLOSE(1)
	CALL CPU_TIME(finish)
	PRINT '("Time = ",f8.3," seconds.")',finish-start
	
END PROGRAM	main

FUNCTION calculateMIL(n, a, b, c, voxel)
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
	
	REAL									:: calculateMIL, h, cv, a, b, c, dr, ru
	REAL, PARAMETER							:: pi = 3.141593
	REAL, DIMENSION(3,1)					:: n, m, v0
	REAL, DIMENSION(4,1)					:: pQ1l4, pQ2l4, pQ3l4, pQ4l4
	
	INTEGER, DIMENSION(nint(a), nint(b), nint(c))	:: voxel
	
	h = 0.0
	cv = 0.0
	
	!write(*,*) SIZE(voxel)
	
	IF (n(1,1) == 1 .AND. n(2,1) == 0 .AND. n(3,1) == 0) THEN
	ELSE
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
		
		! First
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
		!WRITE(*,*) 'l: ', l
		
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
				!WRITE(*,*) 'q', q
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
					x1(ii) = x(d + ii)
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

	END IF
		
	calculateMIL = h / cv
END FUNCTION calculateMIL

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

FUNCTION getSize(fileName) RESULT(dimensions)
	IMPLICIT NONE

	CHARACTER(len=100), INTENT(IN)		:: fileName
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
	CHARACTER(len=20) 						:: dataset, FMT
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
		
	WRITE(FMT,*) a * b
	DO ln = 1, c
		READ(9, '(I1,' // ADJUSTL(FMT) // '(2X, I1))') datas(:,ln)
	END DO
	x = RESHAPE((datas), (/ a * b * c /))
	WRITE(*,*) SIZE(x)
	!WRITE(*,*) x
	CLOSE(9)
END FUNCTION importVTK