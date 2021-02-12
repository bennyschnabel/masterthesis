MODULE calculate_mil
        USE ISO_FORTRAN_ENV
        IMPLICIT NONE

        CONTAINS
                SUBROUTINE test(n, dims, array, mil)
                        IMPLICIT NONE

                        ! Parameter
                        REAL(KIND = real64), PARAMETER :: pi = 3.141593
                        
                        ! External variables
                        INTEGER(KIND = int64), DIMENSION(3), INTENT(IN) :: dims
                        REAL(KIND = real64), DIMENSION(3,1), INTENT(IN) :: n
                        INTEGER(KIND = int64), DIMENSION(dims(1),dims(2),dims(3)), INTENT(IN) :: array
                        REAL(KIND = real64), INTENT(OUT) :: mil

                        !-- Internal variables  
                        INTEGER(KIND = int64) :: i, j, k, l
                        INTEGER(KIND = int64) :: d, lmin, counter, indizesL, indizesU
                        INTEGER(KIND = int64), DIMENSION(3,1) :: ps
                        INTEGER(KIND = int64), DIMENSION(3,1) :: pQ1, pQ2, pQ3, pQ4
                        INTEGER(KIND = int64), DIMENSION(3,1) :: pO1, pO2, pO3, pO4
                        INTEGER(KIND = int64), DIMENSION(3,1) :: pQ12, pQ34, pO12, pO34
                        INTEGER(KIND = int64), DIMENSION(3,1) :: pQ, pO
                        INTEGER(KIND = int64), DIMENSION(3,1) :: SP, EP

                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1Q12, x2Q12, x3Q12
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1Q34, x2Q34, x3Q34
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1O12, x2O12, x3O12
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1O34, x2O34, x3O34
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1Q, x2Q, x3Q
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1O, x2O, x3O
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1, x2, x3
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: x1SE, x2SE, x3SE
                        INTEGER(KIND = int64), DIMENSION(:), ALLOCATABLE :: indizes

                        REAL(KIND = real64) :: h
                        REAL(KIND = real64) :: cv
                        REAL(KIND = real64) :: dr
                        REAL(KIND = real64) :: ru
                        
                        REAL(KIND = real64), DIMENSION(3,1) :: m
                        REAL(KIND = real64), DIMENSION(3,1) :: v0
                        REAL(KIND = real64), DIMENSION(4,1) :: pQ114, pQ214, pQ314, pQ414

                        REAL(KIND = real64), DIMENSION(4,4) :: rot
                        
                        h = 0.0
                        cv = 0.0

                        ! Room diagonal
                        dr = SQRT(REAL(dims(1))**2 + REAL(dims(2))**2 + REAL(dims(3))**2)
                        ! Radius of the sphere (Half room diagonal)
                        ru = 0.5 * dr
                        ! Origin of the sphere (center of the read data)
                        m = RESHAPE((/ REAL(dims(1))/2, REAL(dims(2))/2, REAL(dims(3))/2 /), (/ 3, 1 /))
                        ! Point on the sphere in direction n
                        ps = NINT(m + ru * n)
                        ! First directional vector (direction 1, 0)
                        
                        v0 = RESHAPE((/ 1.0_real64, 0.0_real64, (-n(1,1) / n(3,1)) /), (/ 3, 1 /))
                        ! Normalize v0
                        v0 = (1 / SQRT(SUM(v0**2))) * v0

                        ! Boundary point 1 on the plane Q
                        pQ1 = NINT(ps + ru * v0)
                        pQ114 = RESHAPE((/ REAL(pQ1(1,1)), REAL(pQ1(2,1)), REAL(pQ1(3,1)), 1.0 /), (/ 4,1 /))

                        ! Boundary point 2 on the plane Q
                        CALL rot3Axis(m, n, pi/2, rot)
                        pQ214 = MATMUL(rot, pQ114)
                        pQ2 = RESHAPE(NINT((/ pQ214(1,1), pQ214(2,1), pQ214(3,1) /)), (/ 3,1 /))

                        ! Boundary point 3 on the plane Q
                        CALL rot3Axis(m, n, pi, rot)
                        pQ314 = MATMUL(rot, pQ114)
                        pQ3 = RESHAPE(NINT((/ pQ314(1,1), pQ314(2,1), pQ314(3,1) /)), (/ 3,1 /))

                        ! Boundary point 4 on the plane Q
                        CALL rot3Axis(m, n, (3 * pi) / 2, rot)
                        pQ414 = MATMUL(rot, pQ114)
                        pQ4 = RESHAPE(NINT((/ pQ414(1,1), pQ414(2,1), pQ414(3,1) /)), (/ 3,1 /))

                        ! Boundary point 1 on the plane O
                        pO1 = pQ1 - NINT(dr * n)
                        ! Boundary point 2 on the plane O
                        pO2 = pQ2 - NINT(dr * n)
                        ! Boundary point 3 on the plane O
                        pO3 = pQ3 - NINT(dr * n)
                        ! Boundary point 4 on the plane O
                        pO4 = pQ4 - NINT(dr * n)

                        ! Line between Q1 and Q2 (Q12)
                        d = MAXVAL(ABS(pQ2 - pQ1) + 1)
                        ALLOCATE(x1Q12(d))
                        ALLOCATE(x2Q12(d))
                        ALLOCATE(x3Q12(d))
                        CALL bresenhamAlgorithm(pQ1, pQ2, d, x1Q12, x2Q12, x3Q12)

                        ! Line between Q3 and Q4
                        d = MAXVAL(ABS(pQ4 - pQ3) + 1)
                        ALLOCATE(x1Q34(d))
                        ALLOCATE(x2Q34(d))
                        ALLOCATE(x3Q34(d))
                        CALL bresenhamAlgorithm(pQ3, pQ4, d, x1Q34, x2Q34, x3Q34)

                        ! Line between O1 and O2
                        d = MAXVAL(ABS(pO1 - pO2) + 1)
                        ALLOCATE(x1O12(d))
                        ALLOCATE(x2O12(d))
                        ALLOCATE(x3O12(d))
                        CALL bresenhamAlgorithm(pO1, pO2, d, x1O12, x2O12, x3O12)

                        ! Line between O3 and O4
                        d = MAXVAL(ABS(pO3 - pO4) + 1)
                        ALLOCATE(x1O34(d))
                        ALLOCATE(x2O34(d))
                        ALLOCATE(x3O34(d))
                        CALL bresenhamAlgorithm(pO3, pO4, d, x1O34, x2O34, x3O34)

                        lmin = MIN(SIZE(x1Q12), SIZE(x1Q34), SIZE(x1O12), SIZE(x1O34))

                        DO i = 1, (lmin - 1)
                                pQ12 = RESHAPE((/x1Q12(i), x2Q12(i), x3Q12(i) /), (/ 3, 1 /))
                                pQ34 = RESHAPE((/x1Q34(lmin - i), x2Q34(lmin - i), x3Q34(lmin - i) /), (/ 3, 1 /))
                                pO12 = RESHAPE((/x1O12(i), x2O12(i), x3O12(i) /), (/ 3, 1 /))
                                pO34 = RESHAPE((/x1O34(lmin - i), x2O34(lmin - i), x3O34(lmin - i) /), (/ 3, 1 /))

                                ! Line between Q12 and Q34 (xQ)
                                d = MAXVAL(ABS(pQ12 - pQ34) + 1)
                                ALLOCATE(x1Q(d))
                                ALLOCATE(x2Q(d))
                                ALLOCATE(x3Q(d))
                                CALL bresenhamAlgorithm(pQ12, pQ34, d, x1Q, x2Q, x3Q)

                                ! Line between O12 and O34 (xO)
                                d = MAXVAL(ABS(pO12 - pO34) + 1)
                                ALLOCATE(x1O(d))
                                ALLOCATE(x2O(d))
                                ALLOCATE(x3O(d))
                                CALL bresenhamAlgorithm(pO12, pO34, d, x1O, x2O, x3O)

                                DO j = 1, SIZE(x1Q)
                                        pQ = RESHAPE((/x1Q(j), x2Q(j), x3Q(j) /), (/ 3, 1 /))
                                        pO = RESHAPE((/x1O(j), x2O(j), x3O(j) /), (/ 3, 1 /))

                                        ! Line between Q and O (x)
                                        d = MAXVAL(ABS(pQ - pO) + 1)
                                        ALLOCATE(x1(d))
                                        ALLOCATE(x2(d))
                                        ALLOCATE(x3(d))
                                        CALL bresenhamAlgorithm(pQ, pO, d, x1, x2, x3)

                                        counter = 0
                                        DO k = 1, SIZE(x1)
                                                IF ((x1(k) >= 1 .AND. x1(k) <= dims(1)) .AND. &
                                                        (x2(k) >=1 .AND. x2(k) <= dims(2)) .AND. &
                                                        (x3(k) >=1 .AND. x3(k) <= dims(3))) THEN
                                                        counter = counter + 1
                                                END IF
                                        END DO

                                        IF (counter > 2) THEN
                                                ALLOCATE(indizes(counter))
                                                l = 1
                                                DO k = 1, SIZE(x1)
                                                        IF ((x1(k) >= 1 .AND. x1(k) <= dims(1)) .AND. &
                                                                (x2(k) >=1 .AND. x2(k) <= dims(2)) .AND. &
                                                                (x3(k) >=1 .AND. x3(k) <= dims(3))) THEN
                                                                indizes(l) = k
                                                                l = l + 1
                                                        END IF
                                                END DO
                                                indizesL = indizes(1)
                                                indizesU = indizes(counter)
                                                DEALLOCATE(indizes)

                                                SP = RESHAPE((/ x1(indizesL), x2(indizesL), x3(indizesL) /), (/ 3, 1 /))
                                                EP = RESHAPE((/ x1(indizesU), x2(indizesU), x3(indizesU) /), (/ 3, 1 /))

                                                h = h + SQRT(SUM(REAL(EP - SP)**2))

                                                d = MAXVAL(ABS(SP - EP) + 1)
                                                ALLOCATE(x1SE(d))
                                                ALLOCATE(x2SE(d))
                                                ALLOCATE(x3SE(d))
                                                CALL bresenhamAlgorithm(SP, EP, d, x1SE, x2SE, x3SE)

                                                ! Calculation of the intersection points with the cancellous bone
                                                DO k = 1, (SIZE(x1SE) - 1)
                                                        IF ((array(x1SE(k), x2SE(k), x3SE(k)) == 0) .AND. &
                                                                (array(x1SE(k+1), x2SE(k+1), x3SE(k+1)) == 1)) THEN
                                                                cv = cv + 1.0
                                                        END IF
                                                END DO
                                                !WRITE(*,*) cv

                                                DEALLOCATE(x1SE)
                                                DEALLOCATE(x2SE)
                                                DEALLOCATE(x3SE)                                                
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

                        mil = h / cv
                END SUBROUTINE test

                SUBROUTINE mil_fabric_tensor(fileNameExport, fileNameExportTensor, domainNo, fun1, fun2)
                IMPLICIT NONE

                ! External variables
                CHARACTER(LEN = 100), INTENT(IN) :: fileNameExport
                CHARACTER(LEN = 100), INTENT(IN) :: fileNameExportTensor
                INTEGER(KIND = int64), INTENT(IN) :: domainNo
                INTEGER(KIND = int64), INTENT(IN) :: fun1
                INTEGER(KIND = int64), INTENT(IN) :: fun2

                ! Internal variables
                CHARACTER(LEN = 10) :: dummy
                CHARACTER(LEN = 100) :: formatString

                INTEGER(KIND = int64) :: status, n, ii, jj
                INTEGER(KIND = int64), PARAMETER :: q = 13

                REAL(KIND = real64), PARAMETER :: abserr = 1.0e-09
                REAL(KIND = real64), DIMENSION(3) :: center, evals
                REAL(KIND = real64), DIMENSION(9) :: u
                REAL(KIND = real64), DIMENSION(10) :: v
                REAL(KIND = real64), DIMENSION(3,3) :: eig, evecs, M, Minv, eig2, evecs2, H, D05
                REAL(KIND = real64), DIMENSION(4,4) :: A, T, R
                REAL(KIND = real64), DIMENSION(:), ALLOCATABLE :: theta, phi, mil, x1, x2, x3, d2
                REAL(KIND = real64), DIMENSION(:,:), ALLOCATABLE :: D

                OPEN(UNIT=fun1, FILE=fileNameExport, IOSTAT=status, STATUS='old')
                READ(fun1,*) n

                ALLOCATE(theta(n))
                ALLOCATE(phi(n))
                ALLOCATE(mil(n))
                ALLOCATE(x1(n))
                ALLOCATE(x2(n))
                ALLOCATE(x3(n))

                DO ii = 1, n
                        READ(fun1, '(1f6.4,A1,1f6.4,A1,3f8.4)') theta(ii), dummy, phi(ii), dummy, mil(ii)

                        ! Convert from spherical coordinates to cartesian coordinates
                        x1(ii) = MIL(ii) * SIN(theta(ii)) * COS(phi(ii))
                        x2(ii) = MIL(ii) * SIN(theta(ii)) * SIN(phi(ii))
                        x3(ii) = MIL(ii) * COS(theta(ii))
                END DO

                CLOSE(fun1)

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

                CALL gauss_2(MATMUL(TRANSPOSE(D), D), MATMUL(TRANSPOSE(D), d2), u, 9_int64)

                ! Find the ellipsoid parameters and convert back to the conventional algebraic form
                v(1) = u(1) + u(2) - 1
                v(2) = u(1) - 2 * u(2) - 1
                v(3) = u(2) - 2 * u(1) - 1
                v(4:10) = u(3:9)

                ! Form the algebraic form of the ellipsoid
                A = RESHAPE((/ v(1), v(4), v(5), v(7), v(4), v(2), v(6), v(8), v(5) &
                        , v(6), v(3), v(9), v(7), v(8), v(9), v(10) /), (/ 4, 4/))

                ! Find the center of the ellipsoid
                CALL gauss_2(-A(1:3, 1:3), v(7:9), center, 3_int64)

                ! Form the corresponding translation matrix
                T = 0
                FORALL(ii = 1 : 4) T(ii,ii) = 1
                T(4, 1:3) = center(1:3)

                R = MATMUL(MATMUL(T, A), TRANSPOSE(T))

                eig = R(1:3, 1:3) / (-R(4,4))

                CALL jacobi(eig, evecs, abserr, 3_int64)

                DO ii = 1, 3
                        evals(ii) = eig(ii, ii)
                END DO

                IF (ABS(v(10)) > 1E-6) THEN
                        v = -v / v(10)
                END IF

                ! Calculate MIL tensor M
                M = RESHAPE((/ v(1), v(4), v(5), v(4), v(2), v(6), v(5), v(6), v(3) /), (/ 3, 3/))

                WRITE(*,*) '-------Results-------'
                WRITE(*,*) 'MIL tensor M: '
                DO ii = 1, SIZE(M, 1)
                        WRITE(*,'(20G12.4)') M(ii,:)
                END DO

                ! Calculate fabric tensor F
                CALL inverse(M, Minv, 3_int64)

                ! Important to get a 'symmetric matric' for jacobi algorithm
                eig2 =  nint(Minv * 1000.0) / 1000.0

                CALL jacobi(eig2, evecs2, abserr, 3_int64)

                DO ii = 1, 3
                        DO jj = 1,3
                                D05(ii, jj) = 0
                                IF (ii == jj) THEN
                                        D05(ii, jj) = sqrt(eig2(ii, jj))
                                END IF
                        END DO
                END DO

                H = MATMUL(MATMUL(evecs2, D05), TRANSPOSE(evecs2))

                WRITE(*,*) 'Fabric tensor H: '
                DO ii = 1, SIZE(H, 1)
                        WRITE(*,'(20G12.4)') H(ii,:)
                END DO
                
                OPEN(UNIT = fun2, FILE = fileNameExportTensor, IOSTAT = status, ACCESS = 'append')
                formatString = '(I6, 9(A1, F11.5))'

                WRITE(fun2,formatString) domainNo, ';', H(1,1), ';', H(1,2), ';', H(1,3), ';', H(2,1), ';', H(2,2), &
                        ';', H(3,2), ';', H(1,3), ';', H(2,3), ';', H(3,3)
                !WRITE(fun2,*) 'MIL tensor M: '
                !WRITE(fun2,'(3f12.7)') M
                !WRITE(fun2,*) 'Fabric tensor H: '
                !WRITE(fun2,'(3f12.7)') H
                CLOSE(fun2)

                DEALLOCATE(D)
                DEALLOCATE(d2)
                DEALLOCATE(x1)
                DEALLOCATE(x2)
                DEALLOCATE(x3)
                        
                END SUBROUTINE mil_fabric_tensor

                SUBROUTINE rot3Axis(a, b, phi, rot)

                        IMPLICIT NONE

                        ! External variables
                        REAL(KIND = real64), INTENT(IN) :: phi
                        REAL(KIND = real64), DIMENSION(3,1), INTENT(IN) :: a, b
                        REAL(KIND = real64), DIMENSION(4,4), INTENT(OUT) :: rot
                        
                        ! Internal variables
                        INTEGER(KIND = int32) :: i, j
                        REAL(KIND = real64) :: d
                        REAL(KIND = real64), DIMENSION(4,4) :: Tp, Tn, Rx2p, Rx2n, Rx3p, Rx3n, Rx3

                        d = sqrt(b(1,1)**2 + b(2,1)**2)

                        ! Loop
                        DO i = 1, 4
                                DO j = 1, 4
                                        Tp(i,j) = 0.0
                                        IF (i == j) Tp(i,j) = 1.0
                                END DO
                        END DO
                        Tn = Tp
                        Rx2p = Tp
                        Rx2n = Tp
                        Rx3p = Tp
                        Rx3n = Tp
                        Rx3 = Tp

                        ! Translation matrix positiv
                        Tp(1,4) = a(1,1)
                        Tp(2,4) = a(2,1)
                        Tp(3,4) = a(3,1)
         
                        ! translation matrix negativ
                        Tn(1,4) = -a(1,1)
                        Tn(2,4) = -a(2,1)
                        Tn(3,4) = -a(3,1)

                        ! Rotation matrix x2 positiv
                        Rx2p(1,1) = b(3,1)
                        Rx2p(1,3) = d
                        Rx2p(3,1) = -d
                        Rx2p(3,3) = b(3,1)

                        ! Rotation matrix x2 negativ
                        Rx2n(1,1) = b(3,1)
                        Rx2n(1,3) = -d
                        Rx2n(3,1) = d
                        Rx2n(3,3) = b(3,1)

                        ! Rotation matrix x3 positiv
                        Rx3p(1,1) = (1/d) * b(1,1)
                        Rx3p(1,2) = -(1/d) * b(2,1)
                        Rx3p(2,1) = (1/d) * b(2,1)
                        Rx3p(2,2) = (1/d) * b(1,1)

                        ! Rotation matrix x3 negativ
                        Rx3n(1,1) = (1/d) * b(1,1)
                        Rx3n(1,2) = (1/d) * b(2,1)
                        Rx3n(2,1) = -(1/d) * b(2,1)
                        Rx3n(2,2) = (1/d) * b(1,1)

                        ! Rotation matrix x3
                        Rx3(1,1) = COS(phi)
                        Rx3(1,2) = -SIN(phi)
                        Rx3(2,1) = SIN(phi)
                        Rx3(2,2) = COS(phi)

                        rot = MATMUL(MATMUL(MATMUL(Tp,Rx3p), MATMUL(Rx2p, Rx3)), MATMUL(MATMUL(Rx2n, Rx3n), Tn))
                END SUBROUTINE rot3Axis

                SUBROUTINE bresenhamAlgorithm(P1, P2, d, XR, YR, ZR)
                        IMPLICIT NONE

                        INTEGER(KIND = int64), INTENT(IN) :: d
                        INTEGER(KIND = int64), DIMENSION(3,1), INTENT(IN) :: P1, P2
                        INTEGER(KIND = int64), DIMENSION(d), INTENT(OUT) :: XR, YR, ZR

                        INTEGER(KIND = int64) :: ii, x1, y1, z1, x2, y2, z2, dx, dy, dz
                        INTEGER(KIND = int64) :: ax, ay, az, sx, sy, sz, x, y, z, idx
                        REAL(KIND = real64) :: xd, yd, zd

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
                        sx = SIGN(1_int64, dx)
                        sy = SIGN(1_int64, dy)
                        sz = SIGN(1_int64, dz)
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
                END SUBROUTINE bresenhamAlgorithm

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
                INTEGER(KIND = int64) :: i, j, k, n
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

                SUBROUTINE importVTK(fun, fileName, dims, array)

                        IMPLICIT NONE

                        ! External variables
                        CHARACTER(len=100), INTENT(IN) :: fileName
                        INTEGER(KIND = int64), INTENT(IN) :: fun
                        INTEGER(KIND = int64), DIMENSION(3), INTENT(IN) :: dims
                        INTEGER(KIND = int64), DIMENSION(dims(1),dims(2),dims(3)), INTENT(OUT) :: array

                        ! Internal variables
                        CHARACTER(LEN = 100) :: n2s
                        CHARACTER(LEN = 100) :: dummy
                        INTEGER(KIND = int64) :: i
                        INTEGER(KIND = int64) :: status
                        INTEGER(KIND = int64), DIMENSION(dims(1) * dims(2), dims(3)) :: array_help

                        OPEN(UNIT=fun, FILE=fileName, IOSTAT=status, STATUS='old')

                        ! 1
                        READ(fun, '(A)') dummy
                        ! 2
                        READ(fun, '(A)') dummy
                        ! 3
                        READ(fun, '(A)') dummy
                        ! 4
                        READ(fun, '(A)') dummy
                        ! 5
                        READ(fun, '(A)') dummy
                        ! 6
                        READ(fun, '(A)') dummy
                        ! 7
                        READ(fun, '(A)') dummy
                        ! 8
                        READ(fun, '(A)') dummy
                        ! 9
                        READ(fun, '(A)') dummy
                        ! 10
                        READ(fun, '(A)') dummy

                        !ALLOCATE(array(dims(1), dims(2), dims(3)))

                        WRITE(n2s,*) dims(1) * dims(2)
                        DO i = 1, dims(3)
                                READ(fun, '(I1,' // ADJUSTL(n2s) // '(2X, I1))') array_help(:,i)
                        END DO

                        array = RESHAPE((array_help), (/ dims(1), dims(2), dims(3) /))
        
                        CLOSE(fun)
                END SUBROUTINE importVTK

END MODULE calculate_mil

