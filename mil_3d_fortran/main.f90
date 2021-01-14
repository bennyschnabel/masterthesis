program main
	implicit none

	integer, parameter :: out_unit=20
	integer :: numberOfOrientations, increment, ii
	real, parameter :: pi = 3.141593, r = 1
	real, dimension(3,1) :: p0, p1, n
	real :: theta, phi, thetaR, phiR, calculate_mil_3d, MIL

	! Number of randomly generated orientation (positiv integer, minimum 9)
	numberOfOrientations = 1
	! Distance between two created lines (positiv integer)
	increment = 1
	
	do ii = 1, numberOfOrientations
		call random_number(thetaR)
		!theta = thetaR * pi
		theta = pi / 4
		call random_number(phiR)
		!phi = phiR * 2 * pi
		phi = pi / 2
		! Origin
		p0 = reshape((/0, 0, 0 /), (/ 3, 1 /))
		! Spherical coordinates to Cartesian coordinates
		p1(1,1) = r * sin(theta) * cos(phi)
		p1(2,1) = r * sin(theta) * sin(phi)
		p1(3,1) = r * cos(theta)
		
		n = (1 / sqrt(sum((p1 - p0)**2))) * (p1 - p0)
		write(*,*) 'ii: ',ii , '/ ', numberOfOrientations
		MIL = calculate_mil_3d(n, 73.0, 73.0, 71.0)
		write(*,*) 'ii: ',ii , '/ ', numberOfOrientations, 'theta: ', theta, 'phi: ', phi, 'MIL: ', MIL
	end do
	
end program	main

function calculate_mil_3d(n, a, b, c)
	implicit none
	
	interface
		function rot3axis(a, b, phi) result (rot)
			REAL, dimension(3,1), intent(in)   :: a, b
			real :: phi
			real, dimension(4,4) :: rot
		end function
		
		function bresenham(P1, P2) result(XR1)
			integer, dimension(3,1), intent(in) :: P1, P2
			integer, dimension(:), allocatable :: XR1
		end function
	end interface
	
	integer, dimension(3,1) :: ps, pQ1, pQ2, pQ3, pQ4, pO1, pO2, pO3, pO4
    real :: calculate_mil_3d, h, cv, dr, ru, a, b, c
	real, parameter :: pi = 3.141593
	real, dimension(4, 1) :: pQ14, pQ24, pQ34, pQ44
	real, dimension(3, 1) :: n, m, v0
	integer, dimension(:), allocatable :: x1
	
	h = 0.0
	cv = 0.0
		
	if (n(1,1) == 1 .and. n(2,1) == 0 .and. n(3,1) == 0) then
		print *, '100'
	else if (n(1,1) == 0.0 .and. n(2,1) == 1.0 .and. n(3,1) == 0.0) then
		print *, '010'
	else
		! Room diagonal
		dr = sqrt(real(a)**2 + real(b)**2 + real(c)**2)
		! Radius of the sphere (Half room diagonal)
		ru = 0.5 * dr
		! Origin of the sphere (center of the read data)
		m = reshape((/ a/2, b/2, c/2 /), (/ 3, 1 /))
		write(*,*) 'm: ', m
		! Point on the sphere in direction n
		ps = nint(m + ru * n)
		write(*,*) 'ps: ', ps
		! First directional vector (direction 1, 0)
		v0 = reshape((/ 1.0, 0.0, (-n(1,1) / n(3,1)) /), (/ 3, 1 /))
		v0 = (1 / sqrt(sum(v0**2))) * v0
		write(*,*) 'v0: ', v0
		
		! Boundary point 1 on the plane Q
		pQ1 = nint(ps + ru * v0)
		pQ14 = reshape((/ real(pQ1(1,1)), real(pQ1(2,1)), real(pQ1(3,1)), 1.0 /), (/ 4,1 /))
		write(*,*) 'pQ1: ', pQ1
		
		! Boundary point 2 on the plane Q
		pQ24 = matmul(rot3axis(m, n, pi/2), pQ14)
		pQ2 = reshape(nint((/ pQ24(1,1), pQ24(2,1), pQ24(3,1) /)), (/ 3,1 /))
		write(*,*) 'pQ2: ', pQ2
		
		! Boundary point 3 on the plane Q
		pQ34 = matmul(rot3axis(m, n, pi), pQ14)
		pQ3 = reshape(nint((/ pQ34(1,1), pQ34(2,1), pQ34(3,1) /)), (/ 3,1 /))
		write(*,*) 'pQ3: ', pQ3
		
		! Boundary point 4 on the plane Q
		pQ44 = matmul(rot3axis(m, n, (3 * pi) / 2), pQ14)
		pQ4 = reshape(nint((/ pQ44(1,1), pQ44(2,1), pQ44(3,1) /)), (/ 3,1 /))
		write(*,*) 'pQ4: ', pQ4
		
		pO1 = pQ1 - nint(dr * n)
		write(*,*) 'pO1: ', pO1
		
		pO2 = pQ2 - dr * n
		write(*,*) 'pO2: ', pO2
		
		pO3 = pQ3 - dr * n
		write(*,*) 'pO3: ', pO3
		
		pO4 = pQ4 - dr * n
		write(*,*) 'pO4: ', pO4
		
		x1 = bresenham(pQ1, pQ2)
	end if

    calculate_mil_3d = h / cv
	
end function calculate_mil_3d
 
! Rotation around any axis
                            
function rot3axis(a, b, phi) result (rot)
	implicit none
	
	interface
		function translation_matrix_positiv(t) result(M)
			real, dimension(3,1), intent(in) :: t
			real, dimension(4,4) :: M
		end function
		
		function translation_matrix_negativ(t) result(M)
			real, dimension(3,1), intent(in) :: t
			real, dimension(4,4) :: M
		end function

		function rotation_matrix_x2_positiv(b, d) result(R)
			real, intent(in) :: d
			real, dimension(3,1), intent(in) :: b
			real, dimension(4,4) :: R
		end function
		
		function rotation_matrix_x2_negativ(b, d) result(R)
			real, intent(in) :: d
			real, dimension(3,1), intent(in) :: b
			real, dimension(4,4) :: R
		end function
		
		function rotation_matrix_x3_positiv(b, d) result(R)
			real, intent(in) :: d
			real, dimension(3,1), intent(in) :: b
			real, dimension(4,4) :: R
		end function
		
		function rotation_matrix_x3_negativ(b, d) result(R)
			real, intent(in) :: d
			real, dimension(3,1), intent(in) :: b
			real, dimension(4,4) :: R
		end function
		
		function rotation_matrix_x3(alpha) result(R)
			real, intent(in) :: alpha
			real, dimension(4,4) :: R
		end function
		
	end interface
	
	real, dimension(3,1), intent(in) :: a, b
	real :: phi, d
	real, dimension(4,4) :: rot, Tp, Tn, Rx2p, Rx2n, Rx3p, Rx3n, Rx3
	
	d = sqrt(b(1,1)**2 + b(2,1)**2)
	
	Tp = translation_matrix_positiv(a)
	Tn = translation_matrix_negativ(a)
	Rx2p = rotation_matrix_x2_positiv(b, d)
	Rx2n = rotation_matrix_x2_negativ(b, d)
	Rx3p = rotation_matrix_x3_positiv(b, d)
	Rx3n = rotation_matrix_x3_negativ(b, d)
	Rx3 = rotation_matrix_x3(phi)
	
	rot = matmul(matmul(matmul(Tp,Rx3p), matmul(Rx2p, Rx3)), matmul(matmul(Rx2n, Rx3n), Tn))
end function rot3axis

function translation_matrix_positiv(t) result(M)
	! The translation matrix
	implicit none
	
	integer :: ii, jj
	real, dimension(3, 1), intent(in) :: t
	real, dimension(4, 4) :: M
	
	do ii = 1,4
		do jj = 1,4
			M(ii,jj) = 0.0
			if (ii == jj) M (ii,jj) = 1.0
		end do
	end do
	M(1,4) = t(1,1)
	M(2,4) = t(2,1)
	M(3,4) = t(3,1)
end function translation_matrix_positiv

function translation_matrix_negativ(t) result(M)
	! The translation matrix
	implicit none
	
	integer :: ii, jj
	real, dimension(3,1), intent(in) :: t
	real, dimension(4,4) :: M
	
	do ii = 1,4
		do jj = 1,4
			M(ii,jj) = 0.0
			if (ii == jj) M (ii,jj) = 1.0
		end do
	end do
	M(1,4) = -t(1,1)
	M(2,4) = -t(2,1)
	M(3,4) = -t(3,1)
end function translation_matrix_negativ

function rotation_matrix_x2_positiv(b, d) result(R)
	! The translation matrix
	implicit none
	
	integer :: ii, jj
	real, intent(in) :: d
	real, dimension(3,1), intent(in) :: b
	real, dimension(4,4) :: R
	
	do ii = 1,4
		do jj = 1,4
			R(ii,jj) = 0.0
			if (ii == jj) R (ii,jj) = 1.0
		end do
	end do
	R(1,1) = b(3,1)
	R(1,3) = d
	R(3,1) = -d
	R(3,3) = b(3,1)
end function rotation_matrix_x2_positiv

function rotation_matrix_x2_negativ(b, d) result(R)
	! The translation matrix
	implicit none
	
	integer :: ii, jj
	real, intent(in) :: d
	real, dimension(3,1), intent(in) :: b
	real, dimension(4,4) :: R
	
	do ii = 1,4
		do jj = 1,4
			R(ii,jj) = 0.0
			if (ii == jj) R (ii,jj) = 1.0
		end do
	end do
	R(1,1) = b(3,1)
	R(1,3) = -d
	R(3,1) = d
	R(3,3) = b(3,1)
end function rotation_matrix_x2_negativ

function rotation_matrix_x3_positiv(b, d) result(R)
	! The translation matrix
	implicit none
	
	integer :: ii, jj
	real, intent(in) :: d
	real, dimension(3,1), intent(in) :: b
	real, dimension(4,4) :: R
	
	do ii = 1,4
		do jj = 1,4
			R(ii,jj) = 0.0
			if (ii == jj) R (ii,jj) = 1.0
		end do
	end do
	R(1,1) = (1/d) * b(1,1)
	R(1,2) = - (1/d) * b(2,1)
	R(2,1) = (1/d) * b(2,1)
	R(2,2) = (1/d) * b(1,1)
end function rotation_matrix_x3_positiv

function rotation_matrix_x3_negativ(b, d) result(R)
	! The translation matrix
	implicit none
	
	integer :: ii, jj
	real, intent(in) :: d
	real, dimension(3,1), intent(in) :: b
	real, dimension(4,4) :: R
	
	do ii = 1,4
		do jj = 1,4
			R(ii,jj) = 0.0
			if (ii == jj) R (ii,jj) = 1.0
		end do
	end do
	R(1,1) = (1/d) * b(1,1)
	R(1,2) = (1/d) * b(2,1)
	R(2,1) = - (1/d) * b(2,1)
	R(2,2) = (1/d) * b(1,1)
end function rotation_matrix_x3_negativ

function rotation_matrix_x3(alpha) result(R)
	! The translation matrix
	implicit none
	
	integer :: ii, jj
	real, intent(in) :: alpha
	real, dimension(4,4) :: R
	
	do ii = 1,4
		do jj = 1,4
			R(ii,jj) = 0.0
			if (ii == jj) R (ii,jj) = 1.0
		end do
	end do
	R(1,1) = cos(alpha)
	R(1,2) = -sin(alpha)
	R(2,1) = sin(alpha)
	R(2,2) = cos(alpha)
end function rotation_matrix_x3

function bresenham(P1, P2) result(XR1)
	implicit none
	
	integer :: ii, d, x1, y1, z1, x2, y2, z2, dx, dy, dz, ax, ay, az, sx, sy, sz, x, y, z, idx
	integer, dimension(3,1), intent(in) :: P1, P2
	! Todo
	integer, dimension(:), allocatable :: XR1, XR2, XR3
	
	real :: xd, yd, zd
	
	d = maxval(abs(P2 - P1))
	write(*,*) 'd', d
	
	! allocate memory 
	allocate(XR1(d))
	allocate(XR2(d))
	allocate(XR3(d))
	
	do ii = 1, d
		XR1(ii) = 0
	end do
	write(*,*) 'X1', X1
	XR2 = XR1
	XR3 = XR1
	
	x1 = P1(1,1)
	y1 = P1(2,1)
	z1 = P1(3,1)
	x2 = P2(1,1)
	y2 = P2(2,1)
	z2 = P2(3,1)
	dx = x2 - x1
	write(*,*) 'dx', dx
	dy = y2 - y1
	write(*,*) 'dy', dy
	dz = z2 - z1
	write(*,*) 'dz', dz
	ax = abs(dx) * 2
	write(*,*) 'ax', ax
	ay = abs(dy) * 2
	write(*,*) 'ay', ay
	az = abs(dz) * 2
	write(*,*) 'az', az
	sx = sign(1, dx)
	write(*,*) 'sx', sx
	sy = sign(1, dy)
	write(*,*) 'sy', sy
	sz = sign(1, dz)
	write(*,*) 'sz', sz
	x = x1
	y = y1
	z = z1
	idx = 1
	
	if (ax >= max(ay, az)) then
		! x1 dominant
		write(*,*) 'ax', ax, 'max(ay, az)', max(ay, az)
		yd = ay - (ax / 2)
		zd = az - (ax / 2)
		do while(.true.)
			XR1(idx) = x
			XR2(idx) = y
			XR3(idx) = z
			idx = idx + 1
			if (x == x2) exit
			if (yd >= 0) then
				y = y + sy
				yd = yd - ax
			end if
			if (zd >= 0) then
				z = z + sz
				zd = zd - ax
			end if
			x = x + sx
			yd = yd + ay
			zd = zd + az
		end do
	else if (ay >= max(ax, az)) then
		! x2 dominant
		write(*,*) 'ay'
	else if (az >= max(ax, ay)) then
		! x3 dominant
		write(*,*) 'az'
	end if
	write(*,*) 'XR1', XR1
	
end function bresenham