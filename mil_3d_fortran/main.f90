program main
	implicit none
	integer, parameter :: out_unit=20
	integer :: numberOfOrientations, increment, r, ii
	real, parameter :: pi = 3.141593
	real, dimension(3) :: p0, p1, n
	real :: theta, phi, thetaR, phiR, calculate_mil_3d, MIL
	
	! Number of randomly generated orientation (positiv integer, minimum 9)
	numberOfOrientations = 400
	! Distance between two created lines (positiv integer)
	increment = 5
	
	!
	r = 1
	
	do ii = 1, numberOfOrientations
		call random_number(thetaR)
		theta = thetaR * pi
		call random_number(phiR)
		phi = phiR * 2 * pi
		! Origin
		p0 = (/0, 0, 0 /)
		! Spherical coordinates to Cartesian coordinates
		p1(1) = r * sin(theta) * cos(phi)
		p1(2) = r * sin(theta) * sin(phi)
		p1(3) = r * cos(theta)
		
		n = (1 / sqrt(sum((p1 - p0)**2))) * (p1 - p0)
		MIL = calculate_mil_3d(n, 1, 2, 3)
		!write(*,*) 'ii: ',ii , '/ ', numberOfOrientations, 'theta: ', theta, 'phi: ', phi, 'MIL: ', MIL
	end do
	
	
end program	main

function calculate_mil_3d(n, a, b, c)
	integer :: a, b, c
	integer, dimension(3) :: ps, pQ1
    REAL    :: calculate_mil_3d, h, cv, dr, ru
	real, dimension(3) :: n, pr1, pr2, m, v0
	
	h = 0.0
	cv = 0.0
		
	if (n(1) == 1 .and. n(1) == 0 .and. n(1) == 0) then
		print *, '100'
	else if (n(1) == 0.0 .and. n(1) == 1.0 .and. n(1) == 0.0) then
		print *, '010'
	else
		! Room diagonal
		dr = sqrt(real(a)**2 + real(b)**2 + real(c)**2)
		! Radius of the sphere (Half room diagonal)
		ru = 0.5 * dr
		! Origin of the sphere (center of the read data)
		m = (/ a/2, b/2, c/2 /)
		! Point on the sphere in direction n
		ps = nint(m + ru * n)
		! First directional vector (direction 1, 0)
		v0(1) = 1
		v0(2) = 0
		v0(3) = -n(1) / n(3)
		v0 = (1 / sqrt(sum(v0**2))) * v0
		
		pQ1 = nint(ps + ru * v0)
		write(*,*) pQ1
	end if

    calculate_mil_3d = h / cv
 end function