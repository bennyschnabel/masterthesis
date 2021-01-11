program main
	implicit none
	integer, parameter :: out_unit=20
	integer :: numberOfOrientations, increment, r, ii
	real, parameter :: pi = 3.141593
	real, dimension(3) :: p0, p1, n
	real :: theta, phi, thetaR, phiR, h
	
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
		write(*,*) 'ii: ',ii, 'n: ', n
	end do
	
	call calculate_mil_3d(n)
end program	main