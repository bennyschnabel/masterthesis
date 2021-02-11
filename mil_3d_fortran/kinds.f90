!==============================================================================
!> \file kinds.f90
!> Module with standard precision definitions
!>
!> \author Johannes Gebert
!> \date 04.01.2021
!> \date 19.01.2021
MODULE KINDS

IMPLICIT NONE

INTEGER, PARAMETER :: sik    = 2    ! INTEGER Kind
INTEGER, PARAMETER :: ik     = 8    ! INTEGER Kind
INTEGER, PARAMETER :: rk     = 8    ! Real Kind
INTEGER, PARAMETER :: mcl    = 512  ! Maximal character length
INTEGER, PARAMETER :: scl    = 64   ! Short character length
!-- Mpi-specific kinds
INTEGER, PARAMETER :: mik = 4    ! MPI INTEGER Kind; Compile with corresponding mpi!!

END MODULE KINDS
