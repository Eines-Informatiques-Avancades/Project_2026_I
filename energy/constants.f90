module constants
    ! Module with all the constants the energy computations may need.
    implicit none

    ! System paramters
    integer, parameter :: NPART=500
    double precision, parameter :: TEMP=3.d0, BOX=505.d0, HBOX=0.5d0*BOX
    double precision :: X(NPART), Y(NPART), Z(NPART)

    ! Molecular parameters 
    double precision, parameter :: LBOND=1.d0, BANG=1.954d0
    double precision :: DANG(NPART-3)

    ! Physical constants
    double precision :: PI=4.d0*atan(1.d0)

end module constants