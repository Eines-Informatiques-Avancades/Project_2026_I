module initConfig
    ! Author: AdrianLLJ
    ! Test module to generate initial configurations ()
    use constants
    implicit none

    contains

    subroutine initDihedrals()
        ! Initialises the all the dihedrals to PI, to start from a planar polymer. 
        
        implicit none
        integer :: i

        do i = 1, NPART-3
            DANG(i) = PI
        end do
        
    end subroutine initDihedrals

    subroutine initPolymer()
        ! Initialises the polymer assuming constant bond length, angle and dihedral.

        implicit none
        integer :: imono
        double precision :: bvec1(3), bvec2(3), evec1(3), evec2(3), evec3(3)

        do imono = 1, NPART
            ! Initialise each monomer
            if (imono.eq.1) then
                X(imono) = 0.d0
                Y(imono) = 0.d0
                Z(imono) = 0.d0
            else if (imono.eq.2) then
                X(imono) = LBOND
                Y(imono) = 0.d0
                Z(imono) = 0.d0
            else if (imono.eq.3) then
                X(imono) = LBOND-LBOND*cos(BANG)
                Y(imono) = LBOND*sin(BANG)
                Z(imono) = 0.d0
            else
                ! Define bond vectors
                bvec1(1) = X(imono-2) - X(imono-3)
                bvec1(2) = Y(imono-2) - Y(imono-3)
                bvec1(3) = Z(imono-2) - Z(imono-3)

                bvec2(1) = X(imono-1) - X(imono-2)
                bvec2(2) = Y(imono-1) - Y(imono-2)
                bvec2(3) = Z(imono-1) - Z(imono-2)
                
                ! Create local ortonormal vectors
                evec1(:) = bvec2(:)/sqrt(bvec2(1)**2.d0 + bvec2(2)**2.d0 + bvec2(3)**2.d0)
                call crossProd(bvec1, bvec2, evec2)
                evec2(:) = evec2(:)/sqrt(evec2(1)**2.d0 + evec2(2)**2.d0 + evec2(3)**2.d0)
                call crossProd(evec2, evec1, evec3)

                ! Initialise new monomer's coordinates
                X(imono) = X(imono-1) + LBOND*sphToCart(evec1(1), evec2(1), evec3(1))
                Y(imono) = Y(imono-1) + LBOND*sphToCart(evec1(2), evec2(2), evec3(2))
                Z(imono) = Z(imono-1) + LBOND*sphToCart(evec1(3), evec2(3), evec3(3))
            end if
        end do

    end subroutine initPolymer

    subroutine shiftPolymer(displ)
        ! Translates the entire chain by a constant displacement. 
        ! Be careful not to throw particles out of the box!
        implicit none
        integer :: i 
        double precision, intent(in) :: displ

        do i = 1, NPART
            X(i) = X(i) + displ
            Y(i) = Y(i) + displ
            Z(i) = Z(i) + displ
        end do
    end subroutine shiftPolymer

    subroutine crossProd(a, b, c)
        ! Computes the cross product c = a x b for Cartesian coordinates.
        ! Be careful with the order of the input vectors!

        implicit none
        double precision, intent(in) :: a(3), b(3)
        double precision, intent(out) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine crossProd

    function sphToCart(evec1i, evec2i, evec3i) result(comp)
        ! Converts the component i of a vector and returns its
        ! expression in Cartesian coordinates from the local spherical basis vectors evec1,2,3.
        ! Assumes the same initial dihedral angle for all 4-particle groups

        implicit none
        double precision, intent(in) :: evec1i, evec2i, evec3i
        double precision :: comp

        comp = -cos(BANG)*evec1i + sin(BANG)*sin(DANG(1))*evec2i + sin(BANG)*cos(DANG(1))*evec3i

    end function sphToCart

    subroutine writeXYZ(filename, unitValue, atomSymbol, blankLine)
    ! Writes the configuration of the system in an XYZ format (to visualize in VMD). 
    implicit none
    integer, intent(in) :: unitValue
    integer :: i
    character(len=*), intent(in) :: filename, atomSymbol, blankLine

    open(unitValue, file=filename, status='unknown')
    write(unitValue, '(I3)') NPART
    write(unitValue, *) blankLine

    do i = 1, NPART
        write(unitValue,'(A,3F12.6)') atomSymbol, X(i), Y(i), Z(i)
    end do 

    close(unitValue)
    end subroutine writeXYZ
    

end module initConfig