module energy
    ! ---
    ! This module defines and initialises all the variables related to the energy
    ! as well as all the functions/subroutines for computing the  total energy of the system. 
    ! The current version only considers a C skeleton.
    ! ---
    use constants
    use nonBonded
    use bonded
    implicit none

    double precision :: En = 0.d0

    contains

    subroutine totEnergy(En)
        ! Computes the total energy of the system, considering both the bonded and nonbonded interactions.
        ! Note that the calls to enerBonded and enerNonBonded are internal, thus the *global* variables of 
        ! the bonded/non-bonded energies won't be updated.
        ! If you want, you can call them separatedly in the main program.
        implicit none
        double precision, intent(out) :: En
        double precision :: Eb, Enb

        En = 0.d0
        call enerBonded(Eb)
        call enerNonBond(Enb)
        En = Eb + Enb
    end subroutine totEnergy

end module energy