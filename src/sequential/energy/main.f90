program main
    ! Dummy main for personal tests
    use constants
    use initConfig
    use nonBonded
    use bonded
    use energy

    implicit none
    integer :: i

    call shiftLenJon()
    print*, 'RC, ECUT:', RC, ECUT

    call initDihedrals()
    call initPolymer()
    call shiftPolymer(0.5d0*HBOX)
    call writeXYZ('initPolymer.xyz', 10, 'C', '0')
    
    call enerNonBond(Enb)
    call enerBonded(Eb)
    call totEnergy(En)
    print*, "Enb, Eb, En:", Enb, Eb, En
end program main