program exact_riemann

    use mod_aux
    use mod_constants
    use mod_solver

    implicit none

    integer(4) :: ncells, icell
    real(8)    :: L, dx, cpos, xpos
    real(8)    :: S
    real(8)    :: rhol, ul, pl, cl, Al, Bl
    real(8)    :: rhor, ur, pr, cr, Ar, Br
    real(8)    :: rhos, us, ps
    real(8)    :: du
    real(8)    :: ustar, pstar

    ! Define domain
    ncells = 100
    L = 1.0d0
    cpos = L/2.0d0
    dx = L/dble(ncells)

    ! Define left and right states
    rhol = 1.0d0;
    ul = 0.0d0;
    pl = 1.0d0;
    rhor = 0.125d0;
    ur = 0.0d0;
    pr = 0.1d0;

    ! Compute difference between left and right state velocities
    du = ur-ul

    ! Compute speed of sound for both states
    cl = sqrt(gamma*pl/rhol);
    cr = sqrt(gamma*pr/rhor);

    ! Compute the Riemann invariants
    call compute_constants(rhol,pl,rhor,pr,Al,Bl,Ar,Br)

    ! Start iterative process
    call starpu(rhol,ul,pl,cl,Al,Bl,rhor,ur,pr,cr,Ar,Br,du,ustar,pstar)

    ! Sample the solution and write to a file
    open(UNIT=2,FILE='exact.out',STATUS='unknown')
    do icell = 1,ncells
        xpos = (dble(icell)-0.5d0)*dx
        S = (xpos-cpos)/tout
        call sample_profiles(ustar,pstar,S,rhol,ul,pl,cl,rhor,ur,pr,cr,rhos,us,ps)
        write(2,20) xpos, rhos, us, ps
    end do
    close(2)
20  format(5(F14.6, 2X))

end program exact_riemann