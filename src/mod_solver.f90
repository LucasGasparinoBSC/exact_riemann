module mod_solver

    use mod_functions
    use mod_constants

    contains

    subroutine starpu(rhol,ul,pl,cl,Al,Bl,rhor,ur,pr,cr,Ar,Br,du,ustar,pstar)

        implicit none

        real(8), intent(in)  :: rhol,ul,pl,cl,Al,Bl,rhor,ur,pr,cr,Ar,Br,du
        real(8), intent(out) :: ustar, pstar
        integer(4)           :: iter
        real(8)              :: fl,fld,fr,frd
        real(8)              :: pold, dp
        
        ! Compute pstar guess
        call guess_pstar(rhol,ul,pl,cl,rhor,ur,pr,cr,du,pstar)
        pold = pstar

        ! Iterate using a  Newton method
        write(6,*) '--------------------------------------------------------'
        write(6,*) '     Iteration             Change in pstar'
        write(6,*) '--------------------------------------------------------'
        iterations: do iter = 1,maxIter
            call side_state(pold,rhol,pl,cl,Al,Bl,fl,fld)
            call side_state(pold,rhor,pr,cr,Ar,Br,fr,frd)
            pstar = pold-((fl+fr+du)/(fld+frd))
            dp = 2.0d0*abs((pstar-pold)/(pstar+pold))
            write(6,30) iter, dp
            if (dp .le. eps) then
                exit iterations
            end if
            if (pstar .lt. 0.0d0) then
                pstar = eps
            end if
            pold = pstar
        end do iterations

        ! Compute u in the star region
        ustar = 0.5d0*(ul+ur+fr-fl)
        write(6,*) '--------------------------------------------------------'
        write(6,*) '     Results in the star region'
        write(6,*) ' '
        write(6,*) '     Pressure           Velocity'
        write(6,40) pstar, ustar
        write(6,*) '--------------------------------------------------------'

        ! Formats
30      format(5X,I5,15X,F12.7)
40      format(2(F14.6, 5X))

    end subroutine starpu

end module mod_solver