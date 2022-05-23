module mod_aux
    use mod_constants
    contains
    subroutine compute_constants(rhol,pl,rhor,pr,Al,Bl,Ar,Br)
        implicit none
        real(8), intent(in)  :: rhol,pl,rhor,pr
        real(8), intent(out) :: Al,Bl,Ar,Br
        Al = 2.0d0/((gamma+1.0d0)*rhol)
        Ar = 2.0d0/((gamma+1.0d0)*rhor)
        Bl = ((gamma-1.0d0)/(gamma+1.0d0))*pl
        Br = ((gamma-1.0d0)/(gamma+1.0d0))*pr
    end subroutine compute_constants
    subroutine sample_profiles(ustar,pstar,S,rhol,ul,pl,cl,rhor,ur,pr,cr,rhos,us,ps)
        implicit none
        real(8), intent(in)  :: ustar,pstar,S
        real(8), intent(in)  :: rhol,ul,pl,cl,rhor,ur,pr,cr
        real(8), intent(out) :: rhos, us, ps
        real(8)              :: cml,cmr,pml,pmr,stl,str,shl,shr,sl,sr,cs

        ! Sample point lies to the left of contact discontinuity
        if (S .le. ustar) then
            ! Left rarefaction wave
            if (pstar .le. pl) then
                shl = ul-cl
                ! Sample point is left state
                if (S .le. shl) then
                    rhos = rhol
                    us = ul
                    ps = pl
                else
                    cml = cl*(pstar/pl)**g1
                    stl = ustar-cml
                    ! Sample point is star left state
                    if (S .gt. stl) then
                        rhos = rhol*(pstar/pl)**(1.0d0/gamma)
                        us = ustar
                        ps = pstar
                    ! Sample point inside left fan
                    else
                        us = g5*(cl+g7*ul+S)
                        cs = g5*(cl+g7*(ul-S))
                        rhos = rhol*((cs/cl)**g4)
                        ps = pl*((cs/cl)**g3)
                    end if
                end if
            ! Left shock
            else
                pml = pstar/pl
                sl = ul-cl*sqrt(g2*pml+g1)
                ! Sample pointt is left of data state
                if (S .le. sl) then
                    rhos = rhol
                    us = ul
                    ps = pl
                ! Sample point is star left state
                else
                    rhos = rhol*(pml+g6)/(pml*g6+1.0d0)
                    us = ustar
                    ps = pstar
                end if
            end if
        ! Sample point lies to the right of contact discontinuity
        else
            ! Right shock
            if (pstar .gt. pr) then
                pmr = pstar/pr
                sr = ur+cr*sqrt(g2*pmr+g1)
                ! Sample point is right data state
                if (S .ge. sr) then
                    rhos = rhor
                    us = ur
                    ps = pr
                ! Sample point is star right state
                else
                    rhos = rhor*(pmr+g6)/(pmr*g6+1.0d0)
                    us = ustar
                    ps = pstar
                end if
            ! Right rarefaction
            else
                shr = ur+cr
                ! Sample point is right data state
                if (S .ge. shr) then
                    rhos = rhor
                    us = ur
                    ps = pr
                else
                    cmr = cr*(pstar/pr)**g1
                    str = ustar+cmr
                    ! Sample point is star right state
                    if ( S .le. str ) then
                        rhos = rhor*(pstar/pr)**(1.0d0/gamma)
                        us = ustar
                        ps = pstar
                    ! Sample point is inside right fan
                    else
                        us = g5*(-cr+g7*ur+S)
                        cs = g5*(cr-g7*(ur-S))
                        rhos = rhor*((cs/cr)**g4)
                        ps = pr*((cs/cr)**g3)
                    end if
                end if
            end if
        end if
    end subroutine sample_profiles
end module mod_aux