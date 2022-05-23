module mod_functions

    use mod_constants

    contains

    subroutine side_state(pstar,rhok,pk,ck,Ak,Bk,fk,fkd)

        implicit none

        real(8), intent(in)  :: pstar, rhok, pk, ck, Ak, Bk
        real(8), intent(out) :: fk, fkd
        real(8)              :: qrt

        if (pstar .gt. pk) then
            qrt = sqrt(Ak/(pstar+Bk))
            fk = (pstar-pk)*qrt
            fkd = (1.0d0-0.5d0*((pstar-pk)/(pstar+Bk)))*qrt
        else if (pstar .le. pk) then
            fk = g4*ck*(((pstar/pk)**g1)-1.0d0)
            fkd = (1.0d0/(ck*rhok))*((pstar/pk)**(-g2))
        end if

    end subroutine side_state

    subroutine guess_pstar(rhol,ul,pl,cl,rhor,ur,pr,cr,du,pstar)

        implicit none

        real(8), intent(in)  :: rhol, ul, pl, cl, rhor, ur, pr, cr, du
        real(8), intent(out) :: pstar
        real(8)              :: pmin, pmax, qmax, ppv, cup, quser
        real(8)              :: pq, um, ptl, ptr, gel, ger

        quser = 2.0d0
        cup = 0.25d0*(rhol+rhor)*(cl+cr)
        ppv = 0.5d0*(pl+pr)+0.5d0*(ul-ur)*cup
        ppv = max(0.0d0,ppv)
        pmin = min(pl,pr)
        pmax = max(pl,pr)
        qmax = pmax/pmin

        if (qmax .le. quser .and. (pmin .le. ppv .and. ppv .le. pmax)) then
            pstar = ppv
        else
            if (ppv .lt. pmin) then
                pq = (pl/pr)**g1
                um = (pq*(ul/cl)+(ur/cr)+g4*(pq-1.0d0))/((pq/cl)+(1.0d0/cr))
                ptl = 1.0d0+g7*(ul-um)/cl
                ptr = 1.0d0+g7*(um-ur)/cr
                pstar = 0.5d0*(pl*(ptl**g3)+pr*(ptr**g3))
            else
                gel = sqrt((g5/rhol)/(g6*pl+ppv))
                ger = sqrt((g5/rhor)/(g6*pr+ppv))
                pstar = (gel*pl+ger*pr-du)/(gel+ger)
            end if
        end if

    end subroutine guess_pstar

end module mod_functions