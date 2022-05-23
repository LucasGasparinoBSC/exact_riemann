module mod_constants
    implicit none
    integer(4), parameter :: maxIter=100
    real(8),    parameter :: eps=1.0e-8
    real(8),    parameter :: gamma=1.4d0
    real(8),    parameter :: g1=(gamma-1.0d0)/(2.0d0*gamma)
    real(8),    parameter :: g2=(gamma+1.0d0)/(2.0d0*gamma)
    real(8),    parameter :: g3=(2.0d0*gamma)/(gamma-1.0d0)
    real(8),    parameter :: g4=(2.0d0)/(gamma-1.0d0)
    real(8),    parameter :: g5=(2.0d0)/(gamma+1.0d0)
    real(8),    parameter :: g6=(gamma-1.0d0)/(gamma+1.0d0)
    real(8),    parameter :: g7=(gamma-1.0d0)/(2.0d0)
    real(8),    parameter :: g8=(gamma-1.0d0)
    real(8),    parameter :: tout=0.25d0
end module mod_constants
