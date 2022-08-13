program mirr
    use, intrinsic :: iso_c_binding
    implicit none

    integer(4), parameter :: n = 512
    real(8), parameter :: L = 1, dx = L/n, xl = 0.35*L, xr = L - xl, pi = dacos(-1.0d0)
    integer(4) ifail, i, j, il, ip, in
    real(8) c(n), mir(n), trig(2*n), work(n), r(n, n), mir_st(n), x, imir(n)

    !коэффициенты фунции которую будем отражать
    c = 0
    c(1) = 20
    c(3) = 19
    c(5) = 21
    c(7) = 7
    c(9) = 3

    il = xl/dx
    mir = 1.0d0
    mir(1:il) = 0.0d0
    mir(n - il:n) = 0.0d0

    call c06haf(1, n, c, 'initial', trig, work, ifail)

    !аналитические коэффициенты
    do i = 1, n
        mir_st(i) = 2.0d0/(pi*i)*(dcos(pi*i/L*xl) - dcos(pi*i/L*(1 - xl)))
    end do
    mir_st = mir_st/dsqrt(2.0d0/n)

    !обратное преобразование из аналитики
    imir(:) = mir_st
    call c06haf(1, n, imir, 'subsequent', trig, work, ifail)

    !матрица
    do i = 0, n - 1
        do j = 0, n - 1
            if ((i .eq. 0) .or. (j .eq. 0)) then
                r(i + 1, j + 1) = 0
            elseif (i == j) then
                r(i + 1, i + 1) = (2*i*pi*(-xl + xr) + L*sin((2*i*pi*xl)/L) - L*sin((2*i*pi*xr)/L))/(2.*L*i*pi)               
            else
                r(i + 1, j + 1) = ((-sin(((j - i)*pi*xl)/L) + sin(((j - i)*pi*xr)/L))/(j - i) + &
                                   (sin(((j + i)*pi*xl)/L) - sin(((j + i)*pi*xr)/L))/(j + i))/pi
            end if
        end do
    end do

    open (1, file='test.dat')
    do i = 1, n
        write (1, '(i,5f12.7)') i-1, c(i), mir(i), mir_st(i), imir(i), r(i,10)
    end do
    close (1)

end program mirr
