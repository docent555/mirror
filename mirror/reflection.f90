module reflection

   implicit none

   integer(4) n
   real(8) L, xl, xr, pi

   real(8), allocatable :: sintr_mirr(:), trig(:), work(:), fft_mirr(:), fft_result(:), fft_work(:), r(:, :)

contains

   subroutine reflection_init(nn, LL, xll)

      implicit none

      integer(4) nn, err_alloc, k, i, j
      real(8) LL, xll

      n = nn
      L = LL
      xl = xll
      xr = L - xl

      pi = 2.0d0*dasin(1.0d0)

      allocate (sintr_mirr(n), trig(2*n), work(n), fft_work(2*n), fft_mirr(2*n), fft_result(2*n), r(n, n), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if

      do k = 1, n
         sintr_mirr(k) = 2.0d0/(k*pi)*(cos(k*pi*xl/L) - cos(k*pi*xr/L))
      end do
      sintr_mirr(:) = sintr_mirr/dsqrt(2.0d0/n)

      !make fourier coefficients from sine transform
      fft_mirr = 0.0d0
      fft_mirr(2:n) = -sintr_mirr(1:n - 1)
      fft_mirr(n + 2:2*n) = sintr_mirr(n - 1:1:-1)

      !matrix
      do i = 1, n
         do j = 1, n
            if ((i .eq. n) .or. (j .eq. n)) then
               r(i, j) = 0
            elseif (i == j) then
               r(i, i) = (2*i*pi*(-xl + xr) + L*sin((2*i*pi*xl)/L) - L*sin((2*i*pi*xr)/L))/(2.*L*i*pi)
            else
               r(i, j) = ((-sin(((j - i)*pi*xl)/L) + sin(((j - i)*pi*xr)/L))/(j - i) + &
                          (sin(((j + i)*pi*xl)/L) - sin(((j + i)*pi*xr)/L))/(j + i))/pi
               !r(i, j) = (2*(-(i*cos((i*pi*xl)/L)*sin((j*pi*xl)/L)) + j*cos((j*pi*xl)/L)*sin((i*pi*xl)/L) + &
               !                      i*cos((i*pi*xr)/L)*sin((j*pi*xr)/L) - j*cos((j*pi*xr)/L)*sin((i*pi*xr)/L)))/((j - i)*(j + i)*pi)
            end if
         end do
      end do

   end subroutine reflection_init

   subroutine ifft(fft_re, fft_im)

      implicit none

      real(8), intent(inout) :: fft_re(2*n), fft_im(2*n)
      integer(4) ifail

      call c06gcf(fft_im, 2*n, ifail)
      call c06fcf(fft_re, fft_im, 2*n, fft_work, ifail)
      call c06gcf(fft_im, 2*n, ifail)

   end subroutine ifft

   subroutine maggot(a, b)

      implicit none

      real(8), intent(in) :: a(n)
      real(8), intent(out) :: b(n)
      integer(4) i, j

      b = 0.0d0
      do i = 1, n
         do j = 1, n
            b(i) = b(i) + r(j,i)*a(j)
         end do
      end do

   end subroutine maggot

end module reflection
