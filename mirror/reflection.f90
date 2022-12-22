module reflection

   implicit none

   integer(4) n, n2, n3, n4
   real(8) L, xl, xr, pi

   real(8), allocatable :: mirr_haf(:), trig_haf(:), work_haf(:), work_faf(:), mirr_fcf(:), mirr_faf(:), result_fcf(:), work_fcf(:), r(:, :), a_faf(:), a_fcf(:)

   !private trig_haf, work_haf

contains

   subroutine reflection_init(nn, LL, xll)

      implicit none

      integer(4) nn, err_alloc, k, i, j
      real(8) LL, xll

      n = nn
      n2 = 2*n
      n3 = 3*n
      n4 = 4*n
      L = LL
      xl = xll
      xr = L - xl

      a_fcf = 0.0d0
      mirr_fcf = 0.0d0

      pi = 2.0d0*dasin(1.0d0)

      allocate (mirr_haf(n), trig_haf(2*n), work_haf(n), work_faf(4*n), work_fcf(4*n), mirr_fcf(4*n), mirr_faf(4*n), result_fcf(4*n), r(n, n), &
                a_faf(4*n), a_fcf(4*n), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if

      do k = 1, n
         mirr_haf(k) = 2.0d0/(k*pi)*(cos(k*pi*xl/L) - cos(k*pi*xr/L))
      end do
      mirr_haf(:) = mirr_haf/dsqrt(2.0d0/n)

      call haf2fcf(mirr_haf, mirr_fcf(1:2*n))

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

   subroutine ifft_odd(fft_re, fft_im)

      implicit none

      real(8), intent(inout) :: fft_re(4*n), fft_im(4*n)
      integer(4) ifail

      fft_re = 0.0d0
      call c06gcf(fft_im, n4, ifail)
      call c06fcf(fft_re, fft_im, n4, work_fcf, ifail)
      call c06gcf(fft_im, n4, ifail)

   end subroutine ifft_odd

   subroutine maggot(a, b)

      implicit none

      real(8), intent(in) :: a(n)
      real(8), intent(out) :: b(n)
      integer(4) i, j

      b = 0.0d0
      do i = 1, n
         do j = 1, n
            b(i) = b(i) + r(i, j)*a(j)
         end do
      end do

   end subroutine maggot

   subroutine reflect(a, b)

      implicit none

      integer i, ifail
      real(8), intent(in) :: a(n)
      real(8), intent(inout) :: b(n)

      call haf2fcf(a, a_fcf(1:2*n))

      b(:) = a_fcf(1:n)

   end subroutine reflect

   subroutine haf2fcf(s, f)
      implicit none

      real(8), intent(in) :: s(n)
      real(8), intent(out) :: f(n4)

      !make fourier coefficients from sine transform
      f = 0.0d0
      f(2:n) = -s(1:n - 1)
      f(3*n + 2:n4) = s(n - 1:1:-1)
   end subroutine haf2fcf

   subroutine haf2faf(s, f)
      implicit none

      real(8), intent(in) :: s(n)
      real(8), intent(out) :: f(n4)

      !make fourier coefficients from sine ones
      f = 0.0d0
      f(n4:n3 + 2:-1) = -s(1:n - 1)*dsqrt(2.0d0)

      !open(1, file='test.dat')
      !do i =1,n2
      !   write(1, '(i,f18.7)') i, f(i)
      !enddo
      !stop
   end subroutine haf2faf
   
   subroutine faf2haf(f, s)
      implicit none

      real(8), intent(out) :: s(n)
      real(8), intent(in) :: f(n4)

      !make fourier coefficients from sine ones      
      s(1:n - 1) = -f(n4:n3 + 2:-1)

      !open(1, file='test.dat')
      !do i =1,n2
      !   write(1, '(i,f18.7)') i, f(i)
      !enddo
      !stop
   end subroutine faf2haf

   subroutine fourier2sine(f, s)
      implicit none

      real(8), intent(in) :: f(n4)
      real(8), intent(out) :: s(n)

      !make sine transform  from fourier coefficients
      s = 0.0d0
      s(1:n - 1) = -f(2:n)

   end subroutine fourier2sine

end module reflection
