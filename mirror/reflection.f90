module reflection

   implicit none

   integer(4) n, n2, n3, n4
   real(8) L, xl, xr, pi

   real(8), allocatable :: mirr_haf(:), trig_haf(:), work_haf(:), work_faf(:), &
                           mirr_faf(:), r(:, :), a_faf(:), mirr2_haf(:), trig2_haf(:), work2_haf(:), mirr2(:), mirr(:), &
                           b2(:), b4(:)

   !private trig_haf, work_haf

contains

   subroutine reflection_init(nn, LL, xll)

      implicit none

      integer(4) nn, err_alloc, k, i, j, ifail
      real(8) LL, xll

      n = nn
      n2 = 2*n
      n3 = 3*n
      n4 = 4*n
      L = LL
      xl = xll
      xr = L - xl

      mirr_faf = 0
      mirr2_haf = 0

      pi = 2.0d0*dasin(1.0d0)

      allocate (mirr_haf(n), trig_haf(2*n), work_haf(n), work_faf(4*n), mirr_faf(4*n), r(n, n), &
                a_faf(4*n), trig2_haf(n4), work2_haf(n2), mirr2(2*n), mirr(n), b2(n2), b4(n4), mirr2_haf(n2), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if

      do k = 1, n
         mirr_haf(k) = 2.0d0/(k*pi)*(cos(k*pi*xl/L) - cos(k*pi*xr/L))
      end do
      !mirr_haf(:) = mirr_haf/dsqrt(2.0d0/n)
      mirr_haf(:) = mirr_haf*dsqrt(n/2.0d0)
      mirr2_haf(1:n) = dsqrt(2.0d0)*mirr_haf

      call haf2faf(mirr_haf, mirr_faf)

      !reverse transformation from analytics
      mirr(:) = mirr_haf
      call c06haf(1, n, mirr, 'initial', trig_haf, work_haf, ifail)
      !double n
      mirr2 = 0.0d0
      mirr2(1:n) = dsqrt(2.0d0)*mirr_haf
      call c06haf(1, n2, mirr2, 'initial', trig2_haf, work2_haf, ifail)

      !matrix
      do i = 1, n
         do j = 1, n
            if ((i .eq. n) .or. (j .eq. n)) then
               r(i, j) = 0
            elseif (i == j) then
               r(i, i) = (2*i*pi*(-xl + xr) + L*dsin((2*i*pi*xl)/L) - L*dsin((2*i*pi*xr)/L))/(2.*L*i*pi)
            else
               r(i, j) = ((-dsin(((j - i)*pi*xl)/L) + dsin(((j - i)*pi*xr)/L))/(j - i) + &
                          (dsin(((j + i)*pi*xl)/L) - dsin(((j + i)*pi*xr)/L))/(j + i))/pi
               !r(i, j) = (2*(-(i*cos((i*pi*xl)/L)*sin((j*pi*xl)/L)) + j*cos((j*pi*xl)/L)*sin((i*pi*xl)/L) + &
               !                      i*cos((i*pi*xr)/L)*sin((j*pi*xr)/L) - j*cos((j*pi*xr)/L)*sin((i*pi*xr)/L)))/((j - i)*(j + i)*pi)
            end if
         end do
      end do

      call haf2faf(mirr_haf, mirr_faf)
      mirr_faf = dsqrt(2.0d0)*mirr_faf
      call c06gbf(mirr_faf, 4*n, ifail)
      call c06fbf(mirr_faf, 4*n, work_faf, ifail)

      !open (1, file='test.dat')
      !do i = 1, n
      !   write (1, '(i,3f18.7)') i, mirr2(2*i), mirr(i), mirr_faf(2*i + 1)
      !end do
      !close (1)
      
      !print *, pi, L, xl, xr
      !
      !open (1, file='r.dat')
      !do i = 1, n
      !   do j = 1, n
      !      write (1, '(f18.7,\)') r(i, j)
      !   end do
      !   write (1, '(/)')
      !end do
      !close (1)
      !
      !stop

   end subroutine reflection_init

   !subroutine ifft_odd(fft_re, fft_im)
   !
   !   implicit none
   !
   !   real(8), intent(inout) :: fft_re(4*n), fft_im(4*n)
   !   integer(4) ifail
   !
   !   fft_re = 0.0d0
   !   call c06gcf(fft_im, n4, ifail)
   !   call c06fcf(fft_re, fft_im, n4, work_fcf, ifail)
   !   call c06gcf(fft_im, n4, ifail)
   !
   !end subroutine ifft_odd

   subroutine maggot(a, b)

      implicit none

      real(8), intent(in) :: a(n)
      real(8), intent(out) :: b(n)
      integer(4) i, j

      !call f06pcf ('l', n, 1, r, n, a, 1, 0, b, 1)
      !call F06PAF('n', n, n, 1, r, n, a, 1, 0, b, 1)

      b = 0.0d0
      do i = 1, n
         do j = 1, n
            b(i) = b(i) + r(i, j)*a(j)
         end do
      end do

   end subroutine maggot

   subroutine reflect_faf(a, b)

      implicit none

      integer ifail
      real(8), intent(in) :: a(n)
      real(8), intent(inout) :: b(n)

      call haf2faf(a, b4)
      call c06gbf(b4, n4, ifail)
      call c06fbf(b4, n4, work_faf, ifail)

      b4(:) = mirr_faf*b4
      b4(n2 + 1:) = -b4(n2 + 1:)

      call c06faf(b4, n4, work_faf, ifail) !fourier coeff. of reflected wave

      call faf2haf(b4, b) !sine transform coeff. of reflected wave

      !call c06haf(1, n, b, 'subsequent', trig_haf, work_haf, ifail)
      !
      !open (1, file='test.dat')
      !do i = 1, n
      !   write (1, '(i,f18.7)') i, b(i)
      !end do
      !stop

   end subroutine reflect_faf

   subroutine reflect_haf(a, b)

      implicit none

      integer i, ifail
      real(8), intent(in) :: a(n)
      real(8), intent(inout) :: b(n)

      b2 = 0.0d0
      b2(1:n) = dsqrt(2.0d0)*a

      call c06haf(1, n2, b2, 'subsequent', trig2_haf, work2_haf, ifail)

      b2(:) = mirr2*b2

      call c06haf(1, n2, b2, 'subsequent', trig2_haf, work2_haf, ifail)

      b(1:n - 1) = b2(1:n - 1)/dsqrt(2.0d0)
      b(n) = b2(n2)

      !open (1, file='test.dat')
      !do i = 1, n2
      !   write (1, '(i,f18.7)') i, b2(i)
      !end do
      !stop

   end subroutine reflect_haf

   subroutine haf2faf(s, f)
      implicit none

      real(8), intent(in) :: s(n)
      real(8), intent(out) :: f(n4)

      !make fourier coefficients from sine ones
      f = 0.0d0
      f(n4:n3 + 2:-1) = -s(1:n - 1)

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

end module reflection
