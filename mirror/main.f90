program mirror
   use, intrinsic :: iso_c_binding
   use reflection
   use ifport
   implicit none

   integer(4), parameter :: nn = 1024
   real(8) dx
   integer(4) i, j, il, ifail
   real(8) ak(nn), i_c(nn), mirr(nn), trig4(8*nn), work4(4*nn), x, fsin(nn), fsin2(2*nn), fft_fiction(2*nn), fsintr(nn), fsintr2(2*nn), s(nn)
   integer(c_int) hours1, minutes1, seconds1, hours2, minutes2, seconds2
   real(c_double) start_time1, stop_time1, calc_time1, start_time2, stop_time2, calc_time2

   call reflection_init(nn, 1.0d0, 0.35d0)

   dx = L/n

   do i = 1, n - 1
      fsin(i) = dsin(pi/L*i*dx)
   end do

   s(1) = 0.0d0
   s(2:n) = fsin(1:n - 1)

   !!with fourer
   !fsin2 = 0.0d0
   !fsin2(2:n) = fsin(1:n - 1)
   !fsin2(n + 2:2*n) = -fsin(n - 1:1:-1)
   !!fourier transform
   !call c06fcf(fsin2, fft_fiction, 2*n, work, ifail)
   !!backward fourier
   !call c06gcf(fft_fiction, 2*n, ifail)
   !call c06fcf(fsin2, fft_fiction, 2*n, work, ifail)
   !call c06gcf(fft_fiction, 2*n, ifail)
   !
   !!with sine transorm
   !fsintr(:) = fsin
   !call c06haf(1, n, fsintr, 'initial', trig, work, ifail)
   !!from sine transform to fourier transform
   !fft_result = 0.0d0
   !fft_result(2:n) = -fsintr(1:n - 1)
   !fft_result(n + 2:2*n) = fsintr(n - 1:1:-1)
   !!!backward fourier
   !!call c06gcf(fft_result, 2*n, ifail)
   !!call c06fcf(fsintr2, fft_result, 2*n, work, ifail)
   !!call c06gcf(fft_result, 2*n, ifail)
   !!from fourier to sine
   !fsin(1:n-1) = -fft_result(2:n)
   !fsin(n) = 0
   !call c06haf(1, n, fsin, 'subsequent', trig, work, ifail)

   !coefficients of the function that we will reflect
   ak = 0
   ak(1) = 25
   ak(3) = 17
   ak(5) = 21
   ak(7) = 12
   ak(9) = 7
   ak(11) = 16

   !mirror
   !il = xl/dx
   !mir = 1.0d0
   !mir(1:il) = 0.0d0
   !mir(n - il:n) = 0.0d0

   i_c(:) = ak
   call c06haf(1, n, i_c, 'initial', trig, work, ifail)

   !sintr_mir4(1:n) = sintr_mir/dsqrt(0.5d0/n)
   !sintr_mir4(n+1:4*n) = 0

   !reverse transformation from analytics
   !i_sintr_mir(:) = sintr_mir
   !call c06haf(1, n, i_sintr_mir, 'subsequent', trig, work, ifail)

   !!matrix
   !do i = 1, n
   !   do j = 1, n
   !      if ((i .eq. n) .or. (j .eq. n)) then
   !         r(i, j) = 0
   !      elseif (i == j) then
   !         r(i, i) = (2*i*pi*(-xl + xr) + L*sin((2*i*pi*xl)/L) - L*sin((2*i*pi*xr)/L))/(2.*L*i*pi)
   !      else
   !         r(i, j) = ((-sin(((j - i)*pi*xl)/L) + sin(((j - i)*pi*xr)/L))/(j - i) + &
   !                    (sin(((j + i)*pi*xl)/L) - sin(((j + i)*pi*xr)/L))/(j + i))/pi
   !         !r(i, j) = (2*(-(i*cos((i*pi*xl)/L)*sin((j*pi*xl)/L)) + j*cos((j*pi*xl)/L)*sin((i*pi*xl)/L) + &
   !         !                      i*cos((i*pi*xr)/L)*sin((j*pi*xr)/L) - j*cos((j*pi*xr)/L)*sin((i*pi*xr)/L)))/((j - i)*(j + i)*pi)
   !      end if
   !   end do
   !end do

   !call c06ecf(fft_result, fft_mirr, 2*n, ifail)

   fft_result = 0.0d0   
   call ifft(fft_result, fft_mirr)

   mirr(:) = sintr_mirr
   call c06haf(1, n, mirr, 'initial', trig, work, ifail)

   open (1, file='test.dat')
   do i = 1, n
      !write (1, '(3f12.7)') (i - 1)*dx, fsin2(i), fft_fiction(i)
      !write (1, '(3E18.7)') (i - 1)*dx, fft_result(i), fft_mirr(i)
      write (1, '(2E18.7)') (i - 1)*dx, i_c(i)
   end do
   close (1)

   open (1, file='test2.dat')
   write (1, '(2E18.7)') 0, 0
   do i = 1, n - 1
      !write (1, '(2f12.7)') i*dx, fsin(i)
      write (1, '(2E18.7)') i*dx, mirr(i)
   end do
   close (1)

end program mirror
