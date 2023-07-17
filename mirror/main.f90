program mirror
   use, intrinsic :: iso_c_binding
   use reflection
   use ifport
   implicit none

   integer(4), parameter :: nn = 512
   real(8) dx
   integer(4) i, j, il, ifail
  real(8) aak(nn), ak(nn), a(nn), af(nn), a2(2*nn), x, fsin(nn), s(nn), ffsin(2*nn), ss(2*nn), b(nn), i_sintr_mir(nn), b_fbf(4*nn), c, m(nn)
   integer(c_int) hours1, minutes1, seconds1, hours2, minutes2, seconds2
   real(c_double) start_time1, stop_time1, calc_time1, start_time2, stop_time2, calc_time2

   !call reflection_init(nn, 1.0d0, 0.35d0)
   call reflection_init(nn, 10d0, 2d0)

   dx = L/n

   do i = 1, n
      fsin(i) = dsin(5*pi/L*(i - 1)*dx)
   end do

   s(1:n - 1) = fsin(2:n)
   s(n) = 0

   ffsin(1:n) = fsin
   ffsin(n + 2:2*n) = -fsin(n:2:-1)

   call c06haf(1, n, s, 'initial', trig_haf, work_haf, ifail)

   !call haf2faf(s, ss)
   !ss(:) = ffsin
   !call c06faf(ss, 2*n, work_faf, ifail)
   !call c06gbf(ss, 2*n, ifail)
   !call c06fbf(ss, 2*n, work_faf, ifail)

   !open (1, file='test1.dat')
   !do i = 1, n
   !   write (1, '(i,2f18.7)') i, s(i), fsin(i)
   !end do
   !close (1)
   !
   !open (1, file='test2.dat')
   !do i = 1, 2*n
   !   write (1, '(i,2f18.7)') i, ss(i), ffsin(i)
   !end do
   !close (1)
   !stop

   !coefficients of the function that we will reflect
   ak = 0
   ak(1) = 25
   ak(3) = 17
   ak(5) = 21
   ak(7) = 12
   ak(9) = 7
   ak(11) = 16
   !ak(63) = 1
   !ak(787) = 10
   !ak(2047) = 3

   !mirror
   !il = xl/dx
   !mir = 1.0d0
   !mir(1:il) = 0.0d0
   !mir(n - il:n) = 0.0d0

   a(:) = ak
   !a2 = 0
   !a2(1:n) = dsqrt(2.0d0)*ak
   call c06haf(1, n, a, 'initial', trig_haf, work_haf, ifail)
   !call c06haf(1, n2, a2, 'initial', trig2_haf, work2_haf, ifail)

   !call c06haf(1, n, mirr_haf, 'initial', trig_haf, work_haf, ifail)
   !call c06haf(1, n2, mirr2_haf, 'initial', trig2_haf, work2_haf, ifail)

   call reflect_faf(ak, af)
   call reflect_haf(ak, a)
   call maggot(ak, m)

   !call haf2faf(ak, a_faf)
   !call c06gbf(a_faf, 4*n, ifail)
   !c = sqrt(dble(n)/512.0d0)
   !a_faf = c*a_faf
   !call c06fbf(a_faf, 4*n, work_faf, ifail)
   !
   !call haf2faf(mirr_haf, mirr_faf)
   !mirr_faf = dsqrt(2.0d0)*mirr_faf
   !call c06gbf(mirr_faf, 4*n, ifail)
   !call c06fbf(mirr_faf, 4*n, work_faf, ifail)
   !
   !b_fbf(:) = a_faf*mirr_faf
   !call c06faf(b_fbf, 4*n, work_faf, ifail) ! fourier coeff. of reflected wave
   !
   !call c06gbf(b_fbf, 4*n, ifail)
   !call c06fbf(b_fbf, 4*n, work_faf, ifail)

   !open (1, file='test.dat')
   !do i = 1, n
   !   write (1, '(i,2f18.7)') i, mirr_haf(i), mirr2_haf(2*i)
   !end do
   !close (1)
   
   call c06haf(1, n, m, 'subsequent', trig_haf, work_haf, ifail)
   call c06haf(1, n, a, 'subsequent', trig_haf, work_haf, ifail)
   call c06haf(1, n, af, 'subsequent', trig_haf, work_haf, ifail)
   
   open(1, file = 'test.dat') 
   do i = 1,n
      write(1, '(4e17.8)') (i-1)*dx, m(i), a(i), af(i)
   end do		
   close(1)
   			
  
end program mirror
