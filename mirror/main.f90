program mirr
   use, intrinsic :: iso_c_binding
   use reflection
   use ifport
   implicit none

   integer(4), parameter :: nn = 1024
   real(8) dx
   integer(4) i, j, il
   real(8) ak(nn), i_c(nn), mir(nn), trig4(8*nn), work4(4*nn), x, &
      bk1(nn), bk2(nn), b(nn), d(nn), &
      sintr_mir4(4*nn), i_sintr_mir4(4*nn), s4(4*nn), sk4(4*nn), s1(nn), sk1(nn)
   integer(c_int) hours1, minutes1, seconds1, hours2, minutes2, seconds2
   real(c_double) start_time1, stop_time1, calc_time1, start_time2, stop_time2, calc_time2

   call reflect_init(nn, 0.35d0, 1.0d0)

   dx = L/n

   !coefficients of the function that we will mark
   ak = 0
   ak(1) = 25
   ak(3) = 17
   ak(5) = 21
   ak(7) = 12
   ak(9) = 7
   ak(11) = 16

   il = xl/dx
   mir = 1.0d0
   mir(1:il) = 0.0d0
   mir(n - il:n) = 0.0d0

   i_c(:) = ak
   call c06haf(1, n, i_c, 'initial', trig, work, ifail)         
   
   sintr_mir4(1:n) = sintr_mir/dsqrt(0.5d0/n)
   sintr_mir4(n+1:4*n) = 0
   
            

   !reverse transformation from analytics
   i_sintr_mir(:) = sintr_mir
   call c06haf(1, n, i_sintr_mir, 'subsequent', trig, work, ifail)

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
   !
   !bk = 0
   !!multiply the matrix by the sine transform of the signal
   !do i = 1, n
   !   do j = 1, n
   !      bk(i) = bk(i) + r(i, j)*ak(j)
   !   end do
   !end do
   
   !reflected by maggot
   !start_time1 = dclock()
   !do i=1,100000
   call maggot(ak, bk1)
   !enddo
   !stop_time1 = dclock()
   !
   !calc_time1 = stop_time1 - start_time1
   !
   !hours1 = calc_time1/3600
   !minutes1 = (calc_time1 - hours1*3600)/60
   !seconds1 = calc_time1 - hours1*3600 - minutes1*60
   !
   !!reflected by reflect
   !start_time2 = dclock()
   !do i=1,100000
   call reflect(ak, bk2)
   !enddo
   !stop_time2 = dclock()
   !
   !calc_time2 = stop_time2 - start_time2
   !
   !hours2 = calc_time2/3600
   !minutes2 = (calc_time2 - hours2*3600)/60
   !seconds2 = calc_time2 - hours2*3600 - minutes2*60
   !
   !write (*, '(/)')
   !print *, 'Calcualting maggot took:', hours1, 'h :', minutes1, 'm :', seconds1, 's'
   !print *, 'Calcualting reflect took:', hours2, 'h :', minutes2, 'm :', seconds2, 's'
      
   b(:) = bk1
   call c06haf(1, n, b, 'subsequent', trig, work, ifail)      
   
   d(:) = bk2
   call c06haf(1, n, d, 'subsequent', trig, work, ifail)
   
   !reverse transformation from analytics on n
   sk1(:) = ak   
   
   i_sintr_mir(:) = sintr_mir
   call c06haf(1, n, i_sintr_mir, 'initial', trig, work, ifail)
   
   s1(:) = sk1
   call c06haf(1, n, s1, 'subsequent', trig, work, ifail)
   
   sk1(:) = s1*i_sintr_mir
   call c06haf(1, n, sk1, 'subsequent', trig, work, ifail)  
   
   !reverse transformation from analytics on 4*n
   sk4(1:n) = ak   
   sk4(n+1:4*n) = 0
   
   i_sintr_mir4(:) = sintr_mir4
   call c06haf(1, 4*n, i_sintr_mir4, 'initial', trig4, work4, ifail)
   
   s4(:) = sk4
   call c06haf(1, 4*n, s4, 'subsequent', trig4, work4, ifail)
   
   sk4(:) = s4*i_sintr_mir4
   call c06haf(1, 4*n, sk4, 'subsequent', trig4, work4, ifail)      
   
   !shifting
   do i = n, 2, -1
      b(i) = b(i - 1)
      d(i) = d(i - 1)
      i_sintr_mir(i) = i_sintr_mir(i-1)      
   end do
   b(1) = 0
   d(1) = 0
   i_sintr_mir(1) = 0        
   !do i = 4*n, 2, -1      
   !   i_sintr_mir4(i) = i_sintr_mir4(i-1)         
   !   s4(i)=s4(i-1)
   !end do   
   !i_sintr_mir4(1) = 0

   open (1, file='test1.dat')
   do i = 1, n
      !write (1, '(i,5f12.7)') i - 1, i_c(i), mir(i), sintr_mir(i), i_sintr_mir(i), b(i)      
      write (1, '(4f12.7)') (i - 1)*dx, sk1(i), bk1(i), sk1(i)-bk1(i)
   end do
   close (1)
        
   
   !open (1, file='test4.dat')
   !do i = 1, 4*n      
   !   write (1, '(2f12.7)') (i - 1)*dx/4, i_sintr_mir4(i)   
   !end do
   !close (1)

end program mirr
