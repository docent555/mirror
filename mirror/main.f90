program mirr
   use, intrinsic :: iso_c_binding
   use reflection
   implicit none

   integer(4), parameter :: nn = 1024
   real(8) dx
   integer(4) i, j, il
   real(8) ak(nn), i_c(nn), mir(nn), trig(2*nn), trig4(8*nn), work(nn), work4(4*nn), x, &
      i_sintr_mir(nn), bk(nn), b(nn), sk(2*nn), s(2*nn), d(nn), &
      sintr_mir4(4*nn), i_sintr_mir4(4*nn), s4(4*nn), sk4(4*nn)

   call reflect_init(nn, 0.35d0, 2.0d0)

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
   
   sintr_mir2(1:n) = sintr_mir/dsqrt(1.0d0/n)
   sintr_mir2(n+1:2*n) = 0
   
   sintr_mir4(1:n) = sintr_mir/dsqrt(0.5d0/n)
   sintr_mir4(n+1:4*n) = 0
   
   sintr_mir = sintr_mir/dsqrt(2.0d0/n) ! it seems that analytical transform multiplied by 2/n by default         

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

   call maggot(ak, bk)
   
   !reflected by maggot
   b(:) = bk
   call c06haf(1, n, b, 'subsequent', trig, work, ifail)

   !now with sine transform   
   sk(1:n) = ak   
   sk(n+1:2*n) = 0
   
   !reverse transformation from analytics on 2*n
   i_sintr_mir2(:) = sintr_mir2
   call c06haf(1, 2*n, i_sintr_mir2, 'initial', trig2, work2, ifail)
   
   s(:) = sk
   call c06haf(1, 2*n, s, 'subsequent', trig2, work2, ifail)
   
   sk = s*i_sintr_mir2
   call c06haf(1, 2*n, sk, 'subsequent', trig2, work2, ifail)      
   
   do i=1,n
      d(i)=bk(i)-sk(i)
   enddo
   
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
      i_sintr_mir(i) = i_sintr_mir(i-1)      
   end do
   b(1) = 0
   i_sintr_mir(1) = 0      
   do i = 2*n, 2, -1      
      i_sintr_mir2(i) = i_sintr_mir2(i-1)     
      s(i)=s(i-1)      
   end do   
   i_sintr_mir2(1) = 0
   do i = 4*n, 2, -1      
      i_sintr_mir4(i) = i_sintr_mir4(i-1)         
      s4(i)=s4(i-1)
   end do   
   i_sintr_mir4(1) = 0

   open (1, file='test1.dat')
   do i = 1, n
      !write (1, '(i,5f12.7)') i - 1, i_c(i), mir(i), sintr_mir(i), i_sintr_mir(i), b(i)      
      write (1, '(2f12.7)') (i - 1)*dx, bk(i)
   end do
   close (1)
      
   open (1, file='test2.dat')
   do i = 1, 2*n      
      write (1, '(2f12.7)') (i - 1)*dx/2, i_sintr_mir2(i)      
   end do
   close (1)
   
   open (1, file='test4.dat')
   do i = 1, 4*n      
      write (1, '(2f12.7)') (i - 1)*dx/4, i_sintr_mir4(i)   
   end do
   close (1)

end program mirr
