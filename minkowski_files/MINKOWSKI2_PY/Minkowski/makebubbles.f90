
program bubbles
  
  !       program to produce a single sphere.
  implicit none
  integer:: n,i,j,k, m 
  real:: boxsize, radius, x, y, z, r
  real(kind=8) :: one, zero
  
  m=0
  n=300
  one=1.
  zero=0.
  boxsize=100
  radius=5.0
  open(unit=6,file='sphere.bin',form='BINARY',status='unknown')
  write(6) 12 
  write(6) n,n,n
  write(6) 12
  write(*,*) 'radius=',radius 
  
  ! Write record
  write(6) 4*n*n*n
  
  do i=1,n
     do j=1,n
        do k=1,n
           
           x=(real(i)/real(n)-0.5)*boxsize
           y=(real(j)/real(n)-0.5)*boxsize
           z=(real(k)/real(n)-0.5)*boxsize
           r=sqrt(x**2+y**2+z**2)
           if(r.lt.radius) then
              write(6) one
              m=m+1
           else
              write(6) zero
              m=m+1
           endif
        enddo
     enddo
  enddo
  ! Write record
  write(6) 4*n*n*n
  close(6)
  write(*,*) (m+7)*4

end program bubbles
