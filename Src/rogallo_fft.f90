module rogallo_fft_for_fargo
   complex(8), allocatable, dimension(:) :: trig_table  ! dimension = nphi
end module rogallo_fft_for_fargo

! FFT routines

module cache_data
  integer, parameter :: CACHE = 32768   !4096	! -- cache work buffer size in complex words--
  type scomplex
     real(8), dimension(0:CACHE-1) :: r, i
  end type scomplex
end module cache_data

!====================  cache load/store   =====================

Subroutine load_cache( n, lot, data, workr, worki, inc, jump )
  implicit none
  integer :: n, lot
  complex(8), dimension(0:0) :: data
  real(8), dimension(0:0) :: workr, worki
  integer :: inc, jump 

  integer :: i, j

if ( jump < inc ) then
   do i=0,n-1
      do j=0,lot-1
         workr(i*lot+j) = real(data(i*inc+j*jump))
         worki(i*lot+j) = aimag(data(i*inc+j*jump))
      end do
   end do
else
   do j=0,lot-1
      do i=0,n-1
         workr(i*lot+j) = real(data(i*inc+j*jump)) 
         worki(i*lot+j) = aimag(data(i*inc+j*jump))
      end do
   end do
end if

End Subroutine load_cache

Subroutine load_cache_dri( n, lot, data, workr, worki, inc, jump )
  implicit none
  integer :: n, lot
  real(8), dimension(0:0) :: data
  real(8), dimension(0:0) :: workr, worki
  integer :: inc, jump 

  integer :: i, j

if ( jump < inc ) then
   do i=0,n-1
      do j=0,lot-1
         workr(i*lot+j) = data(2*i*inc    +j*jump)
         worki(i*lot+j) = data(2*i*inc+inc+j*jump)
      end do
   end do
else
   do j=0,lot-1
      do i=0,n-1
         workr(i*lot+j) = data(2*i*inc    +j*jump) 
         worki(i*lot+j) = data(2*i*inc+inc+j*jump)
      end do
   end do
end if

End Subroutine load_cache_dri


Subroutine store_cache(n, lot, workr, worki, data, inc, jump )
  implicit none
  integer :: n, lot
  real(8), dimension(0:0) :: workr, worki
  complex(8), dimension(0:0) :: data
  integer :: inc, jump 

  integer :: i, j

if ( jump < inc ) then
   do i=0,n-1
      do j=0,lot-1
         data(i*inc+j*jump) = cmplx(workr(i*lot+j), worki(i*lot+j), kind=8)
      end do
   end do
else	
   do j=0,lot-1
      do i=0,n-1
         data(i*inc+j*jump) = cmplx(workr(i*lot+j), worki(i*lot+j), kind=8)
      end do
   end do
end if

End Subroutine store_cache

Subroutine store_cache_dri(n, lot, workr, worki, data, inc, jump )
  implicit none
  integer :: n, lot
  real(8), dimension(0:0) :: workr, worki
  real(8), dimension(0:0) :: data
  integer :: inc, jump 

  integer :: i, j

if ( jump < inc ) then
   do i=0,n-1
      do j=0,lot-1
         data(2*i*inc    +j*jump) = workr(i*lot+j)
         data(2*i*inc+inc+j*jump) = worki(i*lot+j)
      end do
   end do
else	
   do j=0,lot-1
      do i=0,n-1
         data(2*i*inc    +j*jump) = workr(i*lot+j)
         data(2*i*inc+inc+j*jump) = worki(i*lot+j)
      end do
   end do
end if

End Subroutine store_cache_dri


!====================  complex --> complex butterfies  =====================

!----------- radix-8 ---------------

Subroutine cache_radix_8(  &
        p   &		! p factor 
       ,m   &		! m factor 
       ,lot & 		! vector length 
       ,srcr, srci	&	! source address 
       ,dstr, dsti	&	! dest. address 
       ,e		&	! trig table address 
       ,isign		)	! forward/inverse direction 
  implicit none
  integer :: p	&		! p factor 
  ,  m		&		! m factor 
  ,  lot  	&		! vector length 
  ,  isign			! forward/inverse direction 
  real(8), dimension(0:0) :: srcr, srci,  &		! source address 
       dstr, dsti		! dest. address 
  complex(8), dimension(0:0) :: e			! trig table address 

  real(8) ::	      c   , er1 , ei1 , er2 , ei2 , er3 , ei3 , er4 , ei4 , er5 , ei5 , er6 , ei6 , er7 , ei7  &
       , rr0 , ir0 , rr1 , ir1 , rr2 , ir2 , rr3 , ir3 , rr4 , ir4 , rr5 , ir5 , rr6 , ir6 , rr7 , ir7  &
       , rt0 , it0 , rt1 , it1 , rt2 , it2 , rt3 , it3 , rt4 , it4 , rt5 , it5 , rt6 , it6 , rt7 , it7 

  integer :: r , v 

  c = sqrt(.5d0) 

  if ( isign > 0 )	then
     r = 0 
     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*8 ) ;				it0 = srci(               v + lot*p*r*8 ) 
        rt1 = srcr( lot*p   +     v + lot*p*r*8 ) ;				it1 = srci( lot*p   +     v + lot*p*r*8 ) 
        rt2 = srcr( lot*p*2 +     v + lot*p*r*8 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*8 ) 
        rt3 = srcr( lot*p*3 +     v + lot*p*r*8 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*8 ) 
        rt4 = srcr( lot*p*4 +     v + lot*p*r*8 ) ;				it4 = srci( lot*p*4 +     v + lot*p*r*8 ) 
        rt5 = srcr( lot*p*5 +     v + lot*p*r*8 ) ;				it5 = srci( lot*p*5 +     v + lot*p*r*8 ) 
        rt6 = srcr( lot*p*6 +     v + lot*p*r*8 ) ;				it6 = srci( lot*p*6 +     v + lot*p*r*8 ) 
        rt7 = srcr( lot*p*7 +     v + lot*p*r*8 ) ;				it7 = srci( lot*p*7 +     v + lot*p*r*8 ) 

        rr0 = rt0 + rt4 ;				ir0 = it0 + it4 
        rr4 = rt0 - rt4 ;				ir4 = it0 - it4 
        rr1 = rt1 + rt5 ;				ir1 = it1 + it5 
        rr5 = rt1 - rt5 ;				ir5 = it1 - it5 
        rr2 = rt2 + rt6 ;				ir2 = it2 + it6 
        rr6 = rt2 - rt6 ;				ir6 = it2 - it6 
        rr3 = rt3 + rt7 ;				ir3 = it3 + it7 
        rr7 = rt3 - rt7 ;				ir7 = it3 - it7 

        rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 
        rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 
        rt1 = rr1 + rr3 ;				it1 = ir1 + ir3 
        rt3 = rr1 - rr3 ;				it3 = ir1 - ir3 
        rt4 = rr4 - ir6 ;				it4 = ir4 + rr6 
        rt6 = rr4 + ir6 ;				it6 = ir4 - rr6 
        rt5 = rr5 + rr7 ;				it5 = ir5 + ir7 
        rt7 = rr5 - rr7 ;				it7 = ir5 - ir7 

        rr5 = rt7 - it5 ;				ir5 = it7 + rt5 
        rr7 = rt7 + it5 ;				ir7 = it7 - rt5 
        rt5 = c * rr5 ;					it5 = c * ir5 
        rt7 = c * rr7 ;					it7 = c * ir7 

        rr0 = rt0 + rt1 ;				ir0 = it0 + it1 
        rr4 = rt0 - rt1 ;				ir4 = it0 - it1 
        rr1 = rt4 + rt5 ;				ir1 = it4 + it5 
        rr5 = rt4 - rt5 ;				ir5 = it4 - it5 
        rr2 = rt2 - it3 ;				ir2 = it2 + rt3 
        rr6 = rt2 + it3 ;				ir6 = it2 - rt3 
        rr3 = rt6 - rt7 ;				ir3 = it6 - it7 
        rr7 = rt6 + rt7 ;				ir7 = it6 + it7 

        dstr(               v + lot*p*r ) = rr0 ;				dsti(               v + lot*p*r ) = ir0 
        dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;				dsti( lot*m*p   +   v + lot*p*r ) = ir1 
        dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 
        dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 
        dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;				dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 
        dstr( lot*m*p*5 +   v + lot*p*r ) = rr5 ;				dsti( lot*m*p*5 +   v + lot*p*r ) = ir5 
        dstr( lot*m*p*6 +   v + lot*p*r ) = rr6 ;				dsti( lot*m*p*6 +   v + lot*p*r ) = ir6 
        dstr( lot*m*p*7 +   v + lot*p*r ) = rr7 ;				dsti( lot*m*p*7 +   v + lot*p*r ) = ir7 
     end do

     do r=1,m-1

        er1 = real(e(p*r)) ;		        ei1 = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;		ei2 = aimag(e(p*r*2)) ;
        er3 = real(e(p*r*3)) ;		ei3 = aimag(e(p*r*3)) ;
        er4 = real(e(p*r*4)) ;		ei4 = aimag(e(p*r*4)) ;
        er5 = real(e(p*r*5)) ;		ei5 = aimag(e(p*r*5)) ;
        er6 = real(e(p*r*6)) ;		ei6 = aimag(e(p*r*6)) ;
        er7 = real(e(p*r*7)) ;		ei7 = aimag(e(p*r*7)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*8 ) ;				it0 = srci(               v + lot*p*r*8 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*8 ) ;				it1 = srci( lot*p   +     v + lot*p*r*8 ) ;
           rt2 = srcr( lot*p*2 +     v + lot*p*r*8 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*8 ) ;
           rt3 = srcr( lot*p*3 +     v + lot*p*r*8 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*8 ) ;
           rt4 = srcr( lot*p*4 +     v + lot*p*r*8 ) ;				it4 = srci( lot*p*4 +     v + lot*p*r*8 ) ;
           rt5 = srcr( lot*p*5 +     v + lot*p*r*8 ) ;				it5 = srci( lot*p*5 +     v + lot*p*r*8 ) ;
           rt6 = srcr( lot*p*6 +     v + lot*p*r*8 ) ;				it6 = srci( lot*p*6 +     v + lot*p*r*8 ) ;
           rt7 = srcr( lot*p*7 +     v + lot*p*r*8 ) ;				it7 = srci( lot*p*7 +     v + lot*p*r*8 ) ;

           rr1 = er1 * rt1 ;				ir1 = ei1 * rt1 ;
           rr2 = er2 * rt2 ;				ir2 = ei2 * rt2 ;
           rr3 = er3 * rt3 ;				ir3 = ei3 * rt3 ;
           rr4 = er4 * rt4 ;				ir4 = ei4 * rt4 ;
           rr5 = er5 * rt5 ;				ir5 = ei5 * rt5 ;
           rr6 = er6 * rt6 ;				ir6 = ei6 * rt6 ;
           rr7 = er7 * rt7 ;				ir7 = ei7 * rt7 ;

           rt1 = rr1 - ei1 * it1 ;			it1 = ir1 + er1 * it1 ;
           rt2 = rr2 - ei2 * it2 ;			it2 = ir2 + er2 * it2 ;
           rt3 = rr3 - ei3 * it3 ;			it3 = ir3 + er3 * it3 ;
           rt4 = rr4 - ei4 * it4 ;			it4 = ir4 + er4 * it4 ;
           rt5 = rr5 - ei5 * it5 ;			it5 = ir5 + er5 * it5 ;
           rt6 = rr6 - ei6 * it6 ;			it6 = ir6 + er6 * it6 ;
           rt7 = rr7 - ei7 * it7 ;			it7 = ir7 + er7 * it7 ;

           rr0 = rt0 + rt4 ;				ir0 = it0 + it4 ;
           rr4 = rt0 - rt4 ;				ir4 = it0 - it4 ;
           rr1 = rt1 + rt5 ;				ir1 = it1 + it5 ;
           rr5 = rt1 - rt5 ;				ir5 = it1 - it5 ;
           rr2 = rt2 + rt6 ;				ir2 = it2 + it6 ;
           rr6 = rt2 - rt6 ;				ir6 = it2 - it6 ;
           rr3 = rt3 + rt7 ;				ir3 = it3 + it7 ;
           rr7 = rt3 - rt7 ;				ir7 = it3 - it7 ;

           rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 ;
           rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 ;
           rt1 = rr1 + rr3 ;				it1 = ir1 + ir3 ;
           rt3 = rr1 - rr3 ;				it3 = ir1 - ir3 ;
           rt4 = rr4 - ir6 ;				it4 = ir4 + rr6 ;
           rt6 = rr4 + ir6 ;				it6 = ir4 - rr6 ;
           rt5 = rr5 + rr7 ;				it5 = ir5 + ir7 ;
           rt7 = rr5 - rr7 ;				it7 = ir5 - ir7 ;

           rr5 = rt7 - it5 ;				ir5 = it7 + rt5 ;
           rr7 = rt7 + it5 ;				ir7 = it7 - rt5 ;
           rt5 = c * rr5 ;					it5 = c * ir5 ;
           rt7 = c * rr7 ;					it7 = c * ir7 ;

           rr0 = rt0 + rt1 ;				ir0 = it0 + it1 ;
           rr4 = rt0 - rt1 ;				ir4 = it0 - it1 ;
           rr1 = rt4 + rt5 ;				ir1 = it4 + it5 ;
           rr5 = rt4 - rt5 ;				ir5 = it4 - it5 ;
           rr2 = rt2 - it3 ;				ir2 = it2 + rt3 ;
           rr6 = rt2 + it3 ;				ir6 = it2 - rt3 ;
           rr3 = rt6 - rt7 ;				ir3 = it6 - it7 ;
           rr7 = rt6 + rt7 ;				ir7 = it6 + it7 ;

           dstr(               v + lot*p*r ) = rr0 ;				dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;				dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
           dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 ;
           dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;				dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 ;
           dstr( lot*m*p*5 +   v + lot*p*r ) = rr5 ;				dsti( lot*m*p*5 +   v + lot*p*r ) = ir5 ;
           dstr( lot*m*p*6 +   v + lot*p*r ) = rr6 ;				dsti( lot*m*p*6 +   v + lot*p*r ) = ir6 ;
           dstr( lot*m*p*7 +   v + lot*p*r ) = rr7 ;				dsti( lot*m*p*7 +   v + lot*p*r ) = ir7 ;
        end do
     end do

  else	
     r = 0 ;
     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*8 ) ;				it0 = srci(               v + lot*p*r*8 ) ;
        rt1 = srcr( lot*p   +     v + lot*p*r*8 ) ;				it1 = srci( lot*p   +     v + lot*p*r*8 ) ;
        rt2 = srcr( lot*p*2 +     v + lot*p*r*8 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*8 ) ;
        rt3 = srcr( lot*p*3 +     v + lot*p*r*8 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*8 ) ;
        rt4 = srcr( lot*p*4 +     v + lot*p*r*8 ) ;				it4 = srci( lot*p*4 +     v + lot*p*r*8 ) ;
        rt5 = srcr( lot*p*5 +     v + lot*p*r*8 ) ;				it5 = srci( lot*p*5 +     v + lot*p*r*8 ) ;
        rt6 = srcr( lot*p*6 +     v + lot*p*r*8 ) ;				it6 = srci( lot*p*6 +     v + lot*p*r*8 ) ;
        rt7 = srcr( lot*p*7 +     v + lot*p*r*8 ) ;				it7 = srci( lot*p*7 +     v + lot*p*r*8 ) ;

        rr0 = rt0 + rt4 ;				ir0 = it0 + it4 ;
        rr4 = rt0 - rt4 ;				ir4 = it0 - it4 ;
        rr1 = rt1 + rt5 ;				ir1 = it1 + it5 ;
        rr5 = rt1 - rt5 ;				ir5 = it1 - it5 ;
        rr2 = rt2 + rt6 ;				ir2 = it2 + it6 ;
        rr6 = rt2 - rt6 ;				ir6 = it2 - it6 ;
        rr3 = rt3 + rt7 ;				ir3 = it3 + it7 ;
        rr7 = rt3 - rt7 ;				ir7 = it3 - it7 ;

        rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 ;
        rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 ;
        rt1 = rr1 + rr3 ;				it1 = ir1 + ir3 ;
        rt3 = rr1 - rr3 ;				it3 = ir1 - ir3 ;
        rt6 = rr4 - ir6 ;				it6 = ir4 + rr6 ;
        rt4 = rr4 + ir6 ;				it4 = ir4 - rr6 ;
        rt5 = rr5 + rr7 ;				it5 = ir5 + ir7 ;
        rt7 = rr5 - rr7 ;				it7 = ir5 - ir7 ;

        rr7 = rt7 - it5 ;				ir7 = it7 + rt5 ;
        rr5 = rt7 + it5 ;				ir5 = it7 - rt5 ;
        rt5 = c * rr5 ;					it5 = c * ir5 ;
        rt7 = c * rr7 ;					it7 = c * ir7 ;

        rr0 = rt0 + rt1 ;				ir0 = it0 + it1 ;
        rr4 = rt0 - rt1 ;				ir4 = it0 - it1 ;
        rr1 = rt4 + rt5 ;				ir1 = it4 + it5 ;
        rr5 = rt4 - rt5 ;				ir5 = it4 - it5 ;
        rr6 = rt2 - it3 ;				ir6 = it2 + rt3 ;
        rr2 = rt2 + it3 ;				ir2 = it2 - rt3 ;
        rr3 = rt6 - rt7 ;				ir3 = it6 - it7 ;
        rr7 = rt6 + rt7 ;				ir7 = it6 + it7 ;

        dstr(               v + lot*p*r ) = rr0 ;				dsti(               v + lot*p*r ) = ir0 ;
        dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;				dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
        dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
        dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 ;
        dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;				dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 ;
        dstr( lot*m*p*5 +   v + lot*p*r ) = rr5 ;				dsti( lot*m*p*5 +   v + lot*p*r ) = ir5 ;
        dstr( lot*m*p*6 +   v + lot*p*r ) = rr6 ;				dsti( lot*m*p*6 +   v + lot*p*r ) = ir6 ;
        dstr( lot*m*p*7 +   v + lot*p*r ) = rr7 ;				dsti( lot*m*p*7 +   v + lot*p*r ) = ir7 ;
     end do

     do r=1,m-1

        er1 = real(e(p*r)) ;		ei1 = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;		ei2 = aimag(e(p*r*2)) ;
        er3 = real(e(p*r*3)) ;		ei3 = aimag(e(p*r*3)) ;
        er4 = real(e(p*r*4)) ;		ei4 = aimag(e(p*r*4)) ;
        er5 = real(e(p*r*5)) ;		ei5 = aimag(e(p*r*5)) ;
        er6 = real(e(p*r*6)) ;		ei6 = aimag(e(p*r*6)) ;
        er7 = real(e(p*r*7)) ;		ei7 = aimag(e(p*r*7)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*8 ) ;				it0 = srci(               v + lot*p*r*8 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*8 ) ;				it1 = srci( lot*p   +     v + lot*p*r*8 ) ;
           rt2 = srcr( lot*p*2 +     v + lot*p*r*8 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*8 ) ;
           rt3 = srcr( lot*p*3 +     v + lot*p*r*8 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*8 ) ;
           rt4 = srcr( lot*p*4 +     v + lot*p*r*8 ) ;				it4 = srci( lot*p*4 +     v + lot*p*r*8 ) ;
           rt5 = srcr( lot*p*5 +     v + lot*p*r*8 ) ;				it5 = srci( lot*p*5 +     v + lot*p*r*8 ) ;
           rt6 = srcr( lot*p*6 +     v + lot*p*r*8 ) ;				it6 = srci( lot*p*6 +     v + lot*p*r*8 ) ;
           rt7 = srcr( lot*p*7 +     v + lot*p*r*8 ) ;				it7 = srci( lot*p*7 +     v + lot*p*r*8 ) ;

           rr1 = er1 * rt1 ;				ir1 = ei1 * rt1 ;
           rr2 = er2 * rt2 ;				ir2 = ei2 * rt2 ;
           rr3 = er3 * rt3 ;				ir3 = ei3 * rt3 ;
           rr4 = er4 * rt4 ;				ir4 = ei4 * rt4 ;
           rr5 = er5 * rt5 ;				ir5 = ei5 * rt5 ;
           rr6 = er6 * rt6 ;				ir6 = ei6 * rt6 ;
           rr7 = er7 * rt7 ;				ir7 = ei7 * rt7 ;

           rt1 = rr1 + ei1 * it1 ;			it1 = -ir1 + er1 * it1 ;
           rt2 = rr2 + ei2 * it2 ;			it2 = -ir2 + er2 * it2 ;
           rt3 = rr3 + ei3 * it3 ;			it3 = -ir3 + er3 * it3 ;
           rt4 = rr4 + ei4 * it4 ;			it4 = -ir4 + er4 * it4 ;
           rt5 = rr5 + ei5 * it5 ;			it5 = -ir5 + er5 * it5 ;
           rt6 = rr6 + ei6 * it6 ;			it6 = -ir6 + er6 * it6 ;
           rt7 = rr7 + ei7 * it7 ;			it7 = -ir7 + er7 * it7 ;

           rr0 = rt0 + rt4 ;				ir0 = it0 + it4 ;
           rr4 = rt0 - rt4 ;				ir4 = it0 - it4 ;
           rr1 = rt1 + rt5 ;				ir1 = it1 + it5 ;
           rr5 = rt1 - rt5 ;				ir5 = it1 - it5 ;
           rr2 = rt2 + rt6 ;				ir2 = it2 + it6 ;
           rr6 = rt2 - rt6 ;				ir6 = it2 - it6 ;
           rr3 = rt3 + rt7 ;				ir3 = it3 + it7 ;
           rr7 = rt3 - rt7 ;				ir7 = it3 - it7 ;

           rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 ;
           rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 ;
           rt1 = rr1 + rr3 ;				it1 = ir1 + ir3 ;
           rt3 = rr1 - rr3 ;				it3 = ir1 - ir3 ;
           rt6 = rr4 - ir6 ;				it6 = ir4 + rr6 ;
           rt4 = rr4 + ir6 ;				it4 = ir4 - rr6 ;
           rt5 = rr5 + rr7 ;				it5 = ir5 + ir7 ;
           rt7 = rr5 - rr7 ;				it7 = ir5 - ir7 ;

           rr7 = rt7 - it5 ;				ir7 = it7 + rt5 ;
           rr5 = rt7 + it5 ;				ir5 = it7 - rt5 ;
           rt5 = c * rr5 ;					it5 = c * ir5 ;
           rt7 = c * rr7 ;					it7 = c * ir7 ;

           rr0 = rt0 + rt1 ;				ir0 = it0 + it1 ;
           rr4 = rt0 - rt1 ;				ir4 = it0 - it1 ;
           rr1 = rt4 + rt5 ;				ir1 = it4 + it5 ;
           rr5 = rt4 - rt5 ;				ir5 = it4 - it5 ;
           rr6 = rt2 - it3 ;				ir6 = it2 + rt3 ;
           rr2 = rt2 + it3 ;				ir2 = it2 - rt3 ;
           rr3 = rt6 - rt7 ;				ir3 = it6 - it7 ;
           rr7 = rt6 + rt7 ;				ir7 = it6 + it7 ;

           dstr(               v + lot*p*r ) = rr0 ;				dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;				dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
           dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 ;
           dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;				dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 ;
           dstr( lot*m*p*5 +   v + lot*p*r ) = rr5 ;				dsti( lot*m*p*5 +   v + lot*p*r ) = ir5 ;
           dstr( lot*m*p*6 +   v + lot*p*r ) = rr6 ;				dsti( lot*m*p*6 +   v + lot*p*r ) = ir6 ;
           dstr( lot*m*p*7 +   v + lot*p*r ) = rr7 ;				dsti( lot*m*p*7 +   v + lot*p*r ) = ir7 ;
        end do
     end do
  end if
End Subroutine cache_radix_8

!----------- radix-3 --------------

Subroutine cache_radix_3(  &
  p   &		! p factor 
       ,m   &		! m factor 
       ,lot & 		! vector length 
       ,srcr, srci	&	! source address 
       ,dstr, dsti	&	! dest. address 
       ,e		&	! trig table address 
       ,isign		)	! forward/inverse direction 
  implicit none
  integer :: p		&	! p factor 
  ,  m			&	! m factor 
  ,  lot  		&	! vector length 
  ,  isign			! forward/inverse direction 
  real(8), dimension(0:0) :: srcr, srci,  &		! source address 
       dstr, dsti		! dest. address 
  complex(8), dimension(0:0) :: e			! trig table address 

  real(8), parameter :: cr = (-.5d0), ci = (sqrt(3.d0)/2.d0)

  real(8) ::  er1, ei1, er2, ei2, rt0, rt1, rt2, it0, it1, it2, rr0, rr1, rr2, ir0, ir1, ir2 ;

  integer :: r , v ;

  if ( isign > 0 )	then
     r = 0 ;
     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*3 ) ;						it0 = srci(               v + lot*p*r*3 ) ;
        rt1 = srcr( lot*p   +     v + lot*p*r*3 ) ;						it1 = srci( lot*p   +     v + lot*p*r*3 ) ;
        rt2 = srcr( lot*p*2 +     v + lot*p*r*3 ) ;						it2 = srci( lot*p*2 +     v + lot*p*r*3 ) ;

        rr1 = rt1 + rt2 ;						ir1 = it1 + it2 ;
        rt2 = rt1 - rt2 ;						it2 = it1 - it2 ;
        rr0 = rt0 + rr1 ;						ir0 = it0 + ir1 ;
        rr1 = rt0 + cr * rr1 ;					ir1 = it0 + cr * ir1 ;
        rr2 = rr1 + ci * it2 ;					ir2 = ir1 - ci * rt2 ;
        rr1 = rr1 - ci * it2 ;					ir1 = ir1 + ci * rt2 ;

        dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
        dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
        dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
     end do

     do r=1,m-1
        er1  = real(e(p*r)) ;							ei1  = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;							ei2 = aimag(e(p*r*2)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*3 ) ;						it0 = srci(               v + lot*p*r*3 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*3 ) ;						it1 = srci( lot*p   +     v + lot*p*r*3 ) ;
           rt2 = srcr( lot*p*2 +     v + lot*p*r*3 ) ;						it2 = srci( lot*p*2 +     v + lot*p*r*3 ) ;

           rr1 = er1 * rt1 ;						ir1 = ei1 * rt1 ;
           rr2 = er2 * rt2 ;						ir2 = ei2 * rt2 ;

           rt1 = rr1 - ei1 * it1 ;					it1 = ir1 + er1 * it1 ;
           rt2 = rr2 - ei2 * it2 ;					it2 = ir2 + er2 * it2 ;

           rr1 = rt1 + rt2 ;						ir1 = it1 + it2 ;
           rt2 = rt1 - rt2 ;						it2 = it1 - it2 ;
           rr0 = rt0 + rr1 ;						ir0 = it0 + ir1 ;
           rr1 = rt0 + cr * rr1 ;					ir1 = it0 + cr * ir1 ;
           rr2 = rr1 + ci * it2 ;					ir2 = ir1 - ci * rt2 ;
           rr1 = rr1 - ci * it2 ;					ir1 = ir1 + ci * rt2 ;

           dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
        end do
     end do

  else	
     r = 0 ;
     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*3 ) ;						it0 = srci(               v + lot*p*r*3 ) ;
        rt1 = srcr( lot*p   +     v + lot*p*r*3 ) ;						it1 = srci( lot*p   +     v + lot*p*r*3 ) ;
        rt2 = srcr( lot*p*2 +     v + lot*p*r*3 ) ;						it2 = srci( lot*p*2 +     v + lot*p*r*3 ) ;

        rr1 = rt1 + rt2 ;						ir1 = it1 + it2 ;
        rt2 = rt1 - rt2 ;						it2 = it1 - it2 ;
        rr0 = rt0 + rr1 ;						ir0 = it0 + ir1 ;
        rr1 = rt0 + cr * rr1 ;					ir1 = it0 + cr * ir1 ;
        rr2 = rr1 - ci * it2 ;					ir2 = ir1 + ci * rt2 ;
        rr1 = rr1 + ci * it2 ;					ir1 = ir1 - ci * rt2 ;

        dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
        dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
        dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
     end do

     do r=1,m-1
        er1  = real(e(p*r)) ;							ei1  = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;							ei2 = aimag(e(p*r*2)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*3 ) ;						it0 = srci(               v + lot*p*r*3 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*3 ) ;						it1 = srci( lot*p   +     v + lot*p*r*3 ) ;
           rt2 = srcr( lot*p*2 +     v + lot*p*r*3 ) ;						it2 = srci( lot*p*2 +     v + lot*p*r*3 ) ;

           rr1 = er1 * rt1 ;						ir1 = ei1 * rt1 ;
           rr2 = er2 * rt2 ;						ir2 = ei2 * rt2 ;

           rt1 = rr1 + ei1 * it1 ;					it1 = -ir1 + er1 * it1 ;
           rt2 = rr2 + ei2 * it2 ;					it2 = -ir2 + er2 * it2 ;

           rr1 = rt1 + rt2 ;						ir1 = it1 + it2 ;
           rt2 = rt1 - rt2 ;						it2 = it1 - it2 ;
           rr0 = rt0 + rr1 ;						ir0 = it0 + ir1 ;
           rr1 = rt0 + cr * rr1 ;					ir1 = it0 + cr * ir1 ;
           rr2 = rr1 - ci * it2 ;					ir2 = ir1 + ci * rt2 ;
           rr1 = rr1 + ci * it2 ;					ir1 = ir1 - ci * rt2 ;

           dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
        end do
     end do
  end if
End Subroutine cache_radix_3


!----------- radix-5 --------------

Subroutine cache_radix_5(  &
  p   &		! p factor 
       ,m   &		! m factor 
       ,lot & 		! vector length 
       ,srcr, srci	&	! source address 
       ,dstr, dsti	&	! dest. address 
       ,e		&	! trig table address 
       ,isign		)	! forward/inverse direction 
  implicit none
  integer :: p	&		! p factor 
  ,  m		&		! m factor 
  ,  lot  	&		! vector length 
  ,  isign			! forward/inverse direction 
  real(8), dimension(0:0) :: srcr, srci,  &		! source address 
       dstr, dsti		! dest. address 
  complex(8), dimension(0:0) :: e			! trig table address 

  real(8) :: er1, ei1, er2, ei2, er3, ei3, er4, ei4, c1, s1, c2, s2, &
       rt0, rt1, rt2, rt3, rt4, it0, it1, it2, it3, it4, &
       ra1, ra2, ra3, ra4,      ia1, ia2, ia3, ia4, &
       rb1, rb2, rb3, rb4,      ib1, ib2, ib3, ib4, &
       rr0, rr1, rr2, rr3, rr4, ir0, ir1, ir2, ir3, ir4 ;

  integer :: r , v ;

  c1 = real(e(p*m)) ;		s1 = aimag(e(p*m)) ;
  c2 = real(e(2*p*m)) ;	s2 = aimag(e(2*p*m)) ;

  if ( isign > 0 )	then

     r = 0 ;

     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*5 ) ;						it0 = srci(               v + lot*p*r*5 ) ;
        rt1 = srcr( lot*p   +     v + lot*p*r*5 ) ;						it1 = srci( lot*p   +     v + lot*p*r*5 ) ;
        rt2 = srcr( lot*p*2 +     v + lot*p*r*5 ) ;						it2 = srci( lot*p*2 +     v + lot*p*r*5 ) ;
        rt3 = srcr( lot*p*3 +     v + lot*p*r*5 ) ;						it3 = srci( lot*p*3 +     v + lot*p*r*5 ) ;
        rt4 = srcr( lot*p*4 +     v + lot*p*r*5 ) ;						it4 = srci( lot*p*4 +     v + lot*p*r*5 ) ;

        ra1 = rt1 + rt4 ;						ia1 = it1 + it4 ;
        ra4 = rt1 - rt4 ;						ia4 = it1 - it4 ;
        ra2 = rt2 + rt3 ;						ia2 = it2 + it3 ;
        ra3 = rt2 - rt3 ;						ia3 = it2 - it3 ;

        rb1 = c1 * ra1 + c2 * ra2 ;				ib1 = c1 * ia1 + c2 * ia2 ;
        rb2 = c2 * ra1 + c1 * ra2 ;				ib2 = c2 * ia1 + c1 * ia2 ;
        rb3 = s1 * ia4 + s2 * ia3 ;				ib3 = s1 * ra4 + s2 * ra3 ;
        rb4 = s1 * ia3 - s2 * ia4 ;				ib4 = s1 * ra3 - s2 * ra4 ;

        rr0 = rt0 + ra1 + ra2 ;					ir0 = it0 + ia1 + ia2 ;
        rr1 = rt0 + rb1 - rb3 ;					ir1 = it0 + ib1 + ib3 ;
        rr2 = rt0 + rb2 + rb4 ;					ir2 = it0 + ib2 - ib4 ;
        rr3 = rt0 + rb2 - rb4 ;					ir3 = it0 + ib2 + ib4 ;
        rr4 = rt0 + rb1 + rb3 ;					ir4 = it0 + ib1 - ib3 ;

        dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
        dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
        dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
        dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;						dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 ;
        dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;						dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 ;		
     end do

     do r=1,m-1
        er1 = real(e(p*r)) ;							ei1 = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;							ei2 = aimag(e(p*r*2)) ;
        er3 = real(e(p*r*3)) ;							ei3 = aimag(e(p*r*3)) ;
        er4 = real(e(p*r*4)) ;							ei4 = aimag(e(p*r*4)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*5 ) ;						it0 = srci(               v + lot*p*r*5 ) ;
           rr1 = srcr( lot*p   +     v + lot*p*r*5 ) ;						ir1 = srci( lot*p   +     v + lot*p*r*5 ) ;
           rr2 = srcr( lot*p*2 +     v + lot*p*r*5 ) ;						ir2 = srci( lot*p*2 +     v + lot*p*r*5 ) ;
           rr3 = srcr( lot*p*3 +     v + lot*p*r*5 ) ;						ir3 = srci( lot*p*3 +     v + lot*p*r*5 ) ;
           rr4 = srcr( lot*p*4 +     v + lot*p*r*5 ) ;						ir4 = srci( lot*p*4 +     v + lot*p*r*5 ) ;

           rt1 = er1 * rr1 - ei1 * ir1 ;			it1 = er1 * ir1 + ei1 * rr1 ;
           rt2 = er2 * rr2 - ei2 * ir2 ;			it2 = er2 * ir2 + ei2 * rr2 ;
           rt3 = er3 * rr3 - ei3 * ir3 ;			it3 = er3 * ir3 + ei3 * rr3 ;
           rt4 = er4 * rr4 - ei4 * ir4 ;			it4 = er4 * ir4 + ei4 * rr4 ;

           ra1 = rt1 + rt4 ;						ia1 = it1 + it4 ;
           ra4 = rt1 - rt4 ;						ia4 = it1 - it4 ;
           ra2 = rt2 + rt3 ;						ia2 = it2 + it3 ;
           ra3 = rt2 - rt3 ;						ia3 = it2 - it3 ;

           rb1 = c1 * ra1 + c2 * ra2 ;				ib1 = c1 * ia1 + c2 * ia2 ;
           rb2 = c2 * ra1 + c1 * ra2 ;				ib2 = c2 * ia1 + c1 * ia2 ;
           rb3 = s1 * ia4 + s2 * ia3 ;				ib3 = s1 * ra4 + s2 * ra3 ;
           rb4 = s1 * ia3 - s2 * ia4 ;				ib4 = s1 * ra3 - s2 * ra4 ;

           rr0 = rt0 + ra1 + ra2 ;					ir0 = it0 + ia1 + ia2 ;
           rr1 = rt0 + rb1 - rb3 ;					ir1 = it0 + ib1 + ib3 ;
           rr2 = rt0 + rb2 + rb4 ;					ir2 = it0 + ib2 - ib4 ;
           rr3 = rt0 + rb2 - rb4 ;					ir3 = it0 + ib2 + ib4 ;
           rr4 = rt0 + rb1 + rb3 ;					ir4 = it0 + ib1 - ib3 ;

           dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
           dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;						dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 ;
           dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;						dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 ;

        end do
     end do

  else	

     r = 0 ;

     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*5 ) ;						it0 = srci(               v + lot*p*r*5 ) ;
        rt1 = srcr( lot*p   +     v + lot*p*r*5 ) ;						it1 = srci( lot*p   +     v + lot*p*r*5 ) ;
        rt2 = srcr( lot*p*2 +     v + lot*p*r*5 ) ;						it2 = srci( lot*p*2 +     v + lot*p*r*5 ) ;
        rt3 = srcr( lot*p*3 +     v + lot*p*r*5 ) ;						it3 = srci( lot*p*3 +     v + lot*p*r*5 ) ;
        rt4 = srcr( lot*p*4 +     v + lot*p*r*5 ) ;						it4 = srci( lot*p*4 +     v + lot*p*r*5 ) ;

        ra1 = rt1 + rt4 ;						ia1 = it1 + it4 ;
        ra4 = rt1 - rt4 ;						ia4 = it1 - it4 ;
        ra2 = rt2 + rt3 ;						ia2 = it2 + it3 ;
        ra3 = rt2 - rt3 ;						ia3 = it2 - it3 ;

        rb1 = c1 * ra1 + c2 * ra2 ;				ib1 = c1 * ia1 + c2 * ia2 ;
        rb2 = c2 * ra1 + c1 * ra2 ;				ib2 = c2 * ia1 + c1 * ia2 ;
        rb3 = s1 * ia4 + s2 * ia3 ;				ib3 = s1 * ra4 + s2 * ra3 ;
        rb4 = s1 * ia3 - s2 * ia4 ;				ib4 = s1 * ra3 - s2 * ra4 ;

        rr0 = rt0 + ra1 + ra2 ;					ir0 = it0 + ia1 + ia2 ;
        rr1 = rt0 + rb1 + rb3 ;					ir1 = it0 + ib1 - ib3 ;
        rr2 = rt0 + rb2 - rb4 ;					ir2 = it0 + ib2 + ib4 ;
        rr3 = rt0 + rb2 + rb4 ;					ir3 = it0 + ib2 - ib4 ;
        rr4 = rt0 + rb1 - rb3 ;					ir4 = it0 + ib1 + ib3 ;

        dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
        dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
        dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
        dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;						dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 ;
        dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;						dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 ;		
     end do

     do r=1,m-1
        er1 = real(e(p*r)) ;							ei1 = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;							ei2 = aimag(e(p*r*2)) ;
        er3 = real(e(p*r*3)) ;							ei3 = aimag(e(p*r*3)) ;
        er4 = real(e(p*r*4)) ;							ei4 = aimag(e(p*r*4)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*5 ) ;						it0 = srci(               v + lot*p*r*5 ) ;
           rr1 = srcr( lot*p   +     v + lot*p*r*5 ) ;						ir1 = srci( lot*p   +     v + lot*p*r*5 ) ;
           rr2 = srcr( lot*p*2 +     v + lot*p*r*5 ) ;						ir2 = srci( lot*p*2 +     v + lot*p*r*5 ) ;
           rr3 = srcr( lot*p*3 +     v + lot*p*r*5 ) ;						ir3 = srci( lot*p*3 +     v + lot*p*r*5 ) ;
           rr4 = srcr( lot*p*4 +     v + lot*p*r*5 ) ;						ir4 = srci( lot*p*4 +     v + lot*p*r*5 ) ;

           rt1 = er1 * rr1 + ei1 * ir1 ;			it1 = er1 * ir1 - ei1 * rr1 ;
           rt2 = er2 * rr2 + ei2 * ir2 ;			it2 = er2 * ir2 - ei2 * rr2 ;
           rt3 = er3 * rr3 + ei3 * ir3 ;			it3 = er3 * ir3 - ei3 * rr3 ;
           rt4 = er4 * rr4 + ei4 * ir4 ;			it4 = er4 * ir4 - ei4 * rr4 ;

           ra1 = rt1 + rt4 ;						ia1 = it1 + it4 ;
           ra4 = rt1 - rt4 ;						ia4 = it1 - it4 ;
           ra2 = rt2 + rt3 ;						ia2 = it2 + it3 ;
           ra3 = rt2 - rt3 ;						ia3 = it2 - it3 ;

           rb1 = c1 * ra1 + c2 * ra2 ;				ib1 = c1 * ia1 + c2 * ia2 ;
           rb2 = c2 * ra1 + c1 * ra2 ;				ib2 = c2 * ia1 + c1 * ia2 ;
           rb3 = s1 * ia4 + s2 * ia3 ;				ib3 = s1 * ra4 + s2 * ra3 ;
           rb4 = s1 * ia3 - s2 * ia4 ;				ib4 = s1 * ra3 - s2 * ra4 ;

           rr0 = rt0 + ra1 + ra2 ;					ir0 = it0 + ia1 + ia2 ;
           rr1 = rt0 + rb1 + rb3 ;					ir1 = it0 + ib1 - ib3 ;
           rr2 = rt0 + rb2 - rb4 ;					ir2 = it0 + ib2 + ib4 ;
           rr3 = rt0 + rb2 + rb4 ;					ir3 = it0 + ib2 - ib4 ;
           rr4 = rt0 + rb1 - rb3 ;					ir4 = it0 + ib1 + ib3 ;

           dstr(               v + lot*p*r ) = rr0 ;						dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;						dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rr2 ;						dsti( lot*m*p*2 +   v + lot*p*r ) = ir2 ;
           dstr( lot*m*p*3 +   v + lot*p*r ) = rr3 ;						dsti( lot*m*p*3 +   v + lot*p*r ) = ir3 ;
           dstr( lot*m*p*4 +   v + lot*p*r ) = rr4 ;						dsti( lot*m*p*4 +   v + lot*p*r ) = ir4 ;


        end do
     end do
  end if
End Subroutine cache_radix_5


!----------- radix-4 ---------------

Subroutine cache_radix_4(  &
  p   &		! p factor 
       ,m   &		! m factor 
       ,lot & 		! vector length 
       ,srcr, srci	&	! source address 
       ,dstr, dsti	&	! dest. address 
       ,e		&	! trig table address 
       ,isign		)	! forward/inverse direction 
  implicit none
  integer :: p	&		! p factor 
  ,  m		&		! m factor 
  ,  lot  	&		! vector length 
  ,  isign			! forward/inverse direction 
  real(8), dimension(0:0) :: srcr, srci,  &		! source address 
       dstr, dsti		! dest. address 
  complex(8), dimension(0:0) :: e			! trig table address 

  real(8) ::	            er1 , ei1 , er2 , ei2 , er3 , ei3  &
       , rr0 , ir0 , rr1 , ir1 , rr2 , ir2 , rr3 , ir3  &
       , rt0 , it0 , rt1 , it1 , rt2 , it2 , rt3 , it3 ;

  integer :: r, v ;

  if ( isign > 0 )	then
     r = 0 ;
     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*4 ) ;				it0 = srci(               v + lot*p*r*4 ) ;
        rt1 = srcr( lot*p   +     v + lot*p*r*4 ) ;				it1 = srci( lot*p   +     v + lot*p*r*4 ) ;
        rt2 = srcr( lot*p*2 +     v + lot*p*r*4 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*4 ) ;
        rt3 = srcr( lot*p*3 +     v + lot*p*r*4 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*4 ) ;

        rr0 = rt0 + rt2 ;				ir0 = it0 + it2 ;
        rr1 = rt0 - rt2 ;				ir1 = it0 - it2 ;
        rr2 = rt1 + rt3 ;				ir2 = it1 + it3 ;
        rr3 = rt1 - rt3 ;				ir3 = it1 - it3 ;

        rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 ;
        rt1 = rr1 - ir3 ;				it1 = ir1 + rr3 ;
        rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 ;
        rt3 = rr1 + ir3 ;				it3 = ir1 - rr3 ;

        dstr(               v + lot*p*r ) = rt0 ;				dsti(               v + lot*p*r ) = it0 ;
        dstr( lot*m*p   +   v + lot*p*r ) = rt1 ;				dsti( lot*m*p   +   v + lot*p*r ) = it1 ;
        dstr( lot*m*p*2 +   v + lot*p*r ) = rt2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = it2 ;
        dstr( lot*m*p*3 +   v + lot*p*r ) = rt3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = it3 ;	
     end do

     do r=1,m-1

        er1 = real(e(p*r)) ;		ei1 = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;		ei2 = aimag(e(p*r*2)) ;
        er3 = real(e(p*r*3)) ;		ei3 = aimag(e(p*r*3)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*4 ) ;				it0 = srci(               v + lot*p*r*4 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*4 ) ;				it1 = srci( lot*p   +     v + lot*p*r*4 ) ;
           rt2 = srcr( lot*p*2 +     v + lot*p*r*4 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*4 ) ;
           rt3 = srcr( lot*p*3 +     v + lot*p*r*4 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*4 ) ;

           rr1 = er1 * rt1 ;				ir1 = ei1 * rt1 ;
           rr2 = er2 * rt2 ;				ir2 = ei2 * rt2 ;
           rr3 = er3 * rt3 ;				ir3 = ei3 * rt3 ;

           rt1 = rr1 - ei1 * it1 ;			it1 = ir1 + er1 * it1 ;
           rt2 = rr2 - ei2 * it2 ;			it2 = ir2 + er2 * it2 ;
           rt3 = rr3 - ei3 * it3 ;			it3 = ir3 + er3 * it3 ;

           rr0 = rt0 + rt2 ;				ir0 = it0 + it2 ;
           rr1 = rt0 - rt2 ;				ir1 = it0 - it2 ;
           rr2 = rt1 + rt3 ;				ir2 = it1 + it3 ;
           rr3 = rt1 - rt3 ;				ir3 = it1 - it3 ;

           rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 ;
           rt1 = rr1 - ir3 ;				it1 = ir1 + rr3 ;
           rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 ;
           rt3 = rr1 + ir3 ;				it3 = ir1 - rr3 ;

           dstr(               v + lot*p*r ) = rt0 ;				dsti(               v + lot*p*r ) = it0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rt1 ;				dsti( lot*m*p   +   v + lot*p*r ) = it1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rt2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = it2 ;
           dstr( lot*m*p*3 +   v + lot*p*r ) = rt3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = it3 ;	
        end do
     end do

  else	
     r = 0 ;
     do v=0,lot*p-1
        rt0 = srcr(               v + lot*p*r*4 ) ;				it0 = srci(               v + lot*p*r*4 ) ;
        rt1 = srcr( lot*p   +     v + lot*p*r*4 ) ;				it1 = srci( lot*p   +     v + lot*p*r*4 ) ;
        rt2 = srcr( lot*p*2 +     v + lot*p*r*4 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*4 ) ;
        rt3 = srcr( lot*p*3 +     v + lot*p*r*4 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*4 ) ;

        rr0 = rt0 + rt2 ;				ir0 = it0 + it2 ;
        rr1 = rt0 - rt2 ;				ir1 = it0 - it2 ;
        rr2 = rt1 + rt3 ;				ir2 = it1 + it3 ;
        rr3 = rt1 - rt3 ;				ir3 = it1 - it3 ;

        rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 ;
        rt1 = rr1 + ir3 ;				it1 = ir1 - rr3 ;
        rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 ;
        rt3 = rr1 - ir3 ;				it3 = ir1 + rr3 ;

        dstr(               v + lot*p*r ) = rt0 ;				dsti(               v + lot*p*r ) = it0 ;
        dstr( lot*m*p   +   v + lot*p*r ) = rt1 ;				dsti( lot*m*p   +   v + lot*p*r ) = it1 ;
        dstr( lot*m*p*2 +   v + lot*p*r ) = rt2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = it2 ;
        dstr( lot*m*p*3 +   v + lot*p*r ) = rt3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = it3 ;	
     end do

     do r=1,m-1

        er1 = real(e(p*r)) ;		ei1 = aimag(e(p*r)) ;
        er2 = real(e(p*r*2)) ;		ei2 = aimag(e(p*r*2)) ;
        er3 = real(e(p*r*3)) ;		ei3 = aimag(e(p*r*3)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*4 ) ;				it0 = srci(               v + lot*p*r*4 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*4 ) ;				it1 = srci( lot*p   +     v + lot*p*r*4 ) ;
           rt2 = srcr( lot*p*2 +     v + lot*p*r*4 ) ;				it2 = srci( lot*p*2 +     v + lot*p*r*4 ) ;
           rt3 = srcr( lot*p*3 +     v + lot*p*r*4 ) ;				it3 = srci( lot*p*3 +     v + lot*p*r*4 ) ;

           rr1 = er1 * rt1 ;				ir1 = ei1 * rt1 ;
           rr2 = er2 * rt2 ;				ir2 = ei2 * rt2 ;
           rr3 = er3 * rt3 ;				ir3 = ei3 * rt3 ;

           rt1 = rr1 + ei1 * it1 ;			it1 = -ir1 + er1 * it1 ;
           rt2 = rr2 + ei2 * it2 ;			it2 = -ir2 + er2 * it2 ;
           rt3 = rr3 + ei3 * it3 ;			it3 = -ir3 + er3 * it3 ;

           rr0 = rt0 + rt2 ;				ir0 = it0 + it2 ;
           rr1 = rt0 - rt2 ;				ir1 = it0 - it2 ;
           rr2 = rt1 + rt3 ;				ir2 = it1 + it3 ;
           rr3 = rt1 - rt3 ;				ir3 = it1 - it3 ;

           rt0 = rr0 + rr2 ;				it0 = ir0 + ir2 ;
           rt1 = rr1 + ir3 ;				it1 = ir1 - rr3 ;
           rt2 = rr0 - rr2 ;				it2 = ir0 - ir2 ;
           rt3 = rr1 - ir3 ;				it3 = ir1 + rr3 ;

           dstr(               v + lot*p*r ) = rt0 ;				dsti(               v + lot*p*r ) = it0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rt1 ;				dsti( lot*m*p   +   v + lot*p*r ) = it1 ;
           dstr( lot*m*p*2 +   v + lot*p*r ) = rt2 ;				dsti( lot*m*p*2 +   v + lot*p*r ) = it2 ;
           dstr( lot*m*p*3 +   v + lot*p*r ) = rt3 ;				dsti( lot*m*p*3 +   v + lot*p*r ) = it3 ;
        end do
     end do
  end if
End Subroutine cache_radix_4

!----------- radix-2 ---------------

Subroutine cache_radix_2(  &
  p               &       ! p factor 
       ,m   		&       ! m factor 
       ,lot  	   	&       ! vector length 
       ,srcr, srci	&	! source address 
       ,dstr, dsti	&	! dest. address 
       ,e		&	! trig table address 
       ,isign		)	! forward/inverse direction 
  implicit none
  integer :: p,	&	! p factor 
  m,		&	! m factor 
  lot,  	&	! vector length 
  isign		        ! forward/inverse direction 
  real(8), dimension(0:0) :: srcr, srci,  &		! source address 
       dstr, dsti   		! dest. address 
  complex(8), dimension(0:0) :: e			! trig table address 

  real(8) ::	            er1 , ei1   &
       , rr0 , ir0 , rr1 , ir1 , rt0 , rt1 , it0 , it1 ;

  integer :: r, v ;
  r = 0 ;
  do v=0,lot*p-1
     rt0 = srcr(               v + lot*p*r*2 ) ;				it0 = srci(               v + lot*p*r*2 ) ;
     rt1 = srcr( lot*p   +     v + lot*p*r*2 ) ;				it1 = srci( lot*p   +     v + lot*p*r*2 ) ;

     rr0 = rt0 + rt1 ;				ir0 = it0 + it1 ;
     rr1 = rt0 - rt1 ;				ir1 = it0 - it1 ;

     dstr(               v + lot*p*r ) = rr0 ;				dsti(               v + lot*p*r ) = ir0 ;
     dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;				dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;		
  end do

  if ( isign > 0 )	then
     do r=1,m-1

        er1 = real(e(p*r)) ;		ei1 = aimag(e(p*r)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*2 ) ;				it0 = srci(               v + lot*p*r*2 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*2 ) ;				it1 = srci( lot*p   +     v + lot*p*r*2 ) ;

           rr1 = er1 * rt1 ;				ir1 = ei1 * rt1 ;
           rt1 = rr1 - ei1 * it1 ;			it1 = ir1 + er1 * it1 ;

           rr0 = rt0 + rt1 ;				ir0 = it0 + it1 ;
           rr1 = rt0 - rt1 ;				ir1 = it0 - it1 ;

           dstr(               v + lot*p*r ) = rr0 ;				dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;				dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;
        end do
     end do


  else	
     do r=1,m-1

        er1 = real(e(p*r)) ;		ei1 = aimag(e(p*r)) ;

        do v=0,lot*p-1
           rt0 = srcr(               v + lot*p*r*2 ) ;				it0 = srci(               v + lot*p*r*2 ) ;
           rt1 = srcr( lot*p   +     v + lot*p*r*2 ) ;				it1 = srci( lot*p   +     v + lot*p*r*2 ) ;

           rr1 = er1 * rt1 ;				ir1 = ei1 * rt1 ;
           rt1 = rr1 + ei1 * it1 ;			it1 = -ir1 + er1 * it1 ;

           rr0 = rt0 + rt1 ;				ir0 = it0 + it1 ;
           rr1 = rt0 - rt1 ;				ir1 = it0 - it1 ;

           dstr(               v + lot*p*r ) = rr0 ;				dsti(               v + lot*p*r ) = ir0 ;
           dstr( lot*m*p   +   v + lot*p*r ) = rr1 ;				dsti( lot*m*p   +   v + lot*p*r ) = ir1 ;		
        end do
     end do
  end if
End Subroutine cache_radix_2

!====================  real --> complex butterfy  ========================

Subroutine real_to_complex (  &
  n,          &		! transform length 
       lot,	 &		! vector length 
       srcr, srci, &		! source address 
       dstr, dsti, &		! destination address 
       e )    		        ! trig table address 
  implicit none
  integer :: n, &				! transform length 
       lot			        ! vector length 
  real(8), dimension(0:0) :: srcr, srci,	 &	! source address 
       dstr, dsti		! destination address 
  complex(8), dimension(0:0) :: e 		! trig table address 

  real(8) ::	er , ei , rr1 , ir1 , rt0 , rt1 , it0 , it1 ;

  integer ::  r , v;

  do v=0,lot-1
     r = 0 ;
     dstr( v + r*lot ) = srcr( v + r*lot ) + srci( v + r*lot );		dsti( v + r*lot ) = srcr( v + r*lot ) - srci( v + r*lot );

     do r=1,n/2
        er = real(e(r-1)) ;						ei = aimag(e(r-1)) ;
        rt0 = srcr( v + r*lot ) + srcr( n*lot   + v - r*lot );		it0 = srci( v + r*lot ) - srci( n*lot   + v - r*lot );						
        rt1 = srcr( v + r*lot ) - srcr( n*lot   + v - r*lot );		it1 = srci( v + r*lot ) + srci( n*lot   + v - r*lot );
        rr1 = er * it1 + ei * rt1;			ir1 = er * rt1 - ei * it1;
        dstr( v + r*lot ) = rt0 - rr1;				dsti( v + r*lot ) =  it0 + ir1;							
        dstr( n*lot   + v - r*lot ) = rt0 + rr1;				dsti( n*lot   + v - r*lot ) = -it0 + ir1;     
     end do
  end do
End Subroutine real_to_complex

!====================  complex --> real butterfy  ========================

Subroutine complex_to_real (  &
  n,		 &		! transform length 
       lot,	 &		! vector length 
       srcr, srci, &		! source address 
       dstr, dsti, &		! destination address 
       e )    		        ! trig table address 
  implicit none
  integer :: n, &				! transform length 
       lot			        ! vector length 
  real(8), dimension(0:0) :: srcr, srci,	 &	! source address 
       dstr, dsti		! destination address 
  complex(8), dimension(0:0) :: e 		! trig table address 

  real(8) ::	er , ei , rr1 , ir1 , rt0 , rt1 , it0 , it1 ;

  integer ::  r , v;

  do v=0,lot-1
     r = 0 ;
     dstr( v + r*lot ) = 2.d0 * ( srcr( v + r*lot ) + srci( v + r*lot ) );
     dsti( v + r*lot ) = 2.d0 * ( srcr( v + r*lot ) - srci( v + r*lot ) );

     do r=1,n/2
        er = real(e(r-1)) ;						ei = aimag(e(r-1)) ;
        rt0 = srcr( v + r*lot ) + srcr( n*lot   + v - r*lot );		it0 = srci( v + r*lot ) - srci( n*lot   + v - r*lot );						
        rt1 = srcr( n*lot   + v - r*lot ) - srcr( v + r*lot );		it1 = srci( v + r*lot ) + srci( n*lot   + v - r*lot );
        rr1 = er * it1 + ei * rt1;			ir1 = er * rt1 - ei * it1;
        dstr( v + r*lot ) = rt0 + rr1;				dsti( v + r*lot ) =  it0 + ir1;							
        dstr( n*lot   + v - r*lot ) = rt0 - rr1;				dsti( n*lot   + v - r*lot ) = -it0 + ir1;     
     end do
  end do
End Subroutine complex_to_real


!====================  exponential <-- cosine butterfy  ========================

Subroutine cosine_to_complex (n, lot, srcr, srci, dstr, dsti, isign)
  implicit none
  integer :: n, lot, isign
  real(8), dimension(0:0) :: srcr, srci, dstr, dsti
  
  real(8) ::	f1, rt0, it0 ;
  complex(8) :: rsv ;

  integer ::  r , v;

  if( isign > 0 ) then
     f1 = 2.d0;
  else
     f1 = 1.d0;
  end if

  do v=0,lot-1
     r = 0;
     dstr( v + r*lot ) = srcr( v + r*2*lot );					dsti( v + r*lot ) = srci( v + r*2*lot );
     r = n/2;
     dstr( v + r*lot ) = f1 * srcr( v + r*2*lot );				dsti( v + r*lot ) = f1 * srci( v + r*2*lot );
     rsv = cmplx(srcr( - lot     + v + r*2*lot ), srci( - lot     + v + r*2*lot ), kind=8)

     do r=1,n/2-1
        rsv = cmplx(real(rsv) + srcr( - lot     + v + r*2*lot ), aimag(rsv) + srci( - lot     + v + r*2*lot ), kind=8)
        rt0 = srcr( lot     + v + r*2*lot ) - srcr( - lot     + v + r*2*lot );
        it0 = srci( lot     + v + r*2*lot ) - srci( - lot     + v + r*2*lot );
        dstr( v + r*lot ) = srcr( v + r*2*lot ) - it0;		dsti( v + r*lot ) = srci( v + r*2*lot ) + rt0;
        dstr( n*lot   + v - r*lot ) = srcr( v + r*2*lot ) + it0;		dsti( n*lot   + v - r*lot ) = srci( v + r*2*lot ) - rt0;
     end do
     r = 0 ;
     dstr( n*lot   + v - r*lot ) = real(rsv);						dsti( n*lot   + v - r*lot ) = aimag(rsv);
     r = n/2;
     srcr( v + r*2*lot ) = real(rsv);						srci( v + r*2*lot ) = aimag(rsv);
  end do
End Subroutine cosine_to_complex


!====================  exponential --> cosine butterfy  ========================

Subroutine complex_to_cosine (n, lot, srcr, srci, dstr, dsti, e, isign)
  implicit none
  integer :: n, lot, isign
  real(8), dimension(0:0) :: srcr, srci, dstr, dsti
  complex(8), dimension(0:0) :: e
  
  real(8), parameter :: f2 = (.5d0), f4 = (.25d0)

  real(8) ::	f3, rt0,it0 , rt1,it1 ;
  complex(8) :: rsv ;

  integer ::  r , v;

  if( isign > 0 ) then
     f3 = f2;
  else		
     f3 = f4;
  end if

  do v=0,lot-1
     r = 0;
     rsv = cmplx(2.d0 * srcr( n*lot   + v - r*lot ), 2.d0 * srci( n*lot   + v - r*lot ), kind=8)
     dstr( v + r*lot ) = f2 * ( srcr( v + r*lot ) + real(rsv) );				dsti( v + r*lot ) = f2 * ( srci( v + r*lot ) + aimag(rsv) );
     dstr( n*lot   + v - r*lot ) = f3 * ( srcr( v + r*lot ) - real(rsv) );
     dsti( n*lot   + v - r*lot ) = f3 * ( srci( v + r*lot ) - aimag(rsv) );
     dstr( n/2*lot + v         ) = f2 * srcr( n/2*lot + v         );
     dsti( n/2*lot + v         ) = f2 * srci( n/2*lot + v         );

     do r=1,n/2-1
        rt0 = f4     * ( srcr( n*lot   + v - r*lot ) + srcr( v + r*lot ) );
        it0 = f4     * ( srci( n*lot   + v - r*lot ) + srci( v + r*lot ) );
        rt1 = real(e(r)) * ( srcr( n*lot   + v - r*lot ) - srcr( v + r*lot ) );
        it1 = real(e(r)) * ( srci( n*lot   + v - r*lot ) - srci( v + r*lot ) );
        dstr( v + r*lot ) = rt0 - rt1;							dsti( v + r*lot ) = it0 - it1;
        dstr( n*lot   + v - r*lot ) = rt0 + rt1;							dsti( n*lot   + v - r*lot ) = it0 + it1;
     end do
  end do
End Subroutine complex_to_cosine

!====================  exponential <-- sine butterfy  ========================

Subroutine sine_to_complex (n, lot, srcr, srci, dstr, dsti )
  implicit none
  integer:: n, lot
  real(8), dimension(0:0) :: srcr, srci, dstr, dsti

  real(8) ::	rt0, it0 ;

  integer ::  r , v;

  do v=0,lot-1
     r = 0;
     dstr( v + r*lot ) = -2.d0 * srcr( lot     + v + r*2*lot );			dsti( v + r*lot ) = -2.d0 * srci( lot     + v + r*2*lot );
     r = n/2;
     dstr( v + r*lot ) = 2.d0 * srcr( - lot     + v + r*2*lot );				dsti( v + r*lot ) = 2.d0 * srci( - lot     + v + r*2*lot );

     do r=1,n/2-1
        rt0 = srcr( lot     + v + r*2*lot ) - srcr( - lot     + v + r*2*lot );
        it0 = srci( lot     + v + r*2*lot ) - srci( - lot     + v + r*2*lot );
        dstr( v + r*lot ) = -srci( v + r*2*lot ) - rt0;		dsti( v + r*lot ) =  srcr( v + r*2*lot ) - it0;
        dstr( n*lot   + v - r*lot ) =  srci( v + r*2*lot ) - rt0;		dsti( n*lot   + v - r*lot ) = -srcr( v + r*2*lot ) - it0;
     end do
  end do
End Subroutine sine_to_complex

!====================  exponential --> sine butterfy  ========================

Subroutine complex_to_sine (n, lot, srcr, srci, dstr, dsti, e)
  implicit none
  integer :: n, int, lot
  real(8), dimension(0:0) :: srcr, srci, dstr, dsti
  complex(8), dimension(0:0) :: e
  
  real(8), parameter :: f2 = (-.25d0), f4 = (.25d0)

  real(8) ::	rt0,it0 , rt1,it1 ;

  integer ::  r , v;

  do v=0,lot-1
     r = 0;
     dstr( v + r*lot ) = 0.;                             		dsti( v + r*lot ) = 0.;
     dstr( n*lot   + v - r*lot ) = 0.;										dsti( n*lot   + v - r*lot ) = 0.;
     dstr( n/2*lot + v         ) = f2 * srcr( n/2*lot + v         );
     dsti( n/2*lot + v         ) = f2 * srci( n/2*lot + v         );

     do r=1,n/2-1
        rt0 = f4     * ( srcr( v + r*lot ) - srcr( n*lot   + v - r*lot ) );
        it0 = f4     * ( srci( v + r*lot ) - srci( n*lot   + v - r*lot ) );
        rt1 = real(e(r)) * ( srcr( v + r*lot ) + srcr( n*lot   + v - r*lot ) );
        it1 = real(e(r)) * ( srci( v + r*lot ) + srci( n*lot   + v - r*lot ) );
        dstr( v + r*lot ) =  rt0 + rt1;							dsti( v + r*lot ) =  it0 + it1;
        dstr( n*lot   + v - r*lot ) = -rt0 + rt1;							dsti( n*lot   + v - r*lot ) = -it0 + it1;
     end do
  end do
End Subroutine complex_to_sine


!****************************************************************************
!real and complex fast fourier transforms  (strip mined)
!*****************************************************************************
!*****************   define trig table   **********************

Subroutine ctrig( n, e )
  implicit none
  integer :: n;
  complex(8) :: e(0:0);

  integer :: i;
  do i=0,n-1
     e(i) = cmplx(cos(2*3.141592653589793d0 * real(i,8) / real(n,8) ), sin(2*3.141592653589793d0 * real(i,8) / real(n,8) ), kind=8)
  end do
End Subroutine ctrig

Subroutine rtrig( n, e )
  implicit none
  integer :: n;
  complex(8) :: e(0:0);

  integer :: i;

  if( (mod(n,2)) /= 0 )  stop "Real transform length is not even."
  call ctrig( n/2, e ) ;

  do i=0,n/4-1
     e(i+n/2) = cmplx(cos(2*3.141592653589793d0 * real(i+1,8) / real(n,8) ), &
          sin(2*3.141592653589793d0 * real(i+1,8) / real(n,8) ), kind=8)
  end do
End Subroutine rtrig

Subroutine ttrig( n, e )
  implicit none
  integer :: n;
  complex(8) :: e(0:0);

  integer :: i;

  call ctrig( n, e ) ;

  do i=1,n/2-1
     e(i+n) = -.125d0 / sin( 3.141592653589793d0 * real(i,8) / real(n,8) );    ! imag part is 0
  end do
End Subroutine ttrig

!*********************   complex transform   ******************************

Subroutine cfft (  &
       n,	&		! transform length 
       lot,	&		! vector length 
       data,	&	        ! source address 
       inc,	&		! source data stride 
       jump,	&		! source vector stride 
       e,	&		! trig table address 
       isign )			! direction of transform 
  use cache_data
  implicit none
  integer :: n,	    & 			! transform length 
       lot,   &			! vector length 
       isign, &			! direction of transform 
       inc,   &	                ! source data stride 
       jump			! source vector stride 
  complex(8), dimension(0:0) :: data		! source address 
  complex(8), dimension(0:0) :: e			! trig table address 

  type(scomplex), target :: work1 , work2
  type(scomplex), pointer :: wp0 , wp1 , wps
  integer ::   p, m, j, vl, vj, lot_ctr ;

  vl = CACHE/n ;
  if( vl < 1 ) stop "Transform exceeds cache in cfft"
  vj = vl * jump ;		j = 0 ;

  lot_ctr = lot
  do while( lot_ctr > 0 )	
     if( lot_ctr < vl )	vl = lot_ctr ;

     wp0 => work1 ;	wp1 => work2 ;
     call load_cache( n, vl, data(j), wp0%r, wp0%i, inc, jump ) ;

     p = 1 ;		m = n ;

     do while( (mod(m,8)) == 0 )	
        m = m/8 ;
        call cache_radix_8( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*8 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,4)) == 0 )	
        m = m/4 ;
        call cache_radix_4( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*4 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,2)) == 0 )	
        m = m/2 ;
        call cache_radix_2( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*2 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,5)) == 0 )	
        m = m/5 ;
        call cache_radix_5( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*5 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,3)) == 0 )	
        m = m/3 ;
        call cache_radix_3( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*3 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     call store_cache( n, vl, wp0%r, wp0%i, data(j), inc, jump ) ;
		j = j + vj ;			
     lot_ctr = lot_ctr - vl ;
  end do

End Subroutine cfft

  !***************************   real transform   ************************

Subroutine rfft (  &
       n,	&		! transform length 
       lot,	&		! vector length 
       data,	&	        ! source address 
       inc,	&		! source data stride 
       jump,	&		! source vector stride 
       e,	&		! trig table address 
       isign )			! direction of transform 
  use cache_data
  implicit none
  integer :: n,	    & 			! transform length 
       lot,   &			! vector length 
       isign, &			! direction of transform 
       inc,   &	                ! source data stride 
       jump			! source vector stride 
  complex(8), dimension(0:0) :: data		! source address 
  complex(8), dimension(0:0) :: e			! trig table address 
  

  type(scomplex), target :: work1 , work2
  type(scomplex), pointer :: wp0 , wp1 , wps
  integer :: p, m, j, vl, vj, lot_ctr ;

  vl = CACHE/(n/2) ;
  if( vl < 1 ) stop "Transform exceeds cache in rfft"
  vj = vl * jump ;		j = 0 ;

  lot_ctr = lot
  do while( lot_ctr > 0 )	
     if( lot_ctr < vl )	vl = lot_ctr ;

     wp0 => work1 ;	wp1 => work2 ;
     call load_cache( n/2, vl, data(j), wp0%r, wp0%i, inc, jump ) ;

     if ( isign > 0 )		then
        call real_to_complex( n/2, vl, wp0%r, wp0%i, wp1%r, wp1%i, e(n/2) ) ;
        wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end if

     p = 1 ;		m = n/2 ;
     do while( (mod(m,8)) == 0 )	
        m = m/8 ;
        call cache_radix_8( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*8 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,4)) == 0 )	
        m = m/4 ;
        call cache_radix_4( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*4 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,2)) == 0 )	
        m = m/2 ;
        call cache_radix_2( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*2 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,5)) == 0 )	
        m = m/5 ;
        call cache_radix_5( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*5 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,3)) == 0 )	
        m = m/3 ;
        call cache_radix_3( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*3 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     if ( isign < 0 )		then
        call complex_to_real( n/2, vl, wp0%r, wp0%i, wp1%r, wp1%i, e(n/2) ) ;	
        wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end if

     call store_cache( n/2, vl, wp0%r, wp0%i, data(j), inc, jump ) ;

     lot_ctr = lot_ctr - vl ;		j = j + vj ;		
  end do

End Subroutine rfft

  !***************************   real transform with disjoint real and imaginary parts  ************************

Subroutine rfft_dri (  &
       n,	&		! transform length 
       lot,	&		! vector length 
       data,	&	        ! source address 
       inc,	&		! source data stride 
       jump,	&		! source vector stride 
       e,	&		! trig table address 
       isign )			! direction of transform 
  use cache_data
  implicit none
  integer :: n,	    & 			! transform length 
       lot,   &			! vector length 
       isign, &			! direction of transform 
       inc,   &	                ! source data stride in real words 
       jump			! source vector stride in real words
  real(8), dimension(0:0) :: data		        ! source address of complex data with disjoint real and imaginary parts
  complex(8), dimension(0:0) :: e			! trig table address 
  

  type(scomplex), target :: work1 , work2
  type(scomplex), pointer :: wp0 , wp1 , wps
  integer :: p, m, j, vl, vj, lot_ctr ;

  vl = CACHE/(n/2) ;
  if( vl < 1 ) stop "Transform exceeds cache in rfft"
  vj = vl * jump ;		j = 0 ;

  lot_ctr = lot
  do while( lot_ctr > 0 )	
     if( lot_ctr < vl )	vl = lot_ctr ;

     wp0 => work1 ;	wp1 => work2 ;
     call load_cache_dri( n/2, vl, data(j), wp0%r, wp0%i, inc, jump ) ;

     if ( isign > 0 )		then
        call real_to_complex( n/2, vl, wp0%r, wp0%i, wp1%r, wp1%i, e(n/2) ) ;
        wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end if

     p = 1 ;		m = n/2 ;
     do while( (mod(m,8)) == 0 )	
        m = m/8 ;
        call cache_radix_8( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*8 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,4)) == 0 )	
        m = m/4 ;
        call cache_radix_4( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*4 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,2)) == 0 )	
        m = m/2 ;
        call cache_radix_2( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*2 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,5)) == 0 )	
        m = m/5 ;
        call cache_radix_5( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*5 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,3)) == 0 )	
        m = m/3 ;
        call cache_radix_3( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*3 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     if ( isign < 0 )		then
        call complex_to_real( n/2, vl, wp0%r, wp0%i, wp1%r, wp1%i, e(n/2) ) ;	
        wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end if

     call store_cache_dri( n/2, vl, wp0%r, wp0%i, data(j), inc, jump ) ;

     lot_ctr = lot_ctr - vl ;		j = j + vj ;		
  end do

End Subroutine rfft_dri


Subroutine rfft2 (  &
  n,	&		! transform length 
  lot,	&		! vector length 
  data,	&	        ! source address 
       inc,	&		! source data stride 
       jump,	&		! source vector stride 
       e,	&		! trig table address 
       isign )			! direction of transform 
  implicit none
  integer :: n,	    & 			! transform length 
       lot,   &			! vector length 
       isign, &			! direction of transform 
       inc,   &	                ! source data stride 
       jump			! source vector stride 
  complex(8), dimension(0:0) :: data		! source address 
  complex(8), dimension(0:0) :: e			! trig table address 
  
  integer :: i ;

  if ( isign > 0 ) then	
     do i=0,lot-1
        data(i*jump) = cmplx(real(data(i*jump)), real(data(n/2*inc+i*jump)), kind=8) ;   ! real part is unchanged
     end do
     call rfft( n, lot, data, inc, jump, e, isign ) ;	
  else				
     call rfft( n, lot, data, inc, jump, e, isign ) ;
     do i=0,lot-1
        data(n/2*inc+i*jump) = aimag(data(i*jump)) ;   ! imag part is set 0
        data(i*jump) 	     = real(data(i*jump)) ;    ! imag part is set 0
     end do
  end if
End Subroutine rfft2

                     !***************************   cosine transform   ************************

Subroutine cfct (  &
  n,	&		! transform length 
  lot,	&		! vector length 
  data,	&	        ! source address 
  inc,	&		! source data stride 
  jump,	&		! source vector stride 
  e,	&		! trig table address 
  isign )			! direction of transform 
  use cache_data
  implicit none
  integer :: n,	    & 			! transform length 
       lot,   &			! vector length 
       isign, &			! direction of transform 
       inc,   &	                ! source data stride 
       jump			! source vector stride 
  complex(8), dimension(0:0) :: data		! source address 
  complex(8), dimension(0:0) :: e			! trig table address 
  

  type(scomplex), target :: work1 , work2
  type(scomplex), pointer :: wp0 , wp1 , wps
  integer :: p, m, j, vl, vj, lot_ctr ;

  vl = CACHE/(n+1) ;
  if( vl < 1 ) stop "Transform exceeds cache in cfct"
  vj = vl * jump ;		j = 0 ;

  lot_ctr = lot
  do while( lot_ctr > 0 )	
     if( lot_ctr < vl )	vl = lot_ctr ;

     wp0 => work1 ;	wp1 => work2 ;
     call load_cache( n+1, vl, data(j), wp0%r, wp0%i, inc, jump ) ;

     call cosine_to_complex( n, vl, wp0%r, wp0%i, wp1%r, wp1%i, isign ) ;
     wps => wp0;		wp0 => wp1;		wp1 => wps;

     p = 1 ;		m = n ;
     do while( (mod(m,8)) == 0 )	
        m = m/8 ;
        call cache_radix_8( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*8 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,4)) == 0 )	
        m = m/4 ;
        call cache_radix_4( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*4 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,2)) == 0 )	
        m = m/2 ;
        call cache_radix_2( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*2 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,5)) == 0 )	
        m = m/5 ;
        call cache_radix_5( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*5 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,3)) == 0 )	
        m = m/3 ;
        call cache_radix_3( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*3 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     call complex_to_cosine( n, vl, wp0%r, wp0%i, wp1%r, wp1%i, e(n), isign ) ;	
     wps => wp0;		wp0 => wp1;		wp1 => wps;

     call store_cache( n+1, vl, wp0%r, wp0%i, data(j), inc, jump ) ;

     lot_ctr = lot_ctr - vl ;		j = j + vj ;	
  end do

End Subroutine cfct



                       !***************************   sine transform   ************************

Subroutine cfst (  &
  n,	&		! transform length 
  lot,	&		! vector length 
  data,	&	        ! source address 
  inc,	&		! source data stride 
  jump,	&		! source vector stride 
  e,	&		! trig table address 
  isign )			! direction of transform 
  use cache_data
  implicit none
  integer :: n,	    & 			! transform length 
       lot,   &			! vector length 
       isign, &			! direction of transform 
       inc,   &	                ! source data stride 
       jump			! source vector stride 
  complex(8), dimension(0:0) :: data		! source address 
  complex(8), dimension(0:0) :: e			! trig table address 
  

  type(scomplex), target :: work1 , work2
  type(scomplex), pointer :: wp0 , wp1 , wps
  integer :: p, m, j, vl, vj, lot_ctr ;

  vl = CACHE/(n+1) ;
  if( vl < 1 ) stop "Transform exceeds cache in cfst"
  vj = vl * jump ;		j = 0 ;

  lot_ctr = lot
  do while( lot_ctr > 0 )	
     if( lot_ctr < vl )	vl = lot_ctr ;

     wp0 => work1 ;	wp1 => work2 ;
     call load_cache( n+1, vl, data(j), wp0%r, wp0%i, inc, jump ) ;

     call sine_to_complex( n, vl, wp0%r, wp0%i, wp1%r, wp1%i ) ;
     wps => wp0;		wp0 => wp1;		wp1 => wps;

     p = 1 ;		m = n ;

     do while( (mod(m,8)) == 0 )	
        m = m/8 ;
        call cache_radix_8( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*8 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,4)) == 0 )	
        m = m/4 ;
        call cache_radix_4( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*4 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,2)) == 0 )	
        m = m/2 ;
        call cache_radix_2( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*2 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,5)) == 0 )	
        m = m/5 ;
        call cache_radix_5( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, isign ) ;
        p = p*5 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     do while( (mod(m,3)) == 0 )	
        m = m/3 ;
        call cache_radix_3( m, p, vl, wp0%r, wp0%i, wp1%r, wp1%i, e, -1 ) ;
        p = p*3 ;	wps => wp0;		wp0 => wp1;		wp1 => wps;		
     end do

     call complex_to_sine( n, vl, wp0%r, wp0%i, wp1%r, wp1%i, e(n) ) ;	
     wps => wp0;		wp0 => wp1;		wp1 => wps;

     call store_cache( n+1, vl, wp0%r, wp0%i, data(j), inc, jump ) ;

     lot_ctr = lot_ctr - vl ;		j = j + vj ;		
  end do

End Subroutine cfst
