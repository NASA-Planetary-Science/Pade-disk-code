!----------------------------------------------------------------------------------85

call brady_livescu_T4_coeffs

end program

!----------------------------------------------------------------------------------85

subroutine brady_livescu_T4_coeffs

! Fourth-order tridiagonal with third-order boundary scheme.

implicit none
real(8) :: alpha03, alpha13
real(8) :: beta0p1
real(8) :: alpha00, alpha01, alpha02
real(8) :: alpha10, alpha11
real(8) :: alpha21, alpha22, alpha23
real(8) :: beta1m1, beta1p1, beta2p1
real(8) :: w0, w1, w2

alpha03 = -2.63678180555208
alpha13 = -0.201231074955226

w0 = -(12.d0*alpha13 - 9.d0)/((96.d0*alpha03 - 8.d0)*alpha13 + &
     48.d0*alpha03+24.d0)
w1 = (60.d0*alpha03 + 9.d0)/((48.d0*alpha03 - 4.d0)*alpha13 + 24.d0*alpha03 + 12.d0)
w2 = -((12.d0*alpha03 + 8.d0)*alpha13 - 9.d0*alpha03-6.d0)/&
     ((24.d0*alpha03 - 2.d0)*alpha13 + 12.d0*alpha03 + 6.d0)

beta0p1 = 2.d0 - 6.d0*alpha03
alpha00 = (4.d0*alpha03 - 5.d0)/2.d0
alpha01 = 3.d0*alpha03 + 2.d0
alpha02 = -(12.d0*alpha03 - 1.d0)/2.d0
alpha10 = -(4.d0*alpha13 + 3.d0)/4.d0

beta1m1 = (3.d0*alpha13 + 1.d0)/4.d0
alpha11 = 0.d0
beta1p1 = (9.d0*alpha13 + 1.d0)/4.d0

alpha12 = 3.d0/4.d0
alpha20 = 0.d0
beta2m1 = 1.d0/4.d0

alpha21 = -3.d0/4.d0
beta2p1 = 1.d0/4.d0

alpha22 = 0.d0
alpha23 = 3.d0/4.d0

print *, ' Brady-Livescu T4 coefficients are:'

print *, ' LHS row 1'
print *, ' beta0p1 = ', beta0p1
print *, ' '

print *, ' LHS row 2'
print *, ' beta1m1 = ', beta1m1
print *, ' beta1p1 = ', beta1p1
print *, ' '

print *, ' LHS interior scheme:'
print *, ' beta2m1 = ', beta2m1
print *, ' beta2p1 = ', beta2p1
print *, ' '

print *, ' RHS row 1'
print *, ' alpha00 = ', alpha00
print *, ' alpha01 = ', alpha01
print *, ' alpha02 = ', alpha02
print *, ' alpha03 = ', alpha03

print *, ' RHS row 2'
print *, ' alpha10 = ', alpha10
print *, ' alpha12 = ', alpha12
print *, ' alpha13 = ', alpha13

print *, ' RHS interior scheme'
print *, ' alpha21 = ', alpha21
print *, ' alpha23 = ', alpha23

print *, ' Weights'
print *, ' w0 = ', w0
print *, ' w1 = ', w1
print *, ' w2 = ', w2

end subroutine brady_livescu_T4_coeffs

!----------------------------------------------------------------------------------85
