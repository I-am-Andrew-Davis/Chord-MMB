C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

c ----------------------------------------------------------------------
c  This file includes canonical flow cases with analytical/known solutions
c  for code verification. 
c ----------------------------------------------------------------------

#include "CONSTANTS.H"
#include "REAL.H"

 
      subroutine derivF(F, Y)
      implicit none
      integer, parameter :: N = 3
      Real_t Y(1:N), F(1:N)
      F(1) = -Y(1)*Y(3)
      F(2) = Y(1)
      F(3) = Y(2)
      return
      end

      
      subroutine RK4(F, Y, H)
      implicit none
      integer, parameter :: N = 3
      integer K
      Real_t H
      Real_t Y(1:N), F(1:N)
      Real_t Yhat(1:N), Ytilda(1:N), Ybar(1:N), Ynxt(1:N)
      Real_t Fhat(1:N), Ftilda(1:N), Fbar(1:N)


      CALL derivF(F, Y)
      do K = 1, N
         Yhat(K) = Y(K) + 0.5*H*F(K)
      enddo

      CALL derivF(Fhat, Yhat)
      do K = 1, N
         Ytilda(K)  = Y(K) + 0.5*H*Fhat(K)
      enddo

      CALL derivF(Ftilda, Ytilda)
      do K = 1, N
         Ybar(K) = Y(K) + H*Ftilda(K)
      enddo

      CALL derivF(Fbar, Ybar)
      do K = 1, N
         Ynxt(K) =  Y(K) + 1.0/six*H*
     &        (F(K) + 2.0*(Fhat(K)+Ftilda(K))+Fbar(K))
         Y(K) = Ynxt(K) 
      enddo
      
      return
      end
