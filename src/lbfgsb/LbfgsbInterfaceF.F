C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

      subroutine lbfgsb_interface(n, m, x, l, u, nbd, f, g, factr,
     &     pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)

      use iso_c_binding, only: C_CHAR, c_null_char
      implicit none

      character (kind=c_char, len=1), dimension (60), intent (inout) :: csave
      character (kind=c_char, len=1), dimension (60), intent (inout) :: task
      integer, intent(inout) :: n, m, iprint
      integer, intent(inout) :: isave(44), iwa(3*n), lsave(4), nbd(n)
      double precision, intent (inout) :: f, factr, pgtol, dsave(29)
      double precision, intent (inout) :: x(n), l(n), u(n), g(n)
      double precision, intent (inout) :: wa(2*m*n + 5*n + 11*m*m + 8*m)
      character (len=60) :: fort_task, fort_csave
      integer :: i

C     I think if I format the c strings correctly then I won't need this
C     conversion. Need to pad the desired values to 60 spaces, no null terminator
c$$$      print*,"Len of task: ",len(task),"fort_version: ",len(fort_task)
      fort_task = " "
      loop_task: do i=1, len(fort_task)
      if ( task (i) == c_null_char ) then
         exit loop_task
      else
         fort_task(i:i) = task(i)
      end if
      end do loop_task

      fort_csave= " "
      loop_csave: do i=1, len(fort_csave)
      if ( csave (i) == c_null_char ) then
         exit loop_csave
      else
         fort_csave(i:i) = csave(i)
      end if
      end do loop_csave

c$$$      print*,"Current task value: ",task,"fort_version: ",fort_task
c$$$      print*,"Current csave value: ",csave,"fort_version: ",fort_csave


      call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa,
     &     iwa, fort_task, iprint, fort_csave, lsave, isave, dsave)

c$$$      print*,"Task after call: ",task,"fort_version: ",fort_task
c$$$      print*,"csave after call: ",csave,"fort_version: ",fort_csave
C     likewise, if the C++ char formatting is perfect, this might not be needed
      do i = 1, len(fort_task)
         task(i) = fort_task(i:i)
         csave(i) = fort_csave(i:i)
      end do

      return
      end
