      program main
            implicit none

            call runge_kutta_nystrom(0.01d0, 10000, 0.0d0)
      end

      subroutine runge_kutta_nystrom(h, N, x0)
            implicit none
            real*8 h, x0
            real*8 k1, k2, k3, k4, K, L
            real*8 x, y, yd, f
            real*8 i, hon2
            integer N

            hon2 = h/2.0d0
            x = x0
            y = 0.0d0
            yd = 1.0d0

            do i = 0, N-1, +1
                  k1 = hon2*f(x, y, yd)

                  K = hon2*(yd + k1/2.0d0)
                  k2 = hon2*f(x + hon2, y + K, yd + k1)

                  k3 = hon2*f(x + hon2, y + K, yd + k2)

                  L = h*(yd + k3)
                  k4 = hon2*f(x + h, y + L, yd + 2.0d0*k3)

                  x = x + h
                  y = y + h*(yd + (k1 + k2 + k3)/3.0d0 )

                  yd = yd + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/3.0d0

                  write(*,*)x, y
            end do
      
      end

      double precision function f(x, y, yd)
            implicit none
            real*8 f2, f3
            real*8 x, y, yd

            f = (f2(x,y,yd)/f3(x,y,yd))**2.0d0
            return
      end

      double precision function f1(x, y, yd)
            implicit none
            real*8 x, y, yd

            f1 = (-3.0d0/(x + 1))*yd
            return
      end

       double precision function f2(x, y, yd)
            implicit none
            real*8 x, y, yd

            f2 = ((-7.0d0/2.0d0)*(x+1)*yd-3.0d0*y/2.0d0)/((x+1)**2.0d0)
            return
      end

       double precision function f3(x, y, yd)
            implicit none
            real*8 x, y, yd

            f3 = ((-7.0d0/2.0d0)*(x+1)*yd)/((x+1)**2.0d0)
            return
      end

