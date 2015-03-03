c---------------------------------------------------------

      program main
            implicit none
            real*8 phi
            real*8 gamma, lambda, epsilon
            real*8 a, b, c
            real*8 u1, u2, x1, x2, x3, x4
            real*8 y1, y2
            real*8 causticY1, causticY2
            real*8 twopi

            open(unit=1, file='critical_line.dat', status='unknown')
            open(unit=2, file='caustic_curve.dat', status='unknown')

            twopi = 2*4*atan(1.0d0)
            epsilon = 0.5d0

            write(*,*)'Enter Gamma: '
            read(*,*)gamma

            do while(phi .le. twopi)
                  lambda = cos(phi)**2.0d0 - sin(phi)**2.0d0
                  a = gamma*epsilon*lambda
                  b = sqrt((gamma**2.0d0)*(lambda**2.0d0-1.0d0)+1.0d0)
                  c = 1 - gamma**2.0d0

                  u1 = (a + b)/c
                  u2 = (a - b)/c
                  
                  if (u1 .gt. 0.0d0) then
                        x1 = sqrt(u1)
                        x2 = -sqrt(u1)
                        write(1,*) x1*cos(phi), x1*sin(phi)
                        write(1,*) x2*cos(phi), x2*sin(phi)

                        y1 = causticY1(gamma, epsilon, x1, lambda)
                        y2 = causticY2(gamma, epsilon, x1, lambda)
                        write(2,*) y1*cos(phi), y2*sin(phi)

                        y1 = causticY1(gamma, epsilon, x2, lambda)
                        y2 = causticY2(gamma, epsilon, x2, lambda)
                        write(2,*) y1*cos(phi), y2*sin(phi)
                  end if

                  if (u2 .gt. 0.0d0) then
                        x3 = sqrt(u2)
                        x4 = sqrt(u2)
                        write(1,*) x3*cos(phi), x3*sin(phi)
                        write(1,*) x4*cos(phi), x4*sin(phi)

                        y1 = causticY1(gamma, epsilon, x3, lambda)
                        y2 = causticY2(gamma, epsilon, x3, lambda)
                        write(2,*) y1*cos(phi), y2*sin(phi)

                        y1 = causticY1(gamma, epsilon, x4, lambda)
                        y2 = causticY2(gamma, epsilon, x4, lambda)
                        write(2,*) y1*cos(phi), y2*sin(phi)
                  endif
                  
                  phi = phi + 0.0001
            end do

            close(unit=1)
            close(unit=2)
      end

c--------------------------------------------------------

      double precision function causticY1(gamma, epsilon, x, lambda)
            implicit none
            real*8 gamma, epsilon, x, lambda
            real*8 cosphi

            cosphi = sqrt((1 + lambda)/2.0d0)
            causticY1 = ((1 + gamma)*x - epsilon/x)*cosphi
            return
      end

      double precision function causticY2(gamma, epsilon, x, lambda)
            implicit none
            real*8 gamma, epsilon, x, lambda
            real*8 sinphi

            sinphi = sqrt((1 - lambda)/2.0d0)
            causticY2 = ((1 - gamma)*x - epsilon/x)*sinphi 
            return
      end



