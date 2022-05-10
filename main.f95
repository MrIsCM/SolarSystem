program sistema_solar
	implicit none

	! Posicion y derivadas
	integer, parameter :: n = 10
	double precision :: x(1:n), y(1:n), vx(1:n), vy(1:n), ax(1:n), ay(1:n), wx(1:n), wy(1:n) 	
	double precision :: m(1:n)									! Masa
	double precision :: t, h, c, G, Ms				 			! Tiempo y parametros de utilidad
	double precision :: f_r, f_t, f_v 							! Factor de escala para r, t y v
	double precision :: E_t, L_t			! Energia mecanica y momento angular total

	integer :: i
	integer, parameter :: iter = 10000

	!###################################################
	! 		Definicion de algunos parametros
	!###################################################
	c = 1.493 * 10.0**11 		! m
	G = 6.67 * 10.0**(-11) 		! m^3* kg^-1* s^-2 
	Ms = 1.99 * 10.0**30 		! kg

	f_r = c 					! m
	f_t = 5022361.106 			! s
	f_v = 29786.7869			! m/s

	t=0
	h=0.001


	!###################################
	! 		Condiciones iniciales
	!###################################
	call condiciones_iniciales(x, y, vx, vy, m, f_r, f_v)

	open(10, file='Datos2/README.txt', status='unknown')
		write(10,*) "Datos simulacion del Sistema Solar - Fortran95"
		write(10,*) "h = 	", h, "Intervalo de tiempo (~h*58 dias)"
		write(10,*) "i = 	", iter, "Numero de iteraciones"

	open(1, file='Datos2/Posiciones.txt', status='unknown')
	!open(2, file='Datos2/Velocidades.txt', status='unknown')
	!open(3, file='Datos2/Aceleraciones.txt', status='unknown')
	open(4, file='Datos2/Constantes.txt', status='unknown')

	do i = 1, iter
		write(1,*) t, x, y
		!write(2,*) t, vx, vy
		!write(3,*) t, ax, ay
		call Ener_mec(x, y, vx, vy, m, E_t)
		call momento_angular(x, y, vx, vy, m, L_t)
		write(4,*) t, E_t, L_t
		call aceleraciones(x, y, m, ax, ay)				! Aceleraciones para las posiciones de la iter anterior
		call vervelet_pos(x, y, vx, vy, ax, ay, h)		! Calculo nuevas posiciones
		call w(vx, vy, ax, ay, wx, wy, h)				! Calculo la f.aux w
		call aceleraciones(x, y, m, ax, ay)				! Calculo las aceleraciones a partir de las nuevas pos
		call vervelet_vel(vx, vy, ax, ay, wx, wy, h)	! Calculo las nuevas velocidades mediante las ac nuevas

		t = t + h
	end do 												! Repite

	close(1)
	!close(2)
	!close(3)
	close(4)
	
end program sistema_solar


subroutine condiciones_iniciales(x, y, vx, vy, m, f_r, f_v)
	implicit none

	integer, parameter :: n = 10
	double precision,intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n), m(1:n)
	double precision, intent(in) :: f_r, f_v

	double precision :: afelio(1:n), perihelio(1:n)			! Datos (posicion)

	integer :: i

	!###############################################
	!
	! 		Condiciones iniciales 
	!
	!###############################################
	! Sistema compuesto por 10 masas:
	! 1- Sol, 2- Mercurio, 3-Venus, 4-Tierra, 5-Marte, 
	! 6-Jupiter, 7-Saturno, 8-Urano, 9-Neptuno, 10-Pluton

	! Posicion inicial:
	! La coord y de todos las masas es 0, suponiendo así que se encuentran en el afelio o perihelio
	! Datos tomados en: https://ergodic.ugr.es/cphys/LECCIONES/ssolar/planetasdatos.html [1]

	! Perihelio y afelio. masas 2 a 10 (x10^6 km)
	! 46.0	107.5	147.1	206.6	740.5	1352.6	2741.3	4444.5	4435.0
	! 69.8	108.9	152.1	249.2	816.6	1514.5	3003.6	4545.7	7304.3
	
	y = 0

	!Sol inicalmente en el origen
	x(1) = 0

	! Perihelio: Tomando los valores de la tabla [1] y reescalando
	perihelio = (/0.0, 46.0, 107.5, 147.1, 206.6, 740.5, 1352.6, 2741.3, 4444.5, 4435.0/)
	perihelio = perihelio*10.0**9/f_r 


	! Afelio
	afelio = (/0.0, 69.8, 108.9, 152.1, 249.2, 816.6, 1514.5, 3003.6, 4545.7, 7304.3/)
	afelio = afelio*10.0**9/f_r 

	! Velocidades orbitales promedio:
	! Orbital Velocity (km/s):  47.9, 35.0, 29.8, 24.1, 13.1, 9.7, 6.8, 5.4, 4.7

	vx = 0
	! Sol inicialmente en reposo
	vy = (/ 0.0, 47.9, 35.0, 29.8, 24.1, 13.1, 9.7, 6.8, 5.4, 4.7/)		!km/s
	! Reescalando
	vy = vy*10**3/f_v 
 


	! Asigno las posiciones iniciales intercalando afelio y perihelio.
	! Supongo el origen tal que el perihelio queda a la izquierda (<0) 
	! y el afelio a la derecha (>0)

	! Dada la posicion incial adjudicada las velocidades serán:
	! vx = 0, para todas las masas y vy, positiva para las masas situadas en el afelio (derecha)
	! y negativa para las masas situadas en el perihelio (izq)
	! (Rotaciones en sentido antihorario)

	do i = 2, n
		if ( mod(i,2) /= 0 ) then
			x(i) = -perihelio(i)
			vy(i) = -vy(i)
		else 
			x(i) = afelio(i)
		end if 
	end do


	! Masas (x10^24 kg ): 1990000 ,0.330, 4.87, 5.97, 0.642, 1899, 568, 86.8, 102, 0.0125
	! Reescalamiento: m' = m/Ms
	m = (/ 1990000.0 ,0.330, 4.87, 5.97, 0.642, 1899.0, 568.0, 86.8, 102.0, 0.0125/)
	m = m/m(1)

end subroutine condiciones_iniciales

subroutine aceleraciones(x, y, m, ax, ay)
	implicit none

	integer, parameter :: n = 10
	double precision, intent(in) :: x(1:n), y(1:n), m(1:n)
	double precision, intent(out) :: ax(1:n), ay(1:n)

	double precision, parameter :: G = 6.67 * 10.0**(-11)

	integer :: i, j

	ax = 0
	ay = 0
	
	do i = 1, n
		do j = 1, n
			if (i /= j) then
				ax(i) = ax(i) - m(j)*(x(i)-x(j))/((x(i)-x(j))**2 + (y(i)-y(j))**2)**(3.0/2.0)
				ay(i) = ay(i) - m(j)*(y(i)-y(j))/((x(i)-x(j))**2 + (y(i)-y(j))**2)**(3.0/2.0)
			end if 
		end do
	end do

end subroutine aceleraciones

subroutine vervelet_pos(x, y, vx, vy, ax, ay, h)
	implicit none

	integer, parameter :: n=10
	double precision,intent(in) :: h, ax(1:n), ay(1:n)
	double precision, intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n)

	integer :: i

	do i = 1, n
		x(i) = x(i) + h*vx(i) + 0.5*h**2*ax(i)
		y(i) = y(i) + h*vy(i) + 0.5*h**2*ay(i)
	end do

end subroutine vervelet_pos


subroutine w(vx, vy, ax, ay, wx, wy, h)
	implicit none

	integer, parameter :: n = 10
	double precision, intent(in) :: h, vx(1:n), vy(1:n), ax(1:n), ay(1:n)
	double precision, intent(inout) :: wx(1:n), wy(1:n)
	
	integer :: i
	
	do i = 1, n
		wx(i) = vx(i) + 0.5*h*ax(i)
		wy(i) = vy(i) + 0.5*h*ay(i)
	end do
end subroutine w

subroutine vervelet_vel(vx, vy, ax, ay, wx, wy, h)
	implicit none

	integer, parameter :: n = 10
	double precision, intent(in) :: h,wx(1:n), wy(1:n), ax(1:n), ay(1:n)
	double precision, intent(out) :: vx(1:n), vy(1:n)

	integer :: i 

	do i = 1, n
		vx(i) = wx(i) + 0.5*h*ax(i)
		vy(i) = wy(i) + 0.5*h*ay(i) 
	end do

end subroutine vervelet_vel

subroutine Ener_mec(x, y, vx, vy, m, E_mec)
	implicit none

	integer, parameter :: n = 10
	double precision, intent(in) :: x(1:10), y(1:10), vx(1:10), vy(1:10), m(1:10)
	double precision, intent(out) :: E_mec
	double precision :: T, V

	integer :: i, j

	T = 0
	V = 0

	do i = 1, n
		T = T + 0.5*m(i)*(vx(i)**2 + vy(i)**2)
	end do
	do i = 1, n
		do j = 1, n
			if (i /= j) then
				V = V - m(i)*m(j)/((x(i)-x(j))**2+(y(i)-y(j))**2)**0.5
			end if
		end do
	end do

	E_mec = T + 0.5*V			! La mitad de E. Potencial ya que la E.Pot de i,j es = a j,i

	
end subroutine Ener_mec

subroutine momento_angular(x, y, vx, vy, m, L)
	implicit none
	double precision,intent(in) :: x(1:10), y(1:10), vx(1:10), vy(1:10), m(1:10)
	double precision,intent(out) ::  L

	integer, parameter :: n = 10
	integer :: i
	L = 0

	do i = 1, n
		L = L + x(i)*m(i)*vy(i) - y(i)*m(i)*vx(i)
	end do

end subroutine momento_angular

subroutine periodo(y, t, period_orbital)
	implicit none
	double precision ,intent(in) :: y, t 
	double precision ,intent(out) :: period_orbital
end subroutine periodo