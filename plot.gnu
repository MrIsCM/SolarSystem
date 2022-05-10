#set output 'Images/Tierra_Sol.png'


stats "Datos2\\Posiciones.txt" nooutput

set xrange [-1.1:1.1]
set yrange [-1.1:1.1]

f(x) = A*sin(f*x+i)+h
A=1
f=1
i=1
h=0.1

fit f(x) "Datos2\\Posiciones.txt" using 1:5 via A, f, i, h

#plot "Datos2\\Posiciones.txt" using 2:12 ps 0.3 title "Sol",\
	#"Datos2\\Posiciones.txt" using 5:15 ps 0.3 title "Tierra", \
	#f(x)
	# "Datos3\\Posiciones.txt" using 3:13 ps 0.3 title "Mercurio", \
	# "Datos3\\Posiciones.txt" using 4:14 ps 0.3 title "Venus", \
	# "Datos3\\Posiciones.txt" using 6:16 ps 0.3 title "Marte", \
	# "Datos3\\Posiciones.txt" using 7:17 ps 0.3 title "Jupiter"
	# "Datos\\Posiciones.txt" using 9:19, \
	# "Datos\\Posiciones.txt" using 10:20,\
	# "Datos\\Posiciones.txt" using 11:21