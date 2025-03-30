% these parameters are not in use here but can be used for other examples
A1=0;
A2=0;
A3=0;

U= inline('sin(pi*x)^2*sin(pi*y)^2/pi^4','x','y','A1','A2','A3');
Ux= inline('sin(2*pi*x)*sin(pi*y)^2/pi^3','x','y','A1','A2','A3');
Uy= inline('sin(2*pi*y)*sin(pi*x)^2/pi^3','x','y','A1','A2','A3');
LaplU= inline('2*(cos(2*pi*x)*sin(pi*y)^2 + cos(2*pi*y)*sin(pi*x)^2)/pi^2','x','y','A1','A2','A3');
F= inline('8*(-cos(2*pi*x)*sin(pi*y)^2 + cos(2*pi*x)*cos(2*pi*y) - cos(2*pi*y)*sin(pi*x)^2)','x','y','A1','A2','A3');