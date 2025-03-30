% these parameters are not in use here but can be used for other examples
A1=0;
A2=0;
A3=0;

U= inline('x^2*(1-x)^2 * y^2*(1-y)^2','x','y','A1','A2','A3');
Ux= inline('(2*x-6*x^2+4*x^3) * y^2*(1-y)^2','x','y','A1','A2','A3');
Uy= inline('(2*y-6*y^2+4*y^3) * x^2*(1-x)^2','x','y','A1','A2','A3');
LaplU= inline('(2-12*x+12*x^2) * y^2*(1-y)^2 + (2-12*y+12*y^2) * x^2*(1-x)^2','x','y','A1','A2','A3');
F= inline('24*( x^2*(1-x)^2 + y^2*(1-y)^2 ) + 2*(2-12*x+12*x^2)*(2-12*y+12*y^2)','x','y','A1','A2','A3');