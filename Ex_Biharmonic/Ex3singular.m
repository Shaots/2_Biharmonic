A1= 10;
A2= 0;
A3= 2;

expx= 'exp(-A1*(x-A2)^2)';
phi =   strcat(expx,' * ','sin(pi*x)^2/pi^2');
phi_x = strcat(expx,' * ','(-2*A1*(x-A2)*sin(pi*x)^2/pi^2 + sin(2*pi*x)/pi)');
phi_2x= strcat(expx,' * ','((4*A1^2*(x-A2)^2 -2*A1)* sin(pi*x)^2/pi^2 - 4*A1*(x-A2)*sin(2*pi*x)/pi + 2*cos(2*pi*x) )');

expr1= ' ( 16*A1^4*(x-A2)^4-48*A1^3*(x-A2)^2+12*A1^2 )* sin(pi*x)^2/pi^2 ';
expr2= ' 4* ( - 8*A1^3*(x-A2)^3 + 12*A1^2*(x-A2) )* sin(2*pi*x)/pi ';
expr3= ' 6* ( 4*A1^2*(x-A2)^2 - 2*A1 )* 2*cos(2*pi*x) ';
expr4= ' 4* (-2*A1*(x-A2))* (-4*pi*sin(2*pi*x)) ';
expr5= ' -8*pi^2* cos(2*pi*x) ';
phi_4x= strcat(expx,' * ','(',expr1,'+',expr2,'+',expr3,'+',expr4,'+',expr5,')');

psi = 'sin(A3*pi*y)^2/(A3^2*pi^2)';
psi_y= 'sin(2*A3*pi*y)/(A3*pi)';
psi_2y= '2*cos(2*A3*pi*y)';
psi_4y= '-8*A3^2*pi^2*cos(2*A3*pi*y)';

ifunc= strcat(phi,' * ',psi);
U= inline(ifunc,'x','y','A1','A2','A3');

ifunc= strcat(phi_x,' * ',psi);
Ux= inline(ifunc,'x','y','A1','A2','A3');

ifunc= strcat(phi,' * ',psi_y);
Uy= inline(ifunc,'x','y','A1','A2','A3');
    
ifunc= strcat(phi_2x,' * ',psi,' + ',phi,' * ',psi_2y);
LaplU= inline(ifunc,'x','y','A1','A2','A3');

ifunc= strcat(phi_4x,' * ',psi,' + ',' 2* ',phi_2x,' * ',psi_2y,' + ',phi,' * ',psi_4y);
F= inline(ifunc,'x','y','A1','A2','A3');