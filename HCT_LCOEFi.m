% Book of Bernadou and Boisserie 
% The FEM Method in Thin Shell Theory: Application to Arch Dam Simulations
% Warning:
%           X, Y, W, Wx, Wy are given in global order
%
function [LCOEF]=HCT_LCOEFi(X,Y,W,Wx,Wy,i)  
                         % value of i is 1,2,3
order= [1,2,3,1,2];
x= [X(order(i:i+2))];
y= [Y(order(i:i+2))];
w= [W(order(i:i+2))];
wx= [Wx(order(i:i+2))];
wy= [Wy(order(i:i+2))];

%% All values are calculated assuming that i=1
%% but due to changing of order from global to local 
%% we get right results for any subtriangle

% 1x9  global degrees of freedom
GDOF= [w(1),w(2),w(3),wx(1),wy(1),wx(2),wy(2),wx(3),wy(3)];
                          
% D - tridiagonal matrix that converts GDOF to LDOF
D= diag([1,1,1, x(3)-x(1),y(2)-y(1), x(1)-x(2),y(3)-y(2), x(2)-x(3),y(1)-y(3)], 0);
D= D + diag([0,0,0, x(2)-x(1),0, x(3)-x(2),0, x(1)-x(3)],  1);
D= D + diag([0,0,0, y(3)-y(1),0, y(1)-y(2),0, y(2)-y(3)], -1);

% 1x9  local degrees of freedom                         
LDOF= GDOF*D;

l(1)= sqrt( (x(2)-x(3))^2 + (y(2)-y(3))^2 );
l(2)= sqrt( (x(3)-x(1))^2 + (y(3)-y(1))^2 );
l(3)= sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 );

et(1)= ( l(3)^2-l(2)^2 )/( l(1)^2 );
et(2)= ( l(1)^2-l(3)^2 )/( l(2)^2 );
et(3)= ( l(2)^2-l(1)^2 )/( l(3)^2 );

% AR is 9x10 matrix of coefficients of shape functions in barycentric representation
AR= zeros(9,10);

AR(1,:)= [ -0.5*(et(2)-et(3)) , 0 , 0 , 1.5*(3+et(2)) , 1.5*(3-et(3)) ,...
            0 , 0 , 0 , 0 , 0 ];
    
AR(2,:)= [ 0.5*(1-2*et(1)-et(3)) , 1 , 0 , -1.5*(1-et(1)) , 1.5*(et(1)+et(3)) ,...
            3 , 3 , 0 , 0 , 3*(1-et(1)) ];
    
AR(3,:)= [ 0.5*(1+2*et(1)+et(2)) , 0 , 1 , -1.5*(et(1)+et(2)) , -1.5*(1+et(1)) ,...
            0 , 0 , 3 , 3 , 3*(1+et(1)) ];
    
AR(4,:)= [ -0.25*(1+et(2)) , 0 , 0 , 0.25*(5+3*et(2)) , 0.5 ,...
            0 , 0 , 0 , 0 , 0 ];
    
AR(5,:)= [ -0.25*(1-et(3)) , 0 , 0 , 0.5 , 0.25*(5-3*et(3)) ,...
            0 , 0 , 0 , 0 , 0 ];
    
AR(6,:)= [ 0.25*(1-et(3)) , 0 , 0 , -0.5 , -0.25*(1-3*et(3)) ,...
            1 , 0 , 0 , 0 , 1 ];

AR(7,:)= [ -0.5*et(1) , 0 , 0 , -0.25*(1-3*et(1)) , 0.25*(1+3*et(1)) ,... 
            0 , 1 , 0 , 0 , 0.5*(1-3*et(1)) ];
    
AR(8,:)= [ 0.5*et(1) , 0 , 0 , 0.25*(1-3*et(1)) , -0.25*(1+3*et(1)) ,...
            0 , 0 , 1 , 0 , 0.5*(1+3*et(1)) ];
    
AR(9,:)= [ 0.25*(1+et(2)) , 0 , 0 , -0.25*(1+3*et(2)) , -0.5 ,...
            0 , 0 , 0 , 1 , 1 ];

% LCOEF -coefficients of polynomial in barycentric coordinates
LCOEF= LDOF*AR;