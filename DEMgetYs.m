    function [YS1,YS2]= DEMgetYs(p,e,t,Vx,Vy,f,Comega,beta)
    timeASSEM= clock;
    nt= size(t,2);
    np= size(p,2);    
    
    i1= 1;
    i2= np;
    j1= np+1;
    j2= 2*np;
    C1= 1+beta;
    C2= (1+1/beta)*Comega^2;
     
    [S,g1x,g1y,g2x,g2y,g3x,g3y]=pdetrg(p,t);
        
    MDEM_POI= sparse(2*np,2*np); 
    
    v1= t(1,:);
    v2= t(2,:);
    v3= t(3,:);
    m= S/12;
    
    YSi_YSi= sparse(v1,v2,m,np,np)+ sparse(v2,v3,m,np,np)+...
             sparse(v3,v1,m,np,np)+ sparse(v1,v1,m,np,np)+...
             sparse(v2,v2,m,np,np)+ sparse(v3,v3,m,np,np);
    YSi_YSi= YSi_YSi + YSi_YSi.';  % diagonal elements doubles
     
    [ii,jj,NZ]= find(YSi_YSi);
    NZ= C1*NZ;
    MDEM_POI= MDEM_POI+ sparse(ii,jj,NZ,2*np,2*np)+ sparse(ii+np,jj+np,NZ,2*np,2*np);
    clear YSi_YSi;    
      
    YS1x_YS1x= sparse(v1,v2,g1x.*g2x.*S,np,np)+...
               sparse(v2,v3,g2x.*g3x.*S,np,np)+...
               sparse(v3,v1,g3x.*g1x.*S,np,np);
    YS1x_YS1x= YS1x_YS1x + YS1x_YS1x.';
    YS1x_YS1x= YS1x_YS1x + sparse(v1,v1,g1x.^2.*S,np,np)+...
                           sparse(v2,v2,g2x.^2.*S,np,np)+...
                           sparse(v3,v3,g3x.^2.*S,np,np);
    
    [ii,jj,NZ]= find(YS1x_YS1x);
    NZ= C2*NZ;
    MDEM_POI= MDEM_POI + sparse(ii,jj,NZ,2*np,2*np);
    clear YS1x_YS1x;   
     
    YS2y_YS2y= sparse(v1,v2,g1y.*g2y.*S,np,np)+...
               sparse(v2,v3,g2y.*g3y.*S,np,np)+...
               sparse(v3,v1,g3y.*g1y.*S,np,np);
    YS2y_YS2y= YS2y_YS2y + YS2y_YS2y.';
    YS2y_YS2y= YS2y_YS2y + sparse(v1,v1,g1y.^2.*S,np,np)+...
                           sparse(v2,v2,g2y.^2.*S,np,np)+...
                           sparse(v3,v3,g3y.^2.*S,np,np);

    [ii,jj,NZ]= find(YS2y_YS2y);
    NZ= C2*NZ;
    MDEM_POI= MDEM_POI + sparse(ii+np,jj+np,NZ,2*np,2*np);
    clear YS2y_YS2y;
    
    YS1x_YS2y= sparse(v1,v2,g1x.*g2y.*S,np,np)+ sparse(v2,v3,g2x.*g3y.*S,np,np)+...
               sparse(v3,v1,g3x.*g1y.*S,np,np)+ sparse(v2,v1,g2x.*g1y.*S,np,np)+...
               sparse(v3,v2,g3x.*g2y.*S,np,np)+ sparse(v1,v3,g1x.*g3y.*S,np,np)+...
               sparse(v1,v1,g1x.*g1y.*S,np,np)+ sparse(v2,v2,g2x.*g2y.*S,np,np)+...
               sparse(v3,v3,g3x.*g3y.*S,np,np);
    
    [ii,jj,NZ]= find(YS1x_YS2y);
    NZ= C2*NZ;
    MDEM_POI= MDEM_POI + sparse(ii,jj+np,NZ,2*np,2*np) + sparse(jj+np,ii,NZ,2*np,2*np);
    clear YS1x_YS2y;

    FVx= zeros(np,1);
    FVy= zeros(np,1);      
    FYS1x= zeros(np,1);
    FYS2y= zeros(np,1);    
    for T=1:nt
        FVx(v1(T))= FVx(v1(T)) + Vx(T).*S(T)/3;
        FVx(v2(T))= FVx(v2(T)) + Vx(T).*S(T)/3;
        FVx(v3(T))= FVx(v3(T)) + Vx(T).*S(T)/3;        
        
        FVy(v1(T))= FVy(v1(T)) + Vy(T).*S(T)/3;
        FVy(v2(T))= FVy(v2(T)) + Vy(T).*S(T)/3;
        FVy(v3(T))= FVy(v3(T)) + Vy(T).*S(T)/3;        
        
        FYS1x(v1(T))= FYS1x(v1(T)) - f(T).*g1x(T).*S(T);
        FYS1x(v2(T))= FYS1x(v2(T)) - f(T).*g2x(T).*S(T);
        FYS1x(v3(T))= FYS1x(v3(T)) - f(T).*g3x(T).*S(T);

        FYS2y(v1(T))= FYS2y(v1(T)) - f(T).*g1y(T).*S(T);
        FYS2y(v2(T))= FYS2y(v2(T)) - f(T).*g2y(T).*S(T);
        FYS2y(v3(T))= FYS2y(v3(T)) - f(T).*g3y(T).*S(T);
    end    
    clear S g1x g2x g3x g1y g2y g3y;
        
    FDEM_POI= zeros(2*np,1);    
    FDEM_POI(i1:i2)= C1* FVx + C2* FYS1x;
    FDEM_POI(j1:j2)= C1* FVy + C2* FYS2y;
    clear FVx FVy FYS1x FYS2y;
    
 %   condM= condest(MDEM_POI);
    fprintf('TIME:   %g sec - assembling\n',etime(clock,timeASSEM));
 %   disp(sprintf('CONDITION NUMBER FOR MATRIX OF SYSTEM %f',condM));    

    X= MDEM_POI\FDEM_POI;
    
    YS1= X(i1:i2,1);
    YS2= X(j1:j2,1);
    clear MDEM_POI FDEM_POI X;