               % values of these parameters are set in files of examples
global A1;
global A2;
global A3;

clc;
clear;
format compact;
diary(strcat('DEMExRep', datestr(datenum(clock),'yyyymmddTHHMMSS'), '.txt'));
beta1= 1;
beta2= 1;
% for unit square i.e any domain inside of unit square can be used
Comega= 1/(sqrt(2)*pi)

%============================================================================
disp('----- Step 1: Selection of the domain - the unit square');
              % 1) create geometry via pdetool 
              % 2) export it by Draw/Export Geometry Description to workspace
              % 3) dl= decsg(gd) - function analyzes model that you draw
              % 4) wgeom(dl,'unsquareg') - write geometry specification 
              % 5) by analogy export boundary conditions b 
              % 6) wbound(b,'file') - write boundary conditions
g='unsquareg';
b='unsquareb';

%============================================================================
% U, Ux,... , F are given as inline functions

disp('----- Step 2: Selection of the exact solution U and F');
[mfile, mpath] = uigetfile('*.m', 'open file with corresponding data');
file= strcat(mpath,mfile)
feval('run', file);
             
disp(sprintf('\ntypes of initial refinement: longest(times) - regular(times):'));
disp(sprintf('\ncases 1) 8-0, 2) 0-4, 3) 10-0'));
tref= input('input a type of refinement: ');

switch(tref)
    case 1
        Nlon= 8;
        Nreg= 0;
    case 2
        Nlon= 0;
        Nreg= 4;
    case 3
        Nlon= 10;
        Nreg= 0;        
    otherwise
        Nlon= 0;
        Nreg= 0;
end        

[p,e,t]=initmesh(g,'Hmax',inf);
for i=1:Nlon
    [p,e,t]=refinemesh(g,p,e,t,'longest');
end    
for i=1:Nreg
    [p,e,t]=refinemesh(g,p,e,t,'regular');
end    

refinement= 0;
RelErr= 100;

% stoping criterion is approx. 4% of the relative error
while ( RelErr > 4 ) 
    np= size(p,2);
    nt= size(t,2);    
    disp(sprintf('\nREFINEMENT %d -> mesh: %d nodes, %d elements',refinement,np,nt));
    
%  A1,A2,A3 are parameters in the exact solution
%  are defined inside of selected example
   
    timeIn= clock;
    for i=1:np    
        x= p(1,i);
        y= p(2,i);
% column representation for values in points
        u(i,1)= U(x,y,A1,A2,A3); 
        ux(i,1)= Ux(x,y,A1,A2,A3);
        uy(i,1)= Uy(x,y,A1,A2,A3);
        Delu(i,1)= LaplU(x,y,A1,A2,A3);
        Fp(i,1)= F(x,y,A1,A2,A3);                                                     
    end
    disp(sprintf('\nTIME:   %g sec - calc. of inline functions',etime(clock,timeIn)));      
    
% approximate solution is an interpolant of u
    V= u;
    GVx= ux;
    GVy= uy;        
    
    for T=1:nt
        node= t(1:3,T);
        x=  p(1,node(1:3));
        y=  p(2,node(1:3));   
        f(T)= (F(x(1),y(1),A1,A2,A3) +...
               F(x(2),y(2),A1,A2,A3) +...
               F(x(3),y(3),A1,A2,A3))/3;    
    end       
                   
% flag= 0 means "no verification"
    timeHCT= clock;                                          
    [DelU,DelV,intDelVPhi]= HCT_DEL_UV(p,t,LaplU,A1,A2,A3,V,GVx,GVy,0);
    disp(sprintf('TIME:   %g sec - running of HCT_DEL_UV',etime(clock,timeHCT)));    
    if ( beta2<0.01 ) beta2= 0.01; end
            timeKs= clock;                                          
            A=1;
            a=[0 0 -1 0]';
            Fks(1,1:nt)= zeros(1,nt);
            Fks(2,1:nt)= f(1:nt);
            [W]=assempde(b,p,e,t,A,a,Fks);
            KS= -W(np+1:2*np);   
            [KSx,KSy]= pdegrad(p,t,KS);
            disp(sprintf('\nTIME:   %g sec - calc. of k*',etime(clock,timeKs)));
  
            timeYs= clock;                                          
            [YS1,YS2]= DEMgetYs(p,e,t,KSx,KSy,-f,Comega,beta2);
            disp(sprintf('TIME:   %g sec - calc. of y*',etime(clock,timeYs)));
    
    figure
    prnf2D( p,e,t, 0,sprintf('Mesh %d',refinement), 'white', 'off' );

%============================================================================
    disp(sprintf('\n----- Step 3: Calculation of loads\n'));
    % calculation of loads to the energy norm of error and
    % to terms of DEM without constants beta_i
    timeLDS= clock;
    [Et,Nt,Dt,Rt1,Rt2]= HCT_INTschP2(p,t,f,DelU,DelV,KS,KSx,KSy,YS1,YS2);
    disp(sprintf('TIME:   %g sec - calc. of loads\n',etime(clock,timeLDS)));    

%============================================================================
    disp(sprintf('\n----- Step 4: Selection of the optimal beta_1 and beta_2\n'));
    Ng= sqrt(sum(Nt));
    Eg= sqrt(sum(Et));
    RelErr= Eg/Ng * 100;
    disp(sprintf('%11.10f - norm of e',Eg));
    disp(sprintf('%11.10f - norm of u',Ng));
    disp(sprintf('%4.2f - relative error\n',RelErr));

    Dg= sqrt(sum(Dt));
    R1g= sqrt(sum(Rt1));
    R2g= sqrt(sum(Rt2));

    disp(sprintf('initial values:   beta1= %6.5f  beta2= %6.5f',beta1,beta2));
    beta1= Comega*(R1g + Comega*R2g)/Dg;
    beta2= Comega* R2g/R1g;
    disp(sprintf('optimal values:   beta1= %6.5f  beta2= %6.5f\n',beta1,beta2));


    Term= [Dg, Comega*R1g, Comega^2*R2g];
    OptEst= Term(1) + Term(2) + Term(3);
    Ieff= OptEst/Eg;
    disp(sprintf('%11.10f - norm of e',Eg));
    disp(sprintf('%11.10f - DEM',OptEst));
    disp(sprintf('%6.5f - eff. index',Ieff));    
    disp(sprintf('\nterms: \n %11.10f \n %11.10f \n %11.10f',Term(1:3)));

    color1= 0;
    color2= 1;
    Eloc= sqrt(Et(:));
    Dloc= sqrt(Dt(:)); 
    EMark= color1* ones(1,nt);
    DMark= color1* ones(1,nt); 
    Emax= max(Eloc);
    Dmax= max(Dloc); 
    tolE= 0.5;
    tolD= 0.5;
    j=0;                      %  [tlist] must be a column
    for T=1:nt 
        if ( Eloc(T) > Emax*tolE  ) EMark(T)= color2; end     
        if ( Dloc(T) > Dmax*tolD  ) 
            DMark(T)= color2;
            j= j+1;
            tlist(j,1)= T;            
        end           
    end 
    
    figure    
    subplot (2,2,1);
    prnf2D( p,e,t, exp(Eloc), 'E', 'summer', 'off' );           
    subplot (2,2,2);
    prnf2D( p,e,t, exp(Dloc), 'M', 'summer', 'off' );           
    subplot (2,2,3);
    prnf2D( p,e,t, EMark, 'E', 'summer', 'off' );
    subplot (2,2,4);
    prnf2D( p,e,t, DMark, 'M', 'summer', 'off' );  
    
%============================================================================
    if ( refinement==0 )
        disp(sprintf('\n----- Step 5: Figures'));        
        figure
        subplot (2,2,1);
        prnf3D( p,t, u, 'u', 'cool');      
        subplot (2,2,2);
        prnf3D( p,t, ux, 'u_x', 'cool');    
        subplot (2,2,3);
        prnf3D( p,t, uy, 'u_y', 'cool');  
        subplot (2,2,4);
        prnf3D( p,t, Fp, 'f', 'cool');    
    
        figure    
        subplot (2,2,1);    
        prnf3D( p,t, Delu, '\Delta U', 'cool');    
        subplot (2,2,2);
        prnf3D( p,t, KS, '\kappa^* \approx \Delta U', 'cool');
        subplot (2,2,3);
        prnf3D( p,t, YS1, 'y*_1', 'cool');
        subplot (2,2,4);
        prnf3D( p,t, YS2, 'y*_2', 'cool');                    
    end        
             
    [p,e,t]= refinemesh(g,p,e,t,tlist,'longest');        
    status= input('press enter to continue or enter <c> for termination ','s');
     
    if ( isempty(status) )
        close all;
        refinement= refinement+1;
    else 
        diary off;
        close all;
        break;    
    end    
end
diary off;
close all;