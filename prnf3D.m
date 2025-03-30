 function prnf3D( p,t, f3D, ftitle, fcolormap )             

 % different format of mesh for calling trisurf
 trisurf(t(1:3,:)', p(1,:)', p(2,:)', f3D);
 title(ftitle,'FontSize',21,'FontWeight','bold');
 colormap(fcolormap);
 brighten(0.3);
 view(10,30);
 axis tight;