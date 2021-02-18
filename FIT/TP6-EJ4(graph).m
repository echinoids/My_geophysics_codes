%Definimos los parámetros a utilizar

clear
clc
clear all

KAPPA=1e-6;	%[m^2/s] difusividad termica
RHO=4000;  	%[kg/m^3] densidad
G=10; 		%[m/s^2] gravedad
TTOP=273;	%[K] temperatura inferior
TBOTTOM=1273; 	%[K] temperatura borde superior
ALPHA=2.5e-5; 	%[1/K] coef. de expansion termica
X=1000e3; 	%[m] ancho de la caja
Y=1000e3;  	%[m] alto de la caja

ETA=[10^21,10^22,10^23,10^24,10^25]
Ra=zeros(1,length(ETA))

for i=1:length(ETA)
 
  [x_temp,y_temp,x_vel,y_vel,Temp,vx_out,vy_out] = Box2Dconvection(KAPPA,ETA(i),RHO,ALPHA,G,TTOP,TBOTTOM,X,Y);
  
  figure(i)
  
  subplot(2,2,1);
  [XX,YY]=meshgrid(x_temp,y_temp);
  vq=griddata(x_temp,y_temp,Temp,XX,YY);
  [c,h]=contourf(XX,YY,vq,50);
  set(h,'LineColor','none');
  xlim([0 1000]);
  ylim([0 1000]);
  title(sprintf('Dist. Temperatura para un ETA=%d :',ETA(i)));
  clb = colorbar;
  ylabel(clb, 'Temp [K]');
  xlabel('X [m]');
  ylabel('Y [m]');
  colormap(jet);
  axis ij;  %Invierto el eje vertical
  
  subplot(2,2,2)
  quiver(x_vel,y_vel,vx_out,vy_out,'color','k');
  xlabel('X [m]');
  ylabel('Y [m]');
  title(sprintf('Campo de velocidades para un ETA=%d :',ETA(i)));
  axis ij;
  
  N=size(Temp,2);
    
  subplot(2,2,[3,4])
  plot(Temp(:,N/2),y_temp);
  xlabel('Temp [K]');
  ylabel('Y [m]');
  title(sprintf('Perfil de Temp. para un ETA=%d :',ETA(i)));
  xlim([TTOP TBOTTOM]);
  ylim([0 1000]);
  axis ij;
  
  Ra(i)=(ALPHA*G*X^3*RHO*(TBOTTOM-TTOP))/(ETA(i)*KAPPA);

  clearvars x_temp y_temp x_vel y_vel Temp vx_out vy_out, XX, YY, vq;
end
