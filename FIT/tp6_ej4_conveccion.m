%====================================================================================================
% ESTUDIO DE CONVECCIÓN TÉRMICA. Caja calentada por debajo 
% Probamos con distintas viscosidades hasta determinar desde cuál sirve teoría de capa límite
% (a partir de cuál predomina la convección por sobre la conducción, y la variación de T ocurre en
% finas capas)
%====================================================================================================
clear all;
close all;

KAPPA=1e-6;	%[m^2/s] difusividad termica
RHO=4000;  	%[kg/m^3] densidad
G=10; 		%[m/s^2] gravedad
TTOP=273;	%[K] temperatura inferior
TBOTTOM=1273; 	%[K] temperatura borde superior
ALPHA=2.5e-5; 	%[1/K] coef. de expansion termica
X=1000e3; 	%[m] ancho de la caja
Y=1000e3;  	%[m] alto de la caja

% ETA=[10^21 10^22 10^23 10^24 10^25]; %[Pa*s]=[N*s/m²] --> viscosidades a probar
ETA=10^21;
%cuando viscosidad es alta no hay deslizamiento--> velocidad casi nula
%viscosidad baja-> deslizamiento de masa--> velocidad aumenta, donde? en bordes

%------------------------------------------------------------------------------------------------------
% Utilizamos función 'Box2Dconvection' brindada por la cátedra
[x_temp,y_temp,x_vel,y_vel,Temp,vx_out,vy_out]=Box2Dconvection(KAPPA,ETA,RHO,ALPHA,G,TTOP,TBOTTOM,X,Y);
% km   ,  km  ,  km , km  , K  ,  m/s ,  m/s
%------------------------------------------------------------------------------------------------------

N=length(x_temp);			% Cantidad de elementos a representar

mm_y= 1e3.*(60*60*24*365.25);           % Para pasar de m/s a mm/years
vx=vx_out.*mm_y;
vy=vy_out.*mm_y;
vel=sqrt(vx.^2+vy.^2);		        % [mm/year] modulo de la velocidad en cada nodo

% CALCULO NUM DE RAYLEIGH
N_Rayleigh=(ALPHA*G*(TBOTTOM-TTOP)*RHO*Y**3)/(KAPPA*ETA)
  %10^25: Ra=100 / 10^24: Ra=1.000 / 10^23: Ra=10.000 / 10^22: Ra=100.000 / 10^21: Ra=1.000.000

%======================================================================================================
% Plots
%======================================================================================================
clf		   
figure(1)
				   
subplot(2,2,1)
colormap(flipud(rainbow))
pcolor(x_temp,y_temp,Temp)
shading interp;                                     % para que 'pcolor' no me plotee bordes de celdas
%contourf(x_temp,y_temp,Temp,'ShowText','off',15)   %'on' si quiero datos en curvas de nivel
colorbar('southoutside')			   
caxis([TTOP TBOTTOM])
title('Temperatura [K]')
xlabel('x [km]')
ylabel('z(x) [km]')
axis ij                                             %'y' es profundidad, no altura

subplot(2,2,2)
pcolor(x_vel,y_vel,vel)
shading interp;
caxis([min(min(vel)) max(max(vel))])
%contourf(x_vel,y_vel,vel,'ShowText','off',80)
caxis([0 51.230])                                   %vmin=vmin(eta=10^25) y vmax=vmax(eta=10^21)
colorbar('southoutside')
title('Modulo de velocidad [mm/year]')
xlabel('x [km]')
ylabel('z(x) [km]')
axis ij

subplot(2,2,3)
plot(Temp(:,N/2),y_temp,'linewidth',2)              %en el centro de la capa para evitar efectos de borde
title('Perfil de Temperatura [K]')
xlabel('Temperatura [K]')
ylabel('z(x)[km]')
axis ij

subplot(2,2,4)
quiver(x_vel,y_vel,vx_out,vy_out)
title('Campo de Velocidad')
xlabel('x [km]')
ylabel('z(x) [km]')
axis ij

print(figure(1),'ej4_21.jpg')

