% Designing rectangular microstrip antenna
clear all
clc

% Input parameters
epsr=input('Specify the substrate dielectric constant : ');
h=input('Specify the substrate height in cm : ');
f=input('What is the required resonance frequency in GHz ? ');

% Calculations
w=(30/(2*f))*(2/(epsr+1))^.5;           % Patch width
epsreff=(epsr+1)/2 + ((epsr-1)/2)*(1+(12*h/w))^-.5;     ...
    % Effective Dielectric constant of the microstrip antenna
deltal=.412*h*(epsreff+.3)*((w/h)+.264)/( (epsreff-.258) * ((w/h)+.8) );
    % Extended incremental length of the patch
lamda=30/(f*sqrt(epsreff));         % Wavelength in microstrip antenna
lamda0=30/f;                        % Wavelength in air
leff=lamda/2;               % Effective length of the patch
l=leff-(2*deltal);          % Actual length of the patch


k0=2*pi/lamda0;             % Propagation constant in air
x=k0*w;
I1=@(theta) (((sin(0.5.*cos(theta).*x)/cos(theta)).^2 ).* ...
    (sin(theta)).^3);
I11=quad(I1,0,pi);          % Numeric Integration
G1=I11/(120*pi^2);          % Conductance of single slot

I2=@(theta1)((sin(x/2.*cos(theta1))./cos(theta1)).^2.*...
    besselj(0,k0.*l.*sin(theta1)).*(sin(theta1)).^3);
I12=quad(I2,0,pi);          % Numeric Integration
G12=I12./(120*pi^2);        % Mutual conductance between 2-slots

Rin=1/(2*(G1+G12));         % Resonant input resistance
y0=(l/pi)*acos(sqrt(50/Rin));   % Distance from slot#1 for matching

D0=(2*pi*w/lamda0)^2/I11;      % Directivity of single slot
g12=G12/G1;
DAF=2/(1+g12);                 % Directivity of array factor AF
D2=D0*DAF;              % Total Directivity for 2-radiating slots
D2dB=10*log10(D2);      % Total Directivity for 2-radiating slots in dB

figure(1)
phi=0:.001:2*pi;
Ephi1=sin(.5*k0*h.*cos(phi)).*cos(.5*k0*leff.*sin(phi))./...
    (.5*k0*h.*cos(phi));        % Normalized electric field in E-plane
polar(phi,Ephi1);
title('E-plane')

figure(2)
seta=0:.001:pi;
Ephi2=sin(seta).*sin(.5*k0*h.*sin(seta)).*(.5*k0*w.*cos(seta))./...
    (.25*k0^2*h*w.*sin(seta).*cos(seta));
            % Normalized electric field in H-plane
polar(seta,Ephi2);
title('H-plane')

% Displaying Results
fprintf('\nPatch width (W) = %4.4f cm \n\n',w)
fprintf('Effective length of the patch (Leff) = %4.4f cm \n\n',leff)
fprintf('Actual length of the patch (L) = %4.4f cm \n\n',l)
fprintf('Conductance of single slot (G1) = %4.4f Siemens \n\n',G1)
fprintf('Mutual conductance between 2-slots (G12) = %4.4e Siemens \n\n',G12)
fprintf('Resonant input resistance (Rin) = %4.4f Ohm \n\n',Rin)
fprintf('Distance from slot#1 for matching (y0) = %4.4f cm \n\n',y0)
fprintf('Total Directivity for 2-radiating slots (D2) = %4.4f = %4.4f dB\n\n',D2,D2dB)
