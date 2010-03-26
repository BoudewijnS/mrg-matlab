vrest = -80;        %mV
%fiberD=16.0;        %um
mysalength=3;       %um  
nodelength=1.0;     %um
space_p1=0.002;     %um
space_p2=0.004;     %um
space_i=0.004;      %um
r=70;               %Ohm*cm
mycm=0.1;           %uF/cm2
mygm=0.001;         %S/cm2
c=2;                %uF/cm2

g_nap = 0.010;      % sodium channel conductivity [S/cm^2]
g_naf = 3;          %S/cm^2
g_k = 0.080;        %S/cm^2
g_kf = 0;           %S/cm^2 (deactivated)
g_l = 0.007;        % leak channel conductivity [S/cm^2] 
e_na = 50.0;        %mV
e_k = -90.0;        %mV
e_l	= -90.0;        %mV

g_p1 = 0.001;       %S/cm^2
g_p2 = 0.0001;      %S/cm^2 
g_i = g_p2;         %S/cm^2

if fiberD == 2.0
    axonD=1.6;          %um 
    nodeD=1.4;          %um
    mysaD=1.4;         %um
    flutD=1.6;         %um
    deltax=373.2;        %um
    flutlength=10;      %um
    nl=30;             %dimensionless
end

if fiberD == 5.7
    axonD=3.4;          %um 
    nodeD=1.9;          %um
    mysaD=1.9;         %um
    flutD=3.4;         %um
    deltax=500;        %um
    flutlength=35;      %um
    nl=80;             %dimensionless
end

if fiberD == 16.0
    axonD=12.7;          %um 
    nodeD=5.5;          %um
    mysaD=5.5;         %um
    flutD=12.7;         %um
    deltax=1500;        %um
    flutlength=60;      %um
    nl=150;             %dimensionless
end

if fiberD == 10.0
    axonD=6.9;          %um 
    nodeD=3.3;          %um
    mysaD=3.3;         %um
    flutD=6.9;         %um
    deltax=1250;        %um
    flutlength=46;      %um
    nl=120;             %dimensionless

end
interlength=(deltax-nodelength-(2*mysalength)-(2*flutlength))/6;


celsius = 36;
q10_1 = 2.2^((celsius-20)/10);
q10_2 = 2.9^((celsius-20)/10);
q10_3 = 3.0^((celsius-36)/10);


[r_mysa,r_pn1,c_mysa,c_mysa_m,g_mysa,g_mysa_m] = paracomp( ...
    mysaD, mysalength, space_p1, fiberD, c, r, g_p1, nl, mycm, mygm);

[r_flut,r_pn2,c_flut,c_flut_m,g_flut,g_flut_m] = paracomp( ...
    flutD, flutlength, space_p2, fiberD, c, r, g_p2, nl, mycm, mygm);

[r_inter,r_px,c_inter,c_inter_m,g_inter,g_inter_m] = paracomp( ...
    axonD, interlength, space_i, fiberD, c, r, g_i, nl, mycm, mygm);

g_nap = calcConductance(nodeD,nodelength,g_nap);
g_k = calcConductance(nodeD,nodelength,g_k);
g_l = calcConductance(nodeD,nodelength,g_l);
g_naf = calcConductance(nodeD,nodelength,g_naf);
r_node = calcResAxial(nodeD,nodelength,r);
c_node = calcCapacity(nodeD,nodelength,c);
r_pn0 = calcResPeriax(nodeD,nodelength,space_p1,r);