function [topo, cc, rLIR, ionic, q10, node] = nk1mpars(M2nodes);
% nk1mpars = all parameters : mirror-symmetric stim & model 
% 1D MRG'02 model "à la Ned" 
% *** see also: nk1pars

% *********************************************************
% ..geom = topology and geometric parameters 
% *********************************************************

% =========================================================
% topological parameters:

% ---------------------------------------------------------
%# total count of Ranvier nodes (RN's)

% ---------------------------------------------------------

topo.m2nodes=M2nodes;
N_nodes = 1 + (M2nodes-1)/2;
topo.nodes=N_nodes;

% =========================================================
% internodal "module"s:

N_modules = N_nodes-1;
%# total count of "module"s
topo.modules=N_modules;

% *** Table Geom#1
% ---------------------------
% compartment: 			type: 

% regular (active) RN's		 0 (*)	
% initial mirror-symmetric RN	 0 (*)
% passive terminal RN		 4 (*)

% MYSA (axonal+periaxonal)	1
% FLUT (")			2
% STIN (")			3
% ---------------------------

ctActvRN = 0; 
% (*)
ctMirRN  = 0;
ctPasRN  = 4; 
 
%	MYSA[pn1]+FLUT[pn2] + 6x STIN[axoninter] + FLUT[pn2]+MYSA[pn1]
% *** see:  texs/nk1m1.PNG

% also EACH "module" is mirror-symettric, hence the half-counts below: 
M_MYSA=1; M_FLUT=1; M_STIN=3;
N_MYSA=2*M_MYSA; N_FLUT=2*M_FLUT; N_STIN=2*M_STIN;

%# total count of compartments in a "module":
K1module = N_MYSA + N_FLUT + N_STIN;
topo.K1module=K1module;

% The index & compartment types of a single "module":
i1module = 1:K1module;
 
ct1module = [ 1 2 3+zeros(1,N_STIN) 2 1 ];

% ---------------------------------------------------------
% Global model sizing:

%  ============================
% count of axonal Potentials & voltage-states:
% == total compartments count

K_vaxon = N_nodes + N_modules*K1module;
topo.vas= K_vaxon;

% ---------------------------
% count of periaxonal Potentials & voltage-states:

K_vperi = N_modules*K1module;
topo.vps=K_vperi;

%  ============================
%# K1node = Ionic gate-states per active RN
K1gates = 4; topo.K1gates = K1gates;

%# total count of gate-states
K_gates = N_nodes*K1gates;
topo.gates=K_gates;

%  ============================
% total count of dynamic model states:

K_xtotal = K_vaxon + K_gates + K_vperi;
topo.xtotal=K_xtotal;

% ---------------------------------------------------------
% assemble the double cable:

% ---------------------------
% list of compartment types:
ctypes = zeros(K_vaxon,1);

%# Start the list with the (mirror-symmetric) RN
% *** see: Table Geom#1 above

K_compart = 1;
ctypes(1) = ctMirRN; 

for jModule = 1:N_modules
% ---------------------------
ctypes(K_compart + i1module) = ct1module; K_compart = K_compart + K1module;

% ---------------------------
if jModule < N_modules
% add an active RN:
K_compart = K_compart + 1;
end %if jModule < N_modules

% ---------------------------
end %jModule 

% ---------------------------
%# End the list with a RN:
K_compart = K_compart + 1; % == K_vaxon

ctypes(K_vaxon) = ctActvRN;
% the terminal RN is passive according to (Foutz and McIntyre, 2010; J. Neural Eng. 7, 066008):
% ctypes(K_vaxon) = ctPasRN; 

% ---------------------------------------------------------

topo.ctypes=ctypes;

% ---------------------------------------------------------
% end of ..geom

% *********************************************************
% ..axial = more electrical parameters (1D)
% ********************************************************* 

% =========================================================
% electrical constants:

% Membrane:
    c=2;              %uF/cm2

% Axoplasm:
    r=70;               %Ohm*cm

% Myelin-related:
    mycm=0.1;           %uF/cm2
    mygm=0.001;         %S/cm2

% passive conductances, S/cm^2  (NEURON tags) 
    g_p1 = 0.001;     	% MYSA (ga)  
    g_p2 = 0.0001;      % FLUT (gf) 
    g_i = g_p2;         % STIN (gi)


% =========================================================
% geometric parameters:

geomKase = 1;

switch geomKase

% ---------------------
case 1,
fiberD = 16.0;
axonD=12.7; node.D=5.5;  flutlength=60;  
deltax=1500; 		% total length of 1 passive-"module" 
nl=150;      		% Number of myelin lamellae:

% ---------------------
case 2, 		%um == trunk/test.m
fiberD = 10.0; 
axonD=6.9; node.D = 3.3; flutlength=46;
deltax=1250;  nl=120;  

% ---------------------
end % switch geomKase

node.nl 	= nl;

% ---------------------------------------------------------
% invariant sizes: 

node.length 	= 1;   %um 
mysalength 	= 3; 

% =========================================================
% more invariant sizes, in [um]
% see: NEURON/MRG3d.000: TABLE 1. Model geometric parameters  

% periaxonal space widths: 

space_p1	= 0.002;  	% @ MYSA 
space_p2	= 0.004; 	% @ FLUT 
space_i		= 0.004;   	% @ STIN

% ---------------------------------------------------------
% derived sizes:

mysaD = node.D;  
flutD = axonD;

interlength=(deltax - node.length-(2*mysalength)-(2*flutlength))/6; 

% =========================================================
%node surface area: [um^2] -> cm2

node.area_cm2 = (pi*node.D *node.length )/1e8;

% -------------------------------------------------------
% >>> nodal electrical quantities: 

%# Overall Axial Resistance, in [GOhm]:
r_node = calcResAxial(node.D, node.length, r);
cc.r_node = r_node; 

% ------------
%# Overall Peri-Axial Resistance, in [GOhm]:
cc.r_pn0 = calcResPeriax(node.D,node.length,space_p1,r);

% -------------------------
 
%# Overall Membrane Capacity, in [pF]:
cc.c_node = calcCapacity(node.D, node.length, c);

% -------------------------------------------------------
% >>> perinodal electrical quantities: 

% -------------------------
% see ..geom : MYSA is c.type #1

cc.p(1) = pcomp(mysaD, mysalength, space_p1, fiberD, c, r, g_p1, nl, mycm, mygm);

r_mysa = cc.p(1).ra; 

rLIR = [r_mysa,r_node,r_mysa];
%# see the nodeEq() function: 
 
% -------------------------
% FLUT is c.type #2
cc.p(2) = pcomp(flutD, flutlength, space_p2, fiberD, c, r, g_p2, nl, mycm, mygm);

% -------------------------
% STIN is c.type #3   
cc.p(3) = pcomp(axonD, interlength, space_i, fiberD, c, r, g_i, nl, mycm, mygm);

% -------------------------------------------------------
% end of ..axial

% =========================================================
% ..ionic = Ranvier-nodes' electrical parameters (0D) 

% ------------------------------------------------------------
% *** almost verbatim from 	NEURON / AXNODE.mod
% related mcintyre1 functions: axnode2(V,m,h,p,s), dpdt(m_alpha,m_beta,m)

% =========================================================

	celsius	= 36	; % [degC]  

ionic.celsius	= celsius;

% ----------------------------

	q10.f1 = 2.2 ^ ((celsius-20)/ 10 );
	q10.f2 = 2.9 ^ ((celsius-20)/ 10 );
	q10.f3 = 3.0 ^ ((celsius-36)/ 10 );

% ------------------------------------------------------------

	gnapbar = 0.01	; % [S/cm^2]  = [mho/cm2]
	gnabar	= 3.0	; % [mho/cm2] 
	gkbar   = 0.08 	; % [mho/cm2] 
	gl	= 0.007 ; % [mho/cm2] 


% -------------------------
% see ..geom : a passive R.N. is c.type #4
cc.p(4) = pcomp(node.D, node.length, space_p1, fiberD, c, r, gl, NaN, NaN, NaN);

% -------------------------------------------------------

topo.lengths = [ mysalength flutlength interlength node.length ];
% ra ~ [GOhm]
% gpas ~ [nS]  
topo.lambdas = topo.lengths ./ sqrt( [cc.p(1:4).ra] .* [cc.p(1:4).gpas] );

topo.taus = [cc.p(1:4).c] ./ [cc.p(1:4).gpas];

% topo, keyboard


% =========================================================

ionic.gnapbar 	= calcConductance(node.D, node.length, gnapbar);
ionic.gnabar	= calcConductance(node.D, node.length, gnabar);
ionic.gkbar   	= calcConductance(node.D, node.length, gkbar);
ionic.gl	= calcConductance(node.D, node.length, gl);

% ----------------------------

	ionic.ena     	= 50.0  ; % [mV] 
	ionic.ek      	= -90.0 ; 
	ionic.el	= -90.0 ; 
	
	ionic.Vtraub	= -80;  ; % [mV]   == Vrest
	 
	ionic.ap.A = 0.01; 
	ionic.ap.B = 27; 
	ionic.ap.C = 10.2; 
	ionic.bp.A = 0.00025; 
	ionic.bp.B = 34; 
	ionic.bp.C = 10; 

	ionic.am.A = 1.86; 
	ionic.am.B = 21.4; 
	ionic.am.C = 10.3; 
	ionic.bm.A = 0.086; 
	ionic.bm.B = 25.7; 
	ionic.bm.C = 9.16; 

	ionic.ah.A = 0.062; 
	ionic.ah.B = 114.0; 
	ionic.ah.C = 11.0; 
	ionic.bh.A = 2.3; 
	ionic.bh.B = 31.8; 
	ionic.bh.C = 13.4; 

	ionic.as.A = 0.3; 
	ionic.as.B = -27; 
	ionic.as.C = -5; 
	ionic.bs.A = 0.03; 
	ionic.bs.B = 10; 
	ionic.bs.C = -1; 

% -------------------------------------------------------
% end of    ..ionic

% *********************************************************
end % function nk1mpars

% =========================================================
% Helper functions:
% *** N.B. these do NOT need to be and hence are NOT nested

% -------------------------------------------------------

function cp = pcomp(diameter,length,space, fiberDia, c, r, g, nl, xc, xg);

%args must be in um, ohm*cm and uF respectively

       cp.ra = calcResAxial(diameter,length,r);

       cp.rp = calcResPeriax(diameter,length,space,r);

% ---------------------------

       cp.c = calcCapacity(diameter, length, c);

       cp.cmy = calcMyelinCap(fiberDia, length, xc, nl);

% ---------------------------

       cp.gpas = calcConductance(diameter, length, g);

       cp.gmy = calcMyelinCond(fiberDia, length, xg, nl);

    end

% -------------------------------------------------------
% Capacity-related:
% -------------------------------------------------------

    function c = calcMyelinCap(fiberDiameter, length, xc, nl)
        xc = xc/(nl*2);
        c = calcCapacity(fiberDiameter, length, xc);
        %c is in uF
    end

% -------------------------------------------------------

function c = calcCapacity(diameter, length, cap)
%cap in uF/cm^2 -> pF/um^2 
  cap = cap*1e-2;   
  
  c = cap*diameter*length*pi;
%Surface in um^2 -> c  in pF
end

% -------------------------------------------------------
% Resist-related:
% -------------------------------------------------------
    function rax = calcResAxial(diameter, length, r)
% sizes in [um]	% r in [Ohm*cm]
% see: calcRes()

       a=pi*diameter^2;
       rax = calcRes(a,length,r);
    end

% -------------------------------------------------------

    function rpx = calcResPeriax(diameter, length, space, r)
% sizes in [um]	% r in [Ohm*cm]
        a = pi*((diameter+(2*space))^2-diameter^2);
        rpx = calcRes(a,length,r);
    end

% -------------------------------------------------------

    function res = calcRes(area, length, r)

        ra = r*1e-5;  
% Ohm*cm  ->  GOhm [1e-9] * um [1e4]

        res = (4*(length)*ra)/area;
    end

% -------------------------------------------------------
% Conductance-related:
% -------------------------------------------------------

    function gmy = calcMyelinCond(fiberDiameter, length, xg, nl)
       xg = xg/(nl*2);
       gmy = calcConductance(fiberDiameter, length, xg);
    end

% -------------------------------------------------------
% calcConductance() is duplicated from nk1ionic.m

function g = calcConductance(diameter, length, gi)
%gi in S/cm^2  -> nS/um^2 
        gi = gi*1e1; 

%Surface in um^2 -> g in nS
        g = gi*diameter*length*pi;
end
