function rp2mir( M2nodes, dur, tStim, kRate, fgCoarse );
% rp2mir() = the mirror-symmetric MRG'02 model & stim-targets 
% see also: rp1mir()
% 

global fgInit dtout 
global t Y topo cc rLIR ionic q10 QQ10 node

global N_nodes K_vaxon K_xtotal K1module i_vaxon i_nodes i_internodes i_vperi 	

global epas ctypes cpas cmy gpas gmy rLKRa rLKRp

global zV K_vperi  		K1gates K_gates	 j_gates i_gates	

% ----------------------------------

persistent c_node Vrest Vtraub 

persistent i_L  i_R 			 
persistent Va dVa Ea 	Vp Ep  		dXa Xa

persistent IC options

% ----------------------------------
% *** N.B. due to "coarse" T.R. the V(t) trajectory may "cuts-through" 
%	the region of initial decay, that follows the end of stim.

if fgCoarse < 0
% very "coarse" temporal resolution: 
dtout = 0.1*tStim; %% 0.25*tStim;

elseif fgCoarse
% "medium" temporal resolution: 
dtout = 0.05*tStim;

else
% "fine" T.R.: 
dtout = 0.01*tStim;
end

% [too] "fine" T.R.: dtout = min( 0.005*dur, 0.05*tStim ); 

% ----------------------------------

%# in order to reuse most of the initialization data many times:

if isempty(fgInit) | fgInit == 0
%  =================================
fgInit = 1;

% ---------------------------------------------------------

%  =========================================================
% p1 = process all the compute-once quantities
%  =========================================================

% compTags = { 'R.N.' 'MYSA' 'FLUT'  'STIN' };

[topo, cc, rLIR, ionic, q10, node] = nk1mpars(M2nodes);

% whos cc; 

%  ===========================

c_node  = cc.c_node;
QQ10 = diag( [q10.f1 q10.f2 q10.f1 q10.f3] );

% ---------------------------------------------------------

Vtraub  = ionic.Vtraub;	 % == -80 [mV] 
Vrest 	= Vtraub;	

epas 	= Vrest;

% =========================================================
% Model sizing:

% ---------------------------
%# total count of active RNs

N_nodes = topo.nodes;  zV = zeros( N_nodes,1 ); 

% ---------------------------------------------------------
% *** N.B. mirror-symmetric model & stim: 

% ---------------------------
% periaxonal Potentials & voltage-states:

K_vaxon = topo.vas; K_vperi=topo.vps;

i_L = (1:K_vaxon-1)'; i_vaxon = (1:K_vaxon)'; i_R = (2:K_vaxon)';

% >>> Va = Y(i_vaxon);
Va = zeros(K_vaxon,1); dVa = Va;
Ep = Va; Ea = Va;

j_vperi = (1:K_vperi)';
i_vperi = K_vaxon	+ j_vperi; 

% >>> Vp = Y(i_vperi);
Vp = zeros(K_vperi,1); dVp = Vp;

% ---------------------------
% gate-states

K1gates = topo.K1gates; K_gates = topo.gates;

j_gates = (1:K_gates)';
i_gates = K_vaxon+K_vperi+j_gates; 

% >>> Xa = Y(i_gates);
Xa = zeros(N_nodes,K1gates);
dXa = Xa;  

% ---------------------------
% Global states index

K_xtotal = topo.xtotal; K1module = topo.K1module; 

i_nodes 	= 1:(K1module+1):K_vaxon;
i_internodes 	= setdiff(i_vaxon,i_nodes); 

% i_RNs 		= i_nodes(1:end-1)

% =========================================================
% assemble the double-cable's compute-once matrix quantities 

% the fields of interest in ea. cc.p(ict), to compose the vectorized c.g.pars structure
% foiTags = { 'c' 'cmy' 'gpas' 'gmy' } 

ctypes = topo.ctypes;
inctypes = ctypes(i_internodes); incc = cc.p(inctypes);

% incc = 1x40 struct array with fields:     ra    rp    c    cmy    gpas    gmy

% size( [ incc.gpas ] ) =  1x40  always! i.e. even if  incc = 40x1 struct 

cpas 	= [ incc.c ]'; 
cmy	= [ incc.cmy  ]';

gpas 	= [ incc.gpas ]'; 
gmy	= [ incc.gmy  ]';

% ---------------------------------------------------------
% only for the internodal compartment types

rLKRa = zeros(K_vperi,3); rLKRp = rLKRa;

for j = 1:K_vperi  
K_compart = i_internodes(j);

% =========================================================

% ---------------------------
L_compart = K_compart - 1; R_compart = K_compart + 1;
jLKR = [L_compart K_compart R_compart];

% ---------------------------
% retrieve the proper electrical parameters from the cc.p(ict) structures,
% 	according to the compartments' type, see nk1axial and nk1geom

for k=1:3
% ---------------------------
ict1 = ctypes( jLKR(k) ); 

if ict1>0
ra = cc.p(ict1).ra; 	rp = cc.p(ict1).rp;
else
ra = cc.r_node; 	rp = cc.r_pn0;
end % if ict1==0

rLKRa(j,k) = ra; 
rLKRp(j,k) = rp;

% ---------------------------
end % k

% =========================================================

end % for j = j_vperi

%  =========================================================
% endof: p1 = process all the compute-once quantities
%  =========================================================

% ------------------------------------------------------- 

IC = odeICs;


%  =========================================================

    options = CVodeSetOptions('RelTol',1.e-3, 'AbsTol',1.e-4);

%  =================================
end % if isempty(fgInit)

%  --------------------------------
% topo, Vrest, fgInit, keyboard

    CVodeInit(@odeCVode,'BDF','Newton',0,IC,options);

t=[0]; Y=[IC];

% ------------------------------------------------------- 
% Call CVode for the desired t points: ODE not solved 'en block' in time.
% Hence the mesh is precise & not redundant

for tout = dtout:dtout:dur

% txt = sprintf( 'Next timestep ends at time=%.2f [ms]', tout ); disp( txt ); 

[status,t1,Y1] = CVode(tout,'Normal'); t=[t,t1]; Y=[Y,Y1];

end % i
% ------------------------------------------------------- 

    CVodeFree;
    Y=Y';

%  =========================================================

% *********************************************************
% Helper functions:
% *** N.B. these are all nested!
% *********************************************************

% -------------------------------------------------------
% odeICs = assembles the  IC's
 
function IC = odeICs;
% ---------------------------

Va = Vrest + zeros(K_vaxon,1);
Vp = zeros(K_vperi,1);

% ---------------------------
[ m0, h0, p0, s0 ] = nk1HH(Vrest,ionic,q10);

y0 = [ m0, h0, p0, s0 ];
Xa = y0(1+zV,:);

% whos Xa; keyboard

% ---------------------------
% package the IC's:
% IC = zeros(topo.xtotal,1);

IC = [Va; Vp; Xa(j_gates)];

%  ==========================
end % of function odeICs()

%  =========================================================
% the ODE-RHS functions:

% -------------------------------------------------------
% odeCVode = the ODE-RHS function, as required by CVode 
   
function [dY,flag,new_data] = odeCVode(t,Y)
dY = odeMRG02(t,Y);
flag = [0]; new_data = [];
end

% -------------------------------------------------------
% odeMRG02 = the actual model's ODE system
    
function dY = odeMRG02(t,Y)

% t, keyboard

%  =========================================================
% extract the components of Y:

% -----------------------------
% voltage diffs V(t):

Va = Y(i_vaxon);  Vp = Y(i_vperi);

% -----------------------------
% gates' states

Xa(j_gates) = Y(i_gates);  

%  =========================================================
% deduced quantities:

% ---------------------------------------------------------
% voltage diffs V(t) --> Potentials \Phi :

% ---------------------------
% Ee(.) 	== 0  ==>   
% 	Ep 	== [ Ee(node) Vp(internodes) Ee(node) ... ]  i.e. same size as Va
%
% 	Ea(i_nodes) 	== Va(i_nodes)
% 	Ea(i_internodes)== Va(i_internodes) + Ep(i_internodes) 

% *** see: texs/nk1mrg02.cdr (p1

% ---------------------------

% Ep = initially all zeros 
Ep(i_internodes)= Vp;

Ea(i_nodes) 	= Va(i_nodes); 
Ea(i_internodes)= Va(i_internodes) + Ep(i_internodes);

%  =========================================================
% vectorized double-cable RHS's:

% ---------------------------------------------------------
% Neighbor potentials:
% E1a == Ea( L k R ); E1p == ...

E1a = [ Ea([2;i_L]), Ea, [Ea(i_R); Vrest] ];  
E1p = [ [0; Ep(i_L)], Ep, [Ep(i_R); 0] ]; 

% ---------------------------------------------------------
% internodal compartments:

[dVa(i_internodes), dVp] = ...
interEq( Va(i_internodes), E1a(i_internodes,:), E1p(i_internodes,:) );

% ---------------------------------------------------------
% nodal compartments

[dVa(i_nodes), dXa ] = nodeEq( Va(i_nodes), E1a(i_nodes,:), t );

if t < tStim 
dVa(1) = kRate;
end

%# d Xa(:,i_RNs) == d [m;h;p;s];

% ---------------------------------------------------------
% distribute into the double cable RHS's
%# assemble all the derivatives into dY:

% whos dVa dVp dXa

dY = [dVa; dVp; dXa(j_gates)];

% ---------------------------------------------------------
end % function odeMRG02()        

%  =========================================================
%  =========================================================

function [dV, dVp] = interEq(V, E1a, E1p);

Ia = axialI(E1a, rLKRa);
Ip = axialI(E1p, rLKRp);

[dV, dVp] = cmyelin(V, Ia, Ip);
end

% -------------------------------------------------------

function [dV, dVp] = cmyelin(V, Iaxonal, Iperiaxonal)
% vectorized implementation: 
% 	Iaxonal, Iperiaxonal, as well as
% 	the cgpars fields (e.g. .gpas) are arrays of the same size as V (and Vp) 

Ipas 	= gpas .* (V - epas);
Imyelin = gmy  .* Vp;

% -----------------------------

Ic 	= - Ipas - Iaxonal;
Imy 	= - Imyelin - Iperiaxonal - Iaxonal;
 
% -----------------------------
       
dV 	= Ic 	./ cpas;
dVp 	= Imy 	./ cmy;

end


%  =========================================================

function [dV,dXa] = nodeEq(V, E1a, t)

% ---------------------------

[I,dXa] = vionode(V, Xa );

% ------------------------------------------------------------

Ia = axialI(E1a, rLIR);

% ------------------------------------------------------------

% Istim(i_nodes,t) ..

Isum = - I - Ia;
dV = Isum / c_node;

%  ======================

% ---------------------------
end

%  =========================================================

function I = axialI(vLKR, rLKR)
% vectorized implementation: 
% 	vLKR & rLKR  are   (nV,3) matrices,  
% 	where  nV = size(V) = size(Vp) 

   I = 2*( (vLKR(:,2)-vLKR(:,1)) ./ (rLKR(:,2)+rLKR(:,1)) + (vLKR(:,2)-vLKR(:,3)) ./ (rLKR(:,2)+rLKR(:,3)) );
end

%  =========================================================
% f1 = vectorized ionode
%  =========================================================

function [I,dXa] = vionode(V, Xa );
% ionode = I_ion & gate-change derivatives at a given voltage V & gate-state

% inspired from  NEURON / AXNODE.mod  & trunk/mcintyre1.m
%  =========================================================

% ------------------------------------------------------------
%# the alpha's & beta's:
%# --> first derivatives of m, h, p, s
% the a's & b's are columns of the same size as V(i_nodes)

[alpha_m,beta_m] = pm_ab(V, ionic.am, ionic.bm);
[alpha_h,beta_h] = h_ab(V, ionic.ah, ionic.bh);
[alpha_p,beta_p] = pm_ab(V, ionic.ap, ionic.bp);
[alpha_s,beta_s] = s_ab(V, ionic.as, ionic.bs);

dXa = dpdt( [alpha_m alpha_h alpha_p alpha_s], [beta_m beta_h beta_p beta_s] );

% ----------------------
%  *** N.B. For V ~ Vrest ==> dXa ~ 0
% Xa, dXa, disp('[alphas, betas]='); 
% disp( [alpha_m alpha_h alpha_p alpha_s, beta_m beta_h beta_p beta_s]); keyboard

% ------------------------------------------------------------
        
I_Naf = ionic.gnabar	* Xa(:,1).^3 .* Xa(:,2)	.*(V-ionic.ena);       
%Fast Sodium current

I_Nap = ionic.gnapbar	* Xa(:,3).^3		.*(V-ionic.ena);          
%Persistent Sodium current

I_Ks = ionic.gkbar	* Xa(:,4)		.*(V-ionic.ek);                 
%Slow Potassium current

I_Lk = ionic.gl		* 			(V-ionic.el);                     
%Leakage current
    
    I = I_Naf + I_Nap + I_Ks + I_Lk;        %Sum of all nodal currents

% ------------------------------------------------------------
end % function vionode

%  =========================================================
% vionode Helper functions:
%  =========================================================

% ------------------------------------------------------------
% dpdt: calculates the derivative of ALL the gate-states according to
%       alpha, beta and the previous values

% vectorized implementation:
% 	the alpha's & beta's and the gate-states 'Y'=Xa  are all   (nV,K1gates=4) matrices
% 	where  nV = size as V(i_nodes)

function dY = dpdt(alphas, betas)
   dY =  ( alphas.*(1-Xa) - betas.*Xa ) * QQ10;
end

% ------------------------------------------------------------
% the equations structure is shared b/n the  'm'  &  'p' gate-state 
% 	with different numerical coefficients 

	function [a,b] = pm_ab(V, ap, bp);
	a = ap.A*ap.C + zV;

	w = V+ap.B; k = find( abs(w) > 0 ); w = w(k);
	a(k) = ap.A * w ./ (1 - exp(-w/ap.C));
% ------------------
	b = bp.A*bp.C + zV;
	
	w = V+bp.B; k = find( abs(w) > 0 ); w = w(k);
	b(k) = - bp.A * w ./ (1 - exp(w/bp.C));

% ------------------

% *** N.B. these a's, b's have yet to be multiplied by the respective q10 !
% ap, keyboard

% ------------------
	end


% ------------------------------------------------------------

	function [a,b] = h_ab(V, ah, bh);
	a = ah.A*ah.C + zV;
	w = V+ah.B; k = find( abs(w) > 0 ); w = w(k);
	a(k) = - ah.A * w ./ (1 - exp(w/ah.C));
% ------------------
	b = bh.A ./ (1 + exp(-(V+bh.B)/bh.C));
	end

% ------------------------------------------------------------
	function [a,b] = s_ab(V, as, bs);
	v2 = V - Vtraub; % convert to Traub convention

	a = as.A ./ (exp((v2+as.B)/as.C) + 1); 
	b = bs.A ./ (exp((v2+bs.B)/bs.C) + 1);
	end
 
%  =========================================================
% endof: f1 = vectorized ionode
%  =========================================================

%  ------------------------------------------------------- 
end % of function nk1mir()
% *********************************************************
    
