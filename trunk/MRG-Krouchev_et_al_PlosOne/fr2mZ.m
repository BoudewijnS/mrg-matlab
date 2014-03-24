% fr2mZ = the spatial profiles (z|tStim) of the membrane voltage & currents
% 1D MRG'02 model "à la Ned"
% see also: fr1m

%  =========================================================
global epas cpas gpas cc 
global j_gates i_gates K1gates ionic QQ10

c_node  = cc.c_node;

FS = 20;

% ---------------------
  
[dn,jtStim] = min( abs(t-tStim) );

jtStim0 = jtStim - 1;

% ---------------------

xTags = num2str( (1:N_nodes)' );  

%  =========================================================

dt = t(jtStim) - t(jtStim0);
Va = Y(jtStim,i_vaxon); Va0 = Y(jtStim0,i_vaxon);  dVa = (Va - Va0)/dt; 

Xa = zeros(N_nodes,K1gates);
Xa(j_gates) = Y(jtStim,i_gates);

% ------------------------------------------------------------
% 	the cgpars fields (e.g. .gpas) are arrays of the same size as V (and Vp)

Ic = zeros(K_vaxon,1);  Imem = zeros(K_vaxon,1);

% ----------------------------------

Imem(i_internodes) = gpas .* (Va(i_internodes)' - epas);
Ic(i_internodes) = cpas .* dVa(i_internodes)';

Ic(i_nodes) = c_node * dVa(i_nodes)';
[Imem(i_nodes),dXa] = vionode1(Va(i_nodes)', Xa, ionic, QQ10); 

Isum = Imem + Ic;

%  =========================================================
% Show the spatial profile of the membrane currents:

figNo = 4;
figure(figNo); clf; hold on; zoom on; 

yLim = [ min(Isum), max(Isum) ];

set(gca, 'FontSize', FS, 'Color', 'none' ); 
set(gca, 'yLim', yLim, 'xTick', i_nodes, 'XTickLabel', xTags );

txt = sprintf( 'Membrane currents at the end of stimulation' )
% title( txt );

hc(3) = plot( i_vaxon, Ic, 'c', 'lineWidth', 4 );
hc(2) = plot( i_vaxon, Imem, 'r', 'lineWidth', 3 );
hc(1) = plot( i_vaxon, Isum, 'k', 'lineWidth', 2 );
xlabel( 'z [axonal compartments]' ); ylabel( '[pA]' ); 

tags = { ['I_\Sigma(z|T_S=' tsTag ')'], 'I_{M}(z|T_S)', 'I_{C}(z|T_S)' };

hl = legend( hc, tags, 4 ); set(hl, 'FontSize', 2*FS)

%  =========================================================
% Show the Ea(z|tStim) spatial profile:

Ea = Va;
Ea(i_internodes) = Ea(i_internodes) + Y(jtStim,i_vperi);

% ----------------------------------

figNo = 2;
figure(figNo); clf; hold on; zoom on; 

set(gca, 'FontSize', FS, 'Color', 'none' ); 
set(gca, 'xTick', i_nodes, 'XTickLabel', xTags );

% txt = sprintf( 'Membrane voltage & Intracellular potential at the end of stimulation' )

txt = [ 'Linear V(t) growth profile, t < T_{STIM} = ' tsTag ]; 
title( txt );
% ---------------------------

hv(2) = plot( i_vaxon, Ea, 'k', 'lineWidth', 2 );
hv(1) = plot( i_vaxon, Va, 'r', 'lineWidth', 2 );
xlabel( 'z [axonal compartments]' ); ylabel( '[mV]' ); 

tags = { 'V(T_{STIM},Z)', 'E_a(T_{STIM},Z)' };
hl = legend( hv, tags ); set(hl, 'FontSize', 2*FS)

%  =========================================================
