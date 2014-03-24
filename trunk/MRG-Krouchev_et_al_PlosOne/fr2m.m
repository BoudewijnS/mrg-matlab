% fr2m = the mirror-symmetric model & stim-targets
% 1D MRG'02 model "à la Ned"
% see also: fr1m
  
echo off; clear all; close all;
dbstop if error; dbstop if warning
randn( 'state', 0 ); rand( 'state', 0 );

format short g

% ------------------------------------------------------------

global dur t Y dtout 
global N_nodes K_vaxon i_vaxon i_nodes i_internodes i_vperi

% K_vaxon = 221 = 21 nodes + 10 * 20 intra-node
% K_vperi = 200
% K_gates = 4 * 21 = 84
% -------
% K_xtotal = 505 (for just half the axon)

%  =========================================================
M2nodes = 41 

V_r	= -80  % mV
V_THR 	= -30  % mV

% ------------------------------------------------------------
ir1Case = 1
% ------------------------------------------------------------

%  ===================================
% [ms]

% ---------------------
% tStim values, [ms]:
%	 	20, [50, 100] 200 and 400 [500] microsec; 1 and 2  [5] ms 
% ir1Case = 	5		4	3		2	1	0

switch ir1Case		% fr3see: vt0redo

case 0, tStim 	= 5, V_THR 	= -59.273 

case 1, tStim 	= 2, V_THR 	= -61.95 

case 2, tStim 	= 1, V_THR 	= -62.378  

case 3, tStim 	= 0.4, V_THR 	= -60.588

case 4, tStim 	= 0.2, V_THR 	= -57.061

case 5, tStim 	= 0.02, V_THR 	= -25.649

end

% ---------------------------

kRate 	= (V_THR - V_r)/tStim;

kStim 	= 1;
dur 	= tStim + 2; 

fgCoarse = tStim < 0.5
rp2mir( M2nodes, dur, tStim, kRate, fgCoarse );

% ------------------------------------------------------------

% N_nodes = length(i_nodes); 

%  =========================================================
if tStim < 1,	tsStr = num2str(1e3*tStim,'%03du');
		tsTag = num2str(1e3*tStim,'%d us');
else,		tsStr = num2str(tStim,'%03dm');
		tsTag = num2str(tStim,'%d ms');
end

%  =========================================================
%  = the spatial profiles (z|tStim) of the membrane voltage & currents

fr2mZ

%  =========================================================
% = draw the AP figure(s):

%  --------------------------------
fig1nk

% txt = sprintf( 'nodes = %d/%d (mirrored)', N_nodes, M2nodes );
txt = ['T_{STIM} = ' tsTag];
title( txt );

%  =========================================================
