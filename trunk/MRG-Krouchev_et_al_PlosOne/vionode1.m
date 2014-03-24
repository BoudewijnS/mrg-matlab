function [I,dXa] = vionode1(V, Xa, ionic, QQ10);
% vionode1 = vectorized ionode 
% ionode : I_ion & gate-change derivatives at a given voltage V & gate-state
%	(standalone version, using the ionic input structure)
% inspired from  NEURON / AXNODE.mod  & trunk/mcintyre1.m
%  =========================================================

nV = length(V); zV = zeros( nV,1 ); 

% Xa( nV,4 )
% Is there a 'force' that limits the Xa() entries \in [0,1] ?

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
	v2 = V - ionic.Vtraub; % convert to Traub convention

	a = as.A ./ (exp((v2+as.B)/as.C) + 1); 
	b = bs.A ./ (exp((v2+bs.B)/bs.C) + 1);
	end

% ------------------------------------------------------------
end % function vionode