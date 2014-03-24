function [miv, hiv, piv, siv, tmv, thv, tpv, tsv ] = nk1HH(V, ionic, q10);
%nk1HH(V) = computes the gate parameters for the various ion-channels 
%  e.g. for the Fast sodium current:
%		miv==m_infty(V) and 	hiv==h_infty(V) 
%		tmv==tau_m(V) and 	thv==tau_h(V) 
%
% *** see also: mrg01HH.M

%  =========================================================

nV = length(V); zV = zeros( size(V) );

% ========================================================= 

% ------------------------------------------------------------
%# the alpha(V)'s & beta(V)'s:
%# --> gate-dynamics parameters m_inf(V), tau_m(V) ...

[a,b] = pm_ab(V, zV, ionic.am, ionic.bm);
	[miv, tmv] = dynab(a,b,q10.f1);

[a,b] = h_ab(V, zV, ionic.ah, ionic.bh);
	[hiv, thv] = dynab(a,b,q10.f2);

% ------------------

[a,b] = pm_ab(V, zV, ionic.ap, ionic.bp);
    	[piv, tpv] = dynab(a,b,q10.f1);

% ------------------

[a,b] = s_ab(V, zV, ionic.as, ionic.bs, ionic.Vtraub);
	[siv, tsv] = dynab(a,b,q10.f3);

% ------------------------------------------------------------

end % function nk1HH

% ========================

function [y_inf, tay] = dynab(alpha, beta, Q10)

tay = 1./(alpha + beta);
y_inf = alpha .* tay;
tay = tay / Q10;

end

    
%  =========================================================
% Helper functions:
% *** N.B. these are no longer nested!

% ------------------------------------------------------------
% the equations structure is shared b/n the  'm'  &  'p' gate-state 
% 	with different numerical coefficients 

	function [a,b] = pm_ab(V, zV, ap, bp);
	a = ap.A*ap.C + zV;

	w = V+ap.B; k = find( abs(w) > 0 ); w = w(k);
	a(k) = ap.A * w ./ (1 - exp(-w/ap.C));
% ------------------
	b = bp.A*bp.C + zV;
	
	w = V+bp.B; k = find( abs(w) > 0 ); w = w(k);
	b(k) = - bp.A * w ./ (1 - exp(w/bp.C));
	end


% ------------------------------------------------------------

	function [a,b] = h_ab(V, zV, ah, bh);
	a = ah.A*ah.C + zV;
	w = V+ah.B; k = find( abs(w) > 0 ); w = w(k);
	a(k) = - ah.A * w ./ (1 - exp(w/ah.C));
% ------------------
	b = bh.A ./ (1 + exp(-(V+bh.B)/bh.C));
	end

% ------------------------------------------------------------
	function [a,b] = s_ab(V, zV, as, bs, Vtraub);
	v2 = V - Vtraub; % convert to Traub convention

	a = as.A ./ (exp((v2+as.B)/as.C) + 1); 
	b = bs.A ./ (exp((v2+bs.B)/bs.C) + 1);
	end