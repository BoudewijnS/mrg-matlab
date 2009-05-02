function solveOde()

    % parameters
    % ----------------------------------------------------------------
    l_n = 1.5e-4;        % node length [cm]
    D = 20e-4;           % fiber diameter with myelin sheath [cm]
    N_layers = 20e4*D;   % represents internodal myelin sheath 
                         % (derived according to velocity check)
    g_nap = 10;         % sodium channel conductivity [mS/cm^2]
    g_naf = 3000;        %
    g_k = 80;
    g_kf = 0;
    g_l = 7;           % leak channel conductivity [mS/cm^2] 
    e_na     = 50.0;
    e_k      = -90.0;
    e_l	= -90.0;

    r = 0.07;           % specific resistivity [kOhm*cm]
    c = 2;             % specific capacity [muF/cm^2]
    mycm=0.1 
	mygm=0.001
    k = 1;               % 37?C  k=3^(0.1*T-3.7)
    V_fem = 50;          % Voltage applied in FEM simulation [V]
    polarity = -1;       % -1/+1 ... neg./pos. active electrode
    % ----------------------------------------------------------------
    mysalength=3.0; 
	nodelength=1.0;
    v_init = -80; %Initial Voltage
    % Parameters from neuron file with fiberD 14 (needs to be corrected)
    fiberD=14.0;
    g=0.739 ;
    axonD=10.4 ;
    nodeD=4.7 ;
    mysaD=4.7 ;
    flutD=10.4 ;
    deltax=1400 ;
    flutlength=56; 
    nl=140;
    interlength=(deltax-nodelength-(2*mysalength)-(2*flutlength))/6;
    
    c_mysa = c*mysaD*pi*mysalength;
    c_flut = c*flutD*pi*flutlength;
    c_inter = c*axonD*pi*interlength;
    c_node = c*nodeD*pi*nodelength;
    
    g_mysa = 0.001.*paraD1./fiberD;
    g_flut = 0.0001.*paraD2./fiberD;
    g_inter = 0.0001.*axonD./fiberD;
    
    
    c_myelin = mycm/(nodelength*2);
    g_myelin = mygm/(nodelength*2);
    
    
    r_mysa = r*4*mysalength/(mysaD^2*pi);
    r_flut = r*4*flutlength/(flutD^2*pi);
    r_inter = r*4*interlength/(axonD^2*pi);
    r_node = r*4*nodelength/(nodeD^2*pi);
    
    %how many nodes to simulate
    N_nodes = 2;
    N_inter = N_nodes-1
    
    i_node = [1,N_nodes];
    i_para_m = [i_node(2)+1,i_node(2)+N_nodes];
    i_para_h = [i_para_m(2)+1,i_para_m(2)+N_nodes];
    i_para_p = [i_para_h(2)+1,i_para_h(2)+N_nodes];
    i_para_s = [i_para_p(2)+1,i_para_p(2)+N_nodes];
    i_mysa = [i_para_s(2)+1,i_para_s(2)+N_inter,i_para_s(2)+N_inter+1,i_para_s(2)+N_nodes*2];
    i_flut = [i_mysa(4)+1,i_mysa(4)+N_inter,i_mysa(4)+N_inter+1,i_mysa(4)+N_inter*2];
    i_para_n = [i_flut(4)+1,i_flut(4)+N_inter,i_flut(4)+N_inter+1,i_flut(4)+N_inter*2];
    i_inter = [i_para_n(4)+1,i_para_n(4)+N_inter,i_para_n(4)+N_inter+1,i_para_n(4)+N_inter*2];
    i_mysa_m = [i_inter(4)+1,i_inter(4)+N_inter,i_inter(4)+N_inter+1,i_inter(4)+N_inter*2];
    i_flut_m = [i_mysa_m(4)+1,i_mysa_m(4)+N_inter,i_mysa_m(4)+N_inter+1,i_mysa_m(4)+N_inter*2];
    i_inter_m = [i_flut_m(4)+1,i_flut_m(4)+N_inter,i_flut_m(4)+N_inter+1,i_flut_m(4)+N_inter*2];
    
    %Dummy Stimulusvoltage
    Ve = zeros((N_nodes-1)*11+1,1);
    
    
    %Dummy IC
    IC = [ones(N_nodes,1)*v_init;zeros((N_nodes)*4,1);ones((N_nodes-1)*2,1) * v_init; ...
        ones((N_nodes-1)*2,1)*v_init;zeros((N_nodes-1)*2,1);ones((N_nodes-1)*6,1)*v_init;...
        zeros(10*N_inter,1)];
        
    [t,Y] = ode15s(@odeMcIntyr, [0,100], IC);

    plot(t,Y(:,1));

    function dY = odeMcIntyr(t,Y)
        dY = zeros(length(Y),1);
        
        %axonnode current
        % what to do at first and last node??
        [Iax,dY(i_para_m(1):i_para_m(2)),dY(i_para_h(1):i_para_h(2)),dY(i_para_p(1):i_para_p(2)),dY(i_para_s(1):i_para_s(2))] = axnode(Y(i_node(1):i_node(2)),Y(i_para_m(1):i_para_m(2)),Y(i_para_h(1):i_para_h(2)),Y(i_para_p(1):i_para_p(2)),Y(i_para_s(1):i_para_s(2)));
        dY(i_node(1):i_node(2)) = cableEq(Iax,Y(i_node(1):i_node(2)),Y(i_mysa(3):i_mysa(4)),Y(i_mysa(1):i_mysa(2)),V_e(i_node(1):i_node(2)),Y(i_mysa_m(3):i_mysa_m(4)),Y(i_mysa_m(1):i_mysa_m(2)),r_node,r_mysa,r_mysa,c_node);
        
        %MYSA current
        Imy = mysa(Y(i_mysa(1):i_mysa(4)));
        
        
        %FLUT with Potassium current
        %Ifl = flut(Y(i_flut(1):i_flut(4)));
        
        %FLUT without Potassium current
        [Ifl,dY(i_para_n(1):i_para_n(4))] = flutPotassium(Y(i_flut(1):i_flut(4)),Y(i_para_n(1):i_para_n(4)));
        
        %internode currents
        Iin = inter(Y(i_inter(1):i_inter(4)));
        
    end

    function dV = cableEq(I,V,V1,V2,Ve,Ve1,Ve2,Ra,Ra1,Ra2,C)
        dV = (-I + (V1-V)./(Ra1./2+Ra./2) + (V2-V)./(Ra2./2+Ra./2) + ...
              (Ve1-Ve)./(Ra1./2+Ra./2) + (Ve2-Ve)./(Ra2./2+Ra./2))./C;
    end

    function [I,dm,dh,dp,ds] = axnode(V,m,h,p,s)
        
        m_alpha = (6.57 .* (V+20.4))./(1-exp(-(V+20.4)./10.3));
        m_beta = (0.304 .* (-(V+25.7)))./(1-exp((V+25.7)./9.16));
        h_alpha = (0.34 .* (-(V+114)))./(1-exp((V+114)./11));
        h_beta = 12.6./(1+exp(-(V+31.8)./13.4));

        p_alpha = (0.0353 .* (V+27))./(1-exp(-(V+27)./10.2));
        p_beta = (0.000883 .* (-(V+34)))./(1-exp((V+34)./10));

        s_alpha = 0.3./(1+exp((V+53)./-5));
        s_beta = 0.03./(1+exp((V+90)./-1));

        %first derivatives of m, h, p, s
        dm = dpdt(m_alpha,m_beta,m);
        dh = dpdt(h_alpha,h_beta,h);
        dp = dpdt(p_alpha,p_beta,p);
        ds = dpdt(s_alpha,s_beta,s);

        I_Naf = g_naf.*m.^3.*h.*(V-e_na);       %Fast Sodium current
        I_Nap = g_nap.*p.^3.*(V-e_na);          %Persistent Sodium cureent
        I_Ks = g_k.*s.*(V-e_k);                 %Slow Potassium current
        I_Lk = passiveCurrent(V,g_l,e_l);       %Leakage current
        I = I_Naf + I_Nap + I_Ks + I_Lk;        %Sum of all nodal currents
    end

    function [I,dn] = flutPotassium(V,n)
        n_alpha = (0.0462 .* (V+83.2))./(1-exp(-(V+83.2)./1.1));
        n_beta = (0.0824 .* (-(V+66))) ./ (1-exp((V+66)./10.5));
        
        dn = dpdt(n_alpha,n_beta,n);
        
        I_Kf = g_kf.*n.^4.*(V-e_k);
        
        I = I_Kf; 
    end


    function I = flut(V)
        g=g_flut;		
		e=v_init;
        I = passiveCurrent(V,g,e);
    end

    function I = mysa(V)
        g=g_mysa;		
		e=v_init;
        I = passiveCurrent(V,g,e);
    end

    function I = inter(V)
        g=g_inter;
		e=v_init;
        I = passiveCurrent(V,g,e);
    end

    function I = passiveCurrent(V,g,e)
        I = g.*(V-e);
    end


    function dp = dpdt(alpha, beta, para)
        dp = alpha.*(1-para)-beta.*para;
    end

end