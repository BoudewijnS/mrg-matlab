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
    g_myelin = 1;
    g_mysa = 1;
    g_flut = 0.1;
    g_inter = 0.1;
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
    
    r_mysa = r*4*mysalength/(mysaD^2*pi);
    r_flut = r*4*flutlength/(flutD^2*pi);
    r_inter = r*4*interlength/(interD^2*pi);
    r_node = r*4*nodelength/(nodeD^2*pi);
    
    %how many nodes to simulate
    N_nodes = 10;
    
    i_node = [1,N_nodes];
    i_para_m = [N_nodes+1,2*N_nodes];
    i_para_h = [2*N_nodes+1,3*N_nodes];
    i_para_p = [3*N_nodes+1,4*N_nodes];
    i_para_s = [4*N_nodes+1,5*N_nodes];
    i_mysa = [5*N_nodes+1,6*N_nodes;6*N_nodes+1,7*N_nodes];
    i_flut = [7*N_nodes+1,8*N_nodes;8*N_nodes+1,9*N_nodes];
    i_para_n = [9*N_nodes+1,10*N_nodes;10*N_nodes+1,11*N_nodes];
    i_inter = [11*N_nodes+1,12*N_nodes;10*N_nodes+1,11*N_nodes];
    
    %Dummy Stimulusvoltage
    Ve = zeros((N_nodes-1)*11+1);
    
    
    %Dummy IC
    IC = [ones(N_nodes,1)*v_init;zeros((N_nodes)*4,1);ones((N_nodes-1)*2,1) * v_init; ...
        ones((N_nodes-1)*2,1)*v_init;zeros((N_nodes-1)*2,1);ones((N_nodes-1)*6,1)*v_init];
    
    [t,Y] = ode15s(@odeMcIntyr, [0,10], IC);



    function dY = odeMcIntyr(t,Y)


        
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
        ds = dsdt(s_alpha,s_beta,s);

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
    %wrong g!!!
    function I = flut(V)
        g=0.0001.*paraD2./fiberD;		
		e=v_init;
        I = passiveCurrent(V,g,e);
    end

    function I = mysa(V)
        g=0.001.*paraD1./fiberD;		
		e=v_init;
        I = passiveCurrent(V,g,e);
    end

    function I = inter(V)
        g=0.0001.*axonD./fiberD;
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