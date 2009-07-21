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
    rhoa=0.7e6;         %Ohm um
    mycm=0.1 ;
	mygm=0.001;
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
    space_p1=0.002;
    space_p2=0.004;
    space_i=0.004;
    deltax=1400 ;
    nl = 140;
    xg = mygm/(nl*2);   %g for membrane
    xc = mycm/(nl*2);   %c for mambrane
    flutlength=56; 
    nl=140;
    interlength=(deltax-nodelength-(2*mysalength)-(2*flutlength))/6;
    
    c_mysa = c*mysaD*pi*mysalength;
    c_flut = c*flutD*pi*flutlength;
    c_inter = c*axonD*pi*interlength;
    c_node = c*nodeD*pi*nodelength;
    
    g_mysa = 0.001.*mysaD./fiberD;
    g_flut = 0.0001.*flutD./fiberD;
    g_inter = 0.0001.*axonD./fiberD;
    
    
    c_myelin = mycm/(nodelength*2);
    g_myelin = mygm/(nodelength*2);
    
    
    r_mysa = r*4*mysalength/(mysaD^2*pi);
    r_flut = r*4*flutlength/(flutD^2*pi);
    r_inter = r*4*interlength/(axonD^2*pi);
    r_node = r*4*nodelength/(nodeD^2*pi);
    
    r_pn0=(rhoa*.01)/(pi*((((nodeD/2)+space_p1)^2)-((nodeD/2)^2)));
	r_pn1=(rhoa*.01)/(pi*((((mysaD/2)+space_p1)^2)-((mysaD/2)^2)));
	r_pn2=(rhoa*.01)/(pi*((((flutD/2)+space_p2)^2)-((flutD/2)^2)));
	r_px=(rhoa*.01)/(pi*((((axonD/2)+space_i)^2)-((axonD/2)^2)));
    
    %how many nodes to simulate
    N_nodes = 2;
    N_inter = N_nodes-1
    
    % index values of the dV Vektor
    i_node = [1,N_nodes];
    i_para_m = [i_node(2)+1,i_node(2)+N_nodes];
    i_para_h = [i_para_m(2)+1,i_para_m(2)+N_nodes];
    i_para_p = [i_para_h(2)+1,i_para_h(2)+N_nodes];
    i_para_s = [i_para_p(2)+1,i_para_p(2)+N_nodes];
    i_mysa = [i_para_s(2)+1,i_para_s(2)+N_inter, ...
        i_para_s(2)+N_inter+1,i_para_s(2)+N_nodes*2];
    i_flut = [i_mysa(4)+1,i_mysa(4)+N_inter, ... 
        i_mysa(4)+N_inter+1,i_mysa(4)+N_inter*2];
    i_para_n = [i_flut(4)+1,i_flut(4)+N_inter, ...
        i_flut(4)+N_inter+1,i_flut(4)+N_inter*2];
    i_inter = [i_para_n(4)+1,i_para_n(4)+N_inter, ... 
        i_para_n(4)+N_inter+1,i_para_n(4)+N_inter*2, ...
        i_para_n(4)+(2*N_inter)+1,i_para_n(4)*(3*N_inter)...
        i_para_n(4)+(3*N_inter)+1,i_para_n(4)*(4*N_inter)...
        i_para_n(4)+(4*N_inter)+1,i_para_n(4)*(5*N_inter)...
        i_para_n(4)+(5*N_inter)+1,i_para_n(4)*(6*N_inter)];
    i_mysa_b = [i_inter(12)+1,i_inter(12)+N_inter, ...
                i_inter(12)+N_inter+1,i_inter(12)+N_inter*2];
    i_flut_b = [i_mysa_b(4)+1,i_mysa_b(4)+N_inter, ... 
        i_mysa_b(4)+N_inter+1,i_mysa_b(4)+N_inter*2];
    i_inter_b = [i_flut_b(4)+1,i_flut_b(4)+N_inter, ... 
        i_flut_b(4)+N_inter+1,i_flut_b(4)+N_inter*2, ...
        i_flut_b(4)+(2*N_inter)+1,i_flut_b(4)*(3*N_inter)...
        i_flut_b(4)+(3*N_inter)+1,i_flut_b(4)*(4*N_inter)...
        i_flut_b(4)+(4*N_inter)+1,i_flut_b(4)*(5*N_inter)...
        i_flut_b(4)+(5*N_inter)+1,i_flut_b(4)*(6*N_inter)];
    %Dummy Stimulusvoltage
    Ve = zeros((N_nodes-1)*11+1,1);
    
    
    %Dummy IC
    IC = [ones(N_nodes,1)*v_init;zeros((N_nodes)*4,1); ...
        ones((N_nodes-1)*2,1) * v_init; ...
        ones((N_nodes-1)*2,1)*v_init;zeros((N_nodes-1)*2,1); ... 
        ones((N_nodes-1)*6,1)*v_init;...
        zeros(10*N_inter,1)];
        
    [t,Y] = ode15s(@odeMcIntyr, [0,100], IC);

    plot(t,Y(:,1));
    
    % odeMcIntyr: calculates the first derivative of all parameters of the
    %             McIntyre nerve model
    function dY = odeMcIntyr(t,Y)
        dY = zeros(length(Y),1);
        node = Y([i_node(1):i_node(2)]);
        para_m = Y([i_para_m(1):i_para_m(2)]);
        para_h = Y([i_para_h(1):i_para_h(2)]);
        para_p = Y([i_para_p(1):i_para_p(2)]);
        para_s = Y([i_para_s(1):i_para_s(2)]);
        mysa_l = Y([i_mysa(1):i_mysa(2)]);
        mysa_r = Y([i_mysa(3):i_mysa(4)]);
        flut_l = Y([i_flut(1):i_flut(2)]);
        flut_r = Y([i_flut(3):i_flut(4)]);
        para_n_l = Y([i_para_n(1):i_para_n(2)]);
        para_n_r = Y([i_para_n(3):i_para_n(4)]);
        inter_1 = Y([i_inter(1):i_inter(2)]);
        inter_2 = Y([i_inter(3):i_inter(4)]);
        inter_3 = Y([i_inter(5):i_inter(6)]);
        inter_4 = Y([i_inter(7):i_inter(8)]);
        inter_5 = Y([i_inter(9):i_inter(10)]);
        inter_6 = Y([i_inter(11):i_inter(12)]);
        
        mysa_l_b = Y([i_mysa_b(1):i_mysa_b(2)]);
        mysa_r_b = Y([i_mysa_b(3):i_mysa_b(4)]);
        flut_l_b = Y([i_flut_b(1):i_flut_b(2)]);
        flut_r_b = Y([i_flut_b(3):i_flut_b(4)]);
        inter_1_b = Y([i_inter_b(1):i_inter_b(2)]);
        inter_2_b = Y([i_inter_b(3):i_inter_b(4)]);
        inter_3_b = Y([i_inter_b(5):i_inter_b(6)]);
        inter_4_b = Y([i_inter_b(7):i_inter_b(8)]);
        inter_5_b = Y([i_inter_b(9):i_inter_b(10)]);
        inter_6_b = Y([i_inter_b(11):i_inter_b(12)]);
        
        
        %axonnode current
        % what to do at first and last node??
        [Iax,dpara_m,dpara_h,dpara_p,dpara_s] = axnode(node,para_m,...
                                                para_h,para_p,para_s);
        dnode = cableEq(Iax,node,[node(1);mysa_l],...
            [mysa_r;node(N_nodes)],V_e(i_node(1):i_node(2)),...
            [V_e(1);mysa_l_m],[mysa_r_m;V_e(N_node)],...
            r_node,r_mysa,r_mysa,c_node);
        
        
        %MYSA current
        Imy_l = mysa(mysa_l,mysa_l_b);
        Imy_r = mysa(mysa_r,mysa_r_b);
        
        dmysa_l = cableEq(Imy_l,mysa_l,node(1:N_inter),flut_l, ...
                          mysa_l_b,V_e(i_node(1):i_node(2)-1), ...
                          flut_l_b,r_mysa,r_node,r_flut,c_mysa);
        dmysa_r = cableEq(Imy_r,mysa_r,node(2:N_nodes),flut_r, ...
                          mysa_r_b,V_e(i_node(1)+1:i_node(2)), ...
                          flut_r_b,r_mysa,r_node,r_flut,c_mysa);
        
        %currents between myelin and node at mysa regions
        Imy_l_b = -Imy_l + passiveCurrent(mysa_l_b,xg,...
            V_e(i_mysa(1):i_mysa(2)));
        Imy_r_b = -Imy_r + passiveCurrent(mysa_r_b,xg,...
            V_e(i_mysa(3):i_mysa(4)));
        
        dmysa_l_b = cableEq(Imy_l_b,mysa_l_b,V_e(1:N_inter),flut_l_b, ...
                            V_e(i_mysa(1):i_mysa(2)),V_e(1:N_inter), ...
                            V_e(i_flut(1):i_flut(2)),r_pn1,r_pn0,r_pn2,xc);
        dmysa_r_b = cableEq(Imy_r_b,mysa_r_b,V_e(2:N_nodes),flut_r_b, ...
                            V_e(i_mysa(3):i_mysa(4)),V_e(2:N_nodes), ...
                            V_e(i_flut(3):i_flut(4)),r_pn1,r_pn0,r_pn2,xc);
                        
        %FLUT currents
        Ifl_l = flut(flut_l,flut_l_b);
        Ifl_r = flut(flut_r,flut,r_b);
        
        %FLUT potassium currents
        [Ifl_lp, dpara_n_l] = flutPotassium(flut_l,para_n_l);
        Ifl_l = Ifl_l + Ifl_lp;
        [Ifl_rp, dpara_n_r] = flutPotassium(flut_r,para_n_r);
        Ifl_r = Ifl_r + Ifl_rp;
        
        dflut_l = cableEq(Ifl_l,flut_l,mysa_l,inter_1, ...
                          flut_l_b,mysa_l_b,inter_1_b, ...
                          r_flut,r_mysa,r_inter,c_flut);
        dflut_r = cableEq(Ifl_l,flut_r,mysa_r,inter_6, ...
                          flut_r_b,mysa_r_b,inter_6_b, ...
                          r_flut,r_mysa,r_inter,c_flut);
        
        %currents between myelin and node at flut regions
        Ifl_l_b = -Ifl_l + passiveCurrent(flut_l_b,xg, ...
                                            V_e(i_flut(1):i_flut(2)));
        Ifl_r_b = -Ifl_r + passiveCurrent(flut_r_b,xg, ...
                                            V_e(i_flut(3):i_flut(4)));
        
        dflut_l_b = cableEq(Ifl_l_b,flut_l_b,mysa_l_b,inter_1_b, ...
                    V_e(i_flut(1):i_flut(2)),V_e(i_mysa(1):i_mysa(2)), ...
                    V_e(i_inter(1):i_inter(2)),r_pn2,r_pn1,r_px,xc);
        dflut_r_b = cableEq(Ifl_r_b,flut_r_b,mysa_r_b,inter_6_b, ...
                    V_e(i_flut(3):i_flut(4)),V_e(i_mysa(3):i_mysa(4)), ...
                    V_e(i_inter(11):i_inter(12)),r_pn2,r_pn1,r_px,xc);              
        
        %internode currents
        Iin_1 = inter(inter_1,inter_1_b);
        %din_1 = cableEq;
        
        dY=[dnode;dpara_m;dpara_h;dPara_p;dpara_s;dmysa_l;dmysa_r];
        
    end

    % assemble: assembles seprate vectors into one for use with the 
    %           ode function.
    function dV = assemble(axon,para_m,para_h,para_p,para_s,mysa,flut, ...
                           para_n,inter,layer2)
                       
       dV = [axon;para_m;para_h;para_p;para_s;mysa;flut;para_n;inter; ...
             layer2];
    end
    
    % cableEq: calculates the cable-equation
    function dV = cableEq(I,V,V1,V2,Ve,Ve1,Ve2,Ra,Ra1,Ra2,C)
        dV = (-I + (V1-V)./(Ra1./2+Ra./2) + (V2-V)./(Ra2./2+Ra./2) + ...
              (Ve1-Ve)./(Ra1./2+Ra./2) + (Ve2-Ve)./(Ra2./2+Ra./2))./C;
    end
    
    % axnode: calculation of the currents of the axon
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
        I_Nap = g_nap.*p.^3.*(V-e_na);          %Persistent Sodium current
        I_Ks = g_k.*s.*(V-e_k);                 %Slow Potassium current
        I_Lk = passiveCurrent(V,g_l,e_l);       %Leakage current
        I = I_Naf + I_Nap + I_Ks + I_Lk;        %Sum of all nodal currents
    end

    % flutPotassim: flut potassium current
    function [I,dn] = flutPotassium(V,n)
        n_alpha = (0.0462 .* (V+83.2))./(1-exp(-(V+83.2)./1.1));
        n_beta = (0.0824 .* (-(V+66))) ./ (1-exp((V+66)./10.5));
        
        dn = dpdt(n_alpha,n_beta,n);
        
        I_Kf = g_kf.*n.^4.*(V-e_k);
        
        I = I_Kf; 
    end

    % flut: simple flut current
    function I = flut(V,e)
        g=g_flut;
        I = passiveCurrent(V,g,e);
    end
    
    % mysa: simple mysa current
    function I = mysa(V,e)
        g=g_mysa;
        I = passiveCurrent(V,g,e);
    end

    % inter: currents of the internode
    function I = inter(V,e)
        g=g_inter;
        I = passiveCurrent(V,g,e);
    end

    % passiveCurrent: calculation of the passive current
    function I = passiveCurrent(V,g,e)
        I = g.*(V-e);
    end

    % dpdt: calculates the derivative of the of the parameter according to
    %       alpha, beta and the previous value
    function dp = dpdt(alpha, beta, para)
        dp = alpha.*(1-para)-beta.*para;
    end

end