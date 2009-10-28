function Y=test(IC)

    vrest = -80;
    fiberD=10.0;        %um
    paralength1=3;      %um  
    nodelength=1.0;
    space_p1=0.002;  
    space_p2=0.004;
    space_i=0.004;
    rhoa=0.7e6;         %Ohm*um
    mycm=0.1;           %uF/cm2
    mygm=0.001;         %S/cm2
    c=2;                %uF/cm2
    
    axnode =10;
       
    g_nap = 0.010;         % sodium channel conductivity [S/cm^2]
    g_naf = 3;        %
    g_k = 0.080;
    g_kf = 0;
    g_l = 0.007;           % leak channel conductivity [S/cm^2] 
    e_na = 50.0;
    e_k = -90.0;
    e_l	= -90.0;

    g=0.690; 
    axonD=6.9; 
    nodeD=3.3; 
    paraD1=3.3; 
    paraD2=6.9; 
    deltax=1150; 
    paralength2=46; 
    nl=120;
    
    mysalength = paralength1;
    flutlength = paralength2;
    
    Rpn0=(rhoa*.01)/(pi*((((nodeD/2)+space_p1)^2)-((nodeD/2)^2)));
    Rpn1=(rhoa*.01)/(pi*((((paraD1/2)+space_p1)^2)-((paraD1/2)^2)));
    Rpn2=(rhoa*.01)/(pi*((((paraD2/2)+space_p2)^2)-((paraD2/2)^2)));
    Rpx=(rhoa*.01)/(pi*((((axonD/2)+space_i)^2)-((axonD/2)^2)));
    interlength=(deltax-nodelength-(2*paralength1)-(2*paralength2))/6;

    %periaxonal resistivity in GOhm
    r_pn0=Rpn0*(nodelength/10000)/1e3;
    r_pn1=Rpn1*(paralength1/10000)/1e3;
    r_pn2=Rpn2*(paralength2/10000)/1e3;
    r_px=Rpx*(interlength/10000)/1e3;

    %surface area in cm2
    sa_node = (pi*nodeD*nodelength)/1e8;
    sa_mysa = (pi*fiberD*paralength1)/1e8;
    sa_flut = (pi*fiberD*paralength2)/1e8;
    sa_inter = (pi*fiberD*interlength)/1e8;

    %resistiviy in Ohm*cm
    r_node = rhoa/10000;
    r_mysa = rhoa*(1/(paraD1/fiberD)^2)/10000;
    r_flut = rhoa*(1/(paraD2/fiberD)^2)/10000;
    r_inter = rhoa*(1/(axonD/fiberD)^2)/10000;

    %resistiviy in GOhm
    r_node = (4*(nodelength/10000)*r_node)/(pi*(nodeD/10000)^2)/1e9;
    r_mysa = (4*(paralength1/10000)*r_mysa)/(pi*(fiberD/10000)^2)/1e9;
    r_flut = (4*(paralength2/10000)*r_flut)/(pi*(fiberD/10000)^2)/1e9;
    r_inter = (4*(interlength/10000)*r_inter)/(pi*(fiberD/10000)^2)/1e9;

    %capacity in uF/cm2
    c_node = c;
    c_mysa = c*paraD1/fiberD;
    c_flut = c*paraD2/fiberD;
    c_inter = c*axonD/fiberD;

    %capacity in pF
    c_node = c_node*sa_node*1e6;
    c_mysa = c_mysa*sa_mysa*1e6;
    c_flut = c_flut*sa_flut*1e6;
    c_inter = c_inter*sa_inter*1e6;

    %passive conductance in S/cm^2
    g_mysa = 0.001*paraD1/fiberD;
    g_flut = 0.0001*paraD2/fiberD;
    g_inter = 0.0001*axonD/fiberD;
    %in nS
    g_mysa = g_mysa*sa_mysa*1e9;
    g_flut = g_flut*sa_flut*1e9;
    g_inter = g_inter*sa_inter*1e9;

    %membrane conductance in S/cm2
    xg = 0.001/(nl*2); 
    %in nS
    g_mysa_m = xg*sa_mysa*1e9;
    g_flut_m = xg*sa_flut*1e9;
    g_inter_m = xg*sa_inter*1e9;

    %membrane capacity in uF/cm2
    xc = 0.1/(nl*2);
    %in pF
    c_mysa_m = xc*sa_mysa*1e6;
    c_flut_m = xc*sa_flut*1e6;
    c_inter_m = xc*sa_inter*1e6;

    %nodal conductances
    g_nap = g_nap*sa_node*1e9;
    g_naf = g_naf*sa_node*1e9;
    g_k = g_k*sa_node*1e9;
    g_l = g_l*sa_node*1e9;
    
    if exist('IC','var') == 0
        IC = [-80,0,0,0,0];
    end
    [t,V] = ode15s(@ode, [0,150], IC);
    %[t1,V1] = CN2(@ode, [0,100] , IC, 2*1e-3);
    figure(1);
   
    plot(t,V(:,1));
    figure(2);
    plot(t,V(:,2:5));
    Y=V;
    
    
    
    function dV = ode(t,V)
        if  t<1
            V(1) = V(1) -20;
        end
       [Iax,dpara_m,dpara_h,dpara_p,dpara_s] = axNode(V(1),V(2),...
                                                V(3),V(4),V(5));
       dv = -Iax./c_node;
       dV = [dv;dpara_m;dpara_h;dpara_p;dpara_s];
        
    end

    function [t,Y] = CN(fun,dur,IC,step)
        t = dur(1):step:dur(2);
        t = t';
        n = length(t);
        Y = zeros(n,length(IC));
        Y(1,:) = IC;
        
        Yt = zeros(1,length(IC));
        
        for i=1:n-1
            Yt = Y(i,:) + step.*fun(t(i),Y(i,:))';
            Y(i+1,:) = Y(i,:)+(step/2).*(fun(t(i),Y(i,:))'+fun(t(i+1),Yt)');
        end
        
    end

    function [t,Y] = CN2(fun,dur,IC,step)
        t = dur(1):step:dur(2);
        t = t';
        n = length(t);
        Y = zeros(n,length(IC));
        Y(1,:) = IC;
        
        Yt = zeros(1,length(IC));
        
        for i=1:n-1
            Yt = Y(i,:) + step.*fun(t(i),Y(i,:)')';
            Y(i+1,:) = Y(i,:)+(step/2).*(fun(t(i),Y(i,:)')'+fun(t(i+1),Yt')');
        end
        
    end
    

    function [I,dm,dh,dp,ds] = axNode(V,m,h,p,s)
        
        m_alpha = (6.57 .* (V+21.4))./(1-exp(-(V+21.4)./10.3));
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

function [I,dm,dh,dp,ds] = axNode2(V,m,h,p,s,celsius)
        
        q10_1 = 2.2^((celsius-20/10));
        q10_2 = 2.9^((celsius-20)/10);
        q10_3 = 3.0^((celsius-36)/10);
        m_alpha = (1.86 .* (V+21.4))./(1-exp(-(V+21.4)./10.3)) * q10_1;
        m_beta = (0.086 .* (-(V+25.7)))./(1-exp((V+25.7)./9.16)) * q10_1;
        h_alpha = (0.062 .* (-(V+114)))./(1-exp((V+114)./11)) * q10_2;
        h_beta = 2.3./(1+exp(-(V+31.8)./13.4)) * q10_2;

        p_alpha = (0.01 .* (V+27))./(1-exp(-(V+27)./10.2)) * q10_1;
        p_beta = (0.00025 .* (-(V+34)))./(1-exp((V+34)./10)) * q10_1;

        s_alpha = 0.3./(1+exp((V-27)./-5))*q10_3;
        s_beta = 0.03./(1+exp((V+10)./-1))*q10_3;

        %first derivatives of m, h, p, s
        dm = dpdt(m_alpha,m_beta,m);
        dh = dpdt(h_alpha,h_beta,h);
        dp = dpdt(p_alpha,p_beta,p);
        ds = dpdt(s_alpha,s_beta,s);

        I_Naf = g_naf.*m.^3.*h.*(V-e_na);       %Fast Sodium current
        I_Nap = g_nap.*p.^3.*(V-e_na);          %Persistent Sodium current
        I_Ks = g_k.*s.*(V-e_k);                 %Slow Potassium current
        I_Lk = g_l*(V-e_l);                     %Leakage current
        I = I_Naf + I_Nap + I_Ks + I_Lk;        %Sum of all nodal currents
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

    function dV = cableEq(I,V,V1,V2,Ve,Ve1,Ve2,Ra,Ra1,Ra2,C)
        dV = (-I + (V1-V)./(Ra1./2+Ra./2) + (V2-V)./(Ra2./2+Ra./2) + ...
              (Ve1-Ve)./(Ra1./2+Ra./2) + (Ve2-Ve)./(Ra2./2+Ra./2))./C;
    end
    
end