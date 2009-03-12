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
    %V_Na = 115;          % Nernst potential for sodium channels [mV]
    %V_l = -0.01;         % Nernst potential for leakage channels [mV]
    e_na     = 50.0;
    e_k      = -90.0;
    e_l	= -90.0;

    r = 0.055;           % specific resistivity [kOhm*cm]
    c = 0.6;             % specific capacity [muF/cm^2]
    g = 1;               % internodal membrane conductivity [mS/cm^2]
    k = 1;               % 37?C  k=3^(0.1*T-3.7)
    V_fem = 50;          % Voltage applied in FEM simulation [V]
    polarity = -1;       % -1/+1 ... neg./pos. active electrode
    % ----------------------------------------------------------------
    l_in = 100*D;
    C_in = l_in * (0.64*D)*pi * c / N_layers;
    C_n = l_n * (0.64*D)*pi * c;
    R_in = 4*r * l_in /((0.64*D)^2*pi); 
    R_n = 4*r * l_n /((0.64*D)^2*pi);  
    g_m = l_in * (0.64*D)*pi * g / N_layers;
    %g_naf = l_n * (0.64*D)*pi * g_naf;
    %g_nap = l_n * (0.64*D)*pi * g_nap;
    %g_k = l_n * (0.64*D)*pi * g_k;
    %g_kf = l_n * (0.64*D)*pi * g_kf;
    %g_l = l_n * (0.64*D)*pi * g_l;

    
        % dY = [ V1 m1 h1 p1 s1 n1 V2 m2 ... ]'

    V=-88.5901;    
    K=koeff(V);
    N=2;
    
    
    % IC inital conditin vector
    IC1 = [V; K(1,1)/(K(1,1)+K(1,2));...
    K(2,1)/(K(2,1)+K(2,2));...
    K(3,1)/(K(3,1)+K(3,2));...
    K(4,1)/(K(4,1)+K(4,2));...
    K(5,1)/(K(5,1)+K(5,2));];
    IC = [];
    for i = 1:N
        IC = [IC;IC1];
    end

    IC
    
    [t,Y] = ode15s(@odeMcIntyr, [0,10], IC);

    %plot(t,Y(:,2),t,Y(:,3),t,Y(:,4),t,Y(:,5))
    plot(t,Y(:,1));





    function dY = odeMcIntyr(t,Y)


        if(t>=1 && t<=2)
            V_ext=-50;
        else
            V_ext=0;
        end

        dY = zeros(6*N,1);
        for i= 1:6:6*N
            %nodal currents
            alpha = zeros(4,1);
            beta = zeros(4,1);
            I = zeros(4,1);

            alpha(1) = (6.57 * (V+20.4))/(1-exp(-(V+20.4)/10.3));
            beta(1) = (0.304 * (-(V+25.7)))/(1-exp((V+25.7)/9.16));
            alpha(2) = (0.34 * (-(V+114)))/(1-exp((V+114)/11));
            beta(2) = 12.6/(1+exp(-(V+31.8)/13.4));

            alpha(3) = (0.0353 * (V+27))/(1-exp(-(V+27)/10.2));
            beta(3) = (0.000883 * (-(V+34)))/(1-exp((V+34)/10));

            alpha(4) = 0.3/(1+exp((V+53)/-5));
            beta(4) = 0.03/(1+exp((V+90)/-1));

            dP = dpdt(alpha, beta, P);
            for j = 1:4
               Y(i+j) = dP(j); 
            end
            
            I(1) = g_naf*P(1)^3*P(2)*(V-e_na);   %INaf
            I(2) = g_nap*P(3)^3*(V-e_na);        %INap
            I(3) = g_k*P(4)*(V-e_k);             %IK
            I(4) = g_l * (V-e_l);                %ILk

            Y(i) = -1*(sum(I)+V_ext) / C_n;
            
            
        end
    end

    function dp = dpdt(alpha, beta, para)
        dp = alpha.*(1-para)-beta.*para;
    end

    function X = koeff(V)
        %compute alpha and beta values according to V
        
        X(1,1) = (6.57 * (V+20.4))/(1-exp(-(V+20.4)/10.3));
        X(1,2) = (0.304 * (-(V+25.7)))/(1-exp((V+25.7)/9.16));
        X(2,1) = (0.34 * (-(V+114)))/(1-exp((V+114)/11));
        X(2,2) = 12.6/(1+exp(-(V+31.8)/13.4));

        X(3,1) = (0.0353 * (V+27))/(1-exp(-(V+27)/10.2));
        X(3,2) = (0.000883 * (-(V+34)))/(1-exp((V+34)/10));

        X(4,1) = 0.3/(1+exp((V+53)/-5));
        X(4,2) = 0.03/(1+exp((V+90)/-1));

        X(5,1) = (0.0462 * (V + 83.2))/(1-exp(-(V+83.2)/1.1));
        X(5,2) = (0.0824 * -(V + 66))/(1-exp((V+66)/10.5));

    end

end