function [t,Y,xs]=solveOde(dur,IC)

%     % parameters
%     % ----------------------------------------------------------------
%     l_n = 1.5e-4;        % node length [cm]
%     D = 20e-4;           % fiber diameter with myelin sheath [cm]
%     g_nap = 0.010;         % sodium channel conductivity [S/cm^2]
%     g_naf = 3;        %
%     g_k = 0.080;
%     g_kf = 0;
%     g_l = 0.007;           % leak channel conductivity [S/cm^2] 
%     e_na = 50.0;
%     e_k = -90.0;
%     e_l	= -90.0;
% 
%     r = 70;           	%Ohm*cm
%     rhoa = 7e5;          %Ohm*um
%     c = 2;             % specific capacity [muF/cm^2]
%     mycm=0.1 ;
% 	mygm=.001;
%     V_fem = 50;          % Voltage applied in FEM simulation [V]
%     polarity = -1;       % -1/+1 ... neg./pos. active electrode
%     % ----------------------------------------------------------------
%     mysalength=3.0e-4; 
% 	nodelength=1.0e-4;
%     vrest = -80; %Initial Voltage
%     % Parameters from neuron file with fiberD 14 (needs to be corrected)
%     fiberD=14.0e-4;
%     
%     axonD=10.4e-4 ;
%     nodeD=4.7e-4 ;
%     mysaD=4.7e-4 ;
%     flutD=10.4e-4 ;
%     space_p1=0.002e-4;
%     space_p2=0.004e-4;
%     space_i=0.004e-4;
%     deltax=1400e-4 ;
%     nl = 140;
%     
%     flutlength=56e-4;
%     interlength=(deltax-nodelength-(2*mysalength)-(2*flutlength))/6;
%     
%     c_mysa = c*mysaD*pi*mysalength*(mysaD/fiberD);
%     c_flut = c*flutD*pi*flutlength*(flutD/fiberD);
%     c_inter = c*axonD*pi*interlength*(axonD/fiberD);
%     c_node = c*nodeD*pi*nodelength;
%     
%     c_mysa_m = mycm/nl*mysaD*pi*mysalength;
%     c_flut_m = mycm/nl*flutD*pi*flutlength;
%     c_inter_m = mycm/nl*axonD*pi*interlength;
%     
%     g_mysa_m = mygm/nl*mysaD*pi*mysalength;
%     g_flut_m = mygm/nl*flutD*pi*flutlength;
%     g_inter_m = mygm/nl*axonD*pi*interlength;
%     
% %     
% %     g_mysa = 0.001.*mysaD./fiberD;
% %     g_flut = 0.0001.*flutD./fiberD;
% %     g_inter = 0.0001.*axonD./fiberD;
%     
%     g_mysa = 0.001*mysaD*pi*mysalength*(mysaD/fiberD);
%     g_flut = 0.0001*flutD*pi*flutlength*(flutD/fiberD);
%     g_inter = 0.0001*axonD*pi*interlength*(axonD/fiberD);
%     g_nap = g_nap*nodeD*pi*nodelength;
%     g_naf = g_naf*nodeD*pi*nodelength;
%     g_k = g_k*nodeD*pi*nodelength;
%     g_l = g_l*nodeD*pi*nodelength;
%     
%     
%     r_mysa = r*4*mysalength/(mysaD^2*pi)*(1/(mysaD/fiberD)^2);
%     r_flut = r*4*flutlength/(flutD^2*pi)*(1/(flutD/fiberD)^2);
%     r_inter = r*4*interlength/(axonD^2*pi)*(1/(axonD/fiberD)^2);
%     r_node = r*4*nodelength/(nodeD^2*pi);
%     
%     r_pn0=(r*0.01*nodelength)/(pi*((((nodeD/2)+space_p1)^2)-((nodeD/2)^2)));
% 	r_pn1=(r*0.01*mysalength)/(pi*((((mysaD/2)+space_p1)^2)-((mysaD/2)^2)));
% 	r_pn2=(r*0.01*flutlength)/(pi*((((flutD/2)+space_p2)^2)-((flutD/2)^2)));
% 	r_px=(r*0.01*interlength)/(pi*((((axonD/2)+space_i)^2)-((axonD/2)^2)));
    
    vrest = -80;
    fiberD=10.0;        %um
    paralength1=3;      %um  
    nodelength=1.0;
    space_p1=0.002;  
    space_p2=0.004;
    space_i=0.004;
    rhoa=7e5;         %Ohm*um
    mycm=0.1;           %uF/cm2
    mygm=0.001;         %S/cm2
    c=2;                %uF/cm2

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
    
    
    %how many nodes to simulate
    N_nodes = 21;
    N_inter = N_nodes-1;
    
    % index values of the dV Vektor
    i_node = [1,N_nodes];
    i_para_m = [i_node(2)+1,i_node(2)+N_nodes];
    i_para_h = [i_para_m(2)+1,i_para_m(2)+N_nodes];
    i_para_p = [i_para_h(2)+1,i_para_h(2)+N_nodes];
    i_para_s = [i_para_p(2)+1,i_para_p(2)+N_nodes];
    i_mysa = [i_para_s(2)+1,i_para_s(2)+N_inter, ...
        i_para_s(2)+N_inter+1,i_para_s(2)+N_inter*2];
    i_flut = [i_mysa(4)+1,i_mysa(4)+N_inter, ... 
        i_mysa(4)+N_inter+1,i_mysa(4)+N_inter*2];
    i_para_n = [i_flut(4)+1,i_flut(4)+N_inter, ...
        i_flut(4)+N_inter+1,i_flut(4)+N_inter*2];
    i_inter = [i_para_n(4)+1,i_para_n(4)+N_inter; ... 
        i_para_n(4)+N_inter+1,i_para_n(4)+N_inter*2; ...
        i_para_n(4)+(2*N_inter)+1,i_para_n(4)+(3*N_inter);...
        i_para_n(4)+(3*N_inter)+1,i_para_n(4)+(4*N_inter);...
        i_para_n(4)+(4*N_inter)+1,i_para_n(4)+(5*N_inter);...
        i_para_n(4)+(5*N_inter)+1,i_para_n(4)+(6*N_inter)];
    i_mysa_b = [i_inter(12)+1,i_inter(12)+N_inter, ...
                i_inter(12)+N_inter+1,i_inter(12)+N_inter*2];
    i_flut_b = [i_mysa_b(4)+1,i_mysa_b(4)+N_inter, ... 
        i_mysa_b(4)+N_inter+1,i_mysa_b(4)+N_inter*2];
    i_inter_b = [i_flut_b(4)+1,i_flut_b(4)+N_inter; ... 
        i_flut_b(4)+N_inter+1,i_flut_b(4)+N_inter*2; ...
        i_flut_b(4)+(2*N_inter)+1,i_flut_b(4)+(3*N_inter);...
        i_flut_b(4)+(3*N_inter)+1,i_flut_b(4)+(4*N_inter);...
        i_flut_b(4)+(4*N_inter)+1,i_flut_b(4)+(5*N_inter);...
        i_flut_b(4)+(5*N_inter)+1,i_flut_b(4)+(6*N_inter)];
    %Dummy Stimulusvoltage
    V_e = zeros(i_inter(6,2),1);
    V_stim = zeros(i_inter(6,2),1);
    q = 500000;
    xe = 14000;
    ye = 5000;
    V_stim(1:N_nodes) = electrode(q,xe,ye,((1:N_nodes)-1)*deltax,zeros(1,N_nodes));
    V_stim(i_mysa(1):i_mysa(2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength/2,zeros(1,N_inter));
    V_stim(i_mysa(3):i_mysa(4)) = electrode(q,xe,ye,((1:N_inter))*deltax-mysalength/2,zeros(1,N_inter));
    V_stim(i_flut(1):i_flut(2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength+flutlength/2,zeros(1,N_inter));
    V_stim(i_flut(3):i_flut(4)) = electrode(q,xe,ye,((1:N_inter))*deltax-mysalength-flutlength/2,zeros(1,N_inter));
    V_stim(i_inter(1,1):i_inter(1,2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength+flutlength+interlength/2,zeros(1,N_inter));
    V_stim(i_inter(2,1):i_inter(2,2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength+flutlength+3*interlength/2,zeros(1,N_inter));
    V_stim(i_inter(3,1):i_inter(3,2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength+flutlength+5*interlength/2,zeros(1,N_inter));
    V_stim(i_inter(4,1):i_inter(4,2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength+flutlength+7*interlength/2,zeros(1,N_inter));
    V_stim(i_inter(5,1):i_inter(5,2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength+flutlength+9*interlength/2,zeros(1,N_inter));
    V_stim(i_inter(6,1):i_inter(6,2)) = electrode(q,xe,ye,((1:N_inter)-1)*deltax+mysalength+flutlength+11*interlength/2,zeros(1,N_inter));
    

    figure(2);
    plot(1:N_nodes,V_stim(1:N_nodes));
    hold off;
    if (exist('IC','var') == 0)
        IC = zeros(i_inter_b(6,2),1);
        IC(1:N_nodes) = -80;
        IC(i_mysa(1):i_flut(4)) = -80;
        IC(i_inter(1,1):i_inter(6,2)) = -80;
        stim = 0;
    end
    if (exist('dur','var')==0)
        dur = 10;
    end
    
    Istim = zeros(N_nodes,1);
    xs=[];
    [t,Y] = ode15s(@odeMcIntyr, [0,dur], IC);

    figure(1);
    for i = 1:N_nodes
        V(i,:) = Y(:,i) - 40*(i-1);
    end
    plot(t,V);
    hold off;
    figure(3);
    %plot(xs(:,4),((xs(:,1)-xs(:,2))./(r_node+r_mysa))./c_node,xs(:,4),((xs(:,1)-xs(:,3))./(r_node+r_mysa))./c_node);
    plot(t,Y(:,i_inter_b(3,1)+11));
    figure(4);
    plot(t,Y(:,11));
    % odeMcIntyr: calculates the first derivative of all parameters of the
    %             McIntyre nerve model
    function dY = odeMcIntyr(t,Y)
        
        if exist('stim','var') ==0
            if mod(t,10) < 1
                %
                %V_e = V_stim;
                %Y(11) = Y(11) -30;
                %V_e(11) = -10;
                %Istim(11) = 0.00010;
            else
                V_e = zeros(i_inter(6,2),1);
                Istim(11) = 0;
            end
        end
        dY = zeros(length(Y),1);
        inter = zeros(N_inter,6);
        inter_b = zeros(N_inter,6);
        
        node = Y(i_node(1):i_node(2));
        para_m = Y(i_para_m(1):i_para_m(2));
        para_h = Y(i_para_h(1):i_para_h(2));
        para_p = Y(i_para_p(1):i_para_p(2));
        para_s = Y(i_para_s(1):i_para_s(2));
        mysa_l = Y(i_mysa(1):i_mysa(2));
        mysa_r = Y(i_mysa(3):i_mysa(4));
        flut_l = Y(i_flut(1):i_flut(2));
        flut_r = Y(i_flut(3):i_flut(4));
        para_n_l = Y(i_para_n(1):i_para_n(2));
        para_n_r = Y(i_para_n(3):i_para_n(4));
        inter(:,1) = Y(i_inter(1,1):i_inter(1,2));
        inter(:,2) = Y(i_inter(2,1):i_inter(2,2));
        inter(:,3) = Y(i_inter(3,1):i_inter(3,2));
        inter(:,4) = Y(i_inter(4,1):i_inter(4,2));
        inter(:,5) = Y(i_inter(5,1):i_inter(5,2));
        inter(:,6) = Y(i_inter(6,1):i_inter(6,2));
       

        mysa_l_b = Y(i_mysa_b(1):i_mysa_b(2));
        mysa_r_b = Y(i_mysa_b(3):i_mysa_b(4));
        flut_l_b = Y(i_flut_b(1):i_flut_b(2));
        flut_r_b = Y(i_flut_b(3):i_flut_b(4));
        inter_b(:,1) = Y(i_inter_b(1,1):i_inter_b(1,2));
        inter_b(:,2) = Y(i_inter_b(2,1):i_inter_b(2,2));
        inter_b(:,3) = Y(i_inter_b(3,1):i_inter_b(3,2));
        inter_b(:,4) = Y(i_inter_b(4,1):i_inter_b(4,2));
        inter_b(:,5) = Y(i_inter_b(5,1):i_inter_b(5,2));
        inter_b(:,6) = Y(i_inter_b(6,1):i_inter_b(6,2)); 
        
        mysa_lv = mysa_l + mysa_l_b;
        mysa_rv = mysa_r + mysa_r_b;
        flut_lv = flut_l + flut_l_b;
        flut_rv = flut_r + flut_r_b;
        interv = inter+inter_b;
        
        xs=[xs;V_e(11),mysa_l_b(11),V_e(i_mysa(1)+11),t];
        
        %axonnode current
        % what to do at first and last node??
        [Iax,dpara_m,dpara_h,dpara_p,dpara_s] = axnode(node,para_m,...
                                                para_h,para_p,para_s);
        dnode = cableEq(Iax-Istim,node,[node(1);mysa_lv],...
            [mysa_rv;node(N_nodes)],V_e(i_node(1):i_node(2)),...
            [V_e(1);V_e(i_mysa(1):i_mysa(2))],[V_e(i_mysa(3):i_mysa(4));V_e(N_nodes)],...
            r_node,r_mysa,r_mysa,c_node);
        
        
        %MYSA current
        Imy_l = mysa(mysa_l);
        Imy_r = mysa(mysa_r);
        
        dmysa_l = cableEq(Imy_l,mysa_lv,node(1:N_inter),flut_lv, ...
                          V_e(i_mysa(1):i_mysa(2)),V_e(i_node(1):i_node(2)-1), ...
                          V_e(i_flut(1):i_flut(2)),r_mysa,r_node,r_flut,c_mysa);
        dmysa_r = cableEq(Imy_r,mysa_rv,node(2:N_nodes),flut_rv, ...
                          V_e(i_mysa(3):i_mysa(4)),V_e(i_node(1)+1:i_node(2)), ...
                          V_e(i_flut(3):i_flut(4)),r_mysa,r_node,r_flut,c_mysa);
        
        %currents between myelin and node at mysa regions
        Imy_l_b = passiveMyelin(mysa_l_b,g_mysa_m,...
            V_e(i_mysa(1):i_mysa(2)));
        Imy_r_b = passiveMyelin(mysa_r_b,g_mysa_m,...
            V_e(i_mysa(3):i_mysa(4)));
        
        dmysa_l_b = extracellular(Imy_l_b,mysa_l_b,V_e(1:N_inter),flut_l_b, ...
                            r_pn1,r_pn0,r_pn2,c_mysa_m,Imy_l,c_mysa);
        
        dmysa_r_b = extracellular(Imy_r_b,mysa_r_b,V_e(2:N_nodes),flut_r_b, ...
                            r_pn1,r_pn0,r_pn2,c_mysa_m,Imy_r,c_mysa);
                        
        %FLUT currents
        Ifl_l = flut(flut_l);
        Ifl_r = flut(flut_r);
        
        %FLUT potassium currents
        [Ifl_lp, dpara_n_l] = flutPotassium(flut_l,para_n_l);
        Ifl_l = Ifl_l + Ifl_lp;
        [Ifl_rp, dpara_n_r] = flutPotassium(flut_r,para_n_r);
        Ifl_r = Ifl_r + Ifl_rp;
        
        dflut_l = cableEq(Ifl_l,flut_lv,mysa_lv,interv(:,1), ...
                          V_e(i_flut(1):i_flut(2)),V_e(i_mysa(1):i_mysa(2)),V_e(i_inter(1,1):i_inter(1,2)), ...
                          r_flut,r_mysa,r_inter,c_flut);
        dflut_r = cableEq(Ifl_r,flut_rv,mysa_rv,interv(:,6), ...
                          V_e(i_flut(3):i_flut(4)),V_e(i_mysa(3):i_mysa(4)),V_e(i_inter(6,1):i_inter(6,2)), ...
                          r_flut,r_mysa,r_inter,c_flut);
        
        %currents between myelin and node at flut regions
        Ifl_l_b = passiveMyelin(flut_l_b,g_flut_m, ...
                                            V_e(i_flut(1):i_flut(2)));
        Ifl_r_b = passiveMyelin(flut_r_b,g_flut_m, ...
                                            V_e(i_flut(3):i_flut(4)));
        
        dflut_l_b = extracellular(Ifl_l_b,flut_l_b,mysa_l_b,inter_b(:,1), ...
                    r_pn2,r_pn1,r_px,c_flut_m,Ifl_l,c_flut);
        dflut_r_b = extracellular(Ifl_r_b,flut_r_b,mysa_r_b,inter_b(:,1), ...
                   r_pn2,r_pn1,r_px,c_flut_m,Ifl_r,c_flut);              
        
        %internode currents
        Iin = internodes(inter);
        inter2 = [flut_lv,interv,flut_rv];
        inter2b = [flut_l_b,inter_b,flut_r_b];
        inter2e = [V_e(i_flut(1):i_flut(2)),...
            V_e(i_inter(1,1):i_inter(1,2)),...
            V_e(i_inter(2,1):i_inter(2,2)),...
            V_e(i_inter(3,1):i_inter(3,2)),...
            V_e(i_inter(4,1):i_inter(4,2)),...
            V_e(i_inter(5,1):i_inter(5,2)),...
            V_e(i_inter(6,1):i_inter(6,2)),...
            V_e(i_flut(3):i_flut(4))];
        
        
        res = ones(6,1)*r_inter;
        res = [r_flut;res;r_flut];
        dinter = zeros(N_inter,6);
        for j = 2:7
            dinter(:,j-1) = cableEq(Iin(:,j-1),inter2(:,j),inter2(:,j-1),...
                        inter2(:,j+1),inter2e(:,j),inter2e(:,j-1),...
                        inter2e(:,j+1),res(j),res(j-1),res(j+1),c_inter);
        end
        
        Iin_b = passiveMyelin(inter_b,g_inter_m,[V_e(i_inter(1,1):i_inter(1,2)),...
            V_e(i_inter(2,1):i_inter(2,2)),V_e(i_inter(3,1):i_inter(3,2)),...
            V_e(i_inter(4,1):i_inter(4,2)),V_e(i_inter(5,1):i_inter(5,2)),...
            V_e(i_inter(6,1):i_inter(6,2))]);
        %Iin_b = Iin_b;
        res = [r_pn2;ones(6,1)*r_px;r_pn2];
        
        dinter_b = zeros(N_inter,6);
        for j = 2:7
            dinter_b(:,j-1) = extracellular(Iin_b(:,j-1),inter2b(:,j),...
                        inter2b(:,j-1),inter2b(:,j+1), ...
                        res(j),res(j-1),res(j+1),c_inter_m,Iin(:,j-1),c_inter);
        end
        
        %finally assemble derivatives into an
        dY=[dnode;dpara_m;dpara_h;dpara_p;dpara_s;dmysa_l;dmysa_r;...
            dflut_l;dflut_r;dpara_n_l;dpara_n_r;dinter(:,1); ...
            dinter(:,2);dinter(:,3);dinter(:,4);dinter(:,5);...
            dinter(:,6);dmysa_l_b;dmysa_r_b;dflut_l_b;dflut_r_b;...
            dinter_b(:,1);dinter_b(:,2);dinter_b(:,3);...
            dinter_b(:,4);dinter_b(:,5);dinter_b(:,6)];
        
    end
    
    % cableEq: calculates the cable-equation
    function dV = cableEq(I,V,V1,V2,Ve,Ve1,Ve2,Ra,Ra1,Ra2,C)
        dV = (-I + (V1-V)./(Ra1./2+Ra./2) + (V2-V)./(Ra2./2+Ra./2) + ...
              (Ve1-Ve)./(Ra1./2+Ra./2) + (Ve2-Ve)./(Ra2./2+Ra./2))./C;
    end

    % extracellular: extracellular currents similar to the Neuron function
    %                with the same name
    function dV = extracellular(I,V,V1,V2,Ra,Ra1,Ra2,C,Imem,Cmem)
       %dV = (-I + Imem+ (V1-V)./(Ra1./2+Ra./2) + (V2-V)./(Ra2./2+Ra./2))./C;
       dV = -I./C + Imem./Cmem + ...
           ((V1-V)./(Ra1./2+Ra./2) + (V2-V)./(Ra2./2+Ra./2))./(C+Cmem);
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
    function I = flut(V)
        g=g_flut;
        e=vrest;
        I = passiveCurrent(V,g,e);
    end
    
    % mysa: simple mysa current
    function I = mysa(V)
        g=g_mysa;
        e=vrest;
        I = passiveCurrent(V,g,e);
    end

    % inter: currents of the internode
    function I = internodes(V)
        g=g_inter;
        e=vrest;
        I = passiveCurrent(V,g,e);
    end

    % passiveCurrent: calculation of the passive current
    function I = passiveCurrent(V,g,e)
        I = g.*(V-e);
    end

    % passiveCurrent: calculation of the passive current
    function I = passiveMyelin(V,g,e)
        I=passiveCurrent(V,g,e);
    end

    % dpdt: calculates the derivative of the of the parameter according to
    %       alpha, beta and the previous value
    function dp = dpdt(alpha, beta, para)
        dp = alpha.*(1-para)-beta.*para;
    end


    function V = electrode(q,xe,ye,x,y )
        d = sqrt((xe-x).^2+(ye-y).^2);
        V = q./(4.*pi.*d);
    end

end