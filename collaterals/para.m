classdef para
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vrest;        %mV
        fiberD;        %um
        mysalength;       %um  
        nodelength;     %um
        space_p1;     %um
        space_p2;     %um
        space_i;      %um
        r;               %Ohm*cm
        mycm;           %uF/cm2
        mygm;         %S/cm2
        c;                %uF/cm2

        g_nap;      % sodium channel conductivity [S/cm^2]
        g_naf;          %S/cm^2
        g_k;        %S/cm^2
        g_kf;           %S/cm^2 (deactivated)
        g_l;        % leak channel conductivity [S/cm^2] 
        e_na;        %mV
        e_k;        %mV
        e_l;        %mV

        g_p1;       %S/cm^2
        g_p2;      %S/cm^2 
        g_i;         %S/cm^2
        axonD;          %um 
        nodeD;          %um
        mysaD;         %um
        flutD;         %um
        deltax;        %um
        flutlength;      %um
        nl;             %dimensionless


       
        interlength;

        celsius;
        q10_1;
        q10_2;
        q10_3;


        r_mysa
        r_pn1
        c_mysa
        c_mysa_m
        g_mysa
        g_mysa_m

        r_flut
        r_pn2
        c_flut
        c_flut_m
        g_flut
        g_flut_m

        r_inter
        r_px
        c_inter
        c_inter_m
        g_inter
        g_inter_m

%         g_nap
%         g_k
%         g_l
%         g_naf
        r_node 
        c_node 
        r_pn0
    end
    
    methods(Static)
        function [rax,rpx,cmem,cmy,gpas,gmy] = ...
            paracomp(diameter,length,space, fiberDia, c, r, g, nl, xc, xg)
           %args must be in um, ohm*cm and uF respectively
           rax = para.calcResAxial(diameter,length,r);
           rpx = para.calcResPeriax(diameter,length,space,r);
           cmem = para.calcCapacity(diameter, length, c);
           cmy = para.calcMyelinCap(fiberDia, length, xc, nl);
           gpas = para.calcConductance(diameter, length, g);
           gmy = para.calcMyelinCond(fiberDia, length, xg, nl);
        end

        function rax = calcResAxial(diameter, length, r)
           %um and Ohm*cm
           a=pi*diameter^2;
           rax = para.calcRes(a,length,r);
        end

        function rpx = calcResPeriax(diameter, length, space, r)
            %um and Ohm*cm
            a = pi*((diameter+(2*space))^2-diameter^2);
            rpx = para.calcRes(a,length,r);
        end

        function res = calcRes(area, length, r)

            ra = r*1e-5;     %GOhm*um
            res = (4*(length)*ra)/area;
        end

        function c = calcMyelinCap(fiberDiameter, length, xc, nl)
            xc = xc/(nl*2);
            c = para.calcCapacity(fiberDiameter, length, xc);
            %c is in uF
        end

        function c = calcCapacity(diameter, length, cap)
           %um and uF/cm^2
           cap = cap*1e-2;   %pF/um^2 
           c = cap*diameter*length*pi;
           %c is in uF
        end

        function gmy = calcMyelinCond(fiberDiameter, length, xg, nl)
           xg = xg/(nl*2);
           gmy = para.calcConductance(fiberDiameter, length, xg);
        end

        function g = calcConductance(diameter, length, gi)
            %in S/cm^2 and um
            gi = gi*1e1; %nS/um^2
            g = gi*diameter*length*pi;
        end
    end
    
    
    methods
        function obj =  para(diameter)
            obj.vrest = -80;        %mV
            obj.mysalength=3;       %um  
            obj.nodelength=1.0;     %um
            obj.space_p1=0.002;     %um
            obj.space_p2=0.004;     %um
            obj.space_i=0.004;      %um
            obj.r=70;               %Ohm*cm
            obj.mycm=0.1;           %uF/cm2
            obj.mygm=0.001;         %S/cm2
            obj.c=2;                %uF/cm2

            obj.g_nap = 0.010;      % sodium channel conductivity [S/cm^2]
            obj.g_naf = 3;          %S/cm^2
            obj.g_k = 0.080;        %S/cm^2
            obj.g_kf = 0;           %S/cm^2 (deactivated)
            obj.g_l = 0.007;        % leak channel conductivity [S/cm^2] 
            obj.e_na = 50.0;        %mV
            obj.e_k = -90.0;        %mV
            obj.e_l	= -90.0;        %mV

            obj.g_p1 = 0.001;       %S/cm^2
            obj.g_p2 = 0.0001;      %S/cm^2 
            obj.g_i = obj.g_p2;         %S/cm^2
            obj.axonD=6.9;          %um 
            obj.nodeD=3.3;          %um
            obj.mysaD=3.3;         %um
            obj.flutD=6.9;         %um
            obj.deltax=1250;        %um
            obj.flutlength=46;      %um
            obj.nl=120;             %dimensionle

            obj.celsius = 36;
            obj.q10_1 = 2.2^((obj.celsius-20)/10);
            obj.q10_2 = 2.9^((obj.celsius-20)/10);
            obj.q10_3 = 3.0^((obj.celsius-36)/10);
            if diameter == 2.0
                obj.axonD=1.6;          %um 
                obj.nodeD=1.4;          %um
                obj.mysaD=1.4;         %um
                obj.flutD=1.6;         %um
                obj.deltax=373.2;        %um
                obj.flutlength=10;      %um
                obj.nl=30;             %dimensionless
            end

            if diameter == 5.7
                obj.axonD=3.4;          %um 
                obj.nodeD=1.9;          %um
                obj.mysaD=1.9;         %um
                obj.flutD=3.4;         %um
                obj.deltax=500;        %um
                obj.flutlength=35;      %um
                obj.nl=80;             %dimensionless
            end

            if diameter == 16.0
                obj.axonD=12.7;          %um 
                obj.nodeD=5.5;          %um
                obj.mysaD=5.5;         %um
                obj.flutD=12.7;         %um
                obj.deltax=1500;        %um
                obj.flutlength=60;      %um
                obj.nl=150;             %dimensionless
            end

            if diameter == 10.0
                obj.axonD=6.9;          %um 
                obj.nodeD=3.3;          %um
                obj.mysaD=3.3;         %um
                obj.flutD=6.9;         %um
                obj.deltax=1150;        %um
                obj.flutlength=46;      %um
                obj.nl=120;             %dimensionless

            end 
            if diameter == 11.5
                obj.axonD=8.1;          %um 
                obj.nodeD=3.7;          %um
                obj.mysaD=3.7;         %um
                obj.flutD=8.1;         %um
                obj.deltax=1250;        %um
                obj.flutlength=50;      %um
                obj.nl=130;             %dimensionless
        
            end
    

            obj.fiberD=diameter;
            obj.interlength=(obj.deltax-obj.nodelength-(2*obj.mysalength)-(2*obj.flutlength))/6;
            [obj.r_mysa,obj.r_pn1,obj.c_mysa,obj.c_mysa_m,obj.g_mysa,obj.g_mysa_m] = obj.paracomp( obj.mysaD, obj.mysalength, obj.space_p1, obj.fiberD, obj.c, obj.r, obj.g_p1, obj.nl, obj.mycm, obj.mygm);

            [obj.r_flut,obj.r_pn2,obj.c_flut,obj.c_flut_m,obj.g_flut,obj.g_flut_m] = obj.paracomp( ...
                obj.flutD, obj.flutlength, obj.space_p2, obj.fiberD, obj.c, obj.r, obj.g_p2, obj.nl, obj.mycm, obj.mygm);

            [obj.r_inter,obj.r_px,obj.c_inter,obj.c_inter_m,obj.g_inter,obj.g_inter_m] = obj.paracomp( ...
                obj.axonD, obj.interlength, obj.space_i, obj.fiberD, obj.c, obj.r, obj.g_i, obj.nl, obj.mycm, obj.mygm);

            obj.g_nap = obj.calcConductance(obj.nodeD,obj.nodelength,obj.g_nap);
            obj.g_k = obj.calcConductance(obj.nodeD,obj.nodelength,obj.g_k);
            obj.g_l = obj.calcConductance(obj.nodeD,obj.nodelength,obj.g_l);
            obj.g_naf = obj.calcConductance(obj.nodeD,obj.nodelength,obj.g_naf);
            obj.r_node = obj.calcResAxial(obj.nodeD,obj.nodelength,obj.r);
            obj.c_node = obj.calcCapacity(obj.nodeD,obj.nodelength,obj.c);
            obj.r_pn0 = obj.calcResPeriax(obj.nodeD,obj.nodelength,obj.space_p1,obj.r);
        end
    end
    
end

