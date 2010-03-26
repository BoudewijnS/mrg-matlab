classdef index
    
    properties
        N_nodes;
        N_inter;
        i_node;
        i_para_m;
        i_para_h;
        i_para_p;
        i_para_s;
        i_mysa;
        i_flut;
        i_para_n;
        i_inter;
        i_mysa_b;
        i_flut_b;
        i_inter_b
    end
    
    methods
        function obj = index(N)
            obj.N_nodes = N;
            obj.N_inter = obj.N_nodes-1;
            % index values of the dV Vektor
            obj.i_node = [1,obj.N_nodes];
            obj.i_para_m = [obj.i_node(2)+1,obj.i_node(2)+obj.N_nodes];
            obj.i_para_h = [obj.i_para_m(2)+1,obj.i_para_m(2)+obj.N_nodes];
            obj.i_para_p = [obj.i_para_h(2)+1,obj.i_para_h(2)+obj.N_nodes];
            obj.i_para_s = [obj.i_para_p(2)+1,obj.i_para_p(2)+obj.N_nodes];
            obj.i_mysa = [obj.i_para_s(2)+1,obj.i_para_s(2)+obj.N_inter, ...
                obj.i_para_s(2)+obj.N_inter+1,obj.i_para_s(2)+obj.N_inter*2];
            obj.i_flut = [obj.i_mysa(4)+1,obj.i_mysa(4)+obj.N_inter, ... 
                obj.i_mysa(4)+obj.N_inter+1,obj.i_mysa(4)+obj.N_inter*2];
            obj.i_para_n = [obj.i_flut(4)+1,obj.i_flut(4)+obj.N_inter, ...
                obj.i_flut(4)+obj.N_inter+1,obj.i_flut(4)+obj.N_inter*2];
            obj.i_inter = [obj.i_para_n(4)+1,obj.i_para_n(4)+obj.N_inter; ... 
                obj.i_para_n(4)+obj.N_inter+1,obj.i_para_n(4)+obj.N_inter*2; ...
                obj.i_para_n(4)+(2*obj.N_inter)+1,obj.i_para_n(4)+(3*obj.N_inter);...
                obj.i_para_n(4)+(3*obj.N_inter)+1,obj.i_para_n(4)+(4*obj.N_inter);...
                obj.i_para_n(4)+(4*obj.N_inter)+1,obj.i_para_n(4)+(5*obj.N_inter);...
                obj.i_para_n(4)+(5*obj.N_inter)+1,obj.i_para_n(4)+(6*obj.N_inter)];
            obj.i_mysa_b = [obj.i_inter(12)+1,obj.i_inter(12)+obj.N_inter, ...
                        obj.i_inter(12)+obj.N_inter+1,obj.i_inter(12)+obj.N_inter*2];
            obj.i_flut_b = [obj.i_mysa_b(4)+1,obj.i_mysa_b(4)+obj.N_inter, ... 
                obj.i_mysa_b(4)+obj.N_inter+1,obj.i_mysa_b(4)+obj.N_inter*2];
            obj.i_inter_b = [obj.i_flut_b(4)+1,obj.i_flut_b(4)+obj.N_inter; ... 
                obj.i_flut_b(4)+obj.N_inter+1,obj.i_flut_b(4)+obj.N_inter*2; ...
                obj.i_flut_b(4)+(2*obj.N_inter)+1,obj.i_flut_b(4)+(3*obj.N_inter);...
                obj.i_flut_b(4)+(3*obj.N_inter)+1,obj.i_flut_b(4)+(4*obj.N_inter);...
                obj.i_flut_b(4)+(4*obj.N_inter)+1,obj.i_flut_b(4)+(5*obj.N_inter);...
                obj.i_flut_b(4)+(5*obj.N_inter)+1,obj.i_flut_b(4)+(6*obj.N_inter)];
        end
    end
end
