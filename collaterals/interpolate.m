function [Ve,N,x_af,af,Vlr] = interpolate(file,fiberD)
        parameters;

        data = importdata(file);
        % to mV
        Ve_pulse = 1e3*data(:,4);
        % to cm
        x = 100*data(:,1);
        y = 100*data(:,2);
        z = 100*data(:,3);
        s(1) = 0;
        for i=1:length(x)-1
            s(i+1) = s(i) + sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2 +...
                (z(i+1)-z(i))^2);
        end
        n = floor(s(length(x)-1)/(deltax*1e-4))+16;
        N = n+1;
        
        half = 0;
        
        if half == 0
            [x_n,x_m,x_f,x_i] = calcX(N);
            x_n=x_n-8*(deltax+1)*1e-4;
            x_m=x_m-8*(deltax+1)*1e-4;
            x_f=x_f-8*(deltax+1)*1e-4;
            x_i=x_i-8*(deltax+1)*1e-4;
           
        end
        if half == 1
            [x_n,x_m,x_f,x_i] = calcX(floor(N/2));
            x_n = x_n + floor(N/2)*(deltax+1)*1e-4;
            x_m = x_m + floor(N/2)*(deltax+1)*1e-4;
            x_f = x_f + floor(N/2)*(deltax+1)*1e-4;
            x_i = x_i + floor(N/2)*(deltax+1)*1e-4;
            N = floor(N/2);
            n = N-1;
        end
        if half == 2
            [x_n,x_m,x_f,x_i] = calcX(floor(N/2));
            N = floor(N/2);
            n = N-1;
        end  
        
        xlr=[x_n(1)-2e-4,x_n(N)+2e-4];
        
        displacment = -1;
        x_n=x_n+displacment*1e-4*deltax/10;
        x_m=x_m+displacment*1e-4*deltax/10;
        x_f=x_f+displacment*1e-4*deltax/10;
        x_i=x_i+displacment*1e-4*deltax/10;
        xlr=xlr+displacment*1e-4*deltax/10;
        
        Vlr = interp1(s,Ve_pulse,xlr,'linear','extrap');
        
        Ve_n = interp1(s,Ve_pulse,x_n,'linear','extrap');
        Ve_m = [interp1(s,Ve_pulse,x_m(:,1),'linear','extrap'); ...
                interp1(s,Ve_pulse,x_m(:,2),'linear','extrap')];
        Ve_f = [interp1(s,Ve_pulse,x_f(:,1),'linear','extrap'); ...
                interp1(s,Ve_pulse,x_f(:,2),'linear','extrap')];
        Ve_i = [interp1(s,Ve_pulse,x_i(:,1),'linear','extrap'); ...
                interp1(s,Ve_pulse,x_i(:,2),'linear','extrap'); ...
                interp1(s,Ve_pulse,x_i(:,3),'linear','extrap'); ...
                interp1(s,Ve_pulse,x_i(:,4),'linear','extrap'); ...
                interp1(s,Ve_pulse,x_i(:,5),'linear','extrap'); ...
                interp1(s,Ve_pulse,x_i(:,6),'linear','extrap')];
        Ve = [Ve_n;zeros(N*4,1);Ve_m;Ve_f;zeros((N-1)*2,1);Ve_i];
        x_af = [];
        af = [];
        expot = [];
        
         for i = 1:n
             
             x_af = [x_af,x_n(i),x_m(i,2),x_f(i,2),x_i(i,6),...
                 x_i(i,5),x_i(i,4),x_i(i,3),x_i(i,2),x_i(i,1),...
                 x_f(i,1),x_m(i,1)];
             expot = [expot,Ve_n(i),Ve_m(i+n),Ve_f(i+n),Ve_i(i+5*n),...
                 Ve_i(i+4*n),Ve_i(i+3*n),Ve_i(i+2*n),Ve_i(i+n),Ve_i(i),...
                 Ve_f(i),Ve_m(i)];
         end
         x_af = [x_af,x_n(N)];
         expot = [expot,Ve_n(N)];
        
         figure(10);
         plot(x_af,expot);
        
         
        function [x_n,x_m1,x_f1,x_i1] = calcX(N)
            x_n(1) = 0.5*1e-4;
            x_m(1,1) = x_n(1) + 0.5*1e-4 + mysalength/2*1e-4;
            x_m(1,2) = x_n(1) + (1+deltax)*1e-4 -mysalength/2*1e-4;
            x_f(1,1) = x_m(1,1) + flutlength/2*1e-4 + mysalength/2*1e-4;
            x_f(1,2) = x_m(1,2) - flutlength*1e-4 - mysalength/2*1e-4;
            x_i = zeros(N-1,6);
            x_i(1,1) = x_f(1,1) + flutlength/2*1e-4 +interlength/2*1e-4;
            for j = 2:6
                x_i(1,j) = x_i(1,1) + (j-1)*interlength*1e-4;
            end
            for i=2:N
                x_n(i) = x_n(i-1) + deltax*1e-4+1e-4;
            end
            N
            for i=2:N-1
                x_m(i,1) = x_m(i-1,1) + (1+deltax)*1e-4;
                x_m(i,2) = x_m(i-1,2) + (1+deltax)*1e-4;
                x_f(i,1) = x_f(i-1,1) + (1+deltax)*1e-4;
                x_f(i,2) = x_f(i-1,2) + (1+deltax)*1e-4;
                for j = 1:6
                    x_i(i,j) = x_i(i-1,j) + (1+deltax)*1e-4;
                end 
            end

            x_m1(:,1) = x_m(:,2);
            x_m1(:,2) = x_m(:,1);
            x_f1(:,1) = x_f(:,2);
            x_f1(:,2) = x_f(:,1);
            x_i1=zeros(N-1,6);
            for i=1:6
                x_i1(:,i) = x_i(:,7-i);
            end
            x_n=x_n';
        end
         
         
    end

