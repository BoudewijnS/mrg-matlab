function draw( t,Y,N_nodes,range )
%draw draws Y in range 

    t1 = [];
    Y1 = [];
    
    i = 1;
    while t(i) <= range(2)
       if t(i) >= range(1)-0.0001
          t1 = [t1;t(i)]; 
          Y1 = [Y1;Y(i,1:N_nodes)];
          i = i+1;
       end
    end
    
    Y2 = [];
    for j = 3:N_nodes
        Y2(j,:) = Y1(:,j) - (j-1)*5; 
    end
    Y3 = [];
    for j = 5:N_nodes-4
        Y3(j,:) = Y1(:,j) - (j-1)*20; 
    end
    %figure(5);
    %subplot(1,2,1);
    plot(t1,Y3,'k');
    %subplot(1,2,2);
    %plot(t1,Y3);
end

