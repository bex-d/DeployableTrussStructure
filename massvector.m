function[M]=massvector(X,iX, mass)      %converts mass from a [1 x number of types of truss] matrix to [1 x number of trusses] matrix

    Materials = zeros(length(iX),1);        
    
        for i=1:length(iX)                  
            j = iX(i,3);                    %j = truss type
            Materials(i) = mass(j);         %creates [1 x number of trusses] matrix with mass of each truss
        end
        
  
    M=zeros(length(X),1);                   %Mass Matrix - empty, preallocated for speed
    
        for i=1:length(iX)                  

            node1 = iX(i,1);                %identifies one end of the truss
            node2 = iX(i,2);                %identifies the other end of the truss

            M(node1) = M(node1) + 0.5*Materials(i); %assigns half the mass of the truss to each node
            M(node2) = M(node2) + 0.5*Materials(i);

        end

end