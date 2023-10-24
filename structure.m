function[X,iX,L0]=structure(N,H,R,R_dash)
%% Dependent Variables & Relations

alpha = 2*pi()/N; %angle, rad                                               
D = 2*R*sin(alpha/2); %Outer horizontal members, m                         
r = R*(1-sin(alpha/2))/(1+sin(alpha/2)); %radius of inner ring, m
d = 2*r*sin(alpha/2); %Inner horizontal members, m

%% Configuration Parameter
r_dash = R_dash*(1-sin(alpha/2))/(1+sin(alpha/2)); %r' is inner radius whilst deploying  
theta = acos(R_dash/R); %deployment angle

%% Matrixes - X & iX
x = zeros(4*N,1);
y = zeros(4*N,1);
z = zeros(4*N,1);
iX = zeros(11*N,3);
L0 = zeros(length(iX),1);

    for i = 1:N  %coordinate definition
        
        m = i+N;        %top coordinate of vertical members outer
        n = i;          %bottom coordinate of vertical members outer
        o = i+(2*N);    %bottom coordinate of vertical members inner
        p = i+(3*N);    %top coordinate of vertical members inner
        
        if mod(i,2) == 1   
            z(m) = H;    
            z(n) = 0;        
            z(o) = (R-r)*sin(theta);
            z(p) = H+(R-r)*sin(theta);
        else 
            z(m) = H + D*sin(theta);
            z(n) = D*sin(theta);
            z(o) = (R-r)*sin(theta) - d*sin(theta);
            z(p) = H+ (R-r)*sin(theta) - d*sin(theta);
        end
        
        x(n) = R_dash * cos(alpha * (n-1));        y(n) = R_dash * sin(alpha * (n-1));                        
        x(m) = R_dash * cos(alpha * (n-1));        y(m) = R_dash * sin(alpha * (n-1));                                    
        x(o) = r_dash * cos(alpha * (n-1));        y(o) = r_dash * sin(alpha * (n-1));                           
        x(p) = r_dash * cos(alpha * (n-1));        y(p) = r_dash * sin(alpha * (n-1));        

        X = [x y z];

    end


    for i = 1:N  %connectivity matrix definition

        t = 1;                                                                  %type 1 - Horizontal outer ring members  LENGTH CONSTANT
        if i<N
            iX(i,:) = [i,i+1,t];                                                    %Outer Ring Bottom Horizontal Members
            iX(i+N,:) = [i+N,i+N+1,t];                                              %Outer Ring Top Horizontal Members
        elseif i==N
            iX(i,:) = [1,N,t];                                                      %Outer Ring Bottom Horizontal Members 
            iX(i+N,:) = [2*N,1+N,t,];                                               %Outer Ring Top Horizontal Members 
        end

        t = 2;                                                                  %type 2 - Horizontal inner ring members  LENGTH CONSTANT
        if i<N
            iX(i+2*N,:) = [i+2*N,i+2*N+1,t];                                        %Inner Ring Bottom Horizontal Members
            iX(i+3*N,:) = [i+3*N,i+3*N+1,t];                                        %Inner Ring Top Horizontal Members
        elseif i==N
            iX(i+2*N,:) = [1+3*N,4*N,t];                                            %Inner Ring Top Horizontal Members
            iX(i+3*N,:) = [1+2*N,3*N,t];                                            %Inner Ring Bottom Horizontal Members 
        end

        t = 3;                                                                  %type 3 - vertical members  LENGTH CONSTANT
        iX(i+4*N,:) = [i,i+N,t];                                                    %VERTICAL OUTSIDE               
        iX(i+5*N,:) = [i+(2*N),i+(3*N),t];                                          %VERTICAL INSIDE 

        t = 4;                                                                  %type 4 - horizontal members  LENGTH CONSTANT
        iX(i+6*N,:) = [i+(2*N),i,t];                                                %HORIZONTAL BOTTOM
        iX(i+7*N,:) = [i+N,i+(3*N),t];                                              %HORIZONTAL TOP

        t = 5;                                                                  %type 5 - OUTSIDE DIAGONAL TRUSSES LENGTH NOT CONSTANT
        if i==1                                                                       
            iX(i+8*N,:) = [i,i+N+1,t];
            iX(i+1+8*N,:) = [i,2*N,t];
        elseif i>1 && (mod(i,2) == 1)                                                
            iX(i+8*N,:) = [i, i+N+1,t];
            iX(i+1+8*N,:) = [i, i+N-1,t]; 
        end

        t = 6;                                                                  %type 6 - INSIDE DIAGONAL TRUSSES   LENGTH NOT CONSTANT
        if i==1                                                                   
            iX(i+9*N,:) = [3*N+i,i+2*N+1,t];                                       
            iX(i+1+9*N,:) = [3*N+i,3*N,t];                                    
        elseif i>1 && (mod(i,2) == 1)                                            
            iX(i+9*N,:) =  [2*N+i-1, i+3*N,t];                                   
            iX(i+1+9*N,:) = [2*N+i+1, i+3*N,t];                                   
        end                                        

        t = 7;                                                                  %type 7 - DIAGONAL TRUSSES BETWEEN RINGS LENGTH NOT CONSTANT
        if mod(i,2)==1                              
            iX(i+10*N,:) = [i,3*N+i,t];                                                   
        else
            iX(i+10*N,:) = [N+i,2*N+i,t];                                       
        end 
        
        t = 8;                                                                  %type 8 - DIAGONAL TRUSSES BETWEEN OUTSIDE AND INSIDE (TOP)
        if i==N
            iX(i+11*N,:) =  [3*N+i,1+N,t];
        elseif mod(i,2)==1   
            iX(i+11*N,:) = [N+i,3*N+i+1,t];
        else             
            iX(i+11*N,:) = [3*N+i,1+N+i,t];
        end
      
        t = 9;                                                                  %type 9 - DIAGONAL TRUSSES BETWEEN OUTSIDE AND INSIDE (BOTTOM)
        if i==N
            iX(i+12*N,:) =  [i,i+1+N,t];
        elseif mod(i,2)==1
            iX(i+12*N,:) = [2*N+i,i+1,t];
        else  
            iX(i+12*N,:) = [i,2*N+1+i,t];
        end
        
    end

    for i = 1:length(iX)
        L0(i) = len(X(iX(i,1),:),X(iX(i,2),:));
    end

end