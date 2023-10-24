function[time] = finishTime(x)
time = 'error'
for i=5:length(x)
   
    a = x(1,i);
    b = x(1,i-1);
    c = x(1,i-2);
    d = x(1,i-3);
    e = x(1,i-4);
    
    if a == b && a==c && a==d && a==e
        time = i-4;
    end
    
end
end