function der = derivative(x, dt)
[rows columns] = size(x);

if(columns==1)
    if(numel(x)>=7)
        der(1,1) = (x(2)-x(1))/dt;
        der(2,1) = (x(3)-x(1))/(2*dt);
        der(rows-1,1) = (x(rows)-x(rows-2))/(2*dt);
        der(rows,1) = (x(rows)-x(rows-1))/dt;
        for i=3:rows-2
            der(i,1) = (-x(i+2)+8*x(i+1)-8*x(i-1)+x(i-2))/(12*dt);
        end
    else
    	der = [0; diff(x)];
    end
    der = sgolayfilt(der,3,7);
else if(rows==1)
        if(numel(x)>=7)
            der(1,1) = (x(2)-x(1))/dt;
            der(1,2) = (x(3)-x(1))/(2*dt);
            der(1,rows-1) = (x(rows)-x(rows-2))/(2*dt);
            der(1,rows) = (x(rows)-x(rows-1))/dt;
            for i=3:rows-2
                der(1,i) = (-x(i+2)+8*x(i+1)-8*x(i-1)+x(i-2))/(12*dt);
            end
        else
            der = [0 diff(x)];
        end
        der = sgolayfilt(der,3,7);
    else
        display('Error: matrix derivatives not allowed')
        der=-1;
    end
end

