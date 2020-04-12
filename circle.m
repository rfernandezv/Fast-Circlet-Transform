function I = circle(r,x0,y0,row,col)
%%% Circle function create a circle with radius r and center
%%% [x0,y0] in th output image I with dimension of row*col.
I = ones(row,col);
for x=1:row
    for y=1:col
        if ((x-x0)^2+(y-y0)^2) <= r^2
            I(x,y) = 0;
        end
    end
end