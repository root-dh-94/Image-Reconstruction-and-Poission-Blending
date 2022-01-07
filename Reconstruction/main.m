clc;
clear;
close all;

imgin = im2double(imread('./target.jpg'));

[imh, imw, nb] = size(imgin);
assert(nb==1);
% the image is grayscale

V = zeros(imh, imw);
V(1:imh*imw) = 1:imh*imw;
% V(y,x) = (y-1)*imw + x
% use V(y,x) to represent the variable index of pixel (x,y)
% Always keep in mind that in matlab indexing starts with 1, not 0

%Initializing i,j,v for sparse call
i = zeros(1,(imh-2)*(imw-2)*5+(2*(imh+imw)-8)*3+4);
j = zeros(1,(imh-2)*(imw-2)*5+(2*(imh+imw)-8)*3+4);
v = zeros(1,(imh-2)*(imw-2)*5+(2*(imh+imw)-8)*3+4);
b = zeros(imh*imw,1);
count = 1;
e = 1;

%looping through inner pixels to assign i,j,v vals-5variables
for y = 2:imh-1
    for x = 2:imw-1
        for der = -1:1            
            if der == 0 
                i(count) = e;
                v(count) = 4;                
                j(count) = V(y,x);
                count = count +1;
            else v(count) = -1;
                v(count+1) = -1;
                i(count) = e;
                i(count+1) = e;
                j(count) = V(y,x+der);
                j(count+1) = V(y+der,x);
                count = count+2;          
            end            
        end
        %Assigning values in b to correspond to equations of inner pixels
        b(e) = 4*imgin(y,x)-imgin(y,x+1)-imgin(y,x-1)-imgin(y+1,x)-imgin(y-1,x);
        e = e+1;
    end
end

%Looping through vertical pixels to assign ijv-3variables
for y = 2:imh-1
    
    
        for der = -1:1
            i(count) = e;
            i(count+1) = e+1;
            if der == 0 
                v(count) = 2;
                v(count+1) = 2;
                j(count) = V(y,1);
                j(count+1) = V(y,imw);
                count = count +2;
            else v(count) = -1;
                v(count+1) = -1;
                j(count) = V(y+der,1);
                j(count+1) = V(y+der,imw);
                count = count+2;
           
            end
            
        end
        b(e) = 2*imgin(y,1)-imgin(y+1,1)-imgin(y-1,1);
        b(e+1) = 2*imgin(y,imw)-imgin(y+1,imw)-imgin(y-1,imw);
        e = e+2;
    
end
%Looping through horizontal pixels to assign ijv-3variables
for x = 2:imw-1
    
    
        for der = -1:1
            i(count) = e;
            i(count+1) = e+1;
            if der == 0 
                v(count) = 2;
                v(count+1) = 2;
                j(count) = V(1,x);
                j(count+1) = V(imh,x);
                count = count +2;
            else v(count) = -1;
                v(count+1) = -1;
                j(count) = V(1,x+der);
                j(count+1) = V(imh,x+der);
                count = count+2;
           
            end
            
        end
        b(e) = 2*imgin(1,x)-imgin(1,x+1)-imgin(1,x-1);
        b(e+1) = 2*imgin(imh,x)-imgin(imh,x+1)-imgin(imh,x-1);
        e = e+2;
    
end

%Initialize control points
i(count) = e;
j(count) = 1;
v(count) = 1;
b(e) = imgin(1,1);

i(count+1) = e+1;
j(count+1) = V(imw,1);
v(count+1) = 1;
b(e+1) = imgin(imw,1);

i(count+2) = e+2;
j(count+2) = V(1,imh);
v(count+2) = 1;
b(e+2) = imgin(1,imh);

i(count+3) = e+3;
j(count+3) = V(imw,imh);
v(count+3) = 1;
b(e+3) = imgin(imw,imh);

A = sparse(i,j,v);

%Solving equation
solution = A\b;
error = sum(abs(A*solution-b));
disp(error)
imgout = reshape(solution,[imh,imw]);

imwrite(imgout,'output.png');
figure(), hold off, imshow(imgout);

