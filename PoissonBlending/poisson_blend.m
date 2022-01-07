function imgout = poisson_blend(im_s, mask_s, im_t)
% -----Input
% im_s     source image (object)
% mask_s   mask for source image (1 meaning inside the selected region)
% im_t     target image (background)
% -----Output
% imgout   the blended image

[imh, imw, nb] = size(im_t);

% V(y,x) = (y-1)*imw + x
% use V(y,x) to represent the variable index of pixel (x,y)
% Always keep in mind that in matlab indexing starts with 1, not 0
V = zeros(imh, imw);
V(1:imh*imw) = 1:imh*imw;
%find number of pixels in mask
k = sum(mask_s(:) ==1);

%identify indeces of inner mask pixels as well as border pixel types
row_inside = [];
column_inside = [];
pixel_inside = [];

%for boundary pixel in x direction(single out of mask pixel neighbour)
row_border_1_x = []; 
column_border_1_x = [];
pixel_border_1_x = [];

%for boundary pixel in y direction(single out of mask pixel neighbour)
row_border_1_y = [];
column_border_1_y = [];
pixel_border_1_y = [];

%for boundary pixel with 2 out of mask pixel as neighbours
row_border_2 = [];
column_border_2 = [];
pixel_border_2 = [];

%for boundary pixel with 3 out of mask pixel as neighbours
row_border_3 = [];
column_border_3 = [];
pixel_border_3 = [];

%Loop through all pixels to determin which are inside mask and which are on
%boundary, along with number of neighbours outside mask
for y = 2:imh-1
    for x = 2:imw-1
        if mask_s(y,x) == 1
            if mask_s(y,x) == mask_s(y-1,x) && mask_s(y,x) == mask_s(y+1,x) && mask_s(y,x) == mask_s(y,x-1) && mask_s(y,x) == mask_s(y,x+1)
                row_inside(end+1) = y;
                column_inside(end+1) = x;
                pixel_inside(end+1) = V(y,x);
                
            elseif (mask_s(y-1,x)==0 && mask_s(y+1,x) == 0) || (mask_s(y,x-1)==0 && mask_s(y,x+1) == 0)
                row_border_3(end+1) = y;
                column_border_3(end+1) = x;
                pixel_border_3(end+1) = V(y,x);
                
            elseif (mask_s(y-1,x)==0 || mask_s(y+1,x) == 0) && (mask_s(y,x-1)==0 || mask_s(y,x+1) == 0)
                row_border_2(end+1) = y;
                column_border_2(end+1) = x;
                pixel_border_2(end+1) = V(y,x);
                
            elseif mask_s(y-1,x)==0 || mask_s(y+1,x) == 0
                row_border_1_y(end+1) = y;
                column_border_1_y(end+1) = x;
                pixel_border_1_y(end+1) = V(y,x);
                
            elseif (mask_s(y,x-1)==0 || mask_s(y,x+1) == 0)
                row_border_1_x(end+1) = y;
                column_border_1_x(end+1) = x;
                pixel_border_1_x(end+1) = V(y,x);
            end
        end
    end
end
                
%calculating number of different types of border pixels               
boundary_dim_x_der = size(column_border_1_x,2);
boundary_dim_y_der = size(column_border_1_y,2);
boundary_dim_xy_der = size(column_border_2,2);
boundary_dim_3_der = size(column_border_3,2);
num_boundary = boundary_dim_x_der + boundary_dim_y_der + boundary_dim_xy_der + boundary_dim_3_der;

%Recalculating number of equations and variables to determine size of i,j,v
sparse_dim =(imh-2)*(imw-2)*5+(2*(imh+imw)-8)*3+4-5*num_boundary+4*(boundary_dim_x_der + boundary_dim_y_der)+3*boundary_dim_xy_der + 2*boundary_dim_3_der;

%Looping though one channel at a time
for channel = 1:nb
    img_s_chan = im_s(:,:,channel);
    imgin = im_t(:,:,channel);
    
%Initializing i,j,v for sparse call   
    i = zeros(1,sparse_dim);
    j = zeros(1,sparse_dim);
    v = zeros(1,sparse_dim);
    b = zeros(imh*imw,1);
    count = 1;
    e = 1;
 
    %looping through inner pixels to assign i,j,v vals-5variables
    for y = 2:imh-1
        for x = 2:imw-1
            if (mask_s(y,x) == mask_s(y-1,x) && mask_s(y,x) == mask_s(y+1,x) && mask_s(y,x) == mask_s(y,x-1) && mask_s(y,x) == mask_s(y,x+1))|| mask_s(y,x)==0            
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
            
                %If inner image pixel but not in mask
                if mask_s(y,x) == 0
                    b(e) = 4*imgin(y,x)-imgin(y,x+1)-imgin(y,x-1)-imgin(y+1,x)-imgin(y-1,x);
                    e = e+1;
                 %else must be a inner mask pixel   
                else b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x);
                    e = e+1;
                end
            end
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
    
    %Looping through all pixels to determine which ones are on mask
    %boundary
    for y = 2:imh-1
        for x = 2:imw-1
            
            %Checking if boundary pixel has just one out of mask neighbour
            %in x direction
            if ismember(y,row_border_1_x) && ismember(x,column_border_1_x) && ismember(V(y,x),pixel_border_1_x)
                
                %if out of mask neigbour is to the left
                if mask_s(y,x) ~= mask_s(y,x-1)
                    for der = -1:1
                        if der == -1
                            i(count) = e;
                            v(count) = -1;
                            j(count) = V(y+der,x);
                            count = count +1;
            
                        elseif der == 0 
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
                    %add constant value to right side of equation,
                    %therefore add it to b
                    b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y,x-1);
                    e = e+1;
                
                %if out of mask neighbour is to the right    
                else
                    for der = -1:1
                        if der == 1
                            i(count) = e;
                            v(count) = -1;
                            j(count) = V(y+der,x);
                            count = count +1;
            
                        elseif der == 0 
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
                    b(e)= 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y,x+1);
                    e = e+1;
                end
                
            %Checking if boundary pixel has just one out of mask neighbour
            %in y direction    
            elseif ismember(y,row_border_1_y) && ismember(x,column_border_1_y) && ismember(V(y,x),pixel_border_1_y)
                if mask_s(y,x) ~= mask_s(y-1,x)
                    for der = -1:1
                        if der == -1
                            i(count) = e;
                            v(count) = -1;
                            j(count) = V(y,x+der);
                            count = count +1;
            
                        elseif der == 0 
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
                    b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y-1,x);
                    e = e+1;
                    
                else
                    for der = -1:1
                        if der == 1
                            i(count) = e;
                            v(count) = -1;
                            j(count) = V(y,x+der);
                            count = count +1;
            
                        elseif der == 0 
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
                    b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y+1,x);
                    e = e+1;
                end
                
            %Checking if boundary pixel has two out of mask neighbour   
            elseif ismember(V(y,x),pixel_border_2)
                
                %if out of mask neigbours are in x+1,y+1 direction
                if mask_s(y+1,x) == 0 && mask_s(y,x+1) == 0
                    for der = -1:0           
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
                    
                    %Add two constants to b
                    b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y+1,x)+imgin(y,x+1);
                    e = e+1;
                
                %if out of mask neigbours are in x-1,y-1 direction
                elseif mask_s(y-1,x) == 0 && mask_s(y,x-1) == 0
                    for der = 0:1            
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
                    b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y-1,x)+imgin(y,x-1);
                    e = e+1;
                
                %if out of mask neigbours are in x+1,y-1 direction    
                elseif mask_s(y-1,x) == 0 && mask_s(y,x+1) == 0
                    for der = -1:1            
                        if der == 0 
                            i(count) = e;
                            v(count) = 4;
                            j(count) = V(y,x);
                            count = count +1;
                            
                        elseif der == -1
                            v(count) = -1;
                            i(count) = e;
                            j(count) = V(y,x+der);
                            count = count+1;
                        
                        else
                            v(count) = -1;
                            i(count) = e;
                            j(count) = V(y+der,x);
                            count = count+1;           
                        end            
                    end
                    b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y-1,x)+imgin(y,x+1);
                    e = e+1;
                
                %if out of mask neigbours are in x-1,y+1 direction    
                elseif mask_s(y+1,x) == 0 && mask_s(y,x-1) == 0
                    for der = -1:1            
                        if der == 0 
                            i(count) = e;
                            v(count) = 4;
                            j(count) = V(y,x);
                            count = count +1;
                            
                        elseif der == 1
                            v(count) = -1;
                            i(count) = e;
                            j(count) = V(y,x+der);
                            count = count+1;
                        
                        else
                            v(count) = -1;
                            i(count) = e;
                            j(count) = V(y+der,x);
                            count = count+1;           
                        end
                    end
                    b(e) = 4*img_s_chan(y,x)-img_s_chan(y,x+1)-img_s_chan(y,x-1)-img_s_chan(y+1,x)-img_s_chan(y-1,x)+imgin(y+1,x)+imgin(y,x-1);
                    e = e+1;
                end
            
            %If 3 neigbouring pixels are out of mask-then set value to a constant equivalent
            %to pixel value in target image(too much computation otherwise)
            elseif ismember(V(y,x),pixel_border_3)
                i(count) = e;
                j(count) = V(y,x);
                v(count) = 1;
                b(e) = imgin(y,x);
                count = count+1;
                e = e+1;
            end
        end
    end
                                                                          
%Initialize control points
    i(count) = e;
    j(count) = 1;
    v(count) = 1;
    b(e) = imgin(1,1);

    i(count+1) = e+1;
    j(count+1) = V(imh,1);
    v(count+1) = 1;
    b(e+1) = imgin(imh,1);

    i(count+2) = e+2;
    j(count+2) = V(1,imw);
    v(count+2) = 1;
    b(e+2) = imgin(1,imw);

    i(count+3) = e+3;
    j(count+3) = V(imh,imw);
    v(count+3) = 1;
    b(e+3) = imgin(imh,imw);

    A = sparse(i,j,v);
    
    %solve equation for particular channel at a time
    solution = A\b;
    imgout(:,:,channel) = reshape(solution,[imh,imw]);
    error = sum(abs(A*solution-b));
    disp(error);
end
end

