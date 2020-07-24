function rotImg = myImgRotation(img ,angle)
% function rotImg = myImgRotation(img ,angle)
% MYIMGROTATION take an image img as an input and rotates by angle degrees
% by the assumption that the background is black (0)
%
% INPUT:
% img       : 2D grayscale or RGB image
% angle     : the degrees that the image will be rotated anti-clockwise
% OUTPUT:
% rotImg    : rotated image by angle degrees
%
% Author    : Nikolaos Giakoumoglou AEM: 9043
% Date      : 20/05/2020

% img dimensions
[N1,N2,N3]= size(img); 

% rotImg dimension
M1=ceil(N1*abs(cosd(angle))+N2*abs(sind(angle)));                      
M2=ceil(N1*abs(sind(angle))+N2*abs(cosd(angle)));                     

% rotated image initialized as black
rotImg=(zeros([M1 M2 N3]));

% original image
xo=ceil(N1/2);                                                            
yo=ceil(N2/2);

% center of final image
xfinal=ceil(M1/2);
yfinal=ceil(M2/2);

% main loop: calculate coordinates of rotated image
for m1=1:M1
    for m2=1:M2                                                       
         x = (m1-xfinal)*cosd(angle)+(m2-yfinal)*sind(angle);                                       
         y = -(m1-xfinal)*sind(angle)+(m2-yfinal)*cosd(angle);                             
         x = round(x)+xo;
         y = round(y)+yo;
         if (x>=1 && y>=1 && x<=size(img,1) &&  y<=size(img,2) ) 
              rotImg(m1,m2,:)=img(x,y,:);  
         end
    end
end

end

