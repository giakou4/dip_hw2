function corners = myDetectHarrisFeatures(I)
% function corners = myDetectHarrisFeatures(I)
% MYDETECTHARRISFEATURES computes the corners of an image I as a nx2 vector
% where n is the total number of corners indentified by specific parameters
% as k, threshold and 9x9 neighbourhood.
%
% INPUT:
% I         : 2D grayscale image in [0, 1]
% 
% OUTPUT:
% corners   : coordinates of the corners in nx2 matrix
%
% Author    : Nikolaos Giakoumoglou AEM: 9043
% Date      : 20/05/2020
k = 0.01;
threshold = 0.1;

[N1, N2] = size(I);

Gx = [-1 0 1; % applying Prewitt edge detector in the horizontal direction
      -1 0 1;
      -1 0 1];
Ix = filter2(Gx,I);
Gy = [1 1 1;% applying Prewitt edge detector in the vertical direction
      0 0 0;
    -1 -1 -1];
Iy = filter2(Gy,I);

Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

%applying gaussian filter on the computed value
h= fspecial('gaussian',[9 9],2); % 9x9 gaussian with sigma=2 arbitarily 
Ix2 = filter2(h,Ix2);
Iy2 = filter2(h,Iy2);
Ixy = filter2(h,Ixy);

R = zeros(N1,N2);
Rmax = -inf; 

for n1=1:N1
    for n2=1:N2
        M = [Ix2(n1,n2) Ixy(n1,n2);Ixy(n1,n2) Iy2(n1,n2)]; 
        R(n1,n2) = det(M)-k*(trace(M))^2;
        if R(n1,n2) > Rmax
            Rmax = R(n1,n2);
        end
    end
end

corners=zeros(3,2);
count = 1;
for n1 = 2:N1-1
    for n2 = 2:N2-1
        if R(n1,n2) > threshold*Rmax && R(n1,n2) > R(n1-1,n2-1) && ...
            R(n1,n2) > R(n1-1,n2) && R(n1,n2) > R(n1-1,n2+1) && ...
                R(n1,n2) > R(n1,n2-1) && R(n1,n2) > R(n1,n2+1) && ...
                R(n1,n2) > R(n1+1,n2-1) && R(n1,n2) > R(n1+1,n2) && ... 
                R(n1,n2) > R(n1+1,n2+1)
            corners(count,1) = n1;
            corners(count,2) = n2;
            count = count+1;
        end
    end
end

end

