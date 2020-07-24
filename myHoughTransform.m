function [H, L, res] = myHoughTransform(img_binary, Drho, Dtheta, n)
% function [H, L, res] = myHoughTransform(img_binary, Drho, Dtheta, n)
% MYHOUGHTRANSFORM implements the Hough transform of a binary image in the
% rho theta plain, given specific steps for rho (Drho) and theta (Dtheta).
% Also it calculates the n strongest lines indicated by maximum of the
% accumulator matrix. Lastly in computes the number of pixels that do not
% belong in any of the n strongest lines
%
% INPUT:
% img_binary    : image as threshhold output of an edge detector
% Drho          : rho resolution (pixels)
% Dtheta        : theta resolution (degrees)
% n             : how many lines to keep
% 
% OUTPUT:
% H             : Hough's transform matrix (accumulator)
% L             : rho and theta of n strongest lines 
% res           : number of image's points that do not belong in the lines
%
% Author        : Nikolaos Giakoumoglou AEM: 9043
% Date          : 20/05/2020

%% Calculate Hough Matrix H aka Accumulator 
img_binary=double(img_binary); % convert binary to double
img_binary = flipud(img_binary);
img_binary = imrotate(img_binary,-90); % the following computes the hough transform of the rotated image ???
[N1, N2] = size(img_binary); 
theta=linspace(-90,0,ceil(90/Dtheta)+1); % -90 to 0 (degrees)
theta=[theta -fliplr(theta(2:end-1))]; % theta from -90 to 90 using flip
ntheta=length(theta);
rhomax=sqrt((N1-1)^2+(N2-1)^2); % rho max
rhomax=ceil(rhomax/Drho); % round rho max 
nrho=2*rhomax-1; % length(rho)
rho=linspace(-rhomax*Drho,rhomax*Drho,nrho); % rho from -rhomax to rhomax with step Drho
[x,y,val]=find(img_binary); % find non negative values
%x=x-1;y=y-1;
H=zeros(nrho,ntheta); % initialize Hough matrix
for k=1:ceil(length(val)/1000) % make loop parallel run from first to last ~1000 elements or less in last loop
    first=(k-1)*1000+1;
    last=min(first+999,length(x));
    x_mat=repmat(x(first:last),1,ntheta); % replicate x ntheta-times in columns (1000)x(ntheta)
    y_mat=repmat(y(first:last),1,ntheta); % replicate y ntheta-times in columns (1000)x(ntheta)
    val_mat=repmat(val(first:last),1,ntheta); % replicate val ntheta-times in columns (1000)x(ntheta)
    theta_mat=repmat(theta,size(x_mat,1),1)*pi/180; % replicate theta ntheta-times in columns (1000)x(ntheta)
    rho_mat=x_mat.*cos(theta_mat)+y_mat.*sin(theta_mat); % calculate rho (1000)x(ntheta)
    slope=(nrho-1)/(rho(end)-rho(1)); % ~1, otherwise array out of bounds
    rho_bin_index=round(slope*(rho_mat-rho(1))+1); % rho intervals (1000)x(ntheta)
    theta_bin_index=repmat(1:ntheta,size(x_mat,1),1); % theta intervals (1000)x(ntheta)
    H=H+full(sparse(rho_bin_index(:),theta_bin_index(:),val_mat(:),nrho,ntheta)); % update Hough matrix
end
%% Calculate L
vectorH = H(:);
vectorH = sort(vectorH);
maxH = vectorH(end-n+1:end);
L = NaN*ones(n,2);
for i=1:n
    [rowH, colH] = find(H==maxH(i));
    L(i,1) = rho(rowH(1));
    L(i,2) = theta(colH(1));
end
%% Calculate res
res = N1*N2;
x0 = 1;
xend = N2;
for i=1:n
    rh = L(i,1); th = L(i,2);
    if th == 0
        res = res - N2;
    else
        y0 = (-cosd(th)/sind(th))*x0 + (rh / sind(th));
        yend = (-cosd(th)/sind(th))*xend + (rh / sind(th));
        linelength = round(sqrt((xend-x0)^2 + (yend-y0)^2));
        res = res - linelength;
    end
end

end

%{
rhomax=round(sqrt(N1^2+N2^2));
acc=zeros(rhomax,180);

for n1=1:N1
    for n2=1:N2
        if(inputimage(n1,n2)==0)
            for d=1:180
            rho=round(n1*cos((d*pi)/180)+n2*sin(d*pi)/180);
            if(rho<rhomax & rho>0) acc(rho,d)=acc(rho,d)+1; end
            end
        end
    end
end
%}