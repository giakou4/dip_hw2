classdef myFunctions
    methods(Static)
%% This function plots the houghlines of Hough Transform given a BW image
        function plotLines(BW,L)
            figure
            imshow(BW); 
            hold on;    
            x0 = 1;
            xend = size(BW,2);
            for idx = 1 : size(L,1)
                r = L(idx,1); th = L(idx,2); % get rho and theta
                if (th == 0) % vertical line
                    line([r r], [1 size(BW,1)], 'Color', 'green','LineWidth',2);
                else
                    y0 = (-cosd(th)/sind(th))*x0 + (r / sind(th)); % theta in degrees
                    yend =(-cosd(th)/sind(th))*xend + (r / sind(th));
                    line([x0 xend], [y0 yend], 'Color', 'green','LineWidth',2);
                end
            end
        end
%% This function plots the Hough Transform in rho-theta plain
        function plotHough(H,L,xdata,ydata)
            figure
            imshow(imadjust(rescale(H)),'XData',xdata,'YData',ydata,'InitialMagnification','fit');
            title('Hough Transform');
            xlabel('\theta'), ylabel('\rho');
            axis on, axis normal, hold on;
            colormap('gray');
            hold on;
            plot(L(:,2),L(:,1),'s','color','green','MarkerSize',10);
        end
%% This functions is like myHoughTransform but call hough (build-in)
        function L = houghTransform(img_binary, Drho, Dtheta, n)
        BW = double(img_binary);
        [H,theta,rho] = hough(BW,'RhoResolution',Drho,'Theta',-90:Dtheta:89);
        P = houghpeaks(H,n);
        lines = houghlines(BW,theta,rho,P);
        L = NaN*ones(length(lines),2);
        L(:,1) = [lines.rho];
        L(:,2) = [lines.theta];
        end
%% This is  a custom function to find edges where nearby pixels are white
        function [edges] = findEdges(BW)
           % BW is binary image, output of edge detector
           % background is white = 1
           edges = NaN*ones(4,2);
           no = 1;
           BW = double(BW);
           [N1, N2] = size(BW);
           for n1=4:N1-4
               for n2=4:N2-4
                   if BW(n1,n2)==0 %check if edge
                       if BW(n1+1,n2)==1 && BW(n1+2,n2)==1 & BW(n1+3,n2)==1 && BW(n1,n2-1) && BW(n1,n2-1) && BW(n1,n2-1) % right upper
                                    numbwhite = 1;
                                    for nn1=(n1-2:n1+2)
                                        for nn2=(n2-2):(n2+2)
                                            if BW(nn1,nn2)==1
                                                numwhite = numbwhite + 1;
                                            end
                                        end
                                    end
                                    if numwhite > 18
                                    else
                                        edges(no,1)=n1; edges(no,2)=n2; no = no + 1;
                                    end
                       elseif BW(n1-1,n2)==1 && BW(n1-2,n2)==1 & BW(n1-3,n2)==1 && BW(n1,n2-1) && BW(n1,n2-1) && BW(n1,n2-1) % left upper
                                    numbwhite = 1;
                                    for nn1=(n1-2:n1+2)
                                        for nn2=(n2-2):(n2+2)
                                            if BW(nn1,nn2)==1
                                                numwhite = numbwhite + 1;
                                            end
                                        end
                                    end
                                    if numwhite > 18
                                    else
                                        edges(no,1)=n1; edges(no,2)=n2; no = no + 1;
                                    end
                       elseif BW(n1-1,n2)==1 && BW(n1-2,n2)==1 & BW(n1-3,n2)==1 && BW(n1,n2+1) && BW(n1,n2+1) && BW(n1,n2+1) % left down
                                    numbwhite = 1;
                                    for nn1=(n1-2:n1+2)
                                        for nn2=(n2-2):(n2+2)
                                            if BW(nn1,nn2)==1
                                                numwhite = numbwhite + 1;
                                            end
                                        end
                                    end
                                    if numwhite > 18
                                    else
                                        edges(no,1)=n1; edges(no,2)=n2; no = no + 1;
                                    end
                       elseif BW(n1+1,n2)==1 && BW(n1+2,n2)==1 & BW(n1+3,n2)==1 && BW(n1,n2+1) && BW(n1,n2+1) && BW(n1,n2+1) % right down
                                    numbwhite = 1;
                                    for nn1=(n1-2:n1+2)
                                        for nn2=(n2-2):(n2+2)
                                            if BW(nn1,nn2)==1
                                                numwhite = numbwhite + 1;
                                            end
                                        end
                                    end
                                    if numwhite > 18
                                    else
                                        edges(no,1)=n1; edges(no,2)=n2; no = no + 1;
                                    end
                       end
                   end
               end
           end  

        end
%% This auxiliary function removes total black rows and column from an image I     
        function [I2] = removeBlackLines(I)
            [N1,N2] = size(I);
           for n1=1:N1
              res(n1) = length(find(mean(abs(I(n1,:)))==0)); % all row black => res(n1) = 1
           end
           I2 = I;
           res = find(res==1);
           I2(res,:) = [];
           
           [N1,N2] = size(I2);
           res = -1;
           for n2=1:N2
              res(n2) = length(find(mean(abs(I(:,n2)))==0)); % all column black => res(n1) = 1
           end
           res = find(res==1);
           I2(:,res) = [];            
        end      
        %% This auxiliary function removes white rows and column from an image I  >240
        function [I2] = removeWhiteLines(I)
           threshold = 240;
           [N1,N2] = size(I);
           for n1=1:N1
              res(n1) = length(find(mean(abs(I(n1,:)))>threshold)); % all row white => res(n1) = 1
           end
           I2 = I;
           res = find(res==1);
           I2(res,:) = [];
           
           [N1,N2] = size(I2);
           res = -1;
           for n2=1:N2
              res(n2) = length(find(mean(abs(I(:,n2)))>threshold)); % all row white => res(n1) = 1
           end
           res = find(res==1);
           I2(:,res) = [];            
        end
%% This auxiliary function makes the black background from a rotated image white
        function [I] = makeBackgroundWhite(I)
            [N1,N2] = size(I);
            for n1=1:N1
                for n2=1:N2
                    if I(n1,n2) == 0
                        I(n1,n2) = 255;
                    end
                end
            end
        end
%% This auxiliary function removes hough lines on borders (up, down, left, right) 
        function lines = removeBadLines(BW, lines)
            if length(lines) == 0 return; end
            k = [];
            for i=1:size(lines,1)
                rho = lines(i,1); theta = lines(i,2);
                if (abs(rho)>0.88*size(BW,2) || abs(rho)<0.12*size(BW,2) ) %&& theta == 0 
                    k = [k i];
                elseif( abs(rho) > 0.88*size(BW,1) || abs(rho) < 0.12*size(BW,1) )  %&&( theta == 90 | theta == -90 ) 
                    k = [k i];
                end
            end
            lines(k,:) = [];
        end             
%% This axiliary function removes duplicate hough lines        
        function lines = removeDuplicateLines(lines)
            if length(lines)==0 return; end
            if length(lines)==1 return; end
            k = [];
            for i=1:size(lines)
                for j=(i+1):size(lines)
                    if lines(i,1)==lines(j,1) && lines(i,2)==lines(j,2)
                        k = [k i];
                    end
                end
            end
            lines(k,:) = [];            
        end
%{        
%% Our solver, first call to find the bests hough lines = no of sub-images
        function [IMAGES_FOUND,NO_IMAGES_FOUND] = recursiveFinder1(BW,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND)
            BW = myFunctions.removeWhiteLines(BW);
            BW = myFunctions.removeBlackLines(BW);
            %BWcunny = edge(imbinarize(BW,0.9),'canny');
            BWcunny = edge(BW,'canny');
            lines = myFunctions.houghTransform(BWcunny, Drho, Dtheta, n);
            myFunctions.plotLines(BW,lines); title('Image Input with Hough lines');
            lines = myFunctions.removeDuplicateLines(lines);
            lines = myFunctions.removeBadLines(BW,lines);
            myFunctions.plotLines(BW,lines);title('Image Input with "good" Hough lines - indicate number of images');
            
            numLines = size(lines,1); 
            
            if (numLines == 0) 
                display('CANNOT FIND ANY SUB IMAGES, SORRY!');
                return;
            end
  
            for i=1:1:numLines % for each BASIC LINE
                    rho = lines(i,1);
                    theta = lines(i,2);
                    if ( theta == 90 | theta == -90 ) %&& ( abs(rho) < 0.95*size(BW,1) && abs(rho) > 0.05*size(BW,1) )
                    subImage1 = imcrop(BW,[0 0 size(BW,2) abs(rho)]); %figure; imshow(subImage1); title('Sub Image 1');
                    subImage2 = imcrop(BW,[0 abs(rho) size(BW,2) size(BW,1)-abs(rho)]);% figure; imshow(subImage2); title('Sub image 2');
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder2(subImage1,Drho,Dtheta,200,IMAGES_FOUND,NO_IMAGES_FOUND);
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder2(subImage2,Drho,Dtheta,200,IMAGES_FOUND,NO_IMAGES_FOUND);           
                    end
                    
                    if theta == 0 %&& (abs(rho)<0.95*size(BW,2) && abs(rho)>0.05*size(BW,2) )
                        subImage1 = imcrop(BW,[0 0 abs(rho) size(BW,2)]);%figure; imshow(subImage1); title('Sub Image 1');
                        subImage2 = imcrop(BW,[abs(rho) 0  size(BW,1)-abs(rho) size(BW,2)]); %figure; imshow(subImage2); title('Sub image 2');
                       [IMAGES_FOUND,NO_IMAGES_FOUND] =  myFunctions.recursiveFinder2(subImage1,Drho,Dtheta,200,IMAGES_FOUND,NO_IMAGES_FOUND);
                       [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder2(subImage2,Drho,Dtheta,200,IMAGES_FOUND,NO_IMAGES_FOUND);              
                    end    
                    
                    if theta~= 90 && theta ~=-90 && theta ~= 0
                        BW=imrotate(BW,90+theta);
                        disp(['ROTATING IMAGE ...']);
                        BW = myFunctions.makeBackgroundWhite(BW);
                        [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder1(BW,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);
                    end                        
             end 
            

        end
%% Our solver, second call that cleans the subImage and finds sub-sub-images lel       
        function [IMAGES_FOUND,NO_IMAGES_FOUND] = recursiveFinder2(BW,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND)
            if size(BW,1)<5 || size(BW,2) <5 return; end % small image, garbage
            if mean(BW(:))>230 return; end % white image, garbage
            
            BW = myFunctions.removeWhiteLines(BW);
            BWcrop = imbinarize(BW,0.93);
            BWcunny = edge(BWcrop,'canny');
            lines = myFunctions.houghTransform(BWcunny, Drho, Dtheta, n);
            %myFunctions.plotLines(BW,lines);
            lines = myFunctions.removeDuplicateLines(lines);
            lines = myFunctions.removeBadLines(BW,lines);
            
             numLines = size(lines,1);
             %if (numLines ==0) return; end
             counter = numLines;
             %myFunctions.plotLines(BW,lines);
             for i=1:1:numLines % for each line, seperate the image to 2 sub images
                    rho = lines(i,1);
                    theta = lines(i,2);
                    %myFunctions.plotLines(BW,[rho theta; 0 0]); title('Line seperation')
                    if ( theta == 90 | theta == -90 ) %%&& ( abs(rho) < 0.95*size(BW,1) && abs(rho) > 0.05*size(BW,1) )
                    counter = counter - 1;
                    subImage1 = imcrop(BW,[0 0 size(BW,2) abs(rho)]); %figure; imshow(subImage1); title('Sub Image 1');
                    subImage2 = imcrop(BW,[0 abs(rho) size(BW,2) size(BW,1)-abs(rho)]); %figure; imshow(subImage2); title('Sub image 2');
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder2(subImage1,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder2(subImage2,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);           
                    end
                    
                    if theta == 0 && (abs(rho)<0.95*size(BW,2) && abs(rho)>0.05*size(BW,2) )
                        counter = counter - 1;
                        subImage1 = imcrop(BW,[0 0 abs(rho) size(BW,2)]);%figure; imshow(subImage1); title('Sub Image 1');
                        subImage2 = imcrop(BW,[abs(rho) 0  size(BW,1)-abs(rho) size(BW,2)]); %figure; imshow(subImage2); title('Sub image 2');
                       [IMAGES_FOUND,NO_IMAGES_FOUND] =  myFunctions.recursiveFinder2(subImage1,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);
                       [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder2(subImage2,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);              
                    end    
             end
             % cannot crop further
             if counter == numLines % did not enter any loop so all lines are around the image - useless
                %figure; imshow(BW); title('IMAGE FOUND!')
                NO_IMAGES_FOUND = NO_IMAGES_FOUND + 1;
                for n1=1:size(BW,1)
                    for n2=1:size(BW,2)          
                        IMAGES_FOUND(n1,n2,NO_IMAGES_FOUND) = BW(n1,n2);
                    end
                end
             end
         
        end
%}        
%% Rotates image for a last time 
function I2 = rotImageFinal(I)
   %figure;   imshow(I/255); title('Before rotation')
   BWcanny = edge(I,'canny');
   L = myFunctions.houghTransform(BWcanny,1,1,4);
   %myFunctions.plotLines(I/255,L);
   theta90pos = 0;
   theta90neg = 0;
   theta0 = 0;
   for i=1:length(L)
       if (L(i,2) == 90 && (L(i,1) > 0.90*norm(size(I)) || L(i,1) < 0.10*norm(size(I)) ))
           theta90pos = theta90pos + 1;
       elseif ( L(i,2) == -90 && (  L(i,1) > 0.90*norm(size(I)) || L(i,1) < 0.10*norm(size(I))   )   )
           theta90neg = theta90neg + 1;
       else %if (L(i,2) == 0 && (L(i,1) > 0.90*norm(size(I)) || L(i,1) < 0.10*norm(size(I)) ))
           theta0 = theta0 + 1;
       end
   end
   if theta90pos>theta90neg && theta90pos>theta0
       I2 = imrotate(I,90);
   elseif theta90neg>theta90pos && theta90neg>theta0
       I2 = imrotate(I,-90);
   else
       I2 = I;
   end
   %figure, imshow(I2/255); title('After rotation');
end
%% Our solver vol4 using myHoughTransform
function [IMAGES_FOUND,NO_IMAGES_FOUND] = recursiveFinder4(BW,edgeFlag,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND)
            if n==0 return; end
            
            % white image, garbage 
            whiteGarbageThreshold = 230;
            if mean(BW(:))>whiteGarbageThreshold return; end 
            
            %almost no white background
            whiteThreshold = 0.1;
            binarizeThreshold = 0.85;
            BWbinarized = imbinarize(BW,binarizeThreshold);
            if (sum(BWbinarized(:))/(size(BW,1)*size(BW,2)))<whiteThreshold 
                NO_IMAGES_FOUND = NO_IMAGES_FOUND + 1;
                for n1=1:size(BW,1)
                    for n2=1:size(BW,2)          
                        IMAGES_FOUND(n1,n2,NO_IMAGES_FOUND) = BW(n1,n2);
                    end
                end
                return;
            end
            
            BW = myFunctions.removeWhiteLines(BW);
            BW = myFunctions.removeBlackLines(BW);
            
            % try to find lines harder - one time use
            if edgeFlag == 1
                BWcunny = edge(imbinarize(BW,0.9),'canny');
            else
            	BWcunny = edge(BW,'canny');
            end
            
            [~,lines,~] = myHoughTransform(BWcunny, Drho, Dtheta, n);
            %myFunctions.plotLines(BW,lines(1,:)); title('Hough lines');
            lines = myFunctions.removeDuplicateLines(lines);
            lines = myFunctions.removeBadLines(BW,lines);
            %myFunctions.plotLines(BW,lines);title('"Good" Hough lines');
            
            numLines = size(lines,1); 
            isEnd = numLines; % see later the use of counter
            for i=1:1:numLines % for each BASIC LINE
                    rho = lines(i,1);
                    theta = lines(i,2);
                    if ( theta == 90 | theta == -90 ) %&& ( abs(rho) < 0.95*size(BW,1) && abs(rho) > 0.05*size(BW,1) )
                    isEnd = isEnd - 1;
                    subImage1 = imcrop(BW,[0 0 size(BW,2) abs(rho)]); %figure; imshow(subImage1); title('Sub Image 1');
                    subImage2 = imcrop(BW,[0 abs(rho) size(BW,2) size(BW,1)-abs(rho)]);% figure; imshow(subImage2); title('Sub image 2');
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder4(subImage1,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder4(subImage2,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);           
                    end
                    
                    if theta == 0 %&& (abs(rho)<0.95*size(BW,2) && abs(rho)>0.05*size(BW,2) )
                        isEnd = isEnd - 1;
                        subImage1 = imcrop(BW,[0 0 abs(rho) size(BW,2)]);%figure; imshow(subImage1); title('Sub Image 1');
                        subImage2 = imcrop(BW,[abs(rho) 0  size(BW,1)-abs(rho) size(BW,2)]); %figure; imshow(subImage2); title('Sub image 2');
                       [IMAGES_FOUND,NO_IMAGES_FOUND] =  myFunctions.recursiveFinder4(subImage1,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);
                       [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder4(subImage2,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);              
                    end    
                    
                    if theta~= 90 && theta ~=-90 && theta ~= 0
                        isEnd = isEnd - 1;
                        BW=imrotate(BW,90+theta);
                        disp(['ROTATING IMAGE ...']);
                        BW = myFunctions.makeBackgroundWhite(BW);
                        [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder4(BW,0,Drho,Dtheta,n-1,IMAGES_FOUND,NO_IMAGES_FOUND);
                    end                        
             end 
                         % cannot crop further
             if isEnd == numLines % did not enter any loop so all lines are around the image - useless, we have our sub image
                %figure; imshow(BW); title('IMAGE FOUND!')
                NO_IMAGES_FOUND = NO_IMAGES_FOUND + 1;
                for n1=1:size(BW,1)
                    for n2=1:size(BW,2)          
                        IMAGES_FOUND(n1,n2,NO_IMAGES_FOUND) = BW(n1,n2);
                    end
                end
             end

end
%% Our solver vol3
function [IMAGES_FOUND,NO_IMAGES_FOUND] = recursiveFinder3(BW,edgeFlag,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND)
            if n==0 return; end
            
            % white image, garbage 
            whiteGarbageThreshold = 230;
            if mean(BW(:))>whiteGarbageThreshold return; end 
            
            %almost no white background
            whiteThreshold = 0.1;
            binarizeThreshold = 0.85;
            BWbinarized = imbinarize(BW,binarizeThreshold);
            if (sum(BWbinarized(:))/(size(BW,1)*size(BW,2)))<whiteThreshold 
                NO_IMAGES_FOUND = NO_IMAGES_FOUND + 1;
                for n1=1:size(BW,1)
                    for n2=1:size(BW,2)          
                        IMAGES_FOUND(n1,n2,NO_IMAGES_FOUND) = BW(n1,n2);
                    end
                end
                return;
            end
            
            BW = myFunctions.removeWhiteLines(BW);
            BW = myFunctions.removeBlackLines(BW);
            
            % try to find lines harder - one time use
            if edgeFlag == 1
                BWcunny = edge(imbinarize(BW,0.9),'canny');
            else
            	BWcunny = edge(BW,'canny');
            end
            
            lines = myFunctions.houghTransform(BWcunny, Drho, Dtheta, n);
            %myFunctions.plotLines(BW,lines(1,:)); title('Hough lines');
            lines = myFunctions.removeDuplicateLines(lines);
            lines = myFunctions.removeBadLines(BW,lines);
            %myFunctions.plotLines(BW,lines);title('"Good" Hough lines');
            
            numLines = size(lines,1); 
            isEnd = numLines; % see later the use of counter
            for i=1:1:numLines % for each BASIC LINE
                    rho = lines(i,1);
                    theta = lines(i,2);
                    if ( theta == 90 | theta == -90 ) %&& ( abs(rho) < 0.95*size(BW,1) && abs(rho) > 0.05*size(BW,1) )
                    isEnd = isEnd - 1;
                    subImage1 = imcrop(BW,[0 0 size(BW,2) abs(rho)]); %figure; imshow(subImage1); title('Sub Image 1');
                    subImage2 = imcrop(BW,[0 abs(rho) size(BW,2) size(BW,1)-abs(rho)]);% figure; imshow(subImage2); title('Sub image 2');
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder3(subImage1,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);
                    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder3(subImage2,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);           
                    end
                    
                    if theta == 0 %&& (abs(rho)<0.95*size(BW,2) && abs(rho)>0.05*size(BW,2) )
                        isEnd = isEnd - 1;
                        subImage1 = imcrop(BW,[0 0 abs(rho) size(BW,2)]);%figure; imshow(subImage1); title('Sub Image 1');
                        subImage2 = imcrop(BW,[abs(rho) 0  size(BW,1)-abs(rho) size(BW,2)]); %figure; imshow(subImage2); title('Sub image 2');
                       [IMAGES_FOUND,NO_IMAGES_FOUND] =  myFunctions.recursiveFinder3(subImage1,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);
                       [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder3(subImage2,0,Drho,Dtheta,1,IMAGES_FOUND,NO_IMAGES_FOUND);              
                    end    
                    
                    if theta~= 90 && theta ~=-90 && theta ~= 0
                        isEnd = isEnd - 1;
                        BW=imrotate(BW,90+theta);
                        disp(['ROTATING IMAGE ...']);
                        BW = myFunctions.makeBackgroundWhite(BW);
                        [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder3(BW,0,Drho,Dtheta,n-1,IMAGES_FOUND,NO_IMAGES_FOUND);
                    end                        
             end 
                         % cannot crop further
             if isEnd == numLines % did not enter any loop so all lines are around the image - useless, we have our sub image
                %figure; imshow(BW); title('IMAGE FOUND!')
                NO_IMAGES_FOUND = NO_IMAGES_FOUND + 1;
                for n1=1:size(BW,1)
                    for n2=1:size(BW,2)          
                        IMAGES_FOUND(n1,n2,NO_IMAGES_FOUND) = BW(n1,n2);
                    end
                end
             end

end

    end
end

