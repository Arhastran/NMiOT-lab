%% Adam Ignaciuk 
%% Circle `detection

folder_name = '/Volumes/DISCWORLD/methlab/Lab 5/';
data1 = imread([folder_name,'image.png']);
%gray = double(rgb2gray(data1));

Bin = imbinarize(data1,0.27);

xd = fspecial("log",5,0.40);
 
C1= conv2(Bin,xd, 'same');

Bin2 = imbinarize(C1);
%figure(Color = 'w');
%imshow(C1)

accumulator = zeros(size(C1));

% BW = edge(C1,'canny');
% [H,T,R] = hough(C1,'RhoResolution',0.5,'Theta',-90:0.5:89);
% imshow(imadjust(rescale(H)),'XData',T,'YData',R)
Si = size(C1);
for r =18:1:25
for i=1+r:Si(1)-r
    for j = 1+r:Si(2)-r
        if C1(i,j) >= 0
            for k = -90:1:90
                a = round(i-r*cos(deg2rad(k)));
                b = round(j-r*sin(deg2rad(k)));
                accumulator(a,b) =+ 1;
            end
        end
    end
end
end

% for x = (1:1:Si(1))
%     for y = (1:1:Si(2))
%         if accumulator(x,y) < = accumulator()
%         end
%        
%     end
% end



figure(Color = 'w');
imshow(accumulator)