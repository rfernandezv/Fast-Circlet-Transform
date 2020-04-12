function pts = ctist(Cl,r0,num,sens)
%%% CTIST implements iterative soft-thresholding for circlet transform for
%%% finding more circle in a given CircleLet (Cl) cell array as input. It
%%% iteratively looks for maximum coefficients with a soft-thresholding
%%% technique. The process is repeated 'num' times and the parameter 'sens'
%%% controls the sensitivity of the method for searching. 'sens' parameter
%%% is usefule while searching for the overlapped circles.
%       Inputs:
%           Cl: Circlelet coefficients given by the fdct2 function. 
%           r0: the radii ranges of the desired circular objects. 
%           num: the number of the desired circular objects that must be
%                localized.
%           sens: sensivity paramets which is utilized for localizing
%                 overlapped cicles. 'sens' should be set between 0 and 1.
%                 For higher values of the 'sens', much overlapped circles
%                 will be missed. For localizing non-overlapped objects the
%                 value of 1 is suitable for 'sens'. Default value is 1.
%       Output:
%           pts: m-by-3 matrix contains the locations and radii of the
%           center of detected circular objects. 'm' equals parametr 'num'.
%
%  written by Omid Sarrafzadeh,
%  Isfahan University of Medical Sciences, Isfahan, Iran
%  Email: o.sarrafzade@gmail.com
%
% If you use the code provided here, please cite the following paper:
% O. Sarrafzadeh, A. Mehri, H. Rabbani, N. Ghane, A. Talebi, "Circlet based
% framework for red blood cells segmentation and counting", in Proc. IEEE
% Workshop on Signal Processing Systems,
% Hangzhou, China, Oct. 14-16, 2015.
%
%  Reference for CT:
%  H. Chauris, I. Karoui, P. Garreau, H. Wackernagel, P. Craneguy, and L.
%  Bertino, "The circlet transform: A robust tool for detecting features
%  with circular shapes," Computers & Geosciences, vol. 37, pp. 331-342,
%  2011.

if nargin<=3
    sens = 1;
end
[row,col] = size(Cl{1,1});
M = zeros(size(Cl,1),size(Cl,2));
pts = zeros(num,3);
for n=1:num
    for x=1:size(Cl,1)
        for y=1:size(Cl,2)
            M(x,y) = max(max(Cl{x,y}));
        end
    end
    [r, c] = find(M == max(max(M)));
    [x0, y0] = find(Cl{r,c} == M(r,c));
    for x=1:size(Cl,1)
        Circ = circle(sens*r0(x),x0,y0,row,col);
        for y=1:size(Cl,2)
            Cl{x,y} = Cl{x,y}.*Circ;
        end
    end
    pts(n,:) = [y0 x0 r0(r)];
end