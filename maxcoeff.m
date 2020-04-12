function pts = maxcoeff(Cl, k)
%%% MAXCOEFF determines the maximum circlelet coefficiesnts in a given
%%% CircleLet (Cl) cell array as input. It looks for maximum coefficient
%%% corresponding to each radius and a predifined frequency, k.
%       Inputs:
%           Cl: Circlelet coefficients given by the fdct2 function.
%           k: predifined frequency to determine maximum coefficiesnts.
%       Output:
%           pts: m-by-2 matrix contains the location of the center of
%           detected circles corresponding to each radius. 'm' equals the
%           number of desired radii, size(Cl,1).
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

Clk = Cl(:,k);
pts = zeros(size(Clk,1), 2);
for r=1:size(Clk,1)
    M = max(max(Clk{r,1}));
    [pts(r,2), pts(r,1)] = find(Clk{r,1} == M);
end