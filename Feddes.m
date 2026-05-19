function Result = Feddes(h)
% clear all;
% close all;
% clc;
% h = 0:-1:-20000;

%Assigning values of feddes function, these defines the 
%*********************************
%           CORN
%*********************************
ha=-15;
hb=-30;
hc=-600;
hd=-8000;

% Making a vector of sink term
%Result = ones(1,length(h));
%finding index in vector h with different conditions
%We are trying to find the nodes where the pressure head
%is defined by certain set of conditions as follows
idx1 = find(h>=ha);
idx2 = find(h>=hb & h<ha);
idx3 = find(h>=hc & h<hb);
idx4 = find(h>=hd & hc>h);
idx5 = find(h<hd);  

%Assigning those nodes with the value of Alpha, the feddes funciton
Result(idx1) = 0;
Result(idx2) = (h(idx2)-ha)/(hb-ha);
Result(idx3) = 1;
Result(idx4) = (hd-h(idx4))/(hd-hc);
Result(idx5) = 0;
% 
% dd=Result;
% %plot(h,Result)