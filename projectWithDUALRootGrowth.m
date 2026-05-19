
% with two-step implicit linearization to solve Richards' equation with
% initial GlenCondition h(z,0)=h-sub-b and BCs -K(dh/dz)+K=q-sub-a and dh/dz=0
% at z=0 and z=L, respectively. Hydraulic properties are calculated by
% calling the vectorized functions "GlenCond.m", "GlenSpecCap.m", and "GlenReten.m".

% Clear MATLAB Environment
clear all;
close all;
clc;

%****************************************************************
%               Soil domain and Time Information
%****************************************************************


m = 90;             % Number of Nodes
L = 180;            % Length of domain [cm]
dz = L/m;           % Space increment [cm]
z = dz*(1:m)'-dz/2; % Column vector of depths [cm]
dt = 0.1;           % Time increment [h]
endtime = 24*84;    % end of simulation[h]

% Vector of times at which pressure heads will be saved.

% -------->> "times" MUST BE MULTIPLE OF dt <<---------
%times = 24*(1:0.5:84);
times = 24*[0.5, 2, 5, 15, 30, 50, 84];
%times = 24*[ 0.5, 2.5, 5.5, 15.5, 30.5, 35.5, 50.5, 80];


%****************************************************************
%               ROOT INFORMAITON
%****************************************************************
Sm = -0.01/24;      % Rootwateruptake: maximum value


LoC = 10;            % Initial rooting depth[cm]
LmC = 120;            % Final rooting depth[cm]
tmaxC = 24*60;       % length of growing season (not the crop seson).
rC = (-2/tmaxC)*log(LoC/(LmC-LoC)); %Constant term for root growth

LoS = 10;            % Initial rooting depth[cm]
LmS = 70;            % Final rooting depth[cm]
tmaxS = 24*30;       % length of growing season (not the crop seson).
rS = (-2/tmaxS)*log(LoS/(LmS-LoS)); %Constant term for root growth

%***************************************************************
%           Initial condition and Flux information
%***************************************************************


qa = 0.5/24;        % Flux at the soil surface [cm/h]
%qa=0;
qb = 0;             % Flux at the bottom [cm/h]
ha = -100;          % Initial pressure head [cm]
h = ha*ones(m,1);   % Initial pressure head distributed in domain 
thetainit = GlenReten(h);% Similarly Initial water content

%*************************************************************
%             TRANSPIRATION DATA
%*************************************************************
T = dt:dt:endtime;
Trans = NaN(length(T),2);                % vector for saving data
Trans(:,1) = T;

path = 'C:\Agron916\Project\ET.xlsx';    % path to transpiration data 
Pdata = xlsread(path);                   % Transpiraiton in mm/hour
TpData = Pdata(:,1)/10;                  % storing transpiration data [cm/h]
Transpir = (TpData(1:(endtime/24)))';    % triming it to endtimes

%Repeating single value of transpiration for 24 times (full day)

Transpir  = repmat(Transpir, [24/dt 1]); %[h] 
Trans_hourly  = Transpir(:);
Trans(:,2) = Trans_hourly;

%*************************************************************
%             PREALLOCATING VECTORS
%*************************************************************

[a,b,c,hstar] = deal(NaN(m,1));  % h* and the coefficients a-sub-i, etc.
[S,Q,bxC,bxS, Sink]= deal(zeros(m,1));                   % vectors for sink terms

[Uptake,RWU,FeddesAlphaS, FeddesAlphaQ]= deal(zeros(m,1));    % Uptake for each run through while loop

% saving results at the times in "times"

[headresults,thetaresults,rho] = deal(NaN(m,length(times)));

[TranspirationC,TranspirationS] = deal(zeros(m,length(times)));
[storage_change, Masserr] = deal(zeros(1,length(times)));
[cumUptakeS, cumUptakeC,RWUInstantS,RWUInstantC,RootgrowthC, RootgrowthS]= deal(zeros(length(T),1));
% Create square matrix A and column vector B with zero entries
A = zeros(m);
B = zeros(m,1);
t = dt;              % Set time for first pass through while loop
M = 1;               % Initialize counter
TC=1;                % counter for transpiration 
SumSink = 0;
sumflux = 0;
counter = 0;
cumulative_RWUS=0;
cumulative_RWUC=0;
%************************************************************
%               Begin stepping in time
%************************************************************

while (t-endtime) <= sqrt(eps)
   
    frC = LoC/(LoC+((LmC-LoC)*exp(-(rC*t))));
    frS = LoS/(LoS+((LmS-LoS)*exp(-(rS*t))));
    LrC=LmC*frC;                 %Computing root depth
    LrS=LmS*frS;
    RootgrowthS(TC) = LrS;
    RootgrowthC(TC) = LrC;
    idxC=find(z<LrC);             %index for Corn root
    idxS=find(z<LrS);             %index for Soybean root
    %S(idx)=-Trans(TC,2)*(1-z(idx)/Lr)*(2/Lr);
   
    Tp = -Trans(TC,2);
   
    bxC(idxC) = (1-z(idxC)/LrC)*(2/LrC);
    bxS(idxS) = (1-z(idxS)/LrS)*(2/LrS);
    %bxS=0;
    %S(idx)=Sm*(1-z(idx)/Lr);
    for j = 1:2
        if j==1
            % Create vectors of K and C values for first step
            K = GlenCond(h);
            C = GlenSpecCap(h);
            AlphaS = (Feddes(h))';
            AlphaQ = (FeddesQ(h))';
            S = AlphaS.*bxC*Tp;
            Q = AlphaQ.*bxS*Tp;
        else
            % Create vectors of K and C values for second step
            K = 0.5*(GlenCond(h)+ GlenCond(hstar));
            C = 0.5*(GlenSpecCap(h)+ GlenSpecCap(hstar));
            AlphaS = (0.5*(Feddes(hstar) + Feddes(h)))';
            AlphaQ = (0.5*(FeddesQ(hstar) + FeddesQ(h)))';
            S = AlphaS.*bxC*Tp;
            Q = AlphaQ.*bxS*Tp;
        end
        % Mass bal => thetafin - thetainit = cum.fluxi(cm/h) - cum.fluxo(cm/h) - cum.sink
        % Calculate coefficients a-sub-i, b-sub-i, c-sub-i        
        a = C/dt;
        b(2:m) = (K(1:m-1) + K(2:m))/(2*dz^2);
        c(1:m-1) = (K(1:m-1) + K(2:m))/(2*dz^2);
        % c(m) = K(m)/dz^2;

        % Solve system of m linear equations        
        A(1,1) = 2*a(1) + c(1);
        A(1,2) = -c(1);
            if (30*24 < t)&& (t <= 33*24) %giving flux values at specified time steps
                B(1) = (2*a(1)-c(1))*h(1) + c(1)*h(2) - 2*c(1)*dz + 2*qa/dz + 2*S(1) +2*Q(1);
                flux = qa;
                if j == 2
                    sumflux = sumflux + flux*dt;
                    counter =counter +1;
                end    
            else
                B(1) = (2*a(1)-c(1))*h(1) + c(1)*h(2) - 2*c(1)*dz +2*S(1) +2*Q(1);%qa ,flux term = 0
              
            end
          
           
        for i = 2:m-1                
            A(i,i-1) = -b(i);
            A(i,i) = 2*a(i) + b(i) + c(i);
            A(i,i+1) = -c(i);
            B(i) = (2*a(i)-b(i)-c(i))*h(i) + c(i)*h(i+1) + b(i)*h(i-1) ...
                   + 2*dz*(b(i)-c(i))+2*S(i) + 2*Q(i);
            
        end

        A(m,m-1) = -b(m);
        A(m,m) = 2*a(m) + b(m);
        B(m) = (2*a(m)-b(m))*h(m) + b(m)*h(m-1) + 2*dz*b(m)-2*qb/dz;
        
        if j==1
            hstar = A\B;
            
        else
            h = A\B;
            RWUC = S*(dt*dz); %S [1/T] changing dimension to [L]
            RWUS = Q*(dt*dz); %Q [1/T] changing dimension to [L]
            RWUInstantC(TC) = sum(S*(dz)); %C [1/T] changing dimension to [L]
            RWUInstantS(TC) = sum(Q*(dz)); %S [1/T] changing dimension to [L]
            
            cumulative_RWUC = cumulative_RWUC + sum(RWUC);%cumulative uptake by Corn
            cumulative_RWUS = cumulative_RWUS + sum(RWUS);%cumulative uptake by Weed
            cumUptakeC(TC) = cumulative_RWUC; %for each run saving cummulative sink for corn
            cumUptakeS(TC) = cumulative_RWUS; %for each run saving cummulative sink for weed
            Sink = RWUC+ RWUS;
            SumSink = SumSink+ sum(Sink);
            P = AlphaS;
            PP = AlphaQ;
                        
        end
    end
        
    
    % Commands executed at times specified in the vector "times"
    if any(abs(times-t) <= sqrt(eps))
        headresults(:,M) = h; % Write pressure heads to "headresults"
        TranspirationC(:,M) = RWUC;%writing cumulative uptake for Corn at each times       
        TranspirationS(:,M) = RWUS;%writing cumulative uptake for Weed at each times
        FeddesAlphaS(:,M)= P;
        FeddesAlphaQ(:,M)= PP;
        thetaresults(:,M) = GlenReten(h); % Do same for water content
        %Calculating change in storage
        storage_change(M) = dz*(sum(thetaresults(:,M)-thetainit));
        %Masserr(M) = SumSink + sumflux -storage_change(M);
        Masserr(M) = abs((SumSink + sumflux -storage_change(M))/(SumSink + sumflux))*100;
        rho(:,M) = dt*K./(C*dz^2);  % Fourier number [-]
              
        M = M + 1;
    end
    
    t = t + dt;                     %Increment time for next pass 
    TC= TC+1;
end



str=strings(length(times),1);       %string for graphical use
for i = 1:length(times)
    str(i)= "t "+ num2str(times(i)/24)+" day";
end

%******************************************************************
%                       PLOTING Grpahs
%******************************************************************

% Plot pressure head versus depth
figure(1);
Pa=plot(headresults,z);
g=legend(Pa,str);
g.Location="best";
set(gca,'YDir','reverse');
xlabel('h (cm)');
ylabel('Depth (cm)');
g.Box="off";

% Plot water contents versus depth
figure(2);
Pb=plot(thetaresults,z);
ga=legend(Pb, str);
ga.Location="best";
set(gca,'YDir','reverse');
xlabel('\theta (cm^{3} cm^{-3})');
ylabel('Depth (cm)');
ga.Box="off";

% % Plot water contents versus depth
% figure(3);
% Pb=plot(TranspirationC,z);
% ga=legend(Pb, str);
% ga.Location="best";
% set(gca,'YDir','reverse');
% xlabel('Corn Rootwater Uptake (cm hr^{-1})');
% ylabel('Depth (cm)');
% ga.Box="off";
% 
% % Plot water contents versus depth
% figure(4);
% Pb=plot(TranspirationS,z);
% ga=legend(Pb, str);
% ga.Location="best";
% set(gca,'YDir','reverse');
% xlabel('weed Rootwater Uptake (cm hr^{-1})');
% ylabel('Depth (cm)');
% ga.Box="off";

% Plot water contents versus depth

forplot =NaN(length(T),2);

forplot(:,1) = -1*cumUptakeS;
forplot(:,2) = -1*cumUptakeC;
stri = strings(2,1);
stri(1,1) = "Weed";
stri(2,1) = "Corn";
figure(5);
pb=plot(T/24,forplot);
hold on 
xline(30);
xline(35);
hold off
ga=legend(pb, stri);
ga.Location="best";
xlabel('Time(h)');
ylabel('Cummulative Sink (cm)');


% %Root growth dynamics for Corn
% figure(6);
% plot(T/24,RootgrowthC)
% xlabel("Days")
% ylabel("Depth(cm)")
% title("Corn Root growth dynamics")
% %Root growth dynamics for weed
% figure(7);
% plot(T/24,RootgrowthS)
% xlabel("Days")
% ylabel("Depth(cm)")
% title("Weed Root growth dynamics")
%Masserr Figure
figure(8)
plot(times/24,Masserr,'o')
xlabel("days")
ylabel("%age error")
title("Mass Blanace Error")


%Plot Etc
figure(10)
plot(T/24,-1*RWUInstantC*24)
xline(30);
xlabel("day")
ylabel("Transpiration (cm/day)")
title("Instantaneous Transpiration Corn")

%Plot Etc
figure(11)
plot(T/24,-1*RWUInstantS*24)
xline(30);
xlabel("day")
ylabel("Transpiration (cm/day)")
title("Instantaneous Transpiration weed")

%Plot Etc 
figure(9)
plot(Trans(:,1)/24,Trans(:,2)*24)
xlabel("days")
ylabel("Tp(cm/day)")
title("Instantaneous Potential Transpiration Rate")
xline(30);