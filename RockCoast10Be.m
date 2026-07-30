% Backward geometric model to explore cosmogenic 10Be concentrations 
%   across an active shore platform as a function of cliff retreat and shore 
%   platform down-wearing scenarios

% Swirad, Z. M. et al. 2020. Cosmogenic exposure dating reveals limited long-term
%   variability in erosion of a rocky coastline. Nature Communcations  11: 3804. 
%   https://doi.org/10.1038/s41467-020-17611-9

% RockCoast10Be.m

% Version 3.0 Pakri
% Modified for the purpose of Baltic Sea sites: 
%  - RSL fall
%  - no tides
%  - lower water denisty (brackish sea)
%  - steady-state shore platform down-wearing

% last update: 2026-07-30 (MATLAB R2023b)
% Zuzanna Swirad (zswirad@igf.edu.pl)
% Institute of Geophysics, Polish Academy of Sciences

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% INPUT FILES
GeomagnScalarRaw = load('geomagnetics.txt'); % Geomagnetic scalar, Sgm; c1: time (yr BP), c2: Sgm: Lifton et al., 2014
ProfileRaw = load('profile.txt'); % c1: distance from the cliff (m); c2: elevation (m a.s.l.)
RSLRaw = load('sealevel.txt'); % Relative sea level, RSL; c1: time (yr BP), c2: RSL (m a.s.l.)
Be_raw = load('measured_raw.txt'); % measured 10Be concentrations; c1: distance from the cliff (m); c2: concentrations (atoms/g); c3: background error (atoms/g)
Be = load('measured_inheritance.txt'); % inheritance-corrected measured 10Be concentrations; c1: distance from the cliff (m); c2: concentrations (atoms/g); c3: total error (atoms/g)

% Interpolate geomagnetic scalar and RSL to 1 year and profile to 1 m spacing
GeomagnScalar = zeros(10001,2);
GeomagnScalar(:,1) = 0:10000;
GeomagnScalar(:,2) = interp1(GeomagnScalarRaw(:,1),GeomagnScalarRaw(:,2),0:10000);
clear GeomagnScalarRaw
Profile = zeros(max(ProfileRaw(:,1))-min(ProfileRaw(:,1))+1,2);
Profile(:,1) = min(ProfileRaw(:,1)):max(ProfileRaw(:,1));
Profile(:,2) = interp1(ProfileRaw(:,1),ProfileRaw(:,2),min(ProfileRaw(:,1)):max(ProfileRaw(:,1)));
clear ProfileRaw
RSL = zeros(max(RSLRaw(:,1))+1,2);
RSL(:,1) = 0:max(RSLRaw(:,1));
RSL(:,2) = interp1(RSLRaw(:,1),RSLRaw(:,2),0:max(RSLRaw(:,1)));
clear RSLRaw

% VARIABLES
ProdRate10Be = 4.009; % 10Be production rate (at/g/yr)
TotalTime = 7000; % total time considered (yr)
Time = 0:TotalTime; % time verctor (yr)
PlatformWidth = 200; % contemporary shore platform width (m)
DistanceCliff = 0:PlatformWidth; % distance from the cliff vector (m)
TidalRange = 0; % tidal range (m)
rho = 1.01; % Baltic seawater density (g/cm3)
lambda = 160; % high-energy neutrons attenuation length (g/cm2): Goesse and Phillips, 2001
lambda2 = 1.3; % attenuation length for particle flux: Dunne et al., 1999
CliffHeight = 20; % (m)
Theta = deg2rad(60); % cliff inclination angle (degrees)
dPhi = deg2rad(180); % subtended azimuth angle (degrees)
CliffWidth = round(CliffHeight/tan(Theta)); % horizontal distance between cliff base and cliff top (0 if vertical)
m_coef = 2.3; % scaling coefficient: Dunne et al., 1999
io = ProdRate10Be*(m_coef+1)/(2*pi); % incidence radiation; calculated from maximal radiation Fmax from Dunne et al., 1999
ThetaPlatform = atan(CliffHeight./(CliffWidth+DistanceCliff)); % cliff inclination across the platform
TopoShield0 = 1 - (io*dPhi/(m_coef+1)*sin(ThetaPlatform).^(m_coef+1)/ProdRate10Be);
TopoShield = zeros(size(DistanceCliff,2),2);
TopoShield(:,1) = DistanceCliff;
TopoShield(:,2) = TopoShield0;
clear TopoShield0

% CLIFF RETREAT SCENARIOS & EXPOSURE AGES
SteadyRetreatRate = [0.05:0.05:1]; % (m/yr)
ChangeRetreatRate = 2:1:10; % at TotalTime retreat rate was X time faster/slower than PresentRetreatRate
PresentRetreatRate = 0.25; % Orviku et al., 2013: used to calculate the rates for acceleration/deceleration scenarios (m/yr)
ScenariosNo = length(SteadyRetreatRate) + 2*length(ChangeRetreatRate);

% Matrices follow such a scheme: top left corner refers to the cliff toe
%   and/or present time; they usually have ScenariosNo/PlatformWidth+1/TotalTime+1 
%   rows/columns:
% TotalTime+1: time BP values (Time vector)
% PlatformWidth+1: distance from the cliff values (DistanceCliff vector)
% ScenariosNo: cliff retreat scenarios

% VISUALISATION
% Re-order the rows (cliff retreat scenarios) so that they follow such a scheme:
% 1. acceleration scenarios such that the rate was X times lower at TotalTime
% 2. steady retreat rate scenarios
% 3. deceleration scenarios such that the rate was X times higher at TotalTime

cc = flipud(jet(ScenariosNo)); % colour scheme for plotting (jet scheme makes acceleration scenarios red and deceleration blue)

RetreatRate = zeros(ScenariosNo,TotalTime+1); % retreat rate through time
for m = 1:length(ChangeRetreatRate) % acceleration scenarios (from the largest change)
    end_rate = PresentRetreatRate/ChangeRetreatRate(m);
    RetreatRate(length(ChangeRetreatRate)+1-m,:) = linspace(PresentRetreatRate,end_rate,TotalTime+1);
end
for m = 1:length(SteadyRetreatRate) % steady retreat scenarios (from the fastest rate)
    RetreatRate(length(ChangeRetreatRate)+length(SteadyRetreatRate)+1-m,:) = SteadyRetreatRate(m);
end
for m = 1:length(ChangeRetreatRate) % deceleration scenarios (from the smallest change)
    end_rate = PresentRetreatRate*ChangeRetreatRate(m);
    RetreatRate(m+length(ChangeRetreatRate)+length(SteadyRetreatRate),:) = linspace(PresentRetreatRate,end_rate,TotalTime+1);
end

CliffPosition = zeros(ScenariosNo,TotalTime+1); % relative cliff position through time
for m = 1:ScenariosNo
    for n = 1:TotalTime
        CliffPosition(m,n+1) = CliffPosition(m,n) + RetreatRate(m,n);
    end
end
CliffPosition = round(CliffPosition,4);

ExpoAges = zeros(ScenariosNo,PlatformWidth+1); % cross-shore exposure ages
for m = 1:ScenariosNo
    CliffPositionScenario = CliffPosition(m,:);
    for w = 1:PlatformWidth+1
        if CliffPositionScenario(end) < (w-1)
            ExpoAges(m,w) = 0;
        else
            ExpoAges(m,w) = find(CliffPositionScenario >= (w-1),1);
        end
    end
end

figure(1)
for n = 1:ScenariosNo
    plot(Time,RetreatRate(n,:),'color',cc(n,:))
    hold on
end
hold off
title('Past cliff retreat rates')
xlabel('Time (yr BP)') 
ylabel('Cliff retreat rate (m/yr)') 

figure(2)
for m = 1:ScenariosNo
    for w = 1:PlatformWidth+1
        if ExpoAges(m,w) == 0
            ExpoAges(m,w) = ExpoAges(m,w)/0;
        end
    end
end
for n=1:ScenariosNo
    plot(DistanceCliff,ExpoAges(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 max(max(ExpoAges))])
title('Cross-shore distribution of exposure ages')
xlabel('Distance from the cliff (m)') 
ylabel('Time (yr BP)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GEOMAGNETIC SCALAR (Sgm)
Sgm = zeros(ScenariosNo,PlatformWidth+1); % cumulative Sgm
for m = 1:ScenariosNo
    for w = 1:PlatformWidth+1
        if ExpoAges(m,w) >-1
            Sgm(m,w) = sum(GeomagnScalar(1:ExpoAges(m,w),2))/ExpoAges(m,w);
        else
            Sgm(m,w) = nan;
        end
    end
end

figure(3)
plot(GeomagnScalar(:,1),GeomagnScalar(:,2),'k')
title('Geomagnetic scalar (Lifton et al., 2014)')
xlabel('Time (yr BP)') 
ylabel('Geomagnetic scalar, Sgm')

figure(4)
for n=1:ScenariosNo
    plot(DistanceCliff,Sgm(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0.9 1.15]) 
title('Cross-shore distribution of total geomagnetic scalar')
xlabel('Distance from the cliff (m)') 
ylabel('Geomagnetic scalar, Sgm') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TOPOGRAPHIC SHIELDING (Stopo)
Stopo = zeros(ScenariosNo,PlatformWidth+1); % cumulative Stopo
for m = 1:ScenariosNo
    DistanceCliffScenario = zeros(TotalTime+1,PlatformWidth+1);
    StopoScenario = zeros(TotalTime+1,PlatformWidth+1);
    for w = 1:PlatformWidth+1
        for n = 1:TotalTime+1
            if ExpoAges(m,w) >= n
                DistanceCliffScenario(n,w) = DistanceCliff(w)-CliffPosition(m,n);
                if DistanceCliffScenario(n,w) <= 0
                    StopoScenario(n,w) = TopoShield(1,2);
                else
                    for w2 = 2:PlatformWidth+1
                        if DistanceCliffScenario(n,w) > TopoShield(w2-1,1) && DistanceCliffScenario(n,w) <= TopoShield(w2,1)
                            StopoScenario(n,w) = (DistanceCliffScenario(n,w)-TopoShield(w2-1,1)) * (TopoShield(w2,2)-TopoShield(w2-1,2)) + TopoShield(w2-1,2);
                        end
                    end
                end
            end
        end
    end
    Stopo(m,:) = sum(StopoScenario,1)./ExpoAges(m,:);
end

figure(5)
plot(TopoShield(:,1),TopoShield(:,2),'k')
title('Present topographic shielding')
xlabel('Distance from the cliff (m)') 
ylabel('Topographic shielding, Stopo')

figure(6)
for n = 1:ScenariosNo
    plot(DistanceCliff,Stopo(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0.5 1])
title('Cross-shore distribution of cumulative topographic shielding')
xlabel('Distance from the cliff (m)') 
ylabel('Topographic shielding, Stopo') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHORE PLATFORM EROSION SCALAR (Ser) && WATER SHIELDING (Sw)
% Assumption: shore platform remains at the same elevation relative to the sea level
% and does not change geometry (width, slope)
% Single platform erosion scenario: steady-state model, down-wearing follows the RSL fall

BestFit = polyfit(Profile(:,1),Profile(:,2),1);
BestFitProfile = BestFit(1)*DistanceCliff + BestFit(2);

figure(7)
plot(Profile(:,1),Profile(:,2),'k')
hold on 
plot(DistanceCliff,BestFitProfile,'r')
hold off
legend('real profile','best-fit profile')
title('Cross-shore profile')
xlabel('Distance from the cliff (m)') 
ylabel('Elevation (m a.s.l.)') 

figure(8)
plot(RSL(:,1),RSL(:,2),'k')
title('Relative sea level')
xlabel('Time (yr BP)') 
ylabel('Elevation (m a.s.l.)')

WaterDepth = -BestFitProfile*100; % in cm
WaterDepth(WaterDepth < 0) = 0;
WaterShield = exp(-rho*WaterDepth/lambda);

figure(9)
plot(DistanceCliff,WaterShield,'k')
title('Present water shielding')
xlabel('Water shielding, Sw') 
ylabel('Elevation (m a.s.l.)')
axis([0 PlatformWidth 0 1])

Elev_noRSL = zeros(ScenariosNo,TotalTime+1); % surface elevation relative to the present without RSL
Down_noRSL = RetreatRate*atan((BestFitProfile(1)-BestFitProfile(end))/PlatformWidth); % shore platform down-wearing rate without RSL
for m = 1:ScenariosNo
    for n = 2:TotalTime+1
        Elev_noRSL(m,n) = Elev_noRSL(m,n-1) + Down_noRSL(m,n-1);
    end
end

Sw = zeros(ScenariosNo,PlatformWidth+1); % cumulative Sw
Elev_RSL = zeros(ScenariosNo,TotalTime+1); % surface elevation relative to the present including RSL
for m = 1:ScenariosNo
    SwScenario = zeros(TotalTime+1,PlatformWidth+1);
    Elev_RSL(m,:) = Elev_noRSL(m,:) + RSL(:,2)';
    for n = 1:TotalTime+1
        Elev_RSL_profile = BestFitProfile + Elev_RSL(m,n);
        for w = 1:PlatformWidth+1
            if w-1 >= CliffPosition(m,n)
                if isnan(Elev_RSL_profile(w)) == 0
                    WaterDepth_RSL = -Elev_RSL_profile(w)*100;
                    if WaterDepth_RSL > 0
                        SwScenario(n,w) = exp(-rho*WaterDepth_RSL/lambda);
                    else
                        SwScenario(n,w) = 1;
                    end
                else
                    SwScenario(n,w) = nan;
                end
            end
        end
    end
    for w = 1:PlatformWidth+1
        if any(isnan(SwScenario(:,w))) == 0
            Sw(m,w) = (sum(SwScenario(:,w))) /ExpoAges(m,w);
        else
            Sw(m,w) = nan;
        end
    end
 end

SerScenario = exp(-Elev_RSL/lambda2);
SerScenario(SerScenario>1) = 1;

Ser = zeros(ScenariosNo,PlatformWidth+1); % cumulative Ser
for m = 1:ScenariosNo
    for n = 1:TotalTime+1
        for w = 1:PlatformWidth+1
            if isnan(Sw(m,w)) == 0
                Ser(m,w) = sum(SerScenario(m,1:ExpoAges(m,w)))/ExpoAges(m,w);
            else
                Ser(m,w) = nan;
            end
        end
    end
end

figure(10)
for n=1:ScenariosNo
    plot(DistanceCliff,Sw(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0.4 1])
title('Cross-shore distribution of cumulative water shielding')
xlabel('Distance from the cliff (m)') 
ylabel('Water shielding, Sw') 

figure(11)
for n = 1:ScenariosNo
    plot(DistanceCliff,Ser(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 1])
title('Cross-shore distribution of cumulative platform erosion scalar')
xlabel('Distance from the cliff (m)') 
ylabel('Platform erosion scalar, Ser') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate modelled concentrations
conc = ExpoAges.*Stopo.*Sgm.*Sw.*ProdRate10Be; 

% Visualise raw concentrations
figure(12)
errorbar(Be_raw(:,1),Be_raw(:,2),Be_raw(:,3),'color','r','LineStyle','none');
hold on
scatter(Be_raw(:,1),Be_raw(:,2),'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
hold on
yline(Be_raw(3,2),'-k')
hold off
title('Cross-shore distribution of measured 10Be concentrations');
xlabel('Distance from the cliff (m)') ;
ylabel('10Be concentrations (atoms/g)', 'Color','k') ;

% Visualise inheritance-corrected concentrations and topography
figure(13)
yyaxis left
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
yyaxis right
plot(Profile(:,1),Profile(:,2),'-b')
plot(DistanceCliff,BestFitProfile,'--b')
yyaxis left;
ax = gca();
ax.YColor = 'k';
title('Cross-shore distribution of measured 10Be concentrations');
xlabel('Distance from the cliff (m)') ;
ylabel('10Be concentrations (atoms/g)', 'Color','k') ;
yyaxis right;
ax = gca();
ax.YColor = 'b';
ylabel('Elevation (m a.s.l.)','Color','b');

% Visualise measured 10Be concentrations and those modelled with all scenarios  
figure(14)
for n=1:ScenariosNo
    plot(DistanceCliff,conc(n,:),'color',cc(n,:))
    hold on
end
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold off
axis([0 PlatformWidth -800 3000])
xlabel('Distance from the cliff (m)') 
ylabel('10Be concentrations (atoms/g)') 

% Find the best-fit scenario using mean squared error
conc_be = conc(:,Be(:,1));

for n = 1: size(conc_be,1)
    conc_msd(n,:) = (Be(:,2) - conc_be(n,:)').^2;
end

figure(15)
for n=1:ScenariosNo
    plot(Be(:,1),conc_msd(n,:),'color',cc(n,:))
    hold on
end
hold off
set(gca,'yscale','log')
xlim([0 PlatformWidth])
title('Mean squared difference')
xlabel('Distance from the cliff (m)') 
ylabel('Mean squared difference, MSD ((atoms/g)^2)') 

mse = mean(conc_msd,2);

figure(16)
for n=1:ScenariosNo
    scatter(n,mse(n),[],cc(n,:),'filled')
    hold on
end
xline(size(ChangeRetreatRate,2)+0.5)
xline(ScenariosNo-size(ChangeRetreatRate,2)+0.5)
hold off
set(gca,'yscale','log')
xlim([1 ScenariosNo])
title('Mean squared error')
xlabel('Scenario (acceleration -> steady [fast to slow] -> deceleration)') 
ylabel('Mean squared error, MSE ((atoms/g)^2)') 

min_value = min(mse);
min_value_diff = (mse-min(mse))/min(mse);
x = find(mse==min_value);
BestFitRetreat = RetreatRate(x,:);
RetreatCum = zeros(size(BestFitRetreat));
RetreatCum(1,1) = BestFitRetreat(1,1);
for n = 2:size(RetreatCum,2)
    RetreatCum(1,n) = RetreatCum(1,n-1) + BestFitRetreat(1,n);
end
PlatformAgeFinder = abs(RetreatCum - PlatformWidth);
min_age_diff = min(PlatformAgeFinder);
PlatformAgeId = find(PlatformAgeFinder==min_age_diff);

disp('present cliff retreat rate (m/yr):')
BestFitRetreat(1,1)
disp('initial cliff retreat rate (m/yr):')
BestFitRetreat(1,end)
if BestFitRetreat(1,1) == BestFitRetreat(1,end)
    disp('steady cliff retreat')
elseif BestFitRetreat(1,1) > BestFitRetreat(1,end)
    disp('acceleration')
else
    disp('deceleration')
end
disp('MSE (atoms/g)^2:')
min_value
disp('Contemporary platform width:')
PlatformWidth
disp('Time of exposure of contemporary platform:')
PlatformAgeId
disp('Uncertainty measure, RMSE (atoms/g):')
rmse = sqrt(min_value)

figure(17)
plot(DistanceCliff,conc(x,:),'k')
hold on
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold off
xlim([0 PlatformWidth])
title('Cross-shore distribution of measured and modelled 10Be concentrations')
xlabel('Distance from the cliff (m)') 
ylabel('10Be concentrations (atoms/g)') 