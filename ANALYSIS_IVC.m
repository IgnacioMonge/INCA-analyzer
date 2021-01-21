%% PV LOOP ANALYZER FOR INCA SYSTEM
%% BY MANUEL IGNACIO MONGE GARCÍA - ignaciomonge@gmail.com
clearvars; close all; clc;

Fs = 250;
dt = 1/Fs;
negative_values = 0;

%% loading interface for LV data from the INCA System
% Vectors are defined in different columns. Please, adjust according to the
% exportated file

[filename,path] = uigetfile('*.csv');
data_INCA= dlmread([path,filename], ',', 5, 0);

    beat_n = data_INCA (:,1)';
    beat_time = data_INCA (:,2)';
    ECG_INCA = data_INCA (:,3)';
    pressure_INCA = data_INCA (:,4)';
    volume_INCA = data_INCA (:,6)';
    

% Constructing vector of interpolation points
DT1 = length(pressure_INCA);
tt1_interp = linspace(0,DT1*dt,DT1);

% Correct volume offset if LVV has negative values to calculate properly
% Emax
 if any(volume_INCA<0)
     negative_values = 1;
     offset = abs(min(volume_INCA))+10; % select an arbitrary addition to avoid 0 values on volume
     volume_INCA = volume_INCA+offset;
 end

% Determining the threshold for R peak finder
 [pks_beats,loc_beats]=findpeaks(volume_INCA,'MinPeakProminence',std(volume_INCA)); 
 % Finding mean cycles duration
 
delta1 = zeros(1,length(loc_beats)-1);
for kk1 = 1:(length(loc_beats)-1)
    delta1(kk1) = (loc_beats(kk1+1)) - loc_beats(kk1);
end
 time_threshold = round(min(rmoutliers(delta1,'gesd')));
 
  %% EKG Analsis
f1 =0.5; 
f2 =45; 
Wn=[f1 f2 ]*2/Fs ; 
N = 3; 
[a,b] = butter(N,Wn);
ecg = filter (a,b,ECG_INCA);
%ecg = sgolayfilt(ECG - mean(ECG),7,21); 
ecg_mean = median(rmoutliers(ecg, 'gesd'));
ecg_std = std(rmoutliers(ecg, 'gesd'));
thresholdx = ecg_mean*1.3;
[rpeaks,locs_Rpeaks] = findpeaks(ecg,'MaxPeakWidth',50,'MinPeakHeight',thresholdx,'MinPeakDistance',time_threshold*0.9);
Rwaves =tt1_interp(locs_Rpeaks);
% figure (1001);plot(ecg,'DisplayName','ecg_tonom');hold on;plot(locs_Rpeaks,rpeaks,'*');hold off

% Finding mean cycles duration
delta2 = zeros(1,length(locs_Rpeaks)-1);
for kk1 = 1:(length(locs_Rpeaks)-1)
    delta2(kk1) = (locs_Rpeaks(kk1+1)) - locs_Rpeaks(kk1);
end
min_duration = round(min(rmoutliers(delta2,'gesd')));

%% Time varying elastance calculation
LV_elastance = pressure_INCA./volume_INCA;


%% Max LV Pressure
max_prominence_P = median(pressure_INCA);
[pks_pressure,locs_pressure]=findpeaks(pressure_INCA,'MinPeakDistance',min_duration*0.8);
systolic_pressure = pks_pressure;
t_systolicpressure = tt1_interp(locs_pressure);

%% Determine EDP and EDV at R wave on the EKG
t_end_diastole = Rwaves;
ED_volume = interp1(tt1_interp,volume_INCA,Rwaves);
ED_pressure = interp1(tt1_interp,pressure_INCA,Rwaves);

% DEFINE THE THRESHOLD FOR SIGNIICANT CHANGES DURING IVC
TF = ischange (ED_volume,'linear','Threshold',max(ED_volume(1:4))*0.8);
point = find(TF == 1);
inflection_point = point(1)-1;

% Determine the Emax and t_Emax
 max_prominence = median(LV_elastance);
 [pks,locs]=findpeaks(LV_elastance,'MinPeakDistance',min_duration*0.8);
 emax = pks;
 time_emax = tt1_interp(locs);
 
 % Determine LV pressure and volume at t_Emax
 LV_ESP = interp1(tt1_interp,pressure_INCA,time_emax);
 LV_ESV = interp1(tt1_interp,volume_INCA,time_emax);
 
%  Correct the offsset of ED and ES values 
 if length(ED_volume) > length(LV_ESP)
     if t_end_diastole(end) > time_emax(end)
         t_end_diastole = t_end_diastole(1:end-1);
         ED_volume = ED_volume(1:end-1);
         ED_pressure = ED_pressure(1:end-1);
     end
 end
 
 % Select on those values with a significant change in EDV
ES_pressure = LV_ESP(inflection_point(end):end);
ES_volume = LV_ESV(inflection_point(end):end);
t_ES_max= time_emax(inflection_point(end):end);
Ees = emax(inflection_point(end):end);


%% LV_Ees calculation: Iterative calculation
[fitresult_ESPVR, gof_ESPVR] = Ees_fit (ES_volume,ES_pressure);
coef = coeffvalues(fitresult_ESPVR);
LV_Ees = coef(1);
V0 = coef(2);

%  Non-linear fit of Emax (ESPVR) J. D. Schipke, Am J Physiol 1988 Vol. 255 Issue 3 Pt 2 Pages H679-84
    yfit2 = polyfit(ES_volume,ES_pressure,2);
    nonlinear_fit = polyval(yfit2,ES_volume);
    Eesnl = sqrt((yfit2(2)^2)-(4*yfit2(1)*yfit2(3)));
  
    
%% ITERATIVE ANALYSIS OF THE PV LOOPS (FOR PRWS, SCI, ETC...)
volumes = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    volume_per_cycle = volume_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    volumes (cycle,1:length(volume_per_cycle)) = volume_per_cycle;
end

pressures = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    pressure_per_cycle = pressure_INCA(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    pressures (cycle,1:length(pressure_per_cycle)) = pressure_per_cycle;
end

timelines = zeros(1,1000);
for cycle = 1:length(t_end_diastole)-1
    time_per_cycle = tt1_interp(tt1_interp<t_end_diastole(cycle+1) & tt1_interp>t_end_diastole(cycle));
    timelines (cycle,1:length(time_per_cycle)) = time_per_cycle;
end


for cycle = 1:length(t_end_diastole)-1
volume_cycles = volumes(cycle,:);
volume_cycles = volume_cycles(volume_cycles~=0);
pressure_cycles = pressures(cycle,:);
pressure_cycles = pressure_cycles(volume_cycles~=0);
timeline_cycle = timelines(cycle,:);
timeline_cycles = timeline_cycle(timeline_cycle~=0)';
EDV_cycles (cycle) = volume_cycles(end);
EDP_cycles (cycle) = pressure_cycles(end);
tED_cycles (cycle) = timeline_cycles(end);
SW_cycles (cycle)= polyarea (volume_cycles,pressure_cycles);
E_cycles = pressure_cycles./volume_cycles;
[~,Emax_cycles_idx] = max(E_cycles);
P_Emax(cycle) = pressure_cycles(Emax_cycles_idx);
V_Emax(cycle) = volume_cycles(Emax_cycles_idx);

figure (33)
hold on;
h = fill(volume_cycles,pressure_cycles,[0.1 0.1 rand(1)]);
alpha = linspace (0,0.02,length(t_end_diastole)-1);
set(h,'facealpha',alpha(cycle));
plot (V_Emax(cycle),P_Emax(cycle),'or');
plot (EDV_cycles(cycle),EDP_cycles(cycle),'og');
end

EDV_ix = EDV_cycles (inflection_point(end):end);
SW_ix = SW_cycles (inflection_point(end):end);    
[fitresult_PRSW, gof_PRSW] = prsw_fit(EDV_ix,SW_ix);
coef_PRSW = coeffvalues(fitresult_PRSW);
PRSW = coef_PRSW(1);
Vw = coef_PRSW(2);

% figure (303)
% plot(fitresult_PRSW, EDV_ix, SW_ix);hold on;
% plot(EDV_ix,SW_ix,'o')   
    
    
% %% END-DIASTOLIC ANALYSIS
% ED_volume_max = ED_volume(inflection_point(end)-1:end);
% ED_pressure_max = ED_pressure (inflection_point(end)-1:end);
% t_ED_max = t_end_diastole (inflection_point(end)-1:end);

% Calculation of EDPVR fit
[fitresult_EDPVR, gof] = createFit (EDV_cycles,EDP_cycles);
% [fitresult_EDPVR, gof] = createFit(ED_volume_max, ED_pressure_max);
coeffvals_EDPVR= coeffvalues(fitresult_EDPVR);
beta_index = coeffvals_EDPVR(2);
 
%% post-hoc analysis if beta-index is negative (usually due to an outlier not detected by the fitting algorithm)
if beta_index < 0
    ED_stiffness = EDV_cycles./EDP_cycles;
    [~,stiffness_outliers] = rmoutliers(ED_stiffness,'mean','ThresholdFactor',1);
    % recalculate EDVPR with the new values
    [fitresult_EDPVR, gof] = createFit(EDV_cycles(~stiffness_outliers), EDP_cycles(~stiffness_outliers));
    coeffvals_EDPVR= coeffvalues(fitresult_EDPVR);
    beta_index = coeffvals_EDPVR(2);
end
    
%% PV-LOOP ANALYSIS
stc_pressure_LV = pressure_INCA(tt1_interp<Rwaves(6));
stc_volume_LV = volume_INCA(tt1_interp<Rwaves(6));
stc_ECG = ECG_INCA(tt1_interp<Rwaves(6));

%% Averaging PV loop
tDE1 = locs_Rpeaks(1:5)';    
t_stc = tt1_interp(tDE1);

delta3 = zeros(1,length(tDE1)-1);
for kk1 = 1:(length(tDE1)-1)
    delta3(kk1) = (tDE1(kk1+1)) - tDE1(kk1);
end

DT1 = mean(delta3);
tt1= linspace(0,DT1*dt,DT1);
pp1 = zeros(length(tDE1)-1,length(tt1));
vv1 = zeros(length(tDE1)-1,length(tt1));
ecg1 =zeros(length(tDE1)-1,length(tt1));
mean_LVpressure = zeros(length(DT1),length(tt1));
mean_LVvolume = zeros(length(DT1),length(tt1));
mean_Pressure = zeros(length(DT1),length(tt1));
mean_ECG = zeros(length(DT1),length(tt1));
kk2 = zeros(length(tDE1)-1,length(tt1));
for kk2 = 1:length(tDE1)-1
    ratio_N1 = ((length(tt1_interp(tDE1(kk2):tDE1(kk2+1)))-1)/DT1);
    T1_interp = tt1*ratio_N1 + tt1_interp(tDE1(kk2));
    pp1(kk2,:)= interp1(tt1_interp(tDE1(kk2):tDE1(kk2+1)),stc_pressure_LV(tDE1(kk2):tDE1(kk2+1)),T1_interp);
    vv1(kk2,:)= interp1(tt1_interp(tDE1(kk2):tDE1(kk2+1)),stc_volume_LV(tDE1(kk2):tDE1(kk2+1)),T1_interp);
    ecg1(kk2,:)= interp1(tt1_interp(tDE1(kk2):tDE1(kk2+1)),stc_ECG(tDE1(kk2):tDE1(kk2+1)),T1_interp);
end

for jj1 = 1:DT1
    mean_LVpressure(jj1) = mean(pp1(:,jj1));
    mean_LVvolume(jj1) = mean(vv1(:,jj1));
    mean_ecg(jj1) = mean(ecg1(:,jj1));
end
mean_Pressure = mean_LVpressure;
mean_Volume = mean_LVvolume;
mean_ECG = mean_ecg;

%% Pressure and Volumes
EDP = mean_Pressure (end);
EDV = mean_Volume (end);
t_enddiastole = tt1(end);
Vd = min(mean_Volume) - (0.5*(max(mean_Volume)-min(mean_Volume)));
Elastance = mean_Pressure./(mean_Volume-Vd);
[Emax,indx_Emax] = max(Elastance);
t_Elastance = tt1(indx_Emax);
ESV = mean_Volume(indx_Emax);
ESP = mean_Pressure(indx_Emax);
SV = EDV - ESV;
Ea = ESP / SV;
t_cycle = length(mean_Pressure)*dt;
heart_rate = 60/(t_cycle);
cardiac_output = (SV * heart_rate)/1000;
LVEF = (SV/EDV)*100;
SW = polyarea (mean_LVvolume,mean_LVpressure);

%% LV dP/dtmax & dP/dtmin
LV_dpdt = gradient(mean_Pressure,tt1);
LV_dpdtmax = max(LV_dpdt);
[~, inddLVP] = unique(LV_dpdt);
t_LVdpdtmax = interp1(LV_dpdt(inddLVP),tt1(inddLVP),LV_dpdtmax);
LV_dpdtmin = min(LV_dpdt);
t_LVdpdtmin = interp1(LV_dpdt(inddLVP),tt1(inddLVP),LV_dpdtmin);
pressure_dpdtmax = interp1(tt1,mean_Pressure,t_LVdpdtmax);
volume_dpdtmax = interp1(tt1,mean_Volume,t_LVdpdtmax);
pressure_dpdtmin = interp1(tt1,mean_Pressure,t_LVdpdtmin);
volume_dpdtmin = interp1(tt1,mean_Volume,t_LVdpdtmin);

%% Single-beat Ees estimation based on Ten Brinke et al. (2010) , Acta physiologica 198(1), 37–46.
% Determine the areas and points of interest
dPdt = LV_dpdt;
dPmax = LV_dpdtmax;
t_dPmax = t_LVdpdtmax; 
dPmin = LV_dpdtmin;
t_dPmin = t_LVdpdtmin;
dPmin85 = dPmin*0.15;
t_dPminlast = tt1(tt1>t_dPmin & tt1<t_dPmin+0.15);
dPdtlast = dPdt(tt1>t_dPmin & tt1<t_dPmin+0.15);
[~, inddP85] = unique(dPdtlast);
t_dPmin85 = interp1(dPdtlast(inddP85),t_dPminlast(inddP85),dPmin85);
poi_sbEes = mean_Pressure(tt1<t_dPmax | (tt1>t_dPmin & tt1<t_dPmin85));
toi_sbEes = tt1(tt1<t_dPmax | (tt1>t_dPmin & tt1<t_dPmin85));
[sbEesresult, gofsbEes]  = sbEes(toi_sbEes, poi_sbEes);
isovolemic_fit = sbEesresult(tt1);
[Pmax_iso,idx_Pmax] = max(isovolemic_fit(tt1>t_dPmax & tt1<t_dPmin));
time_Pmax = tt1(tt1>t_dPmax & tt1<t_dPmin);
t_Pmax = time_Pmax(idx_Pmax);

% Single-beat Ees estimation
Ees_sb = (Pmax_iso - ESP) / (EDV - ESV);

% figure (304) % Single-beat fit plot
% hold on;
% plot(tt1,mean_Pressure);
% plot(toi_sbEes,poi_sbEes,'or');
% plot(tt1,isovolemic_fit);
% plot(t_Pmax,Pmax_iso,'ob');
% ylim ([0 max(Pmax_iso)+5]);
% xlim ([0 toi_sbEes(end)+0.1]);

%% Peak filling rate and Peak ejection rate
LV_dVdt_raw = gradient(mean_Volume,tt1);
LV_dVdt = smoothdata(LV_dVdt_raw,'sgolay','SmoothingFactor',0.05,'Degree',2);
[~,idx_minLVV] = min(mean_Volume);
t_minLVV = tt1(idx_minLVV);
PER = min(LV_dVdt(tt1>0.05 & tt1<t_minLVV));
t_PER = interp1(LV_dVdt,tt1,PER)*1000;
DFT = (tt1(end) - t_LVdpdtmin)*1000;
% Peak Filling rate
[~,idx_minLV] = min(mean_Volume);
t_minvol = tt1(idx_minLV);
minpeak = max(LV_dVdt(tt1>t_minvol)*0.70);
[PFRpks,PFRlocs]=findpeaks(LV_dVdt,'MinPeakHeight',minpeak);
tPFR_candidates = tt1(PFRlocs);
tPFR_candidates2 = tPFR_candidates(tPFR_candidates<tt1(end)*0.8);
PFR_candidates = PFRpks(tPFR_candidates<tt1(end)*0.8);
[~,maxPFRpks_idx] = max(PFR_candidates); % Ensure you take the highest peak
 t_PFR = tPFR_candidates2(maxPFRpks_idx)*1000;
PFR = interp1(tt1,LV_dVdt,t_PFR/1000);

%% Atrial contribution assessment
LV_dVdt_atrial = gradient(mean_Volume ,tt1);
ZeroX = @(v) find(diff(sign(v))); % Returns Zero-Crossing Indices Of dPdt
dVdt_zeros = ZeroX(LV_dVdt);
t_dVdt_end = tt1(dVdt_zeros(end))-dt;
dVdt_end = interp1(tt1,LV_dVdt_atrial,t_dVdt_end);
if dVdt_end > 0
    dVdt_crossing = dVdt_zeros(end-1);
else
        dVdt_crossing = dVdt_zeros(end);
end
t_candidate = tt1(dVdt_crossing);
if t_candidate < t_PFR/1000
    minpeakdVdt = mean(-LV_dVdt(tt1>t_minvol)*0.95);
    [dVdt_peaks,locs_dVdt] = findpeaks(-LV_dVdt_atrial,'MinPeakHeight',minpeakdVdt);
    t_candidates = tt1(locs_dVdt);
    atrial_candidates = t_candidates(t_candidates>t_PFR/1000 & t_candidates<tt1(end)*0.98);
    volume_atrial_candidates = interp1(tt1,mean_Volume,atrial_candidates);
    atrial_candidates2 = atrial_candidates(volume_atrial_candidates<EDV*0.99);
    t_atrialcontraction = atrial_candidates2 (end);
else
    t_atrialcontraction = tt1(dVdt_crossing);
end

atrial_duration = (tt1(end) - t_atrialcontraction)*1000;
atrial_pressure = interp1(tt1,mean_Pressure,t_atrialcontraction);
atrial_volume = interp1(tt1,mean_Volume,t_atrialcontraction);
atrial_contribution = ((EDV-atrial_volume)/EDV) *100;
passive_contribution = (atrial_volume/EDV)*100;
atrial_pressure_gradient = EDP - atrial_pressure;

%% estimation of atrial work (pressure x volume)
atrial_contraction_duration = tt1(tt1>=t_atrialcontraction);
atrial_volume_area = trapz(mean_Volume(tt1>=t_atrialcontraction));
atrial_pressure_area = trapz(abs((mean_Pressure(tt1>=t_atrialcontraction)))*atrial_duration/1000);
atrial_work = atrial_volume_area * atrial_pressure_area;

%% Atrial filling rate
toi_atrial = tt1(tt1>t_atrialcontraction);
dvoi_atrial = LV_dVdt(tt1>t_atrialcontraction);
[PFR_atria,PFRatria_idx] = max(dvoi_atrial);
t_PFRatria = toi_atrial(PFRatria_idx);

%%  Estimation of end isovolumetric contraction and starting ejection
ed = dPmax*0.4;
time_ed = tt1(tt1>t_LVdpdtmax & tt1<t_Elastance);
dP_ed = LV_dpdt(tt1>t_LVdpdtmax & tt1<t_Elastance);
vol_ed = mean_Volume(tt1>t_LVdpdtmax & tt1<t_Elastance);
ela_ed = Elastance(tt1>t_LVdpdtmax & tt1<t_Elastance);
pres_ed = mean_Pressure(tt1>t_LVdpdtmax & tt1<t_Elastance);
[~, ind_dPed] = unique(dP_ed);
t_ted = interp1(dP_ed(ind_dPed),time_ed(ind_dPed),ed);
v_ted = interp1(time_ed,vol_ed,t_ted);
p_ted = interp1(time_ed,pres_ed,t_ted);
e_ted = interp1(time_ed,ela_ed,t_ted);

stiffness = mean_Pressure./mean_Volume;
grad = gradient(stiffness);
grad2 = gradient(grad);
isofix_zeros = ZeroX(grad2);
t_zerofix = tt1(isofix_zeros);
t_isozeros = t_zerofix(t_zerofix>t_ted(1));
t_isostart = t_isozeros(1);
p_isostart = interp1(tt1,mean_Pressure,t_isostart);
v_isostart = interp1(tt1,mean_Volume,t_isostart);
e_isostart = interp1(tt1,Elastance,t_isostart);


relax_dif = abs(t_ted - t_isostart)/3;
if relax_dif == 0
    t_isostart = t_endiso;
else
    if t_ted  < t_isostart
        t_endiso = t_ted +relax_dif;
    else
        t_endiso = t_isostart+relax_dif;
    end
end

iso_Ees = interp1(tt1,Elastance,t_endiso);
pressure_endiso = interp1(tt1,mean_Pressure,t_endiso);
volume_endiso = interp1(tt1,mean_Volume,t_endiso);
pressure_endiso = interp1(tt1,mean_Pressure,t_endiso);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
last_timeline = tt1(tt1 >= t_LVdpdtmin);
last_pressures = mean_Pressure(tt1 >= t_LVdpdtmin);
threshold_tau = round(length(last_timeline)/2);

%% LV suction evaulation: suction pressure (lowest pressure during initial diastole)
suction_pressure = min(last_pressures);
t_suction = interp1(last_pressures,last_timeline,suction_pressure);

% figure (2002) % end of isovolumetric contraction based on PV-loops
%     hold on;    
%     plot(tt1,stiffness,'linewidth',1);
% %     plot(tt1,iso_fit,'-k');
% %     plot(tt1,fit_isoup,'-g');
% %     plot(tt1,fit_isodown,'-m');
%     plot(t_endiso,iso_Ees,'*b');
% %     plot(t_iso,stiffness(maxisofit_idx),'ok');
% %     plot(toi1,Eoi1,'.g','linewidth',2);
% %     plot(toi2,Eoi2,'.r','linewidth',2);
%      ylim ([0 Emax*1.1])
%      xlim ([0 length(Elastance)*dt])
%      yyaxis right;
%      plot(tt1,mean_Pressure);hold on;
%      plot(t_endiso,pressure_endiso,'*k');
%      plot(t_suction,suction_pressure,'ob');

figure (50) % LV dynamic analysis
    axis tight;
    xlabel ('time, s')
    title ('dVdt');
    plot(t_PFR/1000,PFR,'*k');hold on;
%     plot(t_PER/1000,PER,'*g');hold on;
    plot (tt1,LV_dVdt);hold on;
    plot(tt1,LV_dVdt_raw);
    plot(t_PFRatria,PFR_atria,'og');
    yyaxis right;
    plot(tt1,mean_Volume);hold on;
    plot(t_atrialcontraction,atrial_volume,'*r'); hold on;
    plot(t_endiso,volume_endiso,'*k'); hold on;
    plot(t_Elastance,ESV,'*g');

 

%% Determination of the end of relaxation and opening of the mitral valve
% Method 1: use the value of end-diastolic elastance and extrapolate to
% descending limb of pressure (t_relax = EDP);

% end_elastance = Elastance(end);
% elastance_of_interest = Elastance(tt1>t_Elastance & tt1<t_suction);
% time_of_interest = tt1(tt1>t_Elastance & tt1<t_suction);
% volume_of_interest = mean_Volume(tt1>t_Elastance & tt1<t_suction);
% t_endrelax = interp1(elastance_of_interest,time_of_interest,end_elastance);
% volume_endrelax = interp1(tt1,mean_Volume,t_endrelax);
% pressure_endrelax = interp1(tt1,mean_Pressure,t_endrelax);

% Method2: calculate the increase in volume after minimum LV volume
% A) Calculate the rapid filling as the first approximation of the aperture
% of mitral valve using the intersenting tangents
diastolic_volume = transpose(mean_Volume);
diastolic_time = transpose(tt1);
derivative = gradient (diastolic_volume.',diastolic_time);
[~,idx] = max(derivative);
derivative_point = PFR;
tangent_point = interp1(tt1,mean_Volume,t_PFR/1000);
t_tangent = t_PFR/1000;
tangent_dPdt=(diastolic_time-t_tangent)*derivative_point+tangent_point;

tangent_diastolic = (diastolic_time-t_minLVV)*0+min(mean_Volume); 
[intersci,interscy] = polyxpoly(diastolic_time,tangent_dPdt,diastolic_time,tangent_diastolic);
tau_end_volume2 = interp1(tt1,mean_Volume,intersci);
t_tau_end2 = intersci;
tau_end_pressure2 = interp1(tt1,mean_Pressure,t_tau_end2);

% Determine the end using a bifunction fit around the t_tau_end2
toi5 = tt1(tt1>t_tau_end2*0.80 & tt1<t_tau_end2*1.02);
voi5 = mean_Volume(tt1>t_tau_end2*0.80 & tt1<t_tau_end2*1.02);
toi6 = tt1(tt1>t_tau_end2*0.98 & tt1<t_PFR/1000*1.2);
voi6 = mean_Volume(tt1>t_tau_end2*0.98 & tt1<t_PFR/1000*1.2);
[fitup5, gofup5] = createFit4(toi5, voi5);% fit iso 
[fitdown6, gofdown6] = createFit4(toi6, voi6);% fit ejection
fit_relaxup = fitup5(tt1);
fit_relaxdown = fitdown6(tt1);
[intersci_relax,interscy_relax] = polyxpoly(tt1,fit_relaxup,tt1,fit_relaxdown);
intersci_relax_candidates = intersci_relax(intersci_relax > t_minvol & intersci_relax>t_Elastance & intersci_relax<t_PFR/1000+0.025); % fix for slow relaxation rate (t_relax > t_PFR)
t_endrelax1 = intersci_relax_candidates(end);
volume_endrelax1 = interp1(tt1,mean_Volume,t_endrelax1);
pressure_endrelax1 = interp1(tt1,mean_Pressure,t_endrelax1);

end_relax2 = atrial_pressure;
pressure_of_interest2 = mean_Pressure(tt1>t_Elastance & tt1<t_suction);
time_of_interest2 = tt1(tt1>t_Elastance & tt1<t_suction);
volume_of_interest2 = mean_Volume(tt1>t_Elastance & tt1<t_suction);
t_endrelax2 = interp1(pressure_of_interest2 ,time_of_interest2,end_relax2);
volume_endrelax2 = interp1(tt1,mean_Volume,t_endrelax2);
pressure_endrelax2 = interp1(tt1,mean_Pressure,t_endrelax2);

% take the middle point (minimization of the bias)
relax_dif = abs(t_endrelax2-t_endrelax1)/2;
if relax_dif == 0
    t_relax = t_endrelax2;
else
    if t_endrelax2 < t_endrelax1
        t_endrelax = t_endrelax2+relax_dif;
    else
        t_endrelax = t_endrelax1+relax_dif;
    end
end

volume_endrelax = interp1(tt1,mean_Volume,t_endrelax);
pressure_endrelax = interp1(tt1,mean_Pressure,t_endrelax);

% Method 3: use tau_end as end-relaxation time
% t_endrelax = t_tau_end;
% volume_endrelax = tau_end_volume;
% pressure_endrelax = tau_end_pressure;

% figure (51) %% end-relaxation time
% hold on;
%     plot(tt1,mean_Volume);
%     plot(toi5,voi5,'.g');
%     plot(toi6,voi6,'.k')
%     plot(t_tau_end,tau_end_volume,'*k');
% 	plot(tt1,fit_relaxup,'-g');
%     plot(tt1,fit_relaxdown,'-m');
%     plot(t_endrelax,volume_endrelax,'ob');
%     ylim ([min(mean_Volume)-5 max(mean_Volume)+5])
 
%% RELAXATION TIME CONSTANT

% Tau time end-point defined as the next EDP plus a predefined offset (5 mmHg)
tau_end_pressure = last_pressures(end)+5;
t_tau_end = interp1(last_pressures(1:threshold_tau),last_timeline(1:threshold_tau),tau_end_pressure);

% Pressure segment for calculation
tau_timeline = tt1(tt1 >= t_LVdpdtmin & tt1 <= t_tau_end);
tau_pressure_line = mean_Pressure(tt1 >= t_LVdpdtmin & tt1 <= t_tau_end);

 % Check if length of pressure segment is enough for calculation
 if length(tau_pressure_line) < 6
    tau_end_pressure = last_pressures(end);
    t_tau_end = interp1(last_pressures(1:threshold_tau),last_timeline(1:threshold_tau),tau_end_pressure);
    tau_timeline = tt1(tt1 >= t_LVdpdtmin & tt1 <= t_tau_end);
    tau_pressure_line = mean_Pressure(tt1 >= t_LVdpdtmin & tt1 <= t_tau_end);
 end
 
% TAU MIRSKY: tau = time at 50% pressure decrease from pressure at dPdtmin.
pressure_Mirsky = pressure_dpdtmin / 2;
t_end_Mirsky = interp1(last_pressures(1:threshold_tau),last_timeline(1:threshold_tau),pressure_Mirsky);
tau_Mirsky = (t_end_Mirsky - t_LVdpdtmin)*1000;

% TAU WEISS
Weiss_pressure = log(tau_pressure_line);
P_Weiss = robustfit((tau_timeline-tau_timeline(1))*1000,Weiss_pressure);
yfit_Weiss = P_Weiss(2)*(tau_timeline-tau_timeline(1))*1000+P_Weiss(1);
tau_Weiss = 1/-P_Weiss(2);

% TAU GLANTZ
tau_dPdt = LV_dpdt(tt1 >= t_LVdpdtmin & tt1 < t_tau_end);
P_Glantz = robustfit(tau_pressure_line,tau_dPdt);
yfit_Glantz = P_Glantz(2)*tau_pressure_line+P_Glantz(1);
tau_Glantz = 1/-P_Glantz(2)*1000;

% TAU LOGISTIC
Pb = min(mean_Pressure);
Pa = (pressure_dpdtmin*(1+exp(1))) - (Pb*(1+exp(1)));
tau_Logistic = (-Pa/((4*LV_dpdtmin)))*1000;
upstroke_polyfit = polyfit(tau_pressure_line, tau_dPdt,4);
upstroke_fit = polyval(upstroke_polyfit,tau_pressure_line);


figure (1)
subplot(4,1,1)
    axis tight;
    xlabel ('time, s')
    title ('EKG');
    plot(Rwaves,rpeaks,'*k');hold on;
    plot(t_stc,rpeaks(1:5),'*g');hold on;
    plot (tt1_interp,ECG_INCA);hold on
subplot(4,1,2)
    axis tight;
    xlabel ('time, s')
    title ('LV Pressure');
    plot (tt1_interp,pressure_INCA,'r')
    hold on
    plot(tED_cycles,EDV_cycles,'*k');hold on;
    plot(tED_cycles,EDP_cycles,'*g');hold on;
    plot(t_systolicpressure,systolic_pressure,'*b');
subplot (4,1,3)
axis tight;
    xlabel ('time, s')
    title ('LV Volume');
    plot (tt1_interp,volume_INCA)
    hold on
    plot(t_ES_max,ES_volume,'*k');hold on;
    plot(tED_cycles,EDV_cycles,'*g');
subplot(4,1,4)
    axis tight;
    xlabel ('time, s')
    title ('LV elastance');
    plot(t_ES_max,Ees,'*k');hold on;
    plot (tt1_interp,LV_elastance)

figure (2) % PV loop with Emax marks
plot (volume_INCA,pressure_INCA);
hold on;
% plot (ED_volume_max,ED_pressure_max,'*g');
hold on;
plot (mean_LVvolume,mean_LVpressure,'g','linewidth',3)
plot(ESV,ESP,'*m');
plot(EDV,EDP,'*m');
plot(fitresult_ESPVR, ES_volume, ES_pressure);hold on;
% plot(ES_volume,yfit,'r-.','linewidth',3); hold on;
% plot(ES_volume, nonlinear_fit,'g-.','linewidth',3);hold on;
txt_Ees = ['LV Ees: ' num2str(LV_Ees)];
txt_EDPVR = ['EDPVR (beta index): ' num2str(beta_index)];
text(min(volume_INCA)-2,max(pressure_INCA),txt_Ees);
text(min(volume_INCA)-2,max(pressure_INCA)-5,txt_EDPVR);
% edpvr
plot (ES_volume,ES_pressure,'*k'); hold on;
plot(fitresult_EDPVR, EDV_cycles, EDP_cycles);hold on;

figure (3) % PV loop with Emax marks
hold on;
plot (mean_LVvolume,mean_LVpressure,'g','linewidth',2)
plot(ESV,ESP,'ok');
plot(EDV,EDP,'ob');
plot(volume_endiso,pressure_endiso,'Or');
plot(volume_endrelax1,pressure_endrelax1,'*k');
plot(volume_endrelax2,pressure_endrelax2,'*r');
plot(volume_endrelax,pressure_endrelax,'om');
plot(atrial_volume,atrial_pressure,'*b');
plot(v_isostart,p_isostart,'*r');
plot(v_ted,p_ted,'*k');
% legend('LV PV loop', 'end-systole', 'end-diastole', 'end isovolumetric contraction', 'relax1', 'relax2', 'end LV relaxation', 'atrial contraction');
title ('LV Pressue-Volume loop: characteristic points');

figure (4)
subplot (4,1,1)
plot (tt1,mean_Pressure);hold on;
title ('Tau')
xlabel('time, ms')
ylabel('LV pressure, mmHg')
plot (tau_timeline,tau_pressure_line,'r','LineWidth',3);hold on;
plot (t_LVdpdtmin,pressure_dpdtmin,'*k','LineWidth',2);hold on;
plot (t_end_Mirsky,pressure_Mirsky,'*b');hold on;
plot (t_suction,suction_pressure,'*m');hold on;
plot (t_tau_end,tau_end_pressure,'*k','LineWidth',2);hold on;
txt = ['Tau (Mirsky): ' num2str(tau_Mirsky)];
text(0.1,min(tau_pressure_line),txt);
hold on;

subplot (4,1,2)
plot (tau_timeline*100,Weiss_pressure,'*r');hold on;
plot(tau_timeline*100,yfit_Weiss,'r-.');
txt = ['Tau (Weiss): ' num2str(tau_Weiss)];
text((tau_timeline(2)*100),min(Weiss_pressure),txt);
yyaxis right;
plot (tau_timeline*100,tau_pressure_line,'*k'); hold on;

subplot (4,1,3)
plot (tau_pressure_line,tau_dPdt,'*r');hold on;
plot(tau_pressure_line,yfit_Glantz,'k-.');
xlabel('LVP, mmHg')
ylabel('LV dPdt, mmHg/s')
txt = ['Tau (Glantz): ' num2str(tau_Glantz)];
text(tau_pressure_line(end),min(tau_dPdt),txt);

subplot (4,1,4)
plot (tau_pressure_line,tau_dPdt,'*r');hold on;
plot(tau_pressure_line,upstroke_fit,'k-.');hold on;
xlabel('LVP, mmHg')
ylabel('LV dPdt, mmHg/s')
txt = ['Tau (Logistic): ' num2str(tau_Logistic)];
text(tau_pressure_line(end),min(tau_dPdt),txt);


function [fitresult, gof] = createFit(ED_volume_max, ED_pressure_max)
%% Fit: 'EDPVR FIT'.
[xData, yData] = prepareCurveData( ED_volume_max, ED_pressure_max );
[~,outlier] = rmoutliers(ED_pressure_max,'gesd','ThresholdFactor',0.3);
ft = fittype( 'y0+A*exp(B*x)');
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Exclude', outlier );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0 0.0567837103186079 0];
[fitresult, gof] = fit( xData, yData, ft, opts );
 end

 function [fitresult, gof] = createFit2(ES_volume, ES_pressure)
[xData, yData] = prepareCurveData( ES_volume, ES_pressure );
% Set up fittype and options.
ft = fittype( 'poly2' );
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
 end
 
 
function [fitresult, gof] = sbEes(toi_sbEes, poi_sbEes)

%% Fit: '5th order polynomial fit for single-beat Ees estimation based on Pmax'.
[xData, yData] = prepareCurveData( toi_sbEes, poi_sbEes );
% Set up fittype and options.
ft = fittype( 'poly5' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end

function [fitresult, gof] = Ees_fit(ES_volume, ES_pressure)
%CREATEFITS(ES_VOLUME,ES_PRESSURE)
%  Create fits.
%
%  Data for 'untitled fit 1' fit:
%      X Input : ES_volume
%      Y Output: ES_pressure
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( ES_volume, ES_pressure );
[~,outlier] = rmoutliers(ES_pressure,'gesd','ThresholdFactor',0.3);
% Set up fittype and options.
ft = fittype( 'a  * (x - b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Exclude', outlier);
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.Robust = 'Off';
opts.StartPoint = [0 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end