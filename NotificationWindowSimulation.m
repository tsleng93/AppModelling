%%Analagous simulation code to the results for the paper 'The effect of notification window on the
%epidemiological impact of COVID-19 mobile phone apps'
% Written by Trystan Leng  Last update: 1/11/21

clear

%For LFTs
Detection_times = readtable('LFT_curve_summary.csv'); %probability of testing +ve on a given day (in tenths of day)
%For PCRs
%Detection_times = readtable('PCR_curve_summary.csv'); %probability of testing +ve on a given day (in tenths of day)
T = Detection_times.median/sum(Detection_times.median)/0.1; %probability of testing +ve on given day, conditional on that individual having tested positive 



increment = 1000; %number of increments per day
time_horizon = 41; %Time period (in days)
Timeperiod = time_horizon*increment; % Timesteps with an 0.001 increment


%Parameters
alpha = 5.62; beta = 0.98; %Infectiousness profile parameters
symalpha = 5.807; symbeta = 0.948; %Time to symptom onset parameters
Rsym = 1; Rasym = 0.5; %R for symptomatics and asymptomatics
A = 0.3; %Proportion of population asymptomatic
delayinresults = 0; %delay in results (for LFT)
%delayinresults = 2; %delay in results (for PCR)
i_sympt = 10.5*increment; i_notif = 10.5*increment; %isolation period after symptom onset and after notification (from day of last contact)



N = 100000; %index cases 
rs = randsample(1 + (0:300)*100, N, true,T);
for i = 1:Timeperiod   
    G(i) = gampdf((1/increment)*i, alpha, beta); %Populate infectiousness profile
    C(i) = gampdf((1/increment)*i, symalpha, symbeta); %Populate symptoms profile
end

% Construct normalised infectivity profiles
Infectivity_since_infection = [G/sum(G), zeros(1, 10*increment)];
Infectivity_since_infection_sym = Rsym*Infectivity_since_infection; % Symptomatics
Infectivity_since_infection_asym = Rasym*Infectivity_since_infection; % Asymptomatics

% Construct normalised symptoms since infection profile
Symptoms_since_infection = C/sum(C); 


%% ------------------------Generating results for Figure 2B ----------------------------------- Varying Window length
%Here we explore the impact of different notification window lengths,
%assuming 100% active app use (i.e. p = 1). Doing so, we produce an
%analogous figure to Figure 2B, but via explit simulation. To increase accuracy, increase N

for w = 1:11
       
Window = (w-1)*increment;
                    
    tic    
  
        p = 1; %probability of adhering to ping
        
        NumPrim = zeros(1,N); %number of primary (asymptomatic)
        NumPing = zeros(1,N); %number pinged (asymptomatic)
        NumPrimSym = zeros(1,N);  %Primary (symptomatic)
        NumPingSym = zeros(1,N);  %Secondary (symptomatic)
    
        NumAdhere = zeros(1,N); %number adhering
        NumSec = zeros(1,N); %number secondary
        NumSec1 = zeros(1,N); %secondary type 1
        NumSec2 = zeros(1,N); %secondary type 2
        NumSec3 = zeros(1,N); %secondary type 3
      
        NumSec3b = zeros(1,N); %secondary type 3b
        NumSec3c = zeros(1,N); %secondary type 3c
        NumSec3d = zeros(1,N); %secondary type 3d
  
        parfor i = 1:N            

            %r = init;
            r = rs(i);
            
            %Asymptomatics          
            NumPrim(i) = poissrnd(A*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            DayInf = randsample(1:r,NumPrim(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            Pinged =  (r - DayInf) <= Window;  %logical of pinged contacts
            NumPing(i) = sum(Pinged); %Number Pinged


            for j = 1:NumPrim(i)                       

                if Pinged(j)
                    if rand < p
                        %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1;                                
                    
                    NumSecTemp = poissrnd(sum(Infectivity_since_infection_asym(1: delayinresults+r-DayInf(j)))); %before isolation
                    %NumSecTemp2 =poissrnd(sum(Infectivity_since_infection_asym(1+isolationperiod+r-DayInf(j):end)));
                    NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_asym(1+i_notif:end))); %after isolation
                    NumSec3(i) = NumSec3(i) + NumSecTemp+NumSecTemp2;
                    else
                    %does not adhere
                    NumSecTemp = poissrnd(sum(Infectivity_since_infection_asym));                    
                    NumSecTemp2 = 0;
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                   %not pinged -assume at large for entire infectious period
                   NumSecTemp = poissrnd(sum(Infectivity_since_infection_asym));                 
                   NumSecTemp2 = 0;
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end
            
            %Symptomatics
            NumPrimSym(i) = poissrnd((1-A)*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            DayInfSym = randsample(1:r,NumPrimSym(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedSym =  (r - DayInfSym) <= Window;  %logical of pinged contacts
            NumPingSym(i) = sum(PingedSym); %Number Pinged
            
            DaySymptoms = randsample(1:Timeperiod, NumPrimSym(i), true, Symptoms_since_infection);
            
            for j = 1:NumPrimSym(i)                       

                if PingedSym(j) 
                    if rand < p
                    %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1; 
                                        
                    if DaySymptoms(j) < delayinresults+r-DayInfSym(j)
                        
                        %Symptoms before ping
                        NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1: DaySymptoms(j))));                        
                        NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                      
                        NumSec3b(i) = NumSec3b(i) + NumSecTemp + NumSecTemp2;                     
                    elseif DaySymptoms(j) < i_notif+r-DayInfSym(j)
                        %Symptoms after ping but before isolation period is over                 
                        NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3c(i) = NumSec3c(i) + NumSecTemp + NumSecTemp2;
                    else
                        %Symptoms after ping and after isolation period is
                        %over              
                        NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end))) +poissrnd(sum(Infectivity_since_infection_sym(1+i_notif:DaySymptoms(j))));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3d(i) = NumSec3d(i) + NumSecTemp + NumSecTemp2;
                    end
                    else
                    %does not adhere -assume at large until symptom onset
                    NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:DaySymptoms(j))));                     
                    NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                %not pinged - assume at large for entire infectious period
                   NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:DaySymptoms(j)))); 
                   NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end
        end        
            
    Primwindow(w) = mean(NumPrim) + mean(NumPrimSym);
    Pingwindow(w) = mean(NumPing) + mean(NumPingSym);
    Adherewindow(w) = mean(NumAdhere);
    Secwindow(w) = mean(NumSec);
        
    Sec1window(w) = mean(NumSec1);
    Sec2window(w) = mean(NumSec2);
    Sec3window(w) = mean(NumSec3);
    Sec3bwindow(w) = mean(NumSec3b);
    Sec3cwindow(w) = mean(NumSec3c);
    Sec3dwindow(w) = mean(NumSec3d);
       
    toc
    
end


figure;
Rwindow = Secwindow./Primwindow; %Calculate R

plot(0:10, 100*(1 - (Rwindow/Rwindow(1))) , 'Color', 'b', 'LineWidth', 1.5, 'LineStyle', '-'); hold on

ylabel('% reduction in R^*  (with 100% active app use)');
xlabel('Notification window length, w');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);


%% -------------------Generating results regarding lower levels of active app use (Figures 2c, 2d and 2a-f)
%Here we consider the impact of active app use on results, considering 5
%scenarios:

%scenario 1 - 5 day window, equal active app use to 2 day window
%scenario 2 - 5 day window, 80% active app use compared to 2 day window
%scenario 3 - 5 day window, 60% active app use compared to 2 day window
%scenario 4 - 5 day window, 40% active app use compared to 2 day window
%scenario 5 - 2 day window

%These results are then used to generate an analogous Figure to Figure 2c,
%but via explicit stochastic simulation. To increase accuracy, increase N.


for scenario = 1:2
       
       scenario
    
    if scenario ==1
        Window = 5*increment;
    else
        Window = 2*increment;
    end
                
    
     
    tic    
        
    for P = 1:21
        
        P
        
        
        
        p = 0.05*(P-1); %probability of adhering to ping
        
        NumPrim = zeros(1,N); %number of primary (asymptomatic)
        NumPing = zeros(1,N); %number pinged (asymptomatic)
        NumPrimSym = zeros(1,N);  %Primary (symptomatic)
        NumPingSym = zeros(1,N);  %Secondary (symptomatic)
    
        NumAdhere = zeros(1,N); %number adhering
        NumSec = zeros(1,N); %number secondary
        NumSec1 = zeros(1,N); %secondary type 1
        NumSec2 = zeros(1,N); %secondary type 2
        NumSec3 = zeros(1,N); %secondary type 3
      
        NumSec3b = zeros(1,N); %secondary type 3b
        NumSec3c = zeros(1,N); %secondary type 3c
        NumSec3d = zeros(1,N); %secondary type 3d

        

    
        parfor i = 1:N            

            %r = init;
            r = rs(i);
            
            %Asymptomatics          
            NumPrim(i) = poissrnd(A*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            DayInf = randsample(1:r,NumPrim(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            Pinged =  (r - DayInf) <= Window;  %logical of pinged contacts
            NumPing(i) = sum(Pinged); %Number Pinged


            for j = 1:NumPrim(i)                       

                if Pinged(j)
                    if rand < p
                        %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1;                                
                    
                    NumSecTemp = poissrnd(sum(Infectivity_since_infection_asym(1: delayinresults+r-DayInf(j)))); %before isolation
                    %NumSecTemp2 =poissrnd(sum(Infectivity_since_infection_asym(1+isolationperiod+r-DayInf(j):end)));
                    NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_asym(1+i_notif:end))); %after isolation
                    NumSec3(i) = NumSec3(i) + NumSecTemp+NumSecTemp2;
                    else
                    %does not adhere
                    NumSecTemp = poissrnd(sum(Infectivity_since_infection_asym));                    
                    NumSecTemp2 = 0;
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                   %not pinged -assume at large for entire infectious period
                   NumSecTemp = poissrnd(sum(Infectivity_since_infection_asym));                 
                   NumSecTemp2 = 0;
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end
            
            %Symptomatics
            NumPrimSym(i) = poissrnd((1-A)*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            DayInfSym = randsample(1:r,NumPrimSym(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedSym =  (r - DayInfSym) <= Window;  %logical of pinged contacts
            NumPingSym(i) = sum(PingedSym); %Number Pinged
            
            DaySymptoms = randsample(1:Timeperiod, NumPrimSym(i), true, Symptoms_since_infection);
            
            for j = 1:NumPrimSym(i)                       

                if PingedSym(j) 
                    if rand < p
                    %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1; 
                                        
                    if DaySymptoms(j) < delayinresults+r-DayInfSym(j)
                        
                        %Symptoms before ping
                        NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1: DaySymptoms(j))));                        
                        NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                      
                        NumSec3b(i) = NumSec3b(i) + NumSecTemp + NumSecTemp2;                     
                    elseif DaySymptoms(j) < i_notif+r-DayInfSym(j)
                        %Symptoms after ping but before isolation period is over                 
                        NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3c(i) = NumSec3c(i) + NumSecTemp + NumSecTemp2;
                    else
                        %Symptoms after ping and after isolation period is
                        %over              
                        NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end))) +poissrnd(sum(Infectivity_since_infection_sym(1+i_notif:DaySymptoms(j))));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3d(i) = NumSec3d(i) + NumSecTemp + NumSecTemp2;
                    end
                    else
                    %does not adhere -assume at large until symptom onset
                    NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:DaySymptoms(j))));                     
                    NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                %not pinged - assume at large for entire infectious period
                   NumSecTemp = poissrnd(sum(Infectivity_since_infection_sym(1:DaySymptoms(j)))); 
                   NumSecTemp2 = poissrnd(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end
        end        
            
    PrimVec(P) = mean(NumPrim) + mean(NumPrimSym);
    PingVec(P) = mean(NumPing) + mean(NumPingSym);
    AdhereVec(P) = mean(NumAdhere);
    SecVec(P) = mean(NumSec);
        
    Sec1Vec(P) = mean(NumSec1);
    Sec2Vec(P) = mean(NumSec2);
    Sec3Vec(P) = mean(NumSec3);
    Sec3bVec(P) = mean(NumSec3b);
    Sec3cVec(P) = mean(NumSec3c);
    Sec3dVec(P) = mean(NumSec3d);
        
    end
    toc
    
    Secscenario(scenario,:) = SecVec;
    Primscenario(scenario,:) = PrimVec;  
    
end

%Calculate reduction in R
Reduction5 = 100*((Secscenario(1,1)/mean(Primscenario(1,:))) - (Secscenario(1,:)/mean(Primscenario(1,:))))/(Secscenario(1,1)/mean(Primscenario(1,:)));
Reduction2 = 100*((Secscenario(2,1)/mean(Primscenario(2,:))) - (Secscenario(2,:)/mean(Primscenario(2,:))))/(Secscenario(1,1)/mean(Primscenario(2,:)));

figure;
xs = 0:0.05:1;
plot( xs, Reduction5, 'Color', '#3182bd', 'LineWidth', 1.5, 'Marker', '+'); hold on
plot( xs, Reduction2, 'Color', '#e6550d', 'LineWidth', 1.5, 'Marker', 'o'); hold on

ylabel('% reduction in R^*', 'FontSize', 16);
xlabel('Proportion of population who are active app users', 'FontSize', 16);
legend('5-day window', '2-day window');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);
