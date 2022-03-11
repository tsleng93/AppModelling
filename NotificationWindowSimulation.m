%%Analagous simulation code to the results for the paper 'The effect of notification window on the
%epidemiological impact of COVID-19 mobile phone apps'
% Written by Trystan Leng  Last update: 1/11/21

clearvars -except  Rwindowtrue Reduction2true Reduction5true

%For LFTs
Detection_times = readtable('LFT_curve_summary.csv'); %probability of testing +ve on a given day (in tenths of day)
%For PCRs
%Detection_times = readtable('PCR_curve_summary.csv'); %probability of testing +ve on a given day (in tenths of day)
T = Detection_times.median/sum(Detection_times.median)/0.1; %probability of testing +ve on given day, conditional on that individual having tested positive 



increment = 1000; %number of increments per day
time_horizon = 41; %Time period (in days)
Timeperiod = time_horizon*increment; % Timesteps with an 0.001 increment


%Parameters
shape = 5.62; scale = 0.98; %Infectiousness profile parameters
symshape = 5.807; symscale = 0.948; %Time to symptom onset parameters
Rsym = 10; Rasym = 5; %R for symptomatics and asymptomatics
A = 0.3; %Proportion of population asymptomatic
delayinresults = 0; %delay in results (for LFT)
%delayinresults =2*increment; %delay in results (for PCR)
i_sympt = 10.5*increment; i_notif = 10.5*increment; %isolation period after symptom onset and after notification (from day of last contact)


V = 0.7; %proportion vaccinated
c1 = (1-0.63); %susceptibility to infection of vaccinated individuals
c2 = (1-0.63); %transmissibiliy of of infection of vaccinated individuals


Asym_novacc = A*(1-V);
Sym_novacc = (1-A)*(1-V);
Asym_vacc = A*V;
Sym_vacc  = (1-A)*V;

EffPopSize = (1-V) + c1*V;


N = 200000; %index cases 
rs = randsample(1 + (0:300)*100, N, true,T);
%rs = 7001*ones(1,N);
for i = 1:Timeperiod   
    G(i) = gampdf((1/increment)*i, shape, scale); %Populate infectiousness profile
    C(i) = gampdf((1/increment)*i, symshape, symscale); %Populate symptoms profile
end

% Construct normalised infectivity profiles
Infectivity_since_infection = [G/sum(G), zeros(1, 10*increment)];
Infectivity_since_infection_sym = Rsym*Infectivity_since_infection; % Symptomatics
Infectivity_since_infection_asym = Rasym*Infectivity_since_infection; % Asymptomatics

% Construct normalised symptoms since infection profile
Symptoms_since_infection = C/sum(C); 

As = zeros(1,30);
Bs = zeros(1,30);


%HetMarker = 0; %No heterogeneity;
HetMarker = 1; %Heterogeneous contact rates
%HetMarker = 2; %Heterogeneous infectious periods

%% ------------------------Generating results for Figure 2B ----------------------------------- Varying Window length
%Here we explore the impact of different notification window lengths,
%assuming 100% active app use (i.e. p = 1). Doing so, we produce an
%analogous figure to Figure 2B, but via explit simulation. To increase accuracy, increase N


mg = find(G == max(G));

constpot = max(G)/increment;

Primwindowall = zeros(11,N);
Secwindowall = zeros(11,N);


for w = 1:11
    counter = 1;
    w
    
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
        NumSec4 = zeros(1,N); % secondary type 4
      
        NumSec3b = zeros(1,N); %secondary type 3b
        NumSec3c = zeros(1,N); %secondary type 3c
        NumSec3d = zeros(1,N); %secondary type 3d
        
        NumPrimV = zeros(1,N);
        NumPrimSymV = zeros(1,N);
        NumPingV = zeros(1,N);
        NumPingSymV = zeros(1,N);

  
        parfor i = 1:N   
            
        % Construct normalised infectivity profiles
        Infectivity_since_infection = [G/sum(G), zeros(1, 10*increment)];
        Infectivity_since_infection_sym = Rsym*Infectivity_since_infection; % Symptomatics
        Infectivity_since_infection_asym = Rasym*Infectivity_since_infection; % Asymptomatics

            %r = init;
            r = rs(i);
            %HetInf = 1;
            %rand(i);
            
            if HetMarker == 0
                HetInf = 1;
            elseif HetMarker == 1
                HetInf = lognrnd(0,1.1)/1.8314;
            elseif HetMarker == 2
                HetInf = 1;
                indy = rand*max(G);
                G2 = find(G(1:mg) -indy < 0, 1, 'last');
                G3 = find(G(mg:end) -indy < 0, 1);
                Infectivity_since_infection_asym = zeros(1, 51*increment);
                Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
            end
           
            HetInfkeep(i) = HetInf;
            
            %Asymptomatics          
            NumPrim(i) = poissrnd(HetInf*Asym_novacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
                       
            if NumPrim(i) > 0            
            DayInf = randsample(1:r,NumPrim(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            Pinged =  (r - DayInf) <= Window;  %logical of pinged contacts
            NumPing(i) = sum(Pinged); %Number Pinged
            else
            DayInf = [];
            Pinged = [];
            NumPing(i) = 0; %Number Pinged             
            end
            
            %Symptomatics
            NumPrimSym(i) = poissrnd(HetInf*Sym_novacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            if NumPrimSym(i) > 0
            DayInfSym = randsample(1:r,NumPrimSym(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedSym =  (r - DayInfSym) <= Window;  %logical of pinged contacts
            NumPingSym(i) = sum(PingedSym); %Number Pinged
            else
            DayInfSym = [];
            PingedSym =  [];
            NumPingSym(i) = 0; %Number Pinged
            end
                
            DaySymptoms = randsample(1:Timeperiod, NumPrimSym(i), true, Symptoms_since_infection);
            
            %Asymptomatic + Vaccinated
            NumPrimV(i) = poissrnd(c1*HetInf*Asym_vacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            if NumPrimV(i) > 0
            DayInfV = randsample(1:r,NumPrimV(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedV =  (r - DayInfV) <= Window;  %logical of pinged contacts
            NumPingV(i) = sum(PingedV); %Number Pinged
            else
            DayInfV = []; %day primary infected individuals are infected
            PingedV =  [];  %logical of pinged contacts
            NumPingV(i) = 0; %Number Pinged
            end
                       
            %Symptomatic + Vaccinated
            NumPrimSymV(i) = poissrnd(c1*HetInf*Sym_vacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            if NumPrimSymV(i) > 0
            DayInfSymV = randsample(1:r,NumPrimSymV(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedSymV =  (r - DayInfSymV) <= Window;  %logical of pinged contacts
            NumPingSymV(i) = sum(PingedSymV); %Number Pinged
            else
            DayInfSymV = [];
            PingedSymV = [];
            NumPingSymV(i) = 0;
            end
            DaySymptomsV = randsample(1:Timeperiod, NumPrimSymV(i), true, Symptoms_since_infection);
            

            for j = 1:NumPrim(i)
                
                if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
                end
                
                if Pinged(j)
                    if rand < p
                        %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1;
                
                    
                    NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym(1: delayinresults+r-DayInf(j)))); %before isolation
                    %NumSecTemp2 =poissrnd(sum(Infectivity_since_infection_asym(1+isolationperiod+r-DayInf(j):end)));
                    NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym(1+i_notif:end))); %after isolation
                    NumSec3(i) = NumSec3(i) + NumSecTemp+NumSecTemp2;
                    else
                    %does not adhere
                    NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym));                    
                    NumSecTemp2 = 0;
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                   %not pinged -assume at large for entire infectious period
                   NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym));                 
                   NumSecTemp2 = 0;
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end

            
            for j = 1:NumPrimSym(i)
                if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
                end

                if PingedSym(j) 
                    if rand < p
                    %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1; 
                                        
                    if DaySymptoms(j) < delayinresults+r-DayInfSym(j)
                        
                        %Symptoms before ping
                        NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1: DaySymptoms(j))));                        
                        NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                      
                        NumSec3b(i) = NumSec3b(i) + NumSecTemp + NumSecTemp2;                     
                    elseif DaySymptoms(j) < i_notif
                        %Symptoms after ping but before isolation period is over                 
                        NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                   
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3c(i) = NumSec3c(i) + NumSecTemp + NumSecTemp2;
                    else
                        %Symptoms after ping and after isolation period is
                        %over              
                        NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(EffPopSize*HetInf*(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)) + sum(Infectivity_since_infection_sym(1+i_notif:DaySymptoms(j)))));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3d(i) = NumSec3d(i) + NumSecTemp + NumSecTemp2;
                    end
                    else
                    %does not adhere -assume at large until symptom onset
                    NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1:DaySymptoms(j))));                     
                    NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                %not pinged - assume at large for entire infectious period
                   NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1:DaySymptoms(j)))); 
                   NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end
            


             %Asymptomatic + Vaccinated
             for j = 1:NumPrimV(i)  
                                
               %Find individual's 'infectiousness'
               if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
                end
              %Not affected by app  - infectious throughout
              %infectious period 
              NumSecTemp = poissrnd(EffPopSize*HetInf*c2*sum(Infectivity_since_infection_asym));                      
              NumSecTemp2 = 0;
              NumSec4(i) = NumSec4(i) + NumSecTemp + NumSecTemp2; 
              NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
              
             end
            
             
             %Symptomatic + Vaccinated
             for j = 1:NumPrimSymV(i)             
               %Find individual's 'infectiousness'
               if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
               end
               %Not affected by app - infectious until symptom onset
                NumSecTemp = poissrnd(EffPopSize*HetInf*c2*sum(Infectivity_since_infection_sym(1:DaySymptomsV(j))));                     
                NumSecTemp2 = poissrnd(EffPopSize*HetInf*c2*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptomsV(j):end)));
                NumSec4(i) = NumSec4(i) + NumSecTemp + NumSecTemp2;  

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end

            
            
        end        
            
    Primwindow(w) = mean(NumPrim + NumPrimSym + NumPrimSymV +NumPrimV);
   % Primwindow(w) = mean(NumPrim + NumPrimSym);

    Pingwindow(w) = mean(NumPing) + mean(NumPingSym);
    Adherewindow(w) = mean(NumAdhere);
    Secwindow(w) = mean(NumSec);        
    Sec1window(w) = mean(NumSec1);
    Sec2window(w) = mean(NumSec2);
    Sec3window(w) = mean(NumSec3);
    Sec3bwindow(w) = mean(NumSec3b);
    Sec3cwindow(w) = mean(NumSec3c);
    Sec3dwindow(w) = mean(NumSec3d);
    Sec4window(w) = mean(NumSec4);

    
   Primwindowall(w,:) = NumPrim + NumPrimSym + NumPrimSymV + NumPrimV;
   Secwindowall(w,:) = NumSec;

       
    toc
    
end


figure;
Rwindow = Secwindow./Primwindow; %Calculate R

if HetMarker == 0
    Col = [1 0 0];
elseif HetMarker == 1
    Col = [0 0 0];
elseif HetMarker == 2
    Col = [0.49, 0.18, 0.56];
end


plot(0:0.1:10, 100*(1 - (Rwindowtrue/Rwindowtrue(1))) , 'Color', [0.65 0.65 0.65], 'LineWidth', 1.5, 'LineStyle', '-'); hold on %from analytic model
plot(0:10, 100*(1 - (Rwindow/Rwindow(1))) ,  '.', 'MarkerSize', 20, 'Color', Col); hold on


ylabel('% reduction in R^*  (with 100% active app use)');
xlabel('Notification window length, w');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);
if HetMarker == 0
    legend('analytic model', 'simulation');
elseif HetMarker == 1
     legend('analytic model', 'simulation, heterogeneous contact rates');
elseif HetMarker == 2
     legend('analytic model', 'simulation, heterogeneous infectious periods');
end



%}

%
%}
%% -------------------Generating results regarding lower levels of active app use (Figures 2c, 2d and 2a-f)
%Here we consider the impact of active app use on results, considering 5
%scenarios:

%scenario 1 - 5 day window, equal active app use to 2 day window
%scenario 2 - 2 day window

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
   % rng(1);    
    for P = 1:11
        
        P
        
        %counter = 1;
        
        
        p = 0.1*(P-1); %probability of adhering to ping
        
        NumPrim = zeros(1,N); %number of primary (asymptomatic)
        NumPing = zeros(1,N); %number pinged (asymptomatic)
        NumPrimSym = zeros(1,N);  %Primary (symptomatic)
        NumPingSym = zeros(1,N);  %Secondary (symptomatic)
    
        NumAdhere = zeros(1,N); %number adhering
        NumSec = zeros(1,N); %number secondary
        NumSec1 = zeros(1,N); %secondary type 1
        NumSec2 = zeros(1,N); %secondary type 2
        NumSec3 = zeros(1,N); %secondary type 3
        NumSec4 = zeros(1,N); % secondary type 4

      
        NumSec3b = zeros(1,N); %secondary type 3b
        NumSec3c = zeros(1,N); %secondary type 3c
        NumSec3d = zeros(1,N); %secondary type 3d
        
        NumPrimV = zeros(1,N);
        NumPrimSymV = zeros(1,N);
        NumPingV = zeros(1,N);
        NumPingSymV = zeros(1,N);
       
        
    
       parfor i = 1:N            
                

                % Construct normalised infectivity profiles
                Infectivity_since_infection = [G/sum(G), zeros(1, 10*increment)];
                Infectivity_since_infection_sym = Rsym*Infectivity_since_infection; % Symptomatics
                Infectivity_since_infection_asym = Rasym*Infectivity_since_infection; % Asymptomatics
                HetInf = 1;
           

            r = rs(i);
            
             if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
              end
                       
            %Asymptomatics          
            NumPrim(i) = poissrnd(HetInf*Asym_novacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            
            
            if NumPrim(i) > 0            
            DayInf = randsample(1:r,NumPrim(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            Pinged =  (r - DayInf) <= Window;  %logical of pinged contacts
            NumPing(i) = sum(Pinged); %Number Pinged
            else
            DayInf = [];
            Pinged = [];
            NumPing(i) = 0; %Number Pinged             
            end
            
            %Symptomatics
            NumPrimSym(i) = poissrnd(HetInf*Sym_novacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            if NumPrimSym(i) > 0
            DayInfSym = randsample(1:r,NumPrimSym(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedSym =  (r - DayInfSym) <= Window;  %logical of pinged contacts
            NumPingSym(i) = sum(PingedSym); %Number Pinged
            else
            DayInfSym = [];
            PingedSym =  [];
            NumPingSym(i) = 0; %Number Pinged
            end
                
            DaySymptoms = randsample(1:Timeperiod, NumPrimSym(i), true, Symptoms_since_infection);
            
            %Asymptomatic + Vaccinated
            NumPrimV(i) = poissrnd(c1*HetInf*Asym_vacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            if NumPrimV(i) > 0
            DayInfV = randsample(1:r,NumPrimV(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedV =  (r - DayInfV) <= Window;  %logical of pinged contacts
            NumPingV(i) = sum(PingedV); %Number Pinged
            else
            DayInfV = []; %day primary infected individuals are infected
            PingedV =  [];  %logical of pinged contacts
            NumPingV(i) = 0; %Number Pinged
            end
            
            
            %Symptomatic + Vaccinated
            NumPrimSymV(i) = poissrnd(c1*HetInf*Sym_vacc*sum(Infectivity_since_infection_asym(1:r))); %number of primary cases
            if NumPrimSymV(i) > 0
            DayInfSymV = randsample(1:r,NumPrimSymV(i),true,Infectivity_since_infection_asym(1:r)); %day primary infected individuals are infected
            PingedSymV =  (r - DayInfSymV) <= Window;  %logical of pinged contacts
            NumPingSymV(i) = sum(PingedSymV); %Number Pinged
            else
            DayInfSymV = [];
            PingedSymV = [];
            NumPingSymV(i) = 0;
            end
            DaySymptomsV = randsample(1:Timeperiod, NumPrimSymV(i), true, Symptoms_since_infection);

            %Asymptomatics
            for j = 1:NumPrim(i)  
                
               %Find individual's 'infectiousness'
                if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
                end

                if Pinged(j)
                    if rand < p
                        %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1;                                
                    
                    NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym(1: delayinresults+r-DayInf(j)))); %before isolation
                    %NumSecTemp2 =poissrnd(sum(Infectivity_since_infection_asym(1+isolationperiod+r-DayInf(j):end)));
                    NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym(1+i_notif:end))); %after isolation
                    NumSec3(i) = NumSec3(i) + NumSecTemp+NumSecTemp2;
                    else
                    %does not adhere
                    NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym));                    
                    NumSecTemp2 = 0;
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                   %not pinged -assume at large for entire infectious period
                   NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_asym));                 
                   NumSecTemp2 = 0;
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end

            %Symptomatics
            for j = 1:NumPrimSym(i)
                
               %Find individual's 'infectiousness'
                if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
                end

                if PingedSym(j) 
                    if rand < p
                    %adheres and isolates
                    NumAdhere(i) = NumAdhere(i) + 1; 
                                        
                    if DaySymptoms(j) < delayinresults+r-DayInfSym(j)
                        
                        %Symptoms before ping
                        NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1: DaySymptoms(j))));                        
                        NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                      
                        NumSec3b(i) = NumSec3b(i) + NumSecTemp + NumSecTemp2;                     
                   % elseif DaySymptoms(j) < i_notif+r-DayInfSym(j)
                     elseif DaySymptoms(j) < i_notif
                        %Symptoms after ping but before isolation period is over                 
                        NumSecTemp = poissrnd(EffPopSize*sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3c(i) = NumSec3c(i) + NumSecTemp + NumSecTemp2;
                    else
                        %Symptoms after ping and after isolation period is
                        %over              
                        NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1:delayinresults+r-DayInfSym(j))));                                        
                        NumSecTemp2 = poissrnd(EffPopSize*HetInf*(sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)) + sum(Infectivity_since_infection_sym(1+i_notif:DaySymptoms(j)))));                    
                        NumSec3(i) = NumSec3(i) + NumSecTemp + NumSecTemp2;                       
                        NumSec3d(i) = NumSec3d(i) + NumSecTemp + NumSecTemp2;
                    end
                    else
                    %does not adhere -assume at large until symptom onset
                    NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1:DaySymptoms(j))));                     
                    NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                    NumSec2(i) = NumSec2(i) + NumSecTemp + NumSecTemp2;                       
                    end
                else
                %not pinged - assume at large for entire infectious period
                   NumSecTemp = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1:DaySymptoms(j)))); 
                   NumSecTemp2 = poissrnd(EffPopSize*HetInf*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptoms(j):end)));
                   NumSec1(i) = NumSec1(i) + NumSecTemp + NumSecTemp2;
                end

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end
            
            
         %Asymptomatic + Vaccinated
             for j = 1:NumPrimV(i)  
                                
               %Find individual's 'infectiousness'
                if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
                end
              
                  
                    %does not have to isolate - infectious throughout
                    %infectious period 
                    NumSecTemp = poissrnd(EffPopSize*HetInf*c2*sum(Infectivity_since_infection_asym));                      
                    NumSecTemp2 = 0;
                    NumSec4(i) = NumSec4(i) + NumSecTemp + NumSecTemp2; 

               NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
             end
            
             
             %Symptomatic + Vaccinated
             for j = 1:NumPrimSymV(i)
                
               %Find individual's 'infectiousness'
               if HetMarker == 0
                    HetInf = 1;
                elseif HetMarker == 1
                    HetInf = lognrnd(0,1.1)/1.8314;
                elseif HetMarker == 2
                    HetInf = 1;
                    indy = rand*max(G);
                    G2 = find(G(1:mg) -indy < 0, 1, 'last');
                    G3 = find(G(mg:end) -indy < 0, 1);
                    Infectivity_since_infection_asym = zeros(1, 51*increment);
                    Infectivity_since_infection_asym(1, G2:(mg-1+G3)) = Rasym*constpot;
                end
               
               %Not affected by app - infectious until symptom onset
               
                NumSecTemp = poissrnd(EffPopSize*HetInf*c2*sum(Infectivity_since_infection_sym(1:DaySymptomsV(j))));                     
                NumSecTemp2 = poissrnd(EffPopSize*HetInf*c2*sum(Infectivity_since_infection_sym(1+i_sympt+DaySymptomsV(j):end)));
                NumSec4(i) = NumSec4(i) + NumSecTemp + NumSecTemp2;  

                NumSec(i) = NumSec(i) + NumSecTemp + NumSecTemp2;
            end
                        
        end        
            
    PrimVec(P) = mean(NumPrim + NumPrimSym + NumPrimSymV +NumPrimV);
    PingVec(P) = mean(NumPing) + mean(NumPingSym);
    AdhereVec(P) = mean(NumAdhere);
    SecVec(P) = mean(NumSec);
        
    Sec1Vec(P) = mean(NumSec1);
    Sec2Vec(P) = mean(NumSec2);
    Sec3Vec(P) = mean(NumSec3);
    Sec3bVec(P) = mean(NumSec3b);
    Sec3cVec(P) = mean(NumSec3c);
    Sec3dVec(P) = mean(NumSec3d);
    Sec4Vec(P) = mean(NumSec4);
        
    end
    toc
    
    Secscenario(scenario,:) = SecVec;
    Primscenario(scenario,:) = PrimVec;  

    
end
%}


%Calculate reduction in R
if HetMarker == 0
    Col = [1 0 0];
elseif HetMarker == 1
    Col = [0 0 0];
elseif HetMarker == 2
    Col = [0.49, 0.18, 0.56];
end



Reduction5 = 100*((Secscenario(1,1)/Primscenario(1,1)) - (Secscenario(1,:)./Primscenario(1,:)))/(Secscenario(1,1)/mean(Primscenario(1,1)));
Reduction2 = 100*((Secscenario(2,1)/Primscenario(2,1)) - (Secscenario(2,:)./Primscenario(2,:)))/(Secscenario(2,1)/mean(Primscenario(2,1)));

figure;
%xs = 0:0.05:1;
xs = 0:0.1:1;


plot(xs, Reduction5true, 'Color', '#3182bd', 'LineWidth', 1.5); hold on
plot( xs, Reduction5, 'o', 'MarkerSize', 10, 'Color', Col); hold on
plot( xs, Reduction2true, 'k:', 'Color', '#e6550d', 'LineWidth', 2); hold on
plot( xs, Reduction2, '+', 'MarkerSize', 10, 'Color', Col); hold on

ylabel('% reduction in R^*', 'FontSize', 16);
xlabel('Proportion of population who are active app users', 'FontSize', 16);
if HetMarker == 0
    legend('5-day window, analytic', '5-day window, simulation','2-day window, analytic', '2-day window', 'simulation');
elseif HetMarker == 1
    legend('5-day window, analytic', '5-day window, simulation heterogeneous contact rates','2-day window, analytic', '2-day window, simulation heterogeneous contact rates');
elseif HetMarker == 2
    legend('5-day window, analytic', '5-day window, simulation heterogeneous infectious periods','2-day window, analytic', '2-day window, simulation heterogeneous infectious periods');
end
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);
%}