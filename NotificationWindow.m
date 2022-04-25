%Underlying code for the paper 'The effect of notification window on the
%epidemiological impact of COVID-19 mobile contact tracing applications'
% Written by Trystan Leng  Last update: 5/11/21



clearvars 
%% ---------------Figure 2A-------------------------------------------------
%Here we plot the normalised test probability profiles from Hellewell et
%al, which we use to inform the time the base case is detected

%For LFTs
Detection_times = readtable('LFT_curve_summary.csv'); %probability of testing +ve on a given day (in tenths of day)
%For PCRs
%Detection_times = readtable('PCR_curve_summary.csv'); %probability of testing +ve on a given day (in tenths of day)

T = Detection_times.median/sum(Detection_times.median)/0.1; %probability of testing +ve on given day, conditional on that individual having tested positive 
lower = Detection_times.lower_95/sum(Detection_times.lower_95)/0.1;
upper = Detection_times.upper_95/sum(Detection_times.upper_95)/0.1;


figure;
plot(0:0.1:30, T, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); hold on
plot(0:0.1:30, upper, 'Color', '#d7191c', 'LineWidth', 1.5, 'LineStyle', '--'); hold on
plot(0:0.1:30, lower, 'Color', '#2c7bb6', 'LineWidth', 1.5, 'LineStyle', '-.'); hold on
legend('normalised median', 'normalised upper 95% credible interval', 'normalised lower 95% credible interval');
ylabel('Relative probability base case tests positive on that day');
xlabel('Day of base case in their infectiousness profile, d ');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);
%}
%Model parameters
%Parameters
shape = 5.62; scale = 0.98; %Infectiousness profile parameters
symshape = 5.807; symscale = 0.948; %Time to symptom onset parameters
Rsym = 1; Rasym = 0.5; %R for symptomatics and asymptomatics
A = 0.3; %Proportion of population asymptomatic
delayinresults = 0; %delay in results (for LFT)
%delayinresults = 2; %delay in results (for PCR)
i_sympt = 10.5; i_notif = 10.5; %isolation period after symptom onset and after notification (from day of last contact)
V = 0.7; %proportion vaccinated
c1 = 0.37; % susceptibility of vaccinated individuals
c2 = 0.37; %transmissibility of vaccinated individuals



Asym_novacc = A*(1-V);
Sym_novacc = (1-A)*(1-V);
Asym_vacc = A*V;
Sym_vacc  = (1-A)*V;
EffPopSize = (1-V) + c1*V;

Duration = length(T);
ProbInterval = 101;


%% ------------------------Generating results for Figure 2B ----------------------------------- Varying Window length
%Here we explore the impact of different notification window lengths,
%assuming 100% active app use (i.e. p = 1)


parfor windowlength = 1:101
    
    Window = (windowlength - 1)/10;
    
    for INIT = 1:Duration
        %INIT

        init = (INIT-1)/((Duration-1)/30); %day base case has test
                
        
        %Approximations using i_notif and i_sympt = 10.5
        
        fun1 = @(t)(gamcdf(t,shape,scale)+ 1 - gamcdf(t+i_sympt, shape, scale)).*gampdf(t,symshape,symscale);         
        Int1 = integral(fun1, 0, Inf);
        
        fun2 = @(t)(gamcdf(delayinresults+init-t, shape, scale) + 1 - gamcdf(i_notif, shape, scale)).*gampdf(t, shape, scale);
        Int2 = integral(fun2, init-Window, init);  %For subcase 3a

               
        fun3 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(tau,shape,scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau,symshape,symscale);
        fun3max = @(t) init - t + delayinresults; 
        Int3 = integral2(fun3, init-Window, init,  0, fun3max);%For subcase 3b
        

        fun4 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(init+delayinresults-t, shape, scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau, symshape, symscale);   
        fun4min = @(t) init + delayinresults - t;  
        fun4max = i_notif;
        Int4 = integral2(fun4, init-Window, init, fun4min, fun4max);%For subcase 3c
        
        
        fun5 = @(t,tau) (gampdf(t,shape,scale).*(gamcdf(init + delayinresults - t, shape, scale) + (gamcdf(tau, shape, scale) - gamcdf(i_notif, shape, scale)) + (1 - gamcdf(tau+i_sympt, shape, scale)))).*gampdf(tau,symshape,symscale);
        fun5min = i_notif;
        Int5 = integral2(fun5, init-Window, init, fun5min, Inf);%For subcase 3d

        p = 1; %probability of adhering to pings

        
        
        RG = A*Rasym + (1-A)*Rsym*Int1;
        
       %Without vaccination      
        %{
        Secondary1 = Rasym*gamcdf(init-Window,shape,scale)*(A*Rasym + (1-A)*Rsym*Int1);
        Secondary2 = Rasym*(1-p)*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale))*(A*Rasym + (1-A)*Rsym*Int1);
        Secondary3a = Rasym*p*A*Rasym*Int2;
        Secondary3b = Rasym*p*(1-A)*Rsym*Int3;
        Secondary3c = Rasym*p*(1-A)*Rsym*Int4;
        Secondary3d = Rasym*p*(1-A)*Rsym*Int5;
        %}
        
       
       %{
        %Old way - was working fine!
        
        %With vaccination
       
        
        Primary = Rasym*gamcdf(init,shape,scale);
        Secondary1 = Rasym*gamcdf(init-Window,shape,scale)*((Asym_novacc*Rasym + Asym_vacc*R_vasym) + (Sym_novacc*Rsym +Sym_vacc*R_vsym)*Int1);
        Secondary2 = (1-p).*(Rasym*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale))*((Asym_novacc*Rasym + Asym_vacc*R_vasym) + (Sym_novacc*Rsym +Sym_vacc*R_vsym)*Int1)); 
        Secondary2vacc = p.*(Rasym*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale))*( Asym_vacc*R_vasym + Sym_vacc*R_vsym*Int1));
        Secondary2 = Secondary2 + Secondary2vacc;
        Secondary3a = p.*(Rasym*(Asym_novacc*Rasym)*Int2);
        Secondary3b = p.*(Rasym*(Sym_novacc*Rsym)*Int3);
        Secondary3c = p.*(Rasym*(Sym_novacc*Rsym)*Int4);
        Secondary3d = p.*(Rasym*(Sym_novacc*Rsym)*Int5); 
        Secondary3 = Secondary3a + Secondary3b + Secondary3c + Secondary3d;
        Secondary4 = 0;
        %}
       
        %new way with vaccination - same as paper
        
       Primary = Rasym*EffPopSize*gamcdf(init,shape,scale);
       Secondary1 = Rasym*(1-V)*gamcdf(init-Window,shape,scale)*EffPopSize*RG;
       Secondary2 = Rasym*(1-V)*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale))*(1-p)*EffPopSize*RG;
       Secondary3a = p.*((Rasym*(1-V))*EffPopSize*(A*Rasym)*Int2);
       Secondary3b = p.*((Rasym*(1-V))*EffPopSize*((1-A)*Rsym)*Int3);
       Secondary3c = p.*((Rasym*(1-V))*EffPopSize*((1-A)*Rsym)*Int4);
       Secondary3d = p.*((Rasym*(1-V))*EffPopSize*((1-A)*Rsym)*Int5);
       Secondary4 = Rasym*c1*V*gamcdf(init,shape,scale)*EffPopSize*c2*RG;

              
        PrimMat(INIT,windowlength) = Primary;
        PingMat(INIT,windowlength) = (gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale));
        AdhereMat(INIT,windowlength) = p*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale));
        SecMat(INIT,windowlength) = Secondary1 + Secondary2 + Secondary3a + Secondary3b + Secondary3c + Secondary3d + Secondary4;
        Sec1Mat(INIT,windowlength) =  Secondary1;
        Sec2Mat(INIT,windowlength) = Secondary2;
        Sec3Mat(INIT,windowlength) = Secondary3a + Secondary3b + Secondary3c + Secondary3d;
        Sec3aMat(INIT,windowlength) = Secondary3a;
        Sec3bMat(INIT,windowlength) = Secondary3b;
        Sec3cMat(INIT,windowlength) = Secondary3c;
        Sec3dMat(INIT,windowlength) = Secondary3d;
        Sec4Mat(INIT,windowlength) = Secondary4;
    end
    
end




%Plotting Figure 2B

Detectinterp = @(Detect,t)interp1(0:0.1:30, Detect, t);
Secinterp = @(t,w)interp1(0:0.1:30, SecMat(:,w), t);
Priminterp = @(t,w)interp1(0:0.1:30, PrimMat(:,w), t);
DetectSecfun = @(Detect,t,w)Detectinterp(Detect,t).*Secinterp(t,w);
DetectPrimfun = @(Detect,t,w)Detectinterp(Detect,t).*Priminterp(t,w);

DetectSecfun2 = @(t,w)gampdf(t, symshape, symscale).*Secinterp(t,w);
DetectPrimfun2 = @(t,w)gampdf(t, symshape, symscale).*Priminterp(t,w);


for w = 1:101
    IntSec(w) = integral(@(t)DetectSecfun(T,t,w), 0, 30); 
    IntPrim(w) = integral(@(t)DetectPrimfun(T,t,w), 0, 30); 
    
    IntSecup(w) = integral(@(t)DetectSecfun(upper,t,w), 0, 30); 
    IntPrimup(w) = integral(@(t)DetectPrimfun(upper,t,w), 0, 30); 
    
    IntSecdown(w) = integral(@(t)DetectSecfun(lower,t,w), 0, 30); %For subcase 1 and subcase 2
    IntPrimdown(w) = integral(@(t)DetectPrimfun(lower,t,w), 0, 30); %For subcase 1 and subcase 2
    
    
    %IntSec2(w) = integral(@(t)DetectSecfun2(t,w), 0, 30); 
    %IntPrim2(w) = integral(@(t)DetectPrimfun2(t,w), 0, 30); 

end

%Calculate R
Rwindow = IntSec./IntPrim;
Rwindowup = IntSecup./IntPrimup;
Rwindowdown = IntSecdown./IntPrimdown;
above = 100*(1 - (Rwindowup/Rwindowup(1)));
below = 100*(1 - (Rwindowdown/Rwindowdown(1)));


%Rwindow2 = IntSec2./IntPrim2;

figure;
plot(0:0.1:10, 100*(1 - (Rwindow/Rwindow(1))) , 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); hold on
%plot(0:0.1:10, 100*(1 - (Rwindow2/Rwindow2(1))) , 'Color', 'b', 'LineWidth', 1.5, 'LineStyle', '-'); hold on

xs = 0:0.1:10;
p = patch([xs fliplr(xs)], [below above(end:-1:1)], 'k', 'FaceAlpha', '0.25', 'EdgeAlpha', 0, 'HandleVisibility', 'off');
xlim([0 10]);
ylabel('% reduction in R^*  (with 100% active app use)');
xlabel('Notification window length, w');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);
%}

%% -------------------Generating results regarding lower levels of active app use (Figures 2c, 2d and 2a-f)
%Here we consider the impact of active app use on results, considering 5
%scenarios:

%scenario 1 - 5 day window, equal active app use to 2 day window
%scenario 2 - 5 day window, 80% active app use compared to 2 day window
%scenario 3 - 5 day window, 60% active app use compared to 2 day window
%scenario 4 - 5 day window, 40% active app use compared to 2 day window
%scenario 5 - 2 day window

%These results are then used to generate Figures 2c,2d and S3

Duration = 1 + 10*(length(T) - 1); %Finer granularity required for Figure S3 Plots
ProbInterval = 1001;


%Duration = length(T);
clearvars PrimMat SecMat
tic


for scenario = 1:5
    scenario
    if scenario <5
        Window = 5;
        RelAdherence = 1 - (scenario-1)*0.2; %relative probability of adherence for 5-day window
    else
        Window = 2;
        RelAdherence = 1; %relative probably = 1 because 2-day window
    end
    

    parfor INIT = 1:Duration
        init = (INIT-1)/((Duration-1)/30); %day base case has test  
        
        
        
        %Approximations using i_notif and i_sympt = 10.5
        
        fun1 = @(t)(gamcdf(t,shape,scale)+ 1 - gamcdf(t+i_sympt, shape, scale)).*gampdf(t,symshape,symscale);         
        Int1 = integral(fun1, 0, Inf);
        
        fun2 = @(t)(gamcdf(delayinresults+init-t, shape, scale) + 1 - gamcdf(i_notif, shape, scale)).*gampdf(t, shape, scale);
        Int2 = integral(fun2, init-Window, init);  %For subcase 3a

               
        fun3 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(tau,shape,scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau,symshape,symscale);
        fun3max = @(t) init - t + delayinresults; 
        Int3 = integral2(fun3, init-Window, init,  0, fun3max);%For subcase 3b
        

        fun4 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(init+delayinresults-t, shape, scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau, symshape, symscale);   
        fun4min = @(t) init + delayinresults - t;  
        fun4max = i_notif;
        Int4 = integral2(fun4, init-Window, init, fun4min, fun4max);%For subcase 3c
        
        
        fun5 = @(t,tau) (gampdf(t,shape,scale).*(gamcdf(init + delayinresults - t, shape, scale) + (gamcdf(tau, shape, scale) - gamcdf(i_notif, shape, scale)) + (1 - gamcdf(tau+i_sympt, shape, scale)))).*gampdf(tau,symshape,symscale);
        fun5min = i_notif;
        Int5 = integral2(fun5, init-Window, init, fun5min, Inf);%For subcase 3d
             
        RG = A*Rasym + (1-A)*Rsym*Int1;
        

        PrimVec = []; SecVec = []; 
       
        ps = RelAdherence*((1:ProbInterval)-1)/(ProbInterval-1); %probability of adhering to pings
        
        %Old way with vaccination - was working fine!
        %{
        Primary = Rasym*gamcdf(init,shape,scale);
        Secondary1 = Rasym*gamcdf(init-Window,shape,scale)*((Asym_novacc*Rasym + Asym_vacc*R_vasym) + (Sym_novacc*Rsym +Sym_vacc*R_vsym)*Int1);
        Secondary2 = (1-ps).*(Rasym*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale))*((Asym_novacc*Rasym + Asym_vacc*R_vasym) + (Sym_novacc*Rsym +Sym_vacc*R_vsym)*Int1)); 
        Secondary2vacc = ps.*(Rasym*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale))*( Asym_vacc*R_vasym + Sym_vacc*R_vsym*Int1));
        Secondary2 = Secondary2 + Secondary2vacc;
        Secondary3a = ps.*(Rasym*(Asym_novacc*Rasym)*Int2);
        Secondary3b = ps.*(Rasym*(Sym_novacc*Rsym)*Int3);
        Secondary3c = ps.*(Rasym*(Sym_novacc*Rsym)*Int4);
        Secondary3d = ps.*(Rasym*(Sym_novacc*Rsym)*Int5);
        Secondary4 = 0;
        %}
              
        %New way with vaccination - in paper
       
       Primary = Rasym*EffPopSize*gamcdf(init,shape,scale);
       Secondary1 = Rasym*(1-V)*gamcdf(init-Window,shape,scale)*EffPopSize*RG;
       Secondary2 = (1-ps).*Rasym*(1-V)*(gamcdf(init, shape, scale) - gamcdf(init - Window,shape, scale))*EffPopSize*RG;
       Secondary3a = ps.*((Rasym*(1-V))*EffPopSize*(A*Rasym)*Int2);
       Secondary3b = ps.*((Rasym*(1-V))*EffPopSize*((1-A)*Rsym)*Int3);
       Secondary3c = ps.*((Rasym*(1-V))*EffPopSize*((1-A)*Rsym)*Int4);
       Secondary3d = ps.*((Rasym*(1-V))*EffPopSize*((1-A)*Rsym)*Int5);
       Secondary4 = Rasym*c1*V*gamcdf(init,shape,scale)*EffPopSize*c2*RG;
         
        Secondary3 = Secondary3a + Secondary3b + Secondary3c + Secondary3d;       
        PrimMat(INIT,:) = Primary;
        SecMat(INIT,:) = Secondary1 + Secondary2+Secondary3 + Secondary4; 
    end
    
    Secwindow(scenario,:,:) = SecMat;
    Primwindow(scenario,:,:) = PrimMat;    
end


toc

%Figure 2C Plot

figure;

Detectinterp = @(Detect,t)interp1(0:0.1:30, Detect, t);
Secinterp = @(Scenario,t,p)interp1(0:0.01:30, Secwindow(Scenario,:,p), t);
Priminterp = @(Scenario,t)interp1(0:0.01:30, Primwindow(Scenario,:), t);
DetectSecfun = @(Scenario,Detect,t,p)Detectinterp(Detect,t).*Secinterp(Scenario,t,p);
DetectPrimfun = @(Scenario,Detect,t)Detectinterp(Detect,t).*Priminterp(Scenario,t);

IntPrim5 = integral(@(t)DetectPrimfun(1,T,t), 0, 30);
IntPrimup5 =  integral(@(t)DetectPrimfun(1,upper,t), 0, 30);
IntPrimdown5 = integral(@(t)DetectPrimfun(1,lower,t), 0, 30);
IntPrim2 = integral(@(t)DetectPrimfun(5, T,t), 0, 30);
IntPrimup2 =  integral(@(t)DetectPrimfun(5,upper,t), 0, 30);
IntPrimdown2 = integral(@(t)DetectPrimfun(5, lower,t), 0, 30);

for p = 1:ProbInterval
    %Calculating R
    
    R5(p) = integral(@(t)DetectSecfun(1,T,t,p), 0, 30)/IntPrim5;
    R5up(p) = integral(@(t)DetectSecfun(1,upper,t,p), 0, 30)/IntPrimup5;
    R5down(p) = integral(@(t)DetectSecfun(1,lower,t,p), 0, 30)/IntPrimdown5;
        
    R2(p) = integral(@(t)DetectSecfun(5,T,t,p), 0, 30)/IntPrim5;
    R2up(p) = integral(@(t)DetectSecfun(5,upper,t,p), 0, 30)/IntPrimup5;
    R2down(p) = integral(@(t)DetectSecfun(5,lower,t,p), 0, 30)/IntPrimdown5;

end

%Calculate reductions in R
Reduction5 = 100*(R5(1) - R5)/R5(1);
Reduction5 = Reduction5((1:(ProbInterval-1)/10:end));
Reduction5up = 100*(R5up(1) - R5up)/R5up(1);
Reduction5up = Reduction5up((1:(ProbInterval-1)/10:end));
Reduction5down = 100*(R5down(1) - R5down)/R5down(1);
Reduction5down = Reduction5down((1:(ProbInterval-1)/10:end));

Reduction2 = 100*(R2(1) - R2)/R2(1);
Reduction2 = Reduction2((1:(ProbInterval-1)/10:end));
Reduction2up = 100*(R2up(1) - R2up)/R2up(1);
Reduction2up = Reduction2up((1:(ProbInterval-1)/10:end));
Reduction2down = 100*(R2down(1) - R2down)/R2down(1);
Reduction2down = Reduction2down((1:(ProbInterval-1)/10:end));

xs = (0:(ProbInterval-1)/10:(ProbInterval-1))/(ProbInterval-1);
plot(xs, Reduction5, 'Color', '#3182bd', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 10); hold on
p = patch([xs fliplr(xs)], [Reduction5down Reduction5up(end:-1:1)], 'b', 'FaceAlpha', '0.25', 'FaceColor', '#3182bd', 'EdgeAlpha', 0, 'HandleVisibility', 'off');

plot( xs, Reduction2, 'k:', 'Color', '#e6550d', 'LineWidth', 2, 'Marker', '+', 'MarkerSize', 10); hold on
p = patch([xs fliplr(xs)], [Reduction2down Reduction2up(end:-1:1)], 'b', 'FaceAlpha', '0.25', 'FaceColor','#e6550d', 'EdgeAlpha', 0, 'HandleVisibility', 'off');
ylabel('% reduction in R^*', 'FontSize', 16);
xlabel('Proportion of population who are active app users', 'FontSize', 16);
legend('5-day window', '2-day window');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);




figure;
%Figure 2D Plot


for i = 1:ProbInterval
    for j = 1:ProbInterval       
        p2 = (i-1)/1000;       
        R5nonadherence(i,j) = interp1(0:0.001:1, R5, p2*(j-1)/1000); %find R5 at different probabilities at finer granularity via interpolation
    end
end


matrixprob = 100*(R2 - R5nonadherence)/(mean(mean(SecMat(:,1)))/mean(mean(PrimMat))); %divided by R with p = 0

%For LFT
zMin = -40;
zMax = 40;
%For PCR
%zMin = -20;
%zMax = 20;

%CT = [0.105882352941176,0.470588235294118,0.215686274509804;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.462745098039216,0.164705882352941,0.513725490196078];
CT = [0,0.266666666666667,0.105882352941176;0.0784313725490196,0.392156862745098,0.200000000000000;0.105882352941176,0.470588235294118,0.215686274509804;0.168627450980392,0.552941176470588,0.247058823529412;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.800000000000000,0.929411764705882,0.764705882352941;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.858823529411765,0.764705882352941,0.882352941176471;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.505882352941176,0.274509803921569,0.572549019607843;0.462745098039216,0.164705882352941,0.513725490196078;0.400000000000000,0.0823529411764706,0.439215686274510;0.250980392156863,0,0.294117647058824];
ax1 = axes;
a = pcolor(0:0.001:1, 0:0.1:100, matrixprob); hold on   
a.EdgeColor=  'none';   
plot([0 1], [60.4 60.4], 'k--', 'LineWidth', 2); %For LFT
%plot([0 1], [67.6 67.6], 'k--', 'LineWidth', 2); %For PCR

colormap(ax1, CT);
caxis( [zMin zMax]);
xlabeldat = {'Proportion of population who actively use app', 'given 2-day window'};
xlabel(xlabeldat, 'FontSize', 16);
ylabeldat = {'% active app use with 5-day window' 'compared to 2-day window'};
ylabel(ylabeldat, 'FontSize', 16);
h = colorbar;
%set(get(h,'label'),'string','Improvement of 5-day window over 2-day window');
set(get(h,'label'),'string',{'% reduction in R^* with 5-day window - % reduction in R^*', ' with 2-day window'}, 'Fontsize', 14);
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);




%Figure S3
xs = (0:(ProbInterval-1))/(ProbInterval-1);
ys = (0:(Duration-1))/((Duration-1)/30);

%Figure S3A - R surface for 5-day window
 figure;
 
Secsurf = squeeze(Secwindow(1,:,:));
Primsurf = (Primwindow(1,:));
Tempsurf = 100 * (1 - (Secsurf./Primsurf')/(mean(mean(SecMat(:,1)))/mean(mean(PrimMat))));
zMin = 0;
zMax = 80;
ax1 = axes;
a = pcolor(ys, xs, Tempsurf'); hold on   
a.EdgeColor=  'none';
%CT2 = [0.968627450980392,0.984313725490196,1;0.870588235294118,0.921568627450980,0.968627450980392;0.776470588235294,0.858823529411765,0.937254901960784;0.619607843137255,0.792156862745098,0.882352941176471;0.517647058823530,0.745098039215686,0.862745098039216;0.419607843137255,0.682352941176471,0.839215686274510;0.258823529411765,0.572549019607843,0.776470588235294;0.129411764705882,0.443137254901961,0.709803921568628;0.0313725490196078,0.317647058823529,0.611764705882353;0.0313725490196078,0.188235294117647,0.419607843137255];
CT2 = [0.968627450980392,0.984313725490196,1;0.870588235294118,0.921568627450980,0.968627450980392;0.776470588235294,0.858823529411765,0.937254901960784;0.619607843137255,0.792156862745098,0.882352941176471;0.419607843137255,0.682352941176471,0.839215686274510;0.258823529411765,0.572549019607843,0.776470588235294;0.129411764705882,0.443137254901961,0.709803921568628;0.0313725490196078,0.270588235294118,0.580392156862745];
colormap(ax1, CT2);
caxis( [zMin zMax]);
axis([0 100 1 14 -0.1 1.1])
xlim([0 14]);
ylim([0 1]);
ylabel('Proportion of population who actively use app', 'FontSize', 16);
xlabel('Day of base case in their infectiousness profile');
h = colorbar;
set(get(h,'label'),'string','% reduction in R^* with 5-day window');
set(gcf, 'Position', [300, 300, 600, 500])
set(gca, 'FontSize', 16);
   
%Figure S3B - R surface for 2-day window
 figure;
 
Secsurf = squeeze(Secwindow(5,:,:));
Primsurf = (Primwindow(5,:));
Tempsurf = 100 * (1 - (Secsurf./Primsurf')/(mean(mean(SecMat(:,1)))/mean(mean(PrimMat))));
zMin = 0;
zMax = 80;
ax1 = axes;
a = pcolor(ys, xs, Tempsurf'); hold on   
a.EdgeColor=  'none';
%CT2 = [0.968627450980392,0.984313725490196,1;0.870588235294118,0.921568627450980,0.968627450980392;0.776470588235294,0.858823529411765,0.937254901960784;0.619607843137255,0.792156862745098,0.882352941176471;0.517647058823530,0.745098039215686,0.862745098039216;0.419607843137255,0.682352941176471,0.839215686274510;0.258823529411765,0.572549019607843,0.776470588235294;0.129411764705882,0.443137254901961,0.709803921568628;0.0313725490196078,0.317647058823529,0.611764705882353;0.0313725490196078,0.188235294117647,0.419607843137255];
CT2 = [0.968627450980392,0.984313725490196,1;0.870588235294118,0.921568627450980,0.968627450980392;0.776470588235294,0.858823529411765,0.937254901960784;0.619607843137255,0.792156862745098,0.882352941176471;0.419607843137255,0.682352941176471,0.839215686274510;0.258823529411765,0.572549019607843,0.776470588235294;0.129411764705882,0.443137254901961,0.709803921568628;0.0313725490196078,0.270588235294118,0.580392156862745];
colormap(ax1, CT2);
caxis( [zMin zMax]);
axis([0 100 1 14 -0.1 1.1])
xlim([0 14]);
ylim([0 1]);
ylabel('Proportion of population who actively use app', 'FontSize', 16);
xlabel('Day of base case in their infectiousness profile');
h = colorbar;
set(get(h,'label'),'string','% reduction in R^* with 2-day window');   
set(gcf, 'Position', [300, 300, 600, 500])
set(gca, 'FontSize', 16);
   
%Figure S3C - Difference assuming equal adherence
figure;
Secsurf = squeeze(Secwindow(5,:,:) - Secwindow(1,:,:));
Secsurf(abs(Secsurf) < 1e-10) = 0; %remove numerical errors
Primsurf = Primwindow(5,:);
Tempsurf = 100 * ((Secsurf./Primsurf')/(mean(mean(SecMat(:,1)))/mean(mean(PrimMat))));      
zMin = -40;
zMax = 40;
ax1 = axes;
a = pcolor(ys, xs, Tempsurf'); hold on   
a.EdgeColor=  'none';
%CT = [0.105882352941176,0.470588235294118,0.215686274509804;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.462745098039216,0.164705882352941,0.513725490196078];
CT = [0,0.266666666666667,0.105882352941176;0.0784313725490196,0.392156862745098,0.200000000000000;0.105882352941176,0.470588235294118,0.215686274509804;0.168627450980392,0.552941176470588,0.247058823529412;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.800000000000000,0.929411764705882,0.764705882352941;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.858823529411765,0.764705882352941,0.882352941176471;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.505882352941176,0.274509803921569,0.572549019607843;0.462745098039216,0.164705882352941,0.513725490196078;0.400000000000000,0.0823529411764706,0.439215686274510;0.250980392156863,0,0.294117647058824];
colormap(ax1, CT);
caxis( [zMin zMax]);
axis([0 100 1 14 -0.1 1.1])
xlim([0 14]);
ylim([0 1]);
ylabeldat = {'Proportion of population who actively use app', 'given 2-day window'};
ylabel(ylabeldat, 'FontSize', 16);
xlabel('Day of base case in their infectiousness profile');
h = colorbar;
%set(get(h,'label'),'string','Improvement of 5-day window over  2-day window');
set(get(h,'label'),'string',{'% reduction in R^* with 5-day window - % reduction in R^*', ' with 2-day window'}, 'Fontsize', 14);
set(gcf, 'Position', [300, 300, 600, 500])
set(gca, 'FontSize', 16);

%Figure S3D - Difference assuming 80% adherence relative to 2-day adherence
figure;
Secsurf = squeeze(Secwindow(5,:,:) - Secwindow(2,:,:));
Primsurf = Primwindow(5,:);
Tempsurf = 100 * ((Secsurf./Primsurf')/(mean(mean(SecMat(:,1)))/mean(mean(PrimMat))));      
zMin = -40;
zMax = 40;
ax1 = axes;
a = pcolor(ys, xs, Tempsurf'); hold on   
a.EdgeColor=  'none';
%CT = [0.105882352941176,0.470588235294118,0.215686274509804;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.462745098039216,0.164705882352941,0.513725490196078];
CT = [0,0.266666666666667,0.105882352941176;0.0784313725490196,0.392156862745098,0.200000000000000;0.105882352941176,0.470588235294118,0.215686274509804;0.168627450980392,0.552941176470588,0.247058823529412;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.800000000000000,0.929411764705882,0.764705882352941;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.858823529411765,0.764705882352941,0.882352941176471;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.505882352941176,0.274509803921569,0.572549019607843;0.462745098039216,0.164705882352941,0.513725490196078;0.400000000000000,0.0823529411764706,0.439215686274510;0.250980392156863,0,0.294117647058824];
colormap(ax1, CT);
caxis( [zMin zMax]);
axis([0 100 1 14 -0.1 1.1])
xlim([0 14]);
ylim([0 1]);
ylabeldat = {'Proportion of population who actively use app', 'given 2-day window'};
ylabel(ylabeldat, 'FontSize', 16);
xlabel('Day of base case in their infectiousness profile');
h = colorbar;
%set(get(h,'label'),'string','Improvement of 5-day window over  2-day window');
set(get(h,'label'),'string',{'% reduction in R^* with 5-day window - % reduction in R^*', ' with 2-day window'}, 'Fontsize', 14);
set(gcf, 'Position', [300, 300, 600, 500])
set(gca, 'FontSize', 16);

%Figure S3E - Difference assuming 60% adherence relative to 2-day adherence
figure;
Secsurf = squeeze(Secwindow(5,:,:) - Secwindow(3,:,:));
Primsurf = Primwindow(5,:);
Tempsurf = 100 * ((Secsurf./Primsurf')/(mean(mean(SecMat(:,1)))/mean(mean(PrimMat))));      
zMin = -40;
zMax = 40;
ax1 = axes;
a = pcolor(ys, xs, Tempsurf'); hold on   
a.EdgeColor=  'none';
%CT = [0.105882352941176,0.470588235294118,0.215686274509804;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.462745098039216,0.164705882352941,0.513725490196078];
CT = [0,0.266666666666667,0.105882352941176;0.0784313725490196,0.392156862745098,0.200000000000000;0.105882352941176,0.470588235294118,0.215686274509804;0.168627450980392,0.552941176470588,0.247058823529412;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.800000000000000,0.929411764705882,0.764705882352941;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.858823529411765,0.764705882352941,0.882352941176471;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.505882352941176,0.274509803921569,0.572549019607843;0.462745098039216,0.164705882352941,0.513725490196078;0.400000000000000,0.0823529411764706,0.439215686274510;0.250980392156863,0,0.294117647058824];
colormap(ax1, CT);
caxis( [zMin zMax]);
axis([0 100 1 14 -0.1 1.1])
xlim([0 14]);
ylim([0 1]);
ylabeldat = {'Proportion of population who actively use app', 'given 2-day window'};
ylabel(ylabeldat, 'FontSize', 16);
xlabel('Day of base case in their infectiousness profile');
h = colorbar;
%set(get(h,'label'),'string','Improvement of 5-day window over  2-day window');
set(get(h,'label'),'string',{'% reduction in R^* with 5-day window - % reduction in R^*', ' with 2-day window'}, 'Fontsize', 14);

set(gcf, 'Position', [300, 300, 600, 500])
set(gca, 'FontSize', 16);

%Figure S3F - Difference assuming 40% adherence relative to 2-day adherence
figure;
Secsurf = squeeze(Secwindow(5,:,:) - Secwindow(4,:,:));
Primsurf = Primwindow(5,:);
Tempsurf = 100 * ((Secsurf./Primsurf')/(mean(mean(SecMat(:,1)))/mean(mean(PrimMat))));      
zMin = -40;
zMax = 40;
ax1 = axes;
a = pcolor(ys, xs, Tempsurf'); hold on   
a.EdgeColor=  'none';
%CT = [0.105882352941176,0.470588235294118,0.215686274509804;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.462745098039216,0.164705882352941,0.513725490196078];
CT = [0,0.266666666666667,0.105882352941176;0.0784313725490196,0.392156862745098,0.200000000000000;0.105882352941176,0.470588235294118,0.215686274509804;0.168627450980392,0.552941176470588,0.247058823529412;0.352941176470588,0.682352941176471,0.380392156862745;0.650980392156863,0.858823529411765,0.627450980392157;0.800000000000000,0.929411764705882,0.764705882352941;0.850980392156863,0.941176470588235,0.827450980392157;0.905882352941177,0.831372549019608,0.909803921568627;0.858823529411765,0.764705882352941,0.882352941176471;0.760784313725490,0.647058823529412,0.811764705882353;0.600000000000000,0.439215686274510,0.670588235294118;0.505882352941176,0.274509803921569,0.572549019607843;0.462745098039216,0.164705882352941,0.513725490196078;0.400000000000000,0.0823529411764706,0.439215686274510;0.250980392156863,0,0.294117647058824];
colormap(ax1, CT);
caxis( [zMin zMax]);
axis([0 100 1 14 -0.1 1.1])
xlim([0 14]);
ylim([0 1]);
ylabeldat = {'Proportion of population who actively use app', 'given 2-day window'};
ylabel(ylabeldat, 'FontSize', 16);
xlabel('Day of base case in their infectiousness profile');
h = colorbar;
%set(get(h,'label'),'string','Improvement of 5-day window over  2-day window');
set(get(h,'label'),'string',{'% reduction in R^* with 5-day window - % reduction in R^*', ' with 2-day window'}, 'Fontsize', 14);

set(gcf, 'Position', [300, 300, 600, 500])
set(gca, 'FontSize', 16);
%}



%% ----- Results for Figure S3, exploring the impact of regular testing

%In Figure S3A we plot the detection time profiles of an asymptomatic base
%case assuming they take an LFT every r consecutive days, while in Figure
%3B we consider the % active app use required for a 5-day window to be
%optimal (referred to as epsilon in the manuscript)

%Figure S4A - %Detection PDFs given regular testing
Duration = length(T);

ProbDetect =  Detection_times.median;
r = 0;
for reg = [30 70 300]
    r = r+1;

    for i = 1:301

        if i <= reg+1
            LFTreg(r,i) = ProbDetect(i); %i.e. probability they test +ve that day
        else
            temp = floor((i - 1)/(reg));

            temprange = i - (reg).*(1:temp);

            LFTreg(r,i) = ProbDetect(i)*prod(1-(ProbDetect(temprange))); %probably +ve that day and -ve all previous test days
        end

    end
    
end


LFTreg = (LFTreg'./sum(LFTreg'))'/0.1; %normalising and adjusting for timesteps of tenths of a day

figure;
plot(0:0.1:30, LFTreg(3,:), 'k', 'LineWidth', 2); hold on
plot(0:0.1:30, LFTreg(2,:), 'k--', 'LineWidth', 2, 'Color', [0.64, 0.08, 0.18]); hold on
plot(0:0.1:30, LFTreg(1,:), 'k:', 'LineWidth', 2, 'Color', [1.00, 0.07, 0.65]); hold on
ylabel('Relative probability base case tests positive on that day');
xlabel('Day of base case in their infectiousness profile ');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);
legend('baseline (base case takes one-off test)', '7 days between consecutive tests', '3 days between consecutive tests');

%Figure S3B %adherence level required given base case tests positive on day d
for INIT = 1:Duration
    
init = (INIT-1)/((Duration-1)/30); %day base case has test

fun1 = @(t)(gamcdf(t,shape,scale)+ 1 - gamcdf(t+i_sympt, shape, scale)).*gampdf(t,symshape,symscale);    

G= A*Rasym + (1-A)*integral(fun1, 0, Inf, 'RelTol',0,'AbsTol',1e-12);
Gvec(INIT) = G;

Window = 2;

fun1 = @(t)(gamcdf(t,shape,scale)+ 1 - gamcdf(t+i_sympt, shape, scale)).*gampdf(t,symshape,symscale);    
Int1 = integral(fun1, 0, Inf); %For subcase 1 and subcase 2
fun2 = @(t)(gamcdf(delayinresults+init-t, shape, scale) + 1 - gamcdf(i_notif, shape, scale)).*gampdf(t, shape, scale);
Int2 = integral(fun2, init-Window, init);  %For subcase 3a
fun3 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(tau,shape,scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau,symshape,symscale);
fun3max = @(t) init - t + delayinresults; 
Int3 = integral2(fun3, init-Window, init,  delayinresults, fun3max);%For subcase 3b
fun4 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(init+delayinresults-t, shape, scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau, symshape, symscale);   
fun4min = @(t) init + delayinresults - t;  
fun4max = @(t) init + i_notif - t;
Int4 = integral2(fun4, init-Window, init, fun4min, fun4max);%For subcase 3c
fun5 = @(t,tau) (gampdf(t,shape,scale).*(gamcdf(init + delayinresults - t, shape, scale) + (gamcdf(tau, shape, scale) - gamcdf(i_notif, shape, scale)) + (1 - gamcdf(tau+i_sympt, shape, scale)))).*gampdf(tau,symshape,symscale);
fun5min = fun4max;
Int5 = integral2(fun5, init-Window, init, fun5min, Inf);%For subcase 3d
H2 =  A*Rasym*Int2 + (1-A)*(Int3 + Int4 + Int5);
H2vec(INIT) = H2;


Window = 5;


fun1 = @(t)(gamcdf(t,shape,scale)+ 1 - gamcdf(t+i_sympt, shape, scale)).*gampdf(t,symshape,symscale);    
Int1 = integral(fun1, 0, Inf); %For subcase 1 and subcase 2
fun2 = @(t)(gamcdf(delayinresults+init-t, shape, scale) + 1 - gamcdf(i_notif, shape, scale)).*gampdf(t, shape, scale);
Int2 = integral(fun2, init-Window, init);  %For subcase 3a
fun3 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(tau,shape,scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau,symshape,symscale);
fun3max = @(t) init - t + delayinresults; 
Int3 = integral2(fun3, init-Window, init,  delayinresults, fun3max);%For subcase 3b
fun4 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(init+delayinresults-t, shape, scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau, symshape, symscale);   
fun4min = @(t) init + delayinresults - t;  
fun4max = @(t) init + i_notif - t;
Int4 = integral2(fun4, init-Window, init, fun4min, fun4max);%For subcase 3c
fun5 = @(t,tau) (gampdf(t,shape,scale).*(gamcdf(init + delayinresults - t, shape, scale) + (gamcdf(tau, shape, scale) - gamcdf(i_notif, shape, scale)) + (1 - gamcdf(tau+i_sympt, shape, scale)))).*gampdf(tau,symshape,symscale);
fun5min = fun4max;
Int5 = integral2(fun5, init-Window, init, fun5min, Inf);%For subcase 3d
H5 =  A*Rasym*Int2 + (1-A)*(Int3 + Int4 + Int5);
H5vec(INIT) = H2;

EpsOptimal(INIT) = (G*(gamcdf(init, shape, scale) - gamcdf(init-2, shape, scale)) - H2)/(G*(gamcdf(init, shape, scale) - gamcdf(init-5, shape, scale)) - H5);
EpsOptimalAbove(INIT) = (G*(gamcdf(init, shape, scale) - gamcdf(init-2, shape, scale)) - H2);
EpsOptimalBelow(INIT) = (G*(gamcdf(init, shape, scale) - gamcdf(init-5, shape, scale)) - H5);
end


LFTreg = [];
%adherence level required  given base case tests regular 
for INIT = 2:Duration
        
    init = (INIT-1)/((Duration-1)/30);
    INIT
    
    inits = ((1:Duration)-1)/((Duration-1)/30);
        
    reg = INIT;
    
    for i = 1:301

        if i <= reg
            LFTreg(i) = ProbDetect(i); %i.e. probability they test +ve that day
        else
            temp = floor((i - 1)/(reg-1));

            temprange = i - (reg-1).*(1:temp);

            LFTreg(i) = ProbDetect(i)*prod(1-(ProbDetect(temprange))); %probably +ve that day and -ve all previous test days
        end

    end
    
    LFTreg = (LFTreg'./sum(LFTreg'))'/0.1; %normalising and adjusting for timesteps of tenths of a day


   Detectinterp = @(Detect,t)interp1(0:0.1:30, Detect, t);
   epsilonAbove = @(t)interp1(0:0.1:30, EpsOptimalAbove, t);
   epsilonBelow = @(t)interp1(0:0.1:30, EpsOptimalBelow, t);
   
   EpsAbovefun = @(t) epsilonAbove(t).*Detectinterp(LFTreg,t);
   EpsBelowfun = @(t) epsilonBelow(t).*Detectinterp(LFTreg,t);

   EpsOptimal2(INIT) = (integral(@(t)EpsAbovefun(t), 0, 30))./(integral(@(t)EpsBelowfun(t), 0, 30)); 
  
end

%}

EpsOptimal2(1) = 1; %

xs = (0:(Duration-1))/((Duration-1)/30);
figure;
plot(xs, (100*EpsOptimal), 'k', 'LineWidth', 1.5); hold on
xlabeldat = {'Day of base case in infectiousness profile', '/ frequency of testing'};
xlabel(xlabeldat, 'FontSize', 16);
ylabeldat = {'% active app use (relative to 2-day window)', ' required for 5-day window to be optimal'};
ylabel(ylabeldat, 'FontSize', 16);
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);

ylim([0 100]);
xlim([0 14]);
plot(xs, (100*EpsOptimal2), 'r--',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); hold on
p1 = plot(xs(31), (100*EpsOptimal2(31)), 'rd', 'MarkerSize', 10, 'MarkerFaceColor', [1.00 0.07 0.65]);
p2 = plot(xs(71), (100*EpsOptimal2(71)), 'rs', 'MarkerSize', 12, 'MarkerFaceColor', [0.64 0.08 0.18], 'MarkerEdgeColor',[0.64 0.08 0.18]);
legend('day of base case in infectiousness profile', 'days between consecutive tests', '3 days between consecutive tests: 70%', '7 days between consecutive tests: 62%');


%% Results for Figure S5 - Impact of vaccine effectiveness

Duration = 1 + 10*(length(T) - 1);



for scenario = [1 5]
    scenario
    if scenario < 5
        Window = 5;
    else
        Window = 2;
    end
    
   fun1 = @(t)(gamcdf(t,shape,scale)+ 1 - gamcdf(t+i_sympt, shape, scale)).*gampdf(t,symshape,symscale);         
   Int1 = integral(fun1, 0, Inf);
   G= A*Rasym + (1-A)*Int1;

    
    parfor INIT = 1:Duration
        init = (INIT-1)/((Duration-1)/30); %day base case has test  
      
        %Approximations using i_notif and i_sympt = 10.5
        
        p = [0 1];
        
        fun2 = @(t)(gamcdf(delayinresults+init-t, shape, scale) + 1 - gamcdf(i_notif, shape, scale)).*gampdf(t, shape, scale);
        Int2 = integral(fun2, init-Window, init);  %For subcase 3a
              
        fun3 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(tau,shape,scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau,symshape,symscale);
        fun3max = @(t) init - t + delayinresults; 
        Int3 = integral2(fun3, init-Window, init,  0, fun3max);%For subcase 3b
        
        fun4 = @(t, tau) (gampdf(t, shape, scale).*(1 + gamcdf(init+delayinresults-t, shape, scale) - gamcdf(tau+i_sympt, shape, scale))).*gampdf(tau, symshape, symscale);   
        fun4min = @(t) init + delayinresults - t;  
        fun4max = i_notif;
        Int4 = integral2(fun4, init-Window, init, fun4min, fun4max);%For subcase 3c
              
        fun5 = @(t,tau) (gampdf(t,shape,scale).*(gamcdf(init + delayinresults - t, shape, scale) + (gamcdf(tau, shape, scale) - gamcdf(i_notif, shape, scale)) + (1 - gamcdf(tau+i_sympt, shape, scale)))).*gampdf(tau,symshape,symscale);
        fun5min = i_notif;
        Int5 = integral2(fun5, init-Window, init, fun5min, Inf);%For subcase 3d
        
        H =  A*Rasym*Int2 + (1-A)*(Int3 + Int4 + Int5);

        
        Numerator(INIT,:) = ((1-p).*gamcdf(init,shape, scale) + p.*gamcdf(init-Window,shape,scale))*G + p.*H
        Denominator(INIT) = gamcdf(init,shape,scale);
    end
        Numerator_scenario(scenario,:,:) = Numerator;
        Denominator_scenario(scenario,:) = Denominator;
    
end



Detectinterp = @(Detect,t)interp1(0:0.1:30, Detect, t);
Num_interp = @(Scenario,t,p)interp1(0:0.01:30, Numerator_scenario(Scenario,:,p), t);
Denom_interp = @(Scenario,t)interp1(0:0.01:30, Denominator_scenario(Scenario,:), t);
DetectNumfun = @(Scenario,Detect,t,p)Detectinterp(Detect,t).*Num_interp(Scenario,t,p);
DetectDenomfun = @(Scenario,Detect,t)Detectinterp(Detect,t).*Denom_interp(Scenario,t);


%baseline Rs without vaccination
R5_p0 = (1-V)*integral(@(t)DetectNumfun(1,T,t,1),0,30)/integral(@(t)DetectDenomfun(1,T,t),0,30);
R5_p1 = (1-V)*integral(@(t)DetectNumfun(1,T,t,2),0,30)/integral(@(t)DetectDenomfun(1,T,t),0,30);
R2_p0 = (1-V)*integral(@(t)DetectNumfun(5,T,t,1),0,30)/integral(@(t)DetectDenomfun(5,T,t),0,30);
R2_p1 = (1-V)*integral(@(t)DetectNumfun(5,T,t,2),0,30)/integral(@(t)DetectDenomfun(5,T,t),0,30);

%upperCI
R5_p0_up = (1-V)*integral(@(t)DetectNumfun(1,upper,t,1),0,30)/integral(@(t)DetectDenomfun(1,upper,t),0,30);
R5_p1_up = (1-V)*integral(@(t)DetectNumfun(1,upper,t,2),0,30)/integral(@(t)DetectDenomfun(1,upper,t),0,30);
R2_p0_up = (1-V)*integral(@(t)DetectNumfun(5,upper,t,1),0,30)/integral(@(t)DetectDenomfun(5,upper,t),0,30);
R2_p1_up = (1-V)*integral(@(t)DetectNumfun(5,upper,t,2),0,30)/integral(@(t)DetectDenomfun(5,upper,t),0,30);

%lowerCI
R5_p0_down = (1-V)*integral(@(t)DetectNumfun(1,lower,t,1),0,30)/integral(@(t)DetectDenomfun(1,lower,t),0,30);
R5_p1_down = (1-V)*integral(@(t)DetectNumfun(1,lower,t,2),0,30)/integral(@(t)DetectDenomfun(1, lower, t),0,30);
R2_p0_down = (1-V)*integral(@(t)DetectNumfun(5,lower,t,1),0,30)/integral(@(t)DetectDenomfun(5,lower,t),0,30);
R2_p1_down = (1-V)*integral(@(t)DetectNumfun(5,lower,t,2),0,30)/integral(@(t)DetectDenomfun(5,lower,t),0,30);



for j = 1:1001
    
    Effectiveness = 1- (j-1)*0.001;
    
    R5p0vacc(j) = R5_p0 + V*Effectiveness*G;
    R5p1vacc(j) = R5_p1 + V*Effectiveness*G;
    R2p0vacc(j) = R2_p0 + V*Effectiveness*G;
    R2p1vacc(j) = R2_p1 + V*Effectiveness*G;
    
    R5p0vaccup(j) = R5_p0_up + V*Effectiveness*G;
    R5p1vaccup(j) = R5_p1_up + V*Effectiveness*G;
    R2p0vaccup(j) = R2_p0_up + V*Effectiveness*G;
    R2p1vaccup(j) = R2_p1_up + V*Effectiveness*G;
    
    R5p0vaccdown(j) = R5_p0_down + V*Effectiveness*G;
    R5p1vaccdown(j) = R5_p1_down + V*Effectiveness*G;
    R2p0vaccdown(j) = R2_p0_down + V*Effectiveness*G;
    R2p1vaccdown(j) = R2_p1_down + V*Effectiveness*G;
    
end


lim = 20;

ProbInterval = 1001;

Reduction5 = 100*(R5p0vacc - R5p1vacc)./R5p0vacc;
Reduction5 = Reduction5((1:(ProbInterval-1)/lim:end));
%
Reduction5up = 100*(R5p0vaccup - R5p1vaccup)./R5p0vaccup;
Reduction5up = Reduction5up((1:(ProbInterval-1)/lim:end));
Reduction5down = 100*(R5p0vaccdown - R5p1vaccdown)./R5p0vaccdown;
Reduction5down = Reduction5down((1:(ProbInterval-1)/lim:end));
%

Reduction2 = 100*(R2p0vacc - R2p1vacc)./R2p0vacc;
Reduction2 = Reduction2((1:(ProbInterval-1)/lim:end));

Reduction2up = 100*(R2p0vaccup - R2p1vaccup)./R2p0vaccup;
Reduction2up = Reduction2up((1:(ProbInterval-1)/lim:end));
Reduction2down = 100*(R2p0vaccdown - R2p1vaccdown)./R2p0vaccdown;
Reduction2down = Reduction2down((1:(ProbInterval-1)/lim:end));


xs = (0:(ProbInterval-1)/lim:(ProbInterval-1))/(ProbInterval-1);
figure;
plot(xs, Reduction5, 'Color', '#3182bd', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 10); hold on
p = patch([xs fliplr(xs)], [Reduction5down Reduction5up(end:-1:1)], 'b', 'FaceAlpha', '0.25', 'FaceColor', '#3182bd', 'EdgeAlpha', 0, 'HandleVisibility', 'off');

plot( xs, Reduction2, 'k:', 'Color', '#e6550d', 'LineWidth', 2, 'Marker', '+', 'MarkerSize', 10); hold on
p = patch([xs fliplr(xs)], [Reduction2down Reduction2up(end:-1:1)], 'b', 'FaceAlpha', '0.25', 'FaceColor','#e6550d', 'EdgeAlpha', 0, 'HandleVisibility', 'off');
ylabel('% reduction in R^*', 'FontSize', 16);
xlabel('Vaccine effectiveness (1 - c_1*c_2)', 'FontSize', 16);
legend('5-day window', '2-day window');
set(gcf, 'Position', [300, 300, 600, 500]);
set(gca, 'FontSize', 16);
ylim([0 80]);
