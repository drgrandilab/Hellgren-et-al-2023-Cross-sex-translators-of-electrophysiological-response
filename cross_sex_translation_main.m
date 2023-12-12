% Main code of the cross-sex translator of cardiomyocyte electrophysiology
% described in the following manuscript:
% 
% K.T. Hellgren, H. Ni, S. Morotti, E. Grandi. Predictive male-to-female
% translation of cardiac electrophysiological response to drugs.
% JACC Clin Electrophysiol. 2023 Sep 5: S2405-500X(23)00629-1.
% doi: 10.1016/j.jacep.2023.08.016.

close all
clear
clc

color = [0 0 0]; % Black
color_m = [0.968627451 0.305882353 0.839215686]; % Pink
color_f = [0.333333333 0.62745098 0.984313725]; % Blue
%% Selection of translation direction

male2female = 1; % set to 1 for male-to-female translation
%% Selection of pacing frequency

frequency_selection = 1; % set to 0.5, 1, or 3 (Hz)
%% Selection of drug index
% See list at the end of this file (2-99)

testdrugs = 28;
%% Selection of drug concentration

drug_conc_index = 1; % set to 1, 2, 3, or 4x[ETPC]
%% Plotting options

plot_matrix = 1;
plot_fitting = 0;
plot_validation = 1;
flag_plot_application = 1;
%% Load outputs of simulation of control populations

disp('▣ Cross-sex translator of cardiomyocyte electrophysiology:')

% load matrix all_outputs (columns: N outputs, rows: N trials)
% and array 'output_names' (and 'output_units')
if frequency_selection == 0.5
    disp('Frequency: 0.5 Hz')

    % Female outputs:
    load('Populations/APCAT.BCL.2000.CT.1.GT.2.HL.0.Cmax.0.ALL.mat');
    all_outputs_f = output_matrix;

    % Male outputs:
    load('Populations/APCAT.BCL.2000.CT.1.GT.1.HL.0.Cmax.0.ALL.mat');
    all_outputs_m = output_matrix;
elseif frequency_selection == 1
    disp('Frequency: 1 Hz')

    % Female outputs:
    load('Populations/APCAT.BCL.1000.CT.1.GT.2.HL.0.Cmax.0.ALL.mat');
    all_outputs_f = output_matrix;

    % Male outputs:
    load('Populations/APCAT.BCL.1000.CT.1.GT.1.HL.0.Cmax.0.ALL.mat');
    all_outputs_m = output_matrix;
elseif frequency_selection == 2
    disp('Translation @ 2 Hz')

    % Female outputs:
    load('Populations/APCAT.BCL.500.CT.1.GT.2.HL.0.Cmax.0.ALL.mat');
    all_outputs_f = output_matrix;

    % Male outputs:
    load('Populations/APCAT.BCL.500.CT.1.GT.1.HL.0.Cmax.0.ALL.mat');
    all_outputs_m = output_matrix;
elseif frequency_selection == 3
    disp('Frequency: 3 Hz')

    % Female outputs:
    load('Populations/APCAT.BCL.333.CT.1.GT.2.HL.0.Cmax.0.ALL.mat');
    all_outputs_f = output_matrix;

    % Male outputs:
    load('Populations/APCAT.BCL.333.CT.1.GT.1.HL.0.Cmax.0.ALL.mat');
    all_outputs_m = output_matrix;
end

% Strings
output_names = {'UV', 'APpeak', '|MDP|', 'APamp', 'APD90',...
    'APD70', 'APD50', 'APD30', 'CaTmax', 'CaTmin',...
    'CaTamp', 'CaTttp', 'CaTt50', 'CaTtau', 'Na',...
    'CaSRmax', 'CaSRmin', 'CaSRamp','deltaSCR','deltaEAD',...
    'CaTint','|INaLint|','IK1int'};

output_units = {'mV/ms', 'mV', 'mV', 'mV', 'ms',...
    'ms', 'ms', 'ms', 'mM', 'mM',...
    'mM', 'ms', 'ms', 'ms', 'mM',...
    'mM', 'mM', 'mM', 'mM', 'mV',...
    'ms*mM', 'ms*A/F', 'ms*A/F'};

% Convert cytosolic Ca from mM to nM
all_outputs_m(:,9:11) = all_outputs_m(:,9:11)*1e6;
all_outputs_f(:,9:11) = all_outputs_f(:,9:11)*1e6;
output_units(9:11) = {'nM','nM','nM'};
% Convert CaTint from mM*ms to nM*ms
all_outputs_m(:,21) = all_outputs_m(:,21)*1e6;
all_outputs_f(:,21) = all_outputs_f(:,21)*1e6;
output_units(21) = {'ms*nM'};
% Convert INaLint into absolute value
all_outputs_m(:,22) = abs(all_outputs_m(:,22));
all_outputs_f(:,22) = abs(all_outputs_f(:,22));
%% Check basic properties, and define X and Y matrices

% Outputs indexes:
% 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
% 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
% 16) CaSRmax 17) CaSRmin 18) CaSRamp
% 19) delta_SCR 20) delta_EAD 21) CaT_int 22) -INaL_int 23) IK1_int

% Exclude cells based on APD90 value (if less than 50 ms or more that 650 ms)
% and presence of SCR/EAD events
good_trials_m = (all_outputs_m(:,5) > 50) .* (all_outputs_m(:,5) < 650) .* ...
    (all_outputs_m(:,19)==0) .* (all_outputs_m(:,20)==0);
good_trials_f =  (all_outputs_f(:,5) > 50) .* (all_outputs_f(:,5) < 650) .* ...
    (all_outputs_f(:,19)==0) .* (all_outputs_f(:,20)==0);
good_trials = good_trials_m.*good_trials_f;
good_trials_log = logical(good_trials);

% Good count: number of good trials
good_count_total = sum(good_trials_log);

% Output selection for X and Y
% output_selection_X = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]; % all
output_selection_X = [1 3 4 5 7 10 11 12 13 14]; % Input sex
output_selection_Y = output_selection_X; % Output sex

output_names_X = output_names(output_selection_X);
output_units_X = output_units(output_selection_X);

output_names_Y = output_names(output_selection_Y);
output_units_Y = output_units(output_selection_Y);

N_outputs_X = length(output_names_X);
N_outputs_Y = length(output_names_Y);

% Select outputs
if male2female == 1
    disp('Direction: Male to Female')
    all_outputs_X = all_outputs_m(:,output_selection_X);
    all_outputs_Y = all_outputs_f(:,output_selection_Y);

    color_X = color_m; color_Y = color_f;
else
    disp('Direction: Female to Male')
    all_outputs_X = all_outputs_f(:,output_selection_X);
    all_outputs_Y = all_outputs_m(:,output_selection_Y);

    color_X = color_f; color_Y = color_m;
end

% Good_outputs: array with parameters from good trials only
good_outputs_X = all_outputs_X(good_trials_log,:);
good_outputs_Y = all_outputs_Y(good_trials_log,:);
%% Separation Fitting & Validation Groups

Nval = 400; Nfit = good_count_total-Nval;

actual_inputs = good_outputs_X(end-Nval+1:end,:);
actual_outputs = good_outputs_Y(end-Nval+1:end,:);

X = log(good_outputs_X(1:end-Nval,:));
Y = log(good_outputs_Y(1:end-Nval,:));
%% Construction of the regression model
disp('▣ Fitting...')
good_count = Nfit;

% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
    PLS_nipals(X,Y,rank(X));

% Goodness of fit - R^2
% Fraction of the variance in the dependent variable which is explained by the model
% Calculate agreement of values predicted by regression (Yhat = Bpls*X) with original outputs (Y)
% Assessment on log-transformed values
SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
R2each = SSYR./SSYT;
avg_R2_fit = mean(R2each);

% Assessment on (normal) values
R2ord_fit = zeros(1,N_outputs_Y);
R2adj_fit = zeros(1,N_outputs_Y);
rmse_fit = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    mdl = fitlm(exp(Y(:,i)),exp(Yhat(:,i)));
    R2ord_fit(i) = mdl.Rsquared.Ordinary;
    R2adj_fit(i) = mdl.Rsquared.Adjusted;
    rmse_fit(i) = mdl.RMSE;
end
R2ord_fit;
R2adj_fit;

R2_fit = R2adj_fit; % Values plotted in figures
avg_R2calc_fit = mean(R2_fit);

% Residual Standard Deviation
% Standard deviation for a normal distribution, centered on the predicted regression line,
% representing the distribution of actually observed values
oSD = zeros(1,N_outputs_Y);
rSD = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    oSD(i) = std(exp(Y(:,i)));
    rSD(i) = sqrt(sum((exp(Yhat(:,i)) - exp(Y(:,i)) ).^2) / (good_count-2));
end

% Plot regression coefficients
if plot_matrix == 1
    figure; set(gcf,'color','w')
    imagesc(Bpls'); colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'YDir','normal')
    title('Regression coefficients');
    if male2female == 1
        xlabel('Outputs Male');
        ylabel('Outputs Female');
    else
        xlabel('Outputs Female');
        ylabel('Outputs Male');
    end
    set(gca,'XTick',(1:N_outputs_X))
    set(gca,'XTickLabel',output_names_X)
    set(gca,'YTickLabel',output_names_Y)
    set(gca,'YTick',(1:N_outputs_Y))
    rotateXLabels(gca(), 90)
    colorbar
end

% Scatter plots
N_figures = ceil(N_outputs_Y/10);

if plot_fitting == 1
    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(2,5,i),hold on
                % Plot data points
                factor = 1;
                plot(factor*exp(Y(:,dex1)),factor*exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color)
                xlabel(['Actual ',output_names_Y{dex1}])
                ylabel(['Predicted ',output_names_Y{dex1}])
                title(['R^2 = ',num2str(R2_fit(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',12)
                % Plot identity line
                ylim_ind = get(gca,'ylim') ;
                xlim_ind = get(gca,'xlim') ;
                minpoint = min([ylim_ind(1),xlim_ind(1)]);
                maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
                xlim([minpoint, maxpoint])
                ylim([minpoint, maxpoint])
                dex1 = dex1+1;
            end
        end
    end

    % The error terms are assumed to be:
    % 1) Normally distribuited
    % 2) Homoscedatic (same variance at every X)
    % 3) Independent
    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(2,5,i),hold on
                % Plot data points
                factor = 1;
                plot(factor*exp(Yhat(:,dex1)),(exp(Y(:,dex1))-exp(Yhat(:,dex1)))/rSD(dex1),'Marker','o','LineStyle','none','Color',color);
                xlabel(['Predicted ', output_names_Y{dex1}])
                ylabel('Standardized residuals')
                title(['rSD/oSD = ',num2str(rSD(dex1)/oSD(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',12)
                % Plot identity line
                xlim_ind = get(gca,'xlim') ;
                plot([xlim_ind(1), xlim_ind(2)],[0,0],'--k')
                xlim([xlim_ind(1), xlim_ind(2)])

                dex1 = dex1+1;
            end
        end
    end

    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    subplot(2,2,1),bar(R2_fit,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    ylim([0 1])
    rotateXLabels( gca(), 90)
    title('R^2 values')

    subplot(2,2,2),bar(rSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD values')

    subplot(2,2,4),bar(oSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('oSD values')

    subplot(2,2,3),bar(rSD./oSD,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD/oSD values')
end
%% Validation of the regression model

disp('▣ Validation...')
N_validation = Nval;

predicted_outputs = actual_outputs*0;
for i = 1:Nval
    cell_index = i;
    inputs = actual_inputs(cell_index,:);
    x = log(inputs);
    xz = (x-mean(X))./std(X);

    yz = xz*Bpls;
    %    yz = xz;
    y = yz.*std(Y)+mean(Y);
    predicted_outputs(i,:) = exp(y);
end

% with R^2 = 0.75, the model explains approximately 75% of the variability in the predicted variable
R2ord = zeros(1,N_outputs_Y);
R2adj = zeros(1,N_outputs_Y);
rmse_val = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    mdl = fitlm(actual_outputs(:,i),predicted_outputs(:,i));
    R2ord(i) = mdl.Rsquared.Ordinary;
    R2adj(i) = mdl.Rsquared.Adjusted;
    rmse_val(i) = mdl.RMSE;
end
R2ord;
R2adj;
R2 = R2adj;
avg_R2_val = mean(R2);

% Residual Standard Deviation
% Standard deviation for a normal distribution, centered on the predicted regression line,
% representing the distribution of actually observed values
oSD_val = zeros(1,N_outputs_Y);
rSD_val = zeros(1,N_outputs_Y);
for i = 1:N_outputs_Y
    oSD_val(i) = std(actual_outputs(:,i));
    rSD_val(i) = sqrt(sum((predicted_outputs(:,i) - actual_outputs(:,i) ).^2) / (N_validation-2));
end

if plot_validation == 1
    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(2,5,i),hold on
                % Plot data points
                factor = 1;
                %if (output_index == 1) && (i == 6 || i == 7), factor = 1e6; end
                %if (output_index == 2) && (i == 3 || i == 4), factor = 1e6; end
                plot(factor*actual_outputs(:,dex1),factor*predicted_outputs(:,dex1),'Marker','o','LineStyle','none')%,'Color',color)
                xlabel(['Actual ',output_names_Y{dex1}])
                ylabel(['Predicted ',output_names_Y{dex1}])
                title(['R^2 = ',num2str(R2(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',14)
                % Plot identity line
                ylim_ind = get(gca,'ylim') ;
                xlim_ind = get(gca,'xlim') ;
                minpoint = min([ylim_ind(1),xlim_ind(1)]);
                maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
                xlim([minpoint, maxpoint])
                ylim([minpoint, maxpoint])
                dex1 = dex1+1;
            end
        end
    end

    dex1 = 1;
    for figdex = 1:N_figures
        figure
        set(gcf,'color','w','Position',[50,100,1500,750])
        for i = 1:10
            if dex1 <= N_outputs_Y
                subplot(2,5,i),hold on
                % Plot data points
                factor = 1;
                plot(factor*predicted_outputs(:,dex1),(actual_outputs(:,dex1)-predicted_outputs(:,dex1))/rSD_val(dex1),'Marker','o','LineStyle','none');
                xlabel(['Predicted ', output_names_Y{dex1}])
                ylabel('Standardized residuals')
                title(['rSD/oSD = ',num2str(rSD_val(dex1)/oSD_val(dex1),4)])
                set(gca,'box','off','tickdir','out','fontsize',14)
                % Plot identity line
                xlim_ind = get(gca,'xlim') ;
                plot([xlim_ind(1), xlim_ind(2)],[0,0],'--k')
                xlim([xlim_ind(1), xlim_ind(2)])

                dex1 = dex1+1;
            end
        end
    end

    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    subplot(2,2,1),bar(R2)%,'FaceColor',color)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    ylim([0 1])
    rotateXLabels(gca(), 90)
    title('R^2 values')

    % Add residuals also in this population
    subplot(2,2,2),bar(rSD_val)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD values')

    subplot(2,2,4),bar(oSD_val)
    set(gca,'box','off','tickdir','out','fontsize',14)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('oSD values')

    subplot(2,2,3),bar(rSD_val./oSD_val)
    set(gca,'box','off','tickdir','out','fontsize',12)
    set(gca,'XTick',1:N_outputs_Y)
    set(gca,'XTickLabel',output_names_Y)
    set(gca,'XLim',[0 N_outputs_Y+1])
    %ylim([0 1])
    rotateXLabels( gca(), 90)
    title('rSD/oSD values')

end
%% Application of cross-sex translator to predict drug effects

disp('▣ Application to simulated drug:')

% Select drug index
drug_index = testdrugs;

load drug_data
disp(drug_list{drug_index})

if drug_conc_index == 1
    disp('Concentration: 1x[ETPC]')
elseif drug_conc_index == 2
    disp('Concentration: 2x[ETPC]')
elseif drug_conc_index == 3    
    disp('Concentration: 3x[ETPC]')
elseif drug_conc_index == 4
    disp('Concentration: 4x[ETPC]')
end

if frequency_selection == 0.5
    if drug_conc_index == 1
        drug_matrix = drug_1x;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 2
        drug_matrix = drug_2x;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 3
        drug_matrix = drug_3x;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 4
        drug_matrix = drug_4x;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.2000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_f = output_matrix;
    end
elseif frequency_selection == 1
    if drug_conc_index == 1
        drug_matrix = drug_1x;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 2
        drug_matrix = drug_2x;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 3
        drug_matrix = drug_3x;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 4
        drug_matrix = drug_4x;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.1.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.1000.CT.1.GT.2.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_f = output_matrix;
    end
elseif frequency_selection == 2
    if drug_conc_index == 1
        drug_matrix = drug_1x;

        load('Drugs/APCAT.BCL.500.CT.1.GT.1.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.500.CT.1.GT.2.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 2
        drug_matrix = drug_2x;

        load('Drugs/APCAT.BCL.500.CT.1.GT.1.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.500.CT.1.GT.2.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 3
        drug_matrix = drug_3x;

        load('Drugs/APCAT.BCL.500.CT.1.GT.1.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.500.CT.1.GT.2.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 4
        drug_matrix = drug_4x;

        load('Drugs/APCAT.BCL.500.CT.1.GT.1.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.500.CT.1.GT.2.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_f = output_matrix;
    end
elseif frequency_selection == 3
    if drug_conc_index == 1
        drug_matrix = drug_1x;

        load('Drugs/APCAT.BCL.333.CT.1.GT.1.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.333.CT.1.GT.2.HL.0.Cmax.0.ALL.drug1x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 2
        drug_matrix = drug_2x;

        load('Drugs/APCAT.BCL.333.CT.1.GT.1.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.333.CT.1.GT.2.HL.0.Cmax.0.ALL.drug2x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 3
        drug_matrix = drug_3x;

        load('Drugs/APCAT.BCL.333.CT.1.GT.1.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.333.CT.1.GT.2.HL.0.Cmax.0.ALL.drug3x.mat')
        all_outputs_drug_f = output_matrix;
    elseif drug_conc_index == 4
        drug_matrix = drug_4x;

        load('Drugs/APCAT.BCL.333.CT.1.GT.1.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_m = output_matrix;

        load('Drugs/APCAT.BCL.333.CT.1.GT.2.HL.0.Cmax.0.ALL.drug4x.mat')
        all_outputs_drug_f = output_matrix;
    end
end

% Convert cytosolic Ca from mM to nM
all_outputs_drug_m(:,9:11) = all_outputs_drug_m(:,9:11)*1e6;
all_outputs_drug_f(:,9:11) = all_outputs_drug_f(:,9:11)*1e6;
output_units(9:11) = {'nM','nM','nM'};
% Convert CaTint from mM*ms to nM*ms
all_outputs_drug_m(:,21) = all_outputs_drug_m(:,21)*1e6;
all_outputs_drug_f(:,21) = all_outputs_drug_f(:,21)*1e6;
output_units(21) = {'ms*nM'};
% Convert INaLint into absolute value
all_outputs_drug_m(:,22) = abs(all_outputs_drug_m(:,22));
all_outputs_drug_f(:,22) = abs(all_outputs_drug_f(:,22));

if sum(all_outputs_drug_m(drug_index,output_selection_X)) == 0
    disp('Drug induces instabilities in male model!')
    instabilityindex(1) = 1;
else
    instabilityindex(1) = 0;
end

if sum(all_outputs_drug_f(drug_index,output_selection_X)) == 0
    disp('Drug induces instabilities in female model!')
    instabilityindex(2) = 1;
end

if sum(all_outputs_drug_m(drug_index,19)) > 0
    disp('Drug induces SCR in male model!')
end

if sum(all_outputs_drug_m(drug_index,20)) > 0
    disp('Drug induces EAD in male model!')
end

if sum(all_outputs_drug_f(drug_index,19)) > 0
    disp('Drug induces SCR in female model!')
end

if sum(all_outputs_drug_f(drug_index,20)) > 0
    disp('Drug induces EAD in female model!')
end

if flag_plot_application == 1
    figure,hold on,set(gcf,'color','w')%,'Position',[50,100,1500,750])
    b = bar([all_outputs_drug_m(1,output_selection_X);all_outputs_drug_m(drug_index,output_selection_X);...
        all_outputs_drug_f(1,output_selection_X);all_outputs_drug_f(drug_index,output_selection_X)]','FaceColor',[1 1 1]);
    set(b(1),'FaceColor',color_m, EdgeColor='none'); set(b(3),'FaceColor', color_f, EdgeColor='none');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Outputs'),ylabel('Values')
    xticks(1:length(actual_outputs));
    xticklabels(output_names_Y)
    title(drug_list{drug_index})
    legend('Male Ctrl','Male Drug','Female Ctrl','Female Drug')
end

if male2female == 1
    actual_inputs = all_outputs_drug_m(drug_index,output_selection_X);
    actual_outputs = all_outputs_drug_f(drug_index,output_selection_Y);
else
    actual_inputs = all_outputs_drug_f(drug_index,output_selection_X);
    actual_outputs = all_outputs_drug_m(drug_index,output_selection_Y);
end

inputs = actual_inputs;
x = log(inputs);
xz = (x-mean(X))./std(X);

yz = xz*Bpls;
y = yz.*std(Y)+mean(Y);
predicted_outputs = exp(y);

if flag_plot_application == 1
    figure,hold on,set(gcf,'color','w')%,'Position',[50,100,1500,750])
    b = bar([actual_inputs;predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
    set(b(1),'FaceColor',color_X, EdgeColor='none'); set(b(3),'FaceColor',color_Y, EdgeColor='none');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Outputs'),ylabel('Values')
    xticks(1:length(actual_outputs));
    xticklabels(output_names_Y)

    title(drug_list{drug_index})
    if male2female == 1
        legend('Male','M-to-F','Female')
    else
        legend('Female','F-to-M','Male')
    end

end
%% Assess quality control of translation

% Ratio between predicted and actual values
quality_ratio = zeros(1,length(output_selection_Y));
for ii = 1:length(output_selection_Y)
    quality_ratio(ii) = predicted_outputs(ii)/actual_outputs(ii);
end
% Quality index
quality_index = abs(1-quality_ratio);

figure
set(gcf,'color','w')
set(gca,'box','off','tickdir','out','fontsize',14)
xvalues = output_names_Y;
yvalues = drug_list(testdrugs);
h = heatmap(xvalues,yvalues,quality_ratio);
title("Ratio predicted/actual @ "+frequency_selection+"-Hz, " + drug_conc_index + "x");
xlabel('Outputs');
ylabel('Drugs');
colorbar

figure
set(gcf,'color','w')
set(gca,'box','off','tickdir','out','fontsize',14)
xvalues = output_names_Y;
yvalues = drug_list(testdrugs);
h = heatmap(xvalues,yvalues,quality_index);
title("abs(1 - predicted/actual) @ "+frequency_selection+"-Hz, " + drug_conc_index + "x");
xlabel('Outputs');
ylabel('Drugs');
colorbar
%% Drug index list

% Index Drug name
%  1	Control
%  2	Ajmaline
%  3	Amiodarone_1
%  4	Amiodarone_2
%  5	Amitriptyline
%  6	Astemizole
%  7	Bepridil_1
%  8	Bepridil_2
%  9	Bepridil_CiPA
% 10	Ceftriaxone
% 11	Chlorpromazine_1
% 12	Chlorpromazine_2
% 13	Chlorpromazine_CiPA
% 14	Cibenzoline
% 15	Cilostazol
% 16	Cisapride_1
% 17	Cisapride_2
% 18	Cisapride_CiPA
% 19	Clozapine
% 20	Dasatinib
% 21	Desipramine
% 22	Diazepam
% 23	Diltiazem_1
% 24	Diltiazem_2
% 25	Diltiazem_CiPA
% 26	Diphenhydramine
% 27	Disopyramide
% 28	Dofetilide_1
% 29	Dofetilide_2
% 30	Dofetilide_CiPA
% 31	Donepezil
% 32	Droperidol
% 33	Duloxetine
% 34	Flecainide
% 35	Fluvoxamine
% 36	Halofantrine
% 37	Haloperidol_1
% 38	Haloperidol_2
% 39	Ibutilide
% 40	Imipramine
% 41	Lamivudine
% 42	Linezolid
% 43	Loratadine
% 44	Methadone
% 45	Metronidazole
% 46	Mexiletine
% 47	Mexiletine_CiPA
% 48	Mibefradil_1
% 49	Mibefradil_2
% 50	Mitoxantrone
% 51	Moxifloxacin
% 52	Nifedipine_1
% 53	Nifedipine_2
% 54	Nilotinib
% 55	Nitrendipine_1
% 56	Nitrendipine_2
% 57	Ondansetron_CiPA
% 58	Paliperidone
% 59	Paroxetine
% 60	Pentobarbital
% 61	Phenytoin_1
% 62	Phenytoin_2
% 63	Pimozide_1
% 64	Pimozide_2
% 65	Piperacillin
% 66	Prenylamine
% 67	Procainamide
% 68	Propafenone
% 69	Propranolol
% 70	Quetiapine
% 71	Quinidine_1
% 72	Quinidine_2
% 73	Quinidine_CiPA
% 74	Raltegravir
% 75	Ranolazine_CiPA
% 76	Ribavirin
% 77	Risperidone_1
% 78	Risperidone_2
% 79	Saquinavir
% 80	Sertindole_1
% 81	Sertindole_2
% 82	Sitagliptin
% 83	Solifenacin
% 84	Sotalol
% 85	Sotalol_CiPA
% 86	Sparfloxacin
% 87	Sunitinib
% 88	Tedisamil
% 89	Telbivudine
% 90	Terfenadine_1
% 91	Terfenadine_2
% 92	Terfenadine_CiPA
% 93	Terodiline
% 94	Thioridazine_1
% 95	Thioridazine_2
% 96	Verapamil_1
% 97	Verapamil_2
% 98	Verapamil_CiPA
% 99	Voriconazole
