%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to re-create the results of                                        %
% Sebastian Zaunseder, Antti Vehkaoja, Vincent Fleischhauer,              %
% Christoph Hoog Antink                                                   %
%                                                                         %
% Signal-to-noise ratio is more important than sampling rate in           %
% beat-to-beat interval estimation from optical sensors                   %
%                                                                         %
% Biomedical Signal Processing and Control, Volume 74, 2022, 103538,      %
% ISSN 1746-8094,                                                         %
% https://doi.org/10.1016/j.bspc.2022.103538                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add folders to the path to read data (just have "PulseDecompositionAnalysis" in a subfolder
addpath(genpath('PulseDecompositionAnalysis'));%create and add the subfolder(s)
sourceFolder ='PulseDecompositionAnalysis\Comparison\basicModels\Parameter\';%here some parameters are stored
algorithm='GammaGauss2Mode';%decomposition algorihtm (likely to be fixed)
algorithmType='GammaGaussian';%algorithm name (without kernel number)
numberOfKernels=2;%specify number of kernels (must match algorithm)

DO_PLOT = true;     % Flag to generate the plots

f_s_0 = 10000;      % Ground truth sampling rate f^{\uparrow}_{\rm s}
f_s_target = 1000;  % Sampling rate used for up-sampling f_{\rm s_{est}} 
 
n_rand = 100000;    % Maximum number of iterations
th_NSEM = 0.01;     % Threshold on NSEM
 
t_sig = 4;          % Length of the signal (s)

f_s_query = round(logspace(log(5)/log(10),log(50)/log(10),10));
    % array of all f^{\downarrow}_{\rm s}
SNR_query = -3:3:30;
    % array of all SNR

n_iter = []; % array for the actual number of iterations
pulse_examples = cell(4,size(SNR_query,1),size(f_s_query,1)); % cell for example plot of the signals

for class = 1:4 %here 1...4 can be selected to represent the dawber classes --> varied in the loop (or executed 4 times)

    xx=load([sourceFolder  algorithm '\' num2str(class) '.mat']);%load parameters from previosu decomposition
    temporalShift = 0.7;%shift in seconds each template should be shifted to be approx. in the mean of x; the concrete value should not matter, just leave it
    xx.opt_values(2:3:end)=xx.opt_values(2:3:end)+temporalShift;%add shift(seconds) to the mean distribution value to cause shift
    decomposeFuncName = ['decompose' algorithm];
    currentVals=xx.opt_values;
    [pulse_proto,~,~] = decompositionAlgorithm(zeros(1,20001),f_s_0, ...
        'numKernels', numberOfKernels,...
        'kernelTypes', algorithmType,...
        'noOpt', true,...
        'initialValues',currentVals);
    noiseType = 'pink';

    % Arrays to store data for evluation / plotting
    FOM       = zeros(100,100,4);     
    noise_arr = zeros(100,1);   
    snr_arr   = zeros(100,1);     
    fs_arr    = zeros(100,1);   

    noise_indx = 0;    
    for noise_level = SNR_query % Sweep over different SNR

        noise_indx = noise_indx + 1;
        noise_arr(noise_indx) = noise_level;

        fs_indx = 0;
        for f_s_intermediate = f_s_query   % Sweep over different f^{\downarrow}_{\rm s}

            fs_indx = fs_indx + 1;
            fs_arr(fs_indx) = f_s_intermediate;

            delta_accu = zeros(n_rand,2);   % Array storing ground-truth displacements \Delta and estimated displacements \Delta_{\rm est}
            
            for i = 1:n_rand
                delta_pulse = randn(1) / 10;    % Draw random displacement \Delta
                delta_accu(i,1) = delta_pulse;

                sig_0 = zeros(t_sig * f_s_0,1); % Generation of the first pulse
                sig_0(t_sig/2 * f_s_0) = 1;
                sig_0 = conv(sig_0,pulse_proto,'same'); 

                sig_1 = zeros(t_sig * f_s_0,1); % Generation of the second pulse
                sig_1(round((t_sig/2+delta_pulse) * f_s_0)) = 1;
                sig_1 = conv(sig_1,pulse_proto,'same');

                sig_0_intermediate = resample(sig_0,f_s_intermediate,f_s_0); % Downsamplit to f^{\downarrow}_{\rm s}
                sig_1_intermediate = resample(sig_1,f_s_intermediate,f_s_0);

                sig_0_intermediate_n = noise_ppg(sig_0_intermediate',noiseType,noise_level,f_s_intermediate)'; % Coruption with noise at level SNR
                sig_1_intermediate_n = noise_ppg(sig_1_intermediate',noiseType,noise_level,f_s_intermediate)';

                sig_0_resampled = resample(sig_0_intermediate_n,f_s_target,f_s_intermediate); % Upsampling to f_{\rm s_{est}}
                sig_1_resampled = resample(sig_1_intermediate_n,f_s_target,f_s_intermediate);

                delta_est = finddelay(sig_0_resampled,sig_1_resampled,f_s_target*1)/f_s_target; % Estimating the delay \Delta_{\rm est}

                if DO_PLOT
                    if i==1 % Generate one plot for each instance of the sweep
                        t_0 = ((1:size(sig_0,1))'-1)./f_s_0;
                        t_i = ((1:size(sig_0_intermediate,1))'-1)./f_s_intermediate;
                        t_r = ((1:size(sig_1_resampled,1))'-1)./f_s_target;

                        figure(1)
                        clf
                        subplot(2,2,1)
                        plot(t_0,sig_0)
                        hold on
                        plot(t_0,sig_1)
                        title(['Ground Truth, $\Delta$ = ' num2str(delta_pulse*1000) '\,ms'],'interpreter','latex')
                        yticks([])
                        ylabel('a.u.','interpreter','latex')
                        xlabel('Time (s)','interpreter','latex')

                        subplot(2,2,2)
                        plot(t_i,sig_0_intermediate)
                        hold on                
                        plot(t_i,sig_1_intermediate)
                        title(['Downsampled, $f^{\downarrow}_{\rm s}$ = ' num2str(f_s_intermediate) '\,Hz'],'interpreter','latex')
                        yticks([])
                        ylabel('a.u.','interpreter','latex')
                        xlabel('Time (s)','interpreter','latex')

                        subplot(2,2,3)
                        plot(t_i,sig_0_intermediate_n)
                        hold on
                        plot(t_i,sig_1_intermediate_n)
                        title(['Added Noise, SNR = ' num2str(noise_level) '\,dB'],'interpreter','latex')
                        yticks([])
                        ylabel('a.u.','interpreter','latex')
                        xlabel('Time (s)','interpreter','latex')

                        pulse_examples{class,noise_indx,fs_indx} = sig_0_intermediate_n;

                        subplot(2,2,4)
                        plot(t_r,sig_0_resampled)
                        hold on
                        plot(t_r,sig_1_resampled)
                        title(['Resampled, $f_{\rm s_{est}}$ = ' num2str(f_s_target) '\,Hz, $\Delta_{\rm est}$ = ' num2str(delta_est*1000) '\,ms'],'interpreter','latex')
                        yticks([])
                        ylabel('a.u.','interpreter','latex')
                        xlabel('Time (s)','interpreter','latex')
                        pause(.1)
                    end
                end

                delta_accu(i,2) = delta_est;

                NSEM = std(abs(delta_accu(1:i,2) - delta_accu(1:i,1)))/mean(abs(delta_accu(1:i,2) - delta_accu(1:i,1)))/sqrt(i);
                if (i>2)&&(NSEM < th_NSEM) % Condition for termination
                    n_iter = [n_iter;i]; %#ok<AGROW> % Store actual number of iterations
                    break
                end 
            end
            delta_accu = delta_accu(1:i,:);

            FOM(noise_indx,fs_indx,1) = rms(delta_accu(:,1)-delta_accu(:,2))*1000;               % RMSE
            FOM(noise_indx,fs_indx,2) = mean(abs((delta_accu(:,1)-delta_accu(:,2))))*1000;       % Mean absolute Error (not used in the paper)
            FOM(noise_indx,fs_indx,3) = prctile(abs((delta_accu(:,1)-delta_accu(:,2))),5)*1000;  % 5th Percentile of difference (not used in the paper)
            FOM(noise_indx,fs_indx,4) = prctile(abs((delta_accu(:,1)-delta_accu(:,2))),95)*1000; % 95th Percentile of difference (not used in the paper)
        end
        
        snr_arr(noise_indx) = noise_level;

    end

    % Clean and store results
    FOM = FOM(1:noise_indx,1:fs_indx,:);
    snr_arr = snr_arr(1:noise_indx);
    fs_arr = fs_arr(1:fs_indx);
    noise_arr = noise_arr(1:noise_indx);    
    res.FOM = FOM;
    res.snr_arr = snr_arr;
    res.fs_arr = fs_arr;
    res.noise_arr = noise_arr;    
    save(['class_' num2str(class) '.mat'],'res');

    if DO_PLOT
        figure(2) % Plot RMSE over f^{\downarrow}_{\rm s} for class "class" (see Fig. 7 top right in the paper)
        clf
        hold on
        le_str = {};
        indx = 0;
        for i = 1:11
            indx = indx + 1;
            plot(fs_arr,FOM(i,:,1))
            le_str{indx} = ['SNR = ' num2str(snr_arr(i))]; %#ok<SAGROW>
        end
        le = legend(le_str);
        set(le,'Interpreter','Latex')
        xlim([5 50])
        xlabel('$f^{\downarrow}_{\rm s}$\,(in Hz)','Interpreter','Latex')
        ylabel('RMSE BBI\,(in ms)','Interpreter','Latex')
        title(['Class ' num2str(class)],'Interpreter','Latex')
    
        figure(3) % Plot RMSE over SNR for class "class"  (see Fig. 7 top left in the paper)
        clf
        hold on
        le_str = {};
        x_ticks = [];
        indx = 0;
        for i = 1:10
            indx = indx + 1;
            plot(noise_arr,FOM(:,i,1))
            le_str{indx} = ['$f^{\downarrow}_{\rm s}$ = ' num2str(fs_arr(i))];     %#ok<SAGROW>
        end
        le = legend(le_str);
        set(le,'Interpreter','Latex','Location','NorthWest')
        ylim([0 50])
        xlabel('SNR\,(in dB)','Interpreter','Latex')
        ylabel('RMSE BBI\,(in ms)','Interpreter','Latex')
        xticks(noise_arr(1:2:11))
        xticklabels(num2str(snr_arr(1:2:11)))
        title(['Class ' num2str(class)],'Interpreter','Latex')
    end    
end
save('n_iter_pink','n_iter') % Store the actual number of iterations

if DO_PLOT
    figure(4) % Give a visual impression of the four Dawber classes at different SNR-levels (see Fig. 3 the paper)
    clf
    t = tiledlayout(4,4);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    for i=1:4
        for j = [1 5 7 12]
            nexttile
            pulse = pulse_examples{i,j,10};
            t_i = ((1:size(pulse,1))'-1)./50;
            plot(t_i,pulse)
            ylim([-1 1.5])
            yticks([])
            if i==1
                title(['SNR = ' num2str(snr_arr(j)) '\,dB'],'interpreter','latex')
            end
            if i==4
                xlabel('Time (s)','interpreter','latex')
            else
                xticks([])
            end
            if j==1
                ylabel(['Class ' num2str(i) ' (a.u.)'],'interpreter','latex')
            end
        end
    end
end
