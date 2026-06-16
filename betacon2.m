%% IF on CG PC
exp_dir = 'F:\My new data';
tool_dir = 'F:\My new data';
toolboxes = {'arfit' 'cartographie_motrice' 'cartographie_motrice/CM_extra' 'fieldtrip-20190205' 'fieldtrip-20190205/external/fastica/' 'c3d2OpenSim' 'libeep-3.3.177-matlab'};
for tool = toolboxes
    addpath(fullfile(tool_dir,'toolboxes',tool{1}));
end

ft_defaults;

% cd(fullfile(exp_dir,'analysis'));
data_dir = fullfile(exp_dir,'eeg');
emg_dir = fullfile(exp_dir,'eeg');
data_dir = exp_dir;
emg_dir = exp_dir;
matfold = fullfile(exp_dir,'matfiles');
picsfold = fullfile(exp_dir,'pics');

conds = {'UNII' 'RB' 'RILB' 'LIRB'};
subjects = [];
for n_sub = 1
    subjects{n_sub} = ['Subject' num2str(n_sub)];
end

LE_short_all = cell(1, length(subjects));
FD_short_all = cell(1, length(subjects));
CD_short_all = cell(1, length(subjects));

LE_med_all   = cell(1, length(subjects));
FD_med_all   = cell(1, length(subjects));
CD_med_all   = cell(1, length(subjects));

LE_long_all  = cell(1, length(subjects));
FD_long_all  = cell(1, length(subjects));
CD_long_all  = cell(1, length(subjects));


for n_sub = 1:length(subjects)
    for n_cond = 1:4
        subfold = fullfile(data_dir,subjects{n_sub});
        if n_cond == 2
        files = dir(fullfile(data_dir,subjects{n_sub},['*_' conds{n_cond} '_ICA.mat']));
            for n_file = 1:length(files)
                %Load data
                ICAfile = fullfile(files(n_file).folder,files(n_file).name);
                load(ICAfile)
                
                if n_file == 1
                    all_data = data;
                else
                    all_data.trial{1} = cat(2,all_data.trial{1},data.trial{1});
                    all_data.sampleinfo(end) = size(all_data.trial{1},2);
                    all_data.time{1} = (0:all_data.sampleinfo(end)-1)/all_data.fsample;
                    data.bad(1:data.fsample/2) = 1;
                    all_data.bad(end-data.fsample/2:end) = 1;
                    all_data.bad = cat(2,all_data.bad,data.bad);
                    for n_trig = 1:length(all_data.pos_trigs)
                        all_data.pos_trigs{n_trig} = [all_data.pos_trigs{n_trig} data.pos_trigs{n_trig}+length(all_data.trig)];
                    end
                    all_data.trig = cat(2,all_data.trig,data.trig);
                end 
            end
            data = all_data;


            Fs = double(data.fsample);
            picksMEEG = 1:61;
            picksforce = [find(strcmp(data.label,'BIP14')) find(strcmp(data.label,'force_cometa'))];
            picksEMG = [find(strcmp(data.label,'Right')) find(strcmp(data.label,'Left'))];
            picksACC = [find(strcmp(data.label,'acc'))];
            acc = data.trial{1}(picksACC,:);
            bip = double(data.trig == 4);
            [mvt_onset, bip_onset, reaction_time] = CM_movement_onset_bip_acc(bip,acc,Fs,1,1);
            mvt_onset(isnan(mvt_onset)) = [];


            t0 = mvt_onset';
            tt = [-1000:2500];
            t = t0+tt;
            t(end,:) = [];
                    
            Fs = double(all_data.fsample);
            picksMEEG = 1:61;

            LMC = [9 10 14 18 19 41 42 45];
            Burst_data = mean(all_data.trial{1}(LMC,:)); 
            F_h_x = cosine_filter(length(Burst_data),{'high' 'low'},[13 30]/Fs*2,[2 5]/Fs*2);
            Burst_data = real(ifft(fft(Burst_data,[],2).*F_h_x,[],2));
            raw_signal = Burst_data; 

            if size(Burst_data,1) > 1
                [U,S,V] = svd(Burst_data,'econ');
                Burst_data = V(:,1)';
            end

            Burst_data = hilbert(Burst_data);
            F_h_x = cosine_filter(length(Burst_data),{'low'},[0.1]/Fs*2,[0.1]/Fs*2);
            Burst_data_slow_env = real(ifft(fft(abs(Burst_data)).*F_h_x));
            Burst_data = Burst_data./Burst_data_slow_env;

            signal = reshape(Burst_data(1,t),size(t));

            M = max(abs(signal),[],2);
            M = (M-median(M))./diff(prctile(M,[25 75]));
            bad_ep = find(M > 5);
            signal(bad_ep,:) = [];            

            beta_mov = abs(signal);
            mean_beta_mov(n_sub,:) = mean(beta_mov);
            th_mov = prctile(beta_mov(:),75);
            burst_mov = beta_mov >= th_mov; % each row is for a movement and time 3500
            Tburst_mov = [];
            Dburst_mov = [];

            for n_ep = 1:size(burst_mov,1)
                tdeb_mov = find(diff([0 burst_mov(n_ep,:)]) == 1);
                tfin_mov = find(diff([burst_mov(n_ep,:) 0]) == -1);
                Tburst_mov = [Tburst_mov tt(round((tfin_mov+tdeb_mov)/2))];
                Dburst_mov = [Dburst_mov tfin_mov-tdeb_mov+1];
            end

            %% New part
            

            % Step 1: Collect valid bursts (≥ 50 samples)
            valid_bursts = struct('ep', {}, 'start', {}, 'endd', {}, 'duration', {});
            for n_ep = 1:size(burst_mov,1)
                tdeb_mov = find(diff([0 burst_mov(n_ep,:)]) == 1);
                tfin_mov = find(diff([burst_mov(n_ep,:) 0]) == -1);

                burst_start_idx = t(n_ep, tdeb_mov);
                burst_end_idx   = t(n_ep, tfin_mov);

                for b = 1:length(burst_start_idx)
                    idx_range = burst_start_idx(b):burst_end_idx(b);
                    if all(idx_range > 0 & idx_range <= length(raw_signal))
                        dur = length(idx_range);
                        if dur >= 50
                            valid_bursts(end+1).ep = n_ep;
                            valid_bursts(end).start = burst_start_idx(b);
                            valid_bursts(end).endd  = burst_end_idx(b);
                            valid_bursts(end).duration = dur;
                        end
                    end
               end
            end

            % Step 2: Categorize by duration
            durations = [valid_bursts.duration];
            short_idx = find(durations < 100);
            med_idx   = find(durations >= 100 & durations < 150);
            long_idx  = find(durations >= 150);

            % Process Short Bursts
            if ~isempty(short_idx)
                short_durations = durations(short_idx);
                desired_len = min(short_durations);
                numBursts = length(short_idx);
                LE_short = zeros(1, numBursts);
                FD_short = zeros(1, numBursts);
                CD_short = zeros(1, numBursts);
   
                for i = 1:length(short_idx)
                    b = valid_bursts(short_idx(i));
                    burst_segment = raw_signal(b.start:b.endd);
                    burst_segment_resampled = resample(burst_segment, desired_len, length(burst_segment));
                    [lyap, fd, cd] = extract_nonlinear_features(burst_segment_resampled,Fs);
                    LE_short(i) = lyap;
                    FD_short(i) = fd;
                    CD_short(i) = cd;
                end
            end

            % Process Medium Bursts
            if ~isempty(med_idx)
                med_durations = durations(med_idx);
                desired_len = min(med_durations);
                numBursts = length(med_idx);
                LE_med = zeros(1, numBursts);
                FD_med = zeros(1, numBursts);
                CD_med = zeros(1, numBursts);
   
                for i = 1:length(med_idx)
                    b = valid_bursts(med_idx(i));
                    burst_segment = raw_signal(b.start:b.endd);
                    burst_segment_resampled = resample(burst_segment, desired_len, length(burst_segment));
                    [lyap, fd, cd] = extract_nonlinear_features(burst_segment_resampled,Fs);
                    LE_med(i) = lyap;
                    FD_med(i) = fd;
                    CD_med(i) = cd;
                end
            end

            % Process Long Bursts
            if ~isempty(long_idx)
                long_durations = durations(long_idx);
                desired_len = min(long_durations);
                numBursts = length(long_idx);
                LE_long = zeros(1, numBursts);
                FD_long = zeros(1, numBursts);
                CD_long = zeros(1, numBursts);
                for i = 1:length(long_idx)
                    b = valid_bursts(long_idx(i));
                    burst_segment = raw_signal(b.start:b.endd);
                    burst_segment_resampled = resample(burst_segment, desired_len, length(burst_segment));
                    [lyap, fd, cd] = extract_nonlinear_features(burst_segment_resampled,Fs);
                    LE_long(i) = lyap;
                    FD_long(i) = fd;
                    CD_long(i) = cd;
                end
            end


            to_rm_mov = find(Dburst_mov < 50);
            Tburst_mov(to_rm_mov) = [];
            Dburst_mov(to_rm_mov) = [];
            Twin = -900:100:2400;
            Dwin = 200;
                                            
            for n_win = 1:length(Twin)
                keep = find(Tburst_mov > Twin(n_win)-Dwin/2 & Tburst_mov < Twin(n_win)+Dwin/2);
                Nburst_win_mov(n_sub,n_win) = length(keep)/size(beta_mov,1);
                Dburst_win_mov(n_sub,n_win) = mean(Dburst_mov(keep));
                SDDburst_win_mov(n_sub,n_win) = std(Dburst_mov(keep));
            end  

            LE_short_all{n_sub} = exist('LE_short', 'var') * LE_short;
            FD_short_all{n_sub} = exist('FD_short', 'var') * FD_short;
            CD_short_all{n_sub} = exist('CD_short', 'var') * CD_short;

            LE_med_all{n_sub}   = exist('LE_med', 'var') * LE_med;
            FD_med_all{n_sub}   = exist('FD_med', 'var') * FD_med;
            CD_med_all{n_sub}   = exist('CD_med', 'var') * CD_med;

            LE_long_all{n_sub}  = exist('LE_long', 'var') * LE_long;
            FD_long_all{n_sub}  = exist('FD_long', 'var') * FD_long;
            CD_long_all{n_sub}  = exist('CD_long', 'var') * CD_long;

            clear LE_* FD_* CD_*

        end 
    end
end
          




