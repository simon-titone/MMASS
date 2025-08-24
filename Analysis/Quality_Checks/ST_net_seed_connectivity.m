function ST_net_seed_connectivity(source_filename,options_seed,n)
NET_folder = net('path');

% define paths
% ------------
ddx=fileparts(fileparts(source_filename))
[dd,ff,ext] = fileparts(source_filename);
dd2 = [dd filesep 'seed_connectivity'];

if ~isdir(dd2)
    mkdir(dd2);  % Create the output folder if it doesn't exist..
end

% load data
% ---------
load(source_filename,'source');
deformation_to_subj=[ddx filesep 'mr_data' filesep 'y_anatomy_prepro.nii'];
deformation_to_mni=[ddx filesep 'mr_data' filesep 'iy_anatomy_prepro.nii'];
load(n.spath);

% define parameters (channels, ntp, nvox)
% ---------------------------------------
nchan=size(source.sensor_data,1);

Fs = 1/(source.time(2)-source.time(1));
new_Fs = options_seed.fs;
frequencies = [1:1:80];
if not(new_Fs==Fs)
    channel_data = resample(source.sensor_data',new_Fs,Fs)';
    t=resample(source.time',new_Fs,Fs)';
    Ntime  = length(t);
else
    channel_data = source.sensor_data;
    t=source.time;
    Ntime  = length(t);
end

vox_indices = find(source.inside==1);
nvoxels = length(vox_indices);
xdim    = source.dim(1);
ydim    = source.dim(2);
zdim    = source.dim(3);
xyz = source.pos(vox_indices,:)';

% load and project seeds in individual space
% -----------------------------------------
radius=6; % in mm
seed_file=[NET_folder filesep 'template' filesep 'seeds' filesep options_seed.seed_file '.mat'];
load(seed_file,'seed_info');
nrois=length(seed_info);
for i=1:nrois
    seed_info(i).coord_subj=net_project_coord(deformation_to_subj,seed_info(i).coord_mni);
    dist = pdist2(xyz',seed_info(i).coord_subj);
    voxel_list=find(dist<radius);
    if isempty(voxel_list)
        disp(['problem with seed ' num2str(i) '!']);
        [~,voxel_list]=min(dist);
    end
    seed_info(i).seedindx=voxel_list;
end

% mapping
%  -------
% if strcmp(options_seed.map_enable,'on')
    switch options_seed.connectivity_measure
        case {'blpc_spec','blpc_ss'}
            
            winsize  = round(new_Fs*options_seed.window_duration);
            overlap  = round(new_Fs*options_seed.window_overlap);
            [~, F, T] = spectrogram(t, winsize, overlap, frequencies, new_Fs);
            FT_all = zeros(nvoxels,length(F),length(T));
            FT_roi = zeros(nrois,length(F),length(T));
            
            if strcmp(options_seed.connectivity_measure,'blpc_spec')
                % reconstruct source activities
                % ----------------------------
                brain = zeros(nvoxels,Ntime);
                for k = 1:nvoxels
                    brain(k,:) = source.pca_projection(k,(3*k-2):(3*k))*source.imagingkernel((3*k-2):(3*k),:)*channel_data;
                end
                
                [coeff,score] = pca(brain','Centered',false);
                fpc = score(:,1)*coeff(:,1)'; % for regression
                brain = brain - fpc';
                clear fpc coeff score
                
                for k = 1:nvoxels % for normalization
                    tmp = brain(k,:);
                    brain(k,:) = tmp./std(tmp);
                end
                
                % time-frequency decomposition
                % ----------------------------
                %                 for k=1:nvoxels
                %                     disp([num2str(k) ' / ' num2str(nvoxels)]);
                %                     brain_signal = brain(k,:);
                %                     FT_all(k,:,:) = spectrogram(brain_signal, winsize, overlap, frequencies, new_Fs);
                %                 end
                
                for k=1:nvoxels
%                     disp([num2str(k) ' / ' num2str(nvoxels)]);
                    brain_signal = brain(k,:);
                    FT_all(k,:,:) = spectrogram(brain_signal, winsize, overlap, frequencies, new_Fs);
                end
%                 save('FT_all.mat','FT_all');
                
                for k = 1:nrois
                    seedindx = seed_info(k).seedindx;
                    
                    nv      = length(seedindx);
                    mat_sig = zeros(nv,Ntime);
                    for j = 1:nv
                        q = seedindx(j);
                        mat_sig(j,:) = brain(q,:);
                    end
                    coeff = pca(mat_sig');
                    coeffx = inv(coeff');
                    brain_signal   = coeffx(:,1)'*mat_sig; %first PC
                    [~,~,~,p]= spectrogram(brain_signal, winsize, overlap, frequencies, new_Fs);
                    FT_roi(k,:,:) = p;
                    save('FT_roi.mat','FT_roi');
                end
            end
    end
    
    %     % othogonalization and connectivity
    %     % ------------
    %     for k = 1:nrois
    %         corr_map = zeros(nvoxels,length(frequencies));
    %         Si = squeeze(FT_roi(k,:,:))';
    %         for w = 1:nvoxels
    %             if ismember(w,seed_info(k).seedindx) % w belong to the voxel_list corresponding to the k-th seed, ie is in the ROI around the seed
    %                 corr_map(w,:) = nan;
    %             else
    %                 Sj = squeeze(FT_all(w,:,:))';
    %                 corr_map(w,:) = net_blp_corr(Si, Sj,options_seed.orthogonalize);
    %             end
    %         end
    %         corr_map( isnan(corr_map) ) = max(corr_map(:)) - nanstd(corr_map(:));
    %
    %         nbands = length(seed_info(k).frequency);
    %         for zz=1:nbands
    %             band=seed_info(k).frequency{zz};
    %             vect_f=(frequencies>=band(1) & frequencies<=band(2));
    %             seed_map = mean(corr_map(:,vect_f),2)';
    %
    %             Vt.dim      = source.dim;
    %             Vt.pinfo    = [0.000001 ; 0 ; 0];
    %             Vt.dt       = [16 0];
    %             Vt.fname    = [dd2 filesep 'seed_' seed_info(k).label '_(' num2str(band(1)) '-' num2str(band(2)) 'Hz)_maps.nii'];
    %             Vt.mat      = net_pos2transform(source.pos, source.dim);
    %             Vt.n        = [1 1];
    %             image       = seed_map;
    %             seed_image  = zeros(xdim*ydim*zdim,1);
    %             seed_image(vox_indices) = image;
    %             seed_image  = reshape(seed_image, xdim, ydim, zdim);
    %             spm_write_vol(Vt, seed_image);
    %         end
    %     end
    %     net_warp([dd2 filesep 'seed*maps.nii'],deformation_to_mni);
    % end
    
    
    %% matrix
    %  ------
    % if strcmp(options_seed.matrix_enable,'on')
%     switch options_seed.connectivity_measure
%         case {'blpc_spec','blpc_ss'}
%             winsize  = round(new_Fs*options_seed.window_duration);
%             overlap  = round(new_Fs*options_seed.window_overlap);
%             [~, F, T] = spectrogram(t, winsize, overlap, frequencies, new_Fs);
%             FT_roi = zeros(nrois,length(F),length(T));
%             
%             if strcmp(options_seed.connectivity_measure,'blpc_spec')
%                 % reconstruct source activities
%                 % ----------------------------
%                 brain = zeros(nvoxels,Ntime);
%                 for k = 1:nvoxels
%                     brain(k,:) = source.pca_projection(k,(3*k-2):(3*k))*source.imagingkernel((3*k-2):(3*k),:)*channel_data;
%                 end
%                 
%                 [coeff,score] = pca(brain','Centered',false);
%                 fpc = score(:,1)*coeff(:,1)'; % for regression
%                 brain = brain - fpc';
%                 clear fpc coeff score
%                 
%                 for k = 1:nvoxels % for normalization
%                     tmp = brain(k,:);
%                     brain(k,:) = tmp./std(tmp);
%                 end
%                 
%                 % time-frequency decomposition
%                 % ----------------------------
%                 for k = 1:nrois
%                     seedindx = seed_info(k).seedindx;
%                     
%                     nv      = length(seedindx);
%                     mat_sig = zeros(nv,Ntime);
%                     for j = 1:nv
%                         q = seedindx(j);
%                         mat_sig(j,:) = brain(q,:);
%                     end
%                     coeff = pca(mat_sig');
%                     coeffx = inv(coeff');
%                     brain_signal   = coeffx(:,1)'*mat_sig; %first PC
%                     FT_roi(k,:,:)  = spectrogram(brain_signal, winsize, overlap, frequencies, new_Fs);
%                 end
%             end
%     end
%     save('ft_roi', 'FT_roi');
    % othogonalization and connectivity
    % ------------
    %     corr_matrix = zeros(nrois,nrois,length(frequencies));
    %     for k = 1:nrois
    %         Si = squeeze(FT_roi(k,:,:))';
    %         for w = k+1:nrois
    %             Sj = squeeze(FT_roi(w,:,:))';
    %             corr_matrix(k,w,:) = net_blp_corr(Si, Sj,options_seed.orthogonalize);
    %             corr_matrix(w,k,:) = corr_matrix(k,w,:);
    %         end
    %     end
    %     save([dd2 filesep 'matrix_connectivity.mat'],'corr_matrix','F','T','options_seed','seed_info','-v7.3');
end
% end

% last updated 10.10.2019, by JS