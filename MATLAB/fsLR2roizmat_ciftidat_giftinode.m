function [z, r, tp] = fsLR2roizmat_ciftidat_giftinode(cii, roiL, roiR)
    % DESCRIPTION:
    %   Extract node time series from 32k by 32k fsLR saved in CIFTI dtseries. 
    %   Outputs correlation matrix (r), fishers-z-transfromed corr matrix (z),
    %   and node-by-timepoint matrix (tp). 
    %   
    %   Nodes are L & R GIFTI files.
    %
    %   If running in MATLAB 2018+ addpath GIFTIv1_8
    
    % --- Addpaths should be moved outside of this script ----
    % addpath cifti-matlab
    % addpath gifti_v1_8
    % ----------% 
    
    
    % % Load ROIs
    if ischar(roiL) % check if roi is a file name or gifti object
        roiL = gifti(roiL); % load gifti if it is a file name
        roiR = gifti(roiR);
    end
    roiL = roiL.cdata;
    roiR = roiR.cdata;
    
    
    % ---- Load CIFTI data -----
    if ischar(cii) % check if roi is a file name or gifti object
        cii = cifti_read(cii); 
    end

    cL=cifti_struct_dense_extract_surface_data(cii, 'CORTEX_LEFT');
    cR=cifti_struct_dense_extract_surface_data(cii, 'CORTEX_RIGHT');
    
    
    % --- ---------------- ----
    
    if(size(roiL,2)>1)
        % If ROIs are defined by separate columns
        
        % Left Hemisphere    
        clm = zeros(size(roiL,2),size(cL,2)); %Left
        for j = 1:size(roiL,2) 
            clm(j,:) = mean(cL(roiL(:,j)~=0, :));
        end
        
        % Right Hemisphere    
        crm = zeros(size(roiR,2),size(cR,2)); %Right       
        for k = 1:size(roiR,2)
            crm(k,:) = mean(cR(roiR(:,k)~=0, :));
        end
    else
        % If ROIs are defined by a single columns of different number index 
        uroiL = unique(roiL);
        uroiL = uroiL(uroiL~=0);    

        uroiR = unique(roiR);
        uroiR = uroiR(uroiR~=0);    

        % Left Hemisphere
        clm = zeros(length(uroiL),size(cL,2)); %Left
    
        for j = 1:length(uroiL) % loop ROIs
            clm(j,:) = mean(cL(roiL==uroiL(j), :));
        end
        
        % Right Hemisphere    
        crm = zeros(length(uroiR),size(cR,2)); %Right    
        
        for k = 1:length(uroiR) % loop ROIs
            crm(k,:) = mean(cR(roiR==uroiR(k), :));
        end
    end
    
    tp = [clm;crm]; % Combine left and right
    r = corrcoef(tp'); % Correlation matrix
    z = 0.5 * log((1+r)./(1-r)); % Fisher-z transform
    
    
    end % function
    
    