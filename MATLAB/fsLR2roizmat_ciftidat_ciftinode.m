function [z, r, tp] = fsLR2roizmat_ciftidat_ciftinode(cii, roi_cii)
    % DESCRIPTION:
    %   Extract node time series from 32k fsLR saved in CIFTI dtseries. 
    %   Outputs correlation matrix (r), fishers-z-transfromed corr matrix (z),
    %   and node-by-timepoint matrix (tp). 
    %   
    %   cii:        Full path to CIFIT data or CIFTI object
    %
    %   roi_cii:    Full path to nodes in CIFTI format or CIFTI object in
    %   
    %   requires: cifti-matlab
    % ----------% 
    
    
    % --- Load and setup ROIs ---
    if ischar(roi_cii) % check if roi is a file name or gifti object
        roi_cii = cifti_read(roi_cii); 
    end
    
    roiL=cifti_struct_dense_extract_surface_data(roi_cii, 'CORTEX_LEFT');
    roiR=cifti_struct_dense_extract_surface_data(roi_cii, 'CORTEX_RIGHT');
    
    uroiL = unique(roiL);
    uroiL = uroiL(uroiL~=0);   
    
    uroiR = unique(roiR);
    uroiR=uroiR(uroiR~=0);    
    
    % --- Load and setup data ---
    if ischar(cii) % check if roi is a file name or gifti object
        cii = cifti_read(cii); 
    end
    
    % LEFT
    cL=cifti_struct_dense_extract_surface_data(cii, 'CORTEX_LEFT');
    
    
    % RIGHT
    cR=cifti_struct_dense_extract_surface_data(cii, 'CORTEX_RIGHT');
    
    
    % Left Hemisphere
    clm = zeros(length(uroiL),size(cL,2)); 
    
    for j = 1:length(uroiL) % loop ROIs
        clm(j,:) = mean(cL(roiL==uroiL(j), :));
    end
    
    % Right Hemisphere
    crm = zeros(length(uroiR),size(cR,2));
    
    for k = 1:length(uroiR) % loop ROIs
        crm(k,:) = mean(cR(roiR==uroiR(k), :));
    end
    
    tp = [clm;crm]; % Combine left and right
    r = corrcoef(tp'); % Correlation matrix
    z = 0.5 * log((1+r)./(1-r)); % Fisher-z transform
    
    end % function
    
    