
% ------------------------------------------------------>
# 1. MATLAB ket - pla diff pca
% ------------------------------------------------------>

% load image
GBCdiff=gmrimages('Ket-Pla_HCS40_hpss_res-mVWMWB1d__CAB-NP_P718_BSNIP_parcels_gbc_mFz.ptseries.nii');

% change to data file
ndiff = size(GBCdiff.data,2)
dlmwrite('GBC_ket-pla_data.dat', GBCdiff.data,'\t');	

% read in datafile
GBCdiffdat=GBCdiff.data';

% PCA

% Principal component scores are the representations of X in the principal component space. Rows of score correspond to observations, and columns correspond to components.
[coeffdiff,scorediff,latentdiff,tsquareddiff,explaineddiff,mudiff] = pca(GBCdiffdat);
dlmwrite('ket-pla_pca_explainedvar.dat',explaineddiff,'\t');	
dlmwrite('ket-pla_pca_scores.dat',scorediff,'\t');	
dlmwrite('ket-pla_pca_image_coeff.txt',coeffdiff,'\t');	

% save the nii image
diffpcacoeffs = GBCdiff.zeroframes(ndiff-1);
diffpcacoeffs.data = coeffdiff;
diffpcacoeffs.img_SaveNIfTI('ket-pla_pca_image_');


% permutation analysis

nPCs = size(explaineddiff,1);

%% Shuffle for significance permutation testing
k=5000 % No. of shuffles
explaineddiffallshuffs = zeros(nPCs,k); % Initiate matrix for storing shuffled explained variance

A = GBCdiffdat; % Original data matrix
[M,N] = size(A);
rowIndex = repmat((1:M)',[1 N]); % Preserve row indices

for shuff=1:k
	[~,randomizedColIndex] = sort(rand(M,N),2); % Get randomized column indices by sorting random array
	newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex); % Linear indexing to create ordering for shuffled matrix B
	B = A(newLinearIndex); % B is matrix of shuffled (within subject) data
	[coeffdiffB,scorediffB,latentdiffB,tsquareddiffB,explaineddiffB,mudiffB] = pca(B); % PCA on shuffled data
	explaineddiffallshuffs(:,shuff)=explaineddiffB; % Store explainedvar
end

% Compute 0.95 chance
diffcrit95=zeros(1,39);

for pc = 1:39
	sorted=sort(explaineddiffallshuffs(pc,:));
	diffcrit95(pc)=sorted(4750);
end

diffcrit95_t = transpose(diffcrit95)

dlmwrite('ket-pla_pca_explainedvar_5kshuffle_pcrit95.dat',diffcrit95_t,'\t');



