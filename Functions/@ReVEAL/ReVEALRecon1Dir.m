function ReVEALRecon1Dir(obj,dir)
%REVAMPEST Summary of this function goes here
%   Detailed explanation goes here

if isempty(dir)
    dir = 'x';
end

%==========================================================================
fprintf('\n')
% star_display('Beginning ReVAMP Est',0)

% Set fftw planner to estimate best fft algorithm
fftw('planner','patient');

%==========================================================================
% Do some input checks
if ~isequal(size(obj.data.Yb),size(obj.data.Yx))
    error('Size of obj.data.Yb and obj.data.Yx are not equal')
end
if ~isequal(size(obj.data.sampX),size(obj.data.sampB))
    error('Sampling pattern sizes are not consistent')
end
if ~isequal(size(obj.data.sensMaps),size(obj.data.Yb))
    error('Sensitivity map size is not consistent with data size')
end
if ~isequal(size(squeeze(obj.data.Yb(:,:,1,:))),size(obj.data.sampB))
    error('Sampling pattern and data sizes are not consistent')
end
if ~isequal(size(obj.data.x0),size(squeeze(obj.data.Yb(:,:,1,:))))
    error('x0 size is not consistent with data size')
end
if isempty(obj.data.maxwellCorrX)
    warning('No Maxwell Correction was supplied.  This may affect Results');
end

%==========================================================================
% Set input options

% fftshifts
obj.data.Yb = fftshift(fftshift(obj.data.Yb,1),2);
obj.data.sampB = fftshift(fftshift(obj.data.sampB,1),2);

% Choose which direction to reconstruct
if strcmpi(dir,'x')
    Yv = fftshift(fftshift(obj.data.Yx,1),2);
    sampV = fftshift(fftshift(obj.data.sampX,1),2);
    maxwellCorr =  fftshift(fftshift(obj.data.maxwellCorrX,1),2);
elseif strcmpi(dir,'y')
    Yv = fftshift(fftshift(obj.data.Yy,1),2);
    sampV = fftshift(fftshift(obj.data.sampY,1),2);
    maxwellCorr =  fftshift(fftshift(obj.data.maxwellCorrY,1),2);
elseif strcmpi(dir,'z')
    Yv = fftshift(fftshift(obj.data.Yz,1),2);
    sampV = fftshift(fftshift(obj.data.sampZ,1),2);
    maxwellCorr =  fftshift(fftshift(obj.data.maxwellCorrZ,1),2);
end

obj.data.sensMaps = fftshift(fftshift(obj.data.sensMaps,1),2);
obj.data.x0 =  fftshift(fftshift(obj.data.x0,1),2);
obj.options.ReVEALOpts.gamma = fftshift(fftshift(obj.options.ReVEALOpts.gamma,1),2);
obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);

% Normalize The columns of A to be unit norm
R = numel(obj.data.sampB)/length(find(obj.data.sampB ==1));
obj.data.sensMaps = obj.data.sensMaps*sqrt(R);
obj.data.Yb = obj.data.Yb*sqrt(R);
Yv = Yv*sqrt(R);
obj.options.ReVEALOpts.wvar = obj.options.ReVEALOpts.wvar*R;
obj.options.ReVEALOpts.sigmaSq = obj.options.ReVEALOpts.sigmaSq*R;
fprintf(sprintf('R = %s\n',num2str(R)))

% Copy variables
GAMPoptB = obj.options.GAMPOptB;
if strcmpi(dir,'x')
    GAMPoptV = obj.options.GAMPOptX;
elseif strcmpi(dir,'y')
    GAMPoptV = obj.options.GAMPOptY;
elseif strcmpi(dir,'z')
    GAMPoptV = obj.options.GAMPOptZ;
end

% Remove any warm start
GAMPoptB.xvar0 = [];
GAMPoptB.shat0 = [];
GAMPoptB.svar0 = [];
GAMPoptB.xhatPrev0 = [];
GAMPoptB.scaleFac = 1;

% Remove any warm start
GAMPoptV.xvar0 = [];
GAMPoptV.shat0 = [];
GAMPoptV.svar0 = [];
GAMPoptV.xhatPrev0 = [];
GAMPoptV.scaleFac = 1;

ReVEALOpts = obj.options.ReVEALOpts;
SparseTrans = obj.options.SparseTrans;
lambda0 = obj.options.ReVEALOpts.lambda0;
wvar = obj.options.ReVEALOpts.wvar;
gamma = obj.options.ReVEALOpts.gamma;
sigmaSq = obj.options.ReVEALOpts.sigmaSq;
L1 = obj.options.ReVEALOpts.L1;

% Set initializations
GAMPoptB.xhat0 = obj.data.x0(:);
GAMPoptV.xhat0 = obj.data.x0(:);

% Downsample Data
dataB = downsample_data(obj.data.Yb,obj.data.sampB);
dataV = downsample_data(Yv,sampV);

%==========================================================================
% Create Input Estimation Class
if strcmpi(obj.options.ReVEALOpts.compute,'gpu')
    gpu = 1;
else
    gpu = 0;
end
pMRI_b = pMRI_Op_2D_t(obj.data.sensMaps,obj.data.sampB,'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);
pMRI_v = pMRI_Op_2D_t(obj.data.sensMaps,sampV,'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);

% Check the norm of the columns of A
% tmp = [1;zeros(numel(obj.data.sampB)-1,1)];
% for ind = 1:10
%     tmp = circshift(tmp,round(100000*rand(1,1)));
%     a = pMRI_b.mult(tmp(:));
%     norm(a(:))
%     mean(a(:))
% end
% clear a tmp

% Create nd-DWT Linear Transform
W = ndDWTLinTrans(SparseTrans,size(obj.data.sampB),'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);

% Concatonate All The Linear Transform Operator Together
Op_b = LinTransConcat({pMRI_b;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
Op_v = LinTransConcat({pMRI_v;W},[1,1],obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);

%==========================================================================
% Create Output Estimation classes
MeasEstimOut_b = CAwgnEstimOut(dataB,wvar,obj.options.ReVEALOpts.MAP);
MeasEstimOut_v = CAwgnEstimOut(dataV,wvar,obj.options.ReVEALOpts.MAP);

% Set the Regularization per band
lambda = setLambda(size(obj.data.x0),lambda0);

% Choose L1 or SNIPE Estimator
if L1
    AnaEstimOut1 = CplxLaplaceEstimOut(lambda);
else
    % GM outpe
    [ lambda_b, ~, ~, sigma1_b, sigma2_b] = estBGParams(reshape(obj.data.x0,size(obj.data.sampB)));
    lambda_v = lambda_b;
    sigma1_v = sigma1_b;
    sigma2_v = sigma2_b;
    AnaEstimOut1  = CGaussMixEstimOut(zeros(8*numel(obj.data.x0),1),sigma1_b,sigma2_b,lambda_b);
    
    sigma1_v1 = reshape(sigma1_v,[size(obj.data.sampB),8]);
    sigma2_v1 = reshape(sigma2_v,[size(obj.data.sampB),8]);
    lambda_v1 = reshape(lambda_v,[size(obj.data.sampB),8]);

    sigma1_v1 = squeeze(sigma1_v1(1,1,1,1:8));
    sigma2_v1 = squeeze(sigma2_v1(1,1,1,1:8));
    lambda_v1 = squeeze(lambda_v1(1,1,1,1:8));

    fprintf(sprintf('sigma1\t\t sigma2\t\t lambda\n'))
    for ii = 1:8
        fprintf(sprintf('%f\t%f\t%f\n',sigma1_v1(ii),sigma2_v1(ii),lambda_v1(ii)));
    end
end

% Create Output Estimation class
EstimOut_b = EstimOutConcat({MeasEstimOut_b;AnaEstimOut1},[pMRI_b.M,W.M],...
    obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
EstimOut_v = EstimOutConcat({MeasEstimOut_v;AnaEstimOut1},[pMRI_v.M,W.M],...
    obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);

clear AnaEstimOut1 MeasEstimOut_b pMRI_b pMRI_v W sampV dataB dataV lambda Yv ...
      MeasEstimOut_v

% Create Input Estimation class
EstimIn_b = NullEstimIn(0,1);
EstimIn_v = NullEstimIn(0,1);

%==========================================================================
% Initial Reconstruction 
reveal_start = tic;
if L1
    xhatGAMP_b = gampEst(EstimIn_b,EstimOut_b,Op_b,GAMPoptB);
    xhatGAMP_v = gampEst(EstimIn_v,EstimOut_v,Op_v,GAMPoptV);
else
    [xhatb,maps,xhatGAMP_b,paramb]= pMRIEMGM(ifftshift(ifftshift(obj.data.Yb,1),2)/sqrt(R),...
        ifftshift(ifftshift(obj.data.sampB,1),2),wvar/R);
    [xhatv,~,xhatGAMP_v,paramv] = pMRIEMGM(ifftshift(ifftshift(Yv,1),2)/sqrt(R),...
        ifftshift(ifftshift(sampV,1),2),wvar/R,ifftshift(ifftshift(maps,1),2));
		
        % Learn weights for xb
    AnaEstimOut1 = CGaussMixEstimOut(zeros(8*numel(obj.data.sampB),1),paramb.sigma1,paramb.sigma2,paramb.lambda);
    EstimOut_b = EstimOutConcat({MeasEstimOut_b;AnaEstimOut1},[pMRI_b.M,W.M],...
            obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
        
        % Learn weights for xv
    AnaEstimOut1 = CGaussMixEstimOut(zeros(8*numel(obj.data.sampB),1),paramv.sigma1,paramv.sigma2,paramv.lambda);
    EstimOut_v = EstimOutConcat({MeasEstimOut_v;AnaEstimOut1},[pMRI_v.M,W.M],...
            obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
            
    xhatGAMP_b = gampEst(EstimIn_b,EstimOut_b,Op_b,GAMPoptB);
    xhatGAMP_v = gampEst(EstimIn_v,EstimOut_v,Op_v,GAMPoptV);
end

%% ========================================================================
%  Estimate the Maxwell Correction from the data and apply it to the
%  ssenitvity maps

% Discard  Imaginary parts of the variance
fprintf('It = %d\t|dx_b|/|x_b| = %0.4f ----- |dx_v|/|x_v| = %0.4f\n', 1,1,1)

% Save Previous Results
xhatb_prev = xhatGAMP_b.xhat;
xhatv_prev = xhatGAMP_v.xhat;

% Add Maxwell Correction
if ~isempty(maxwellCorr)
    maxwell_phase = maxwellCorr;
    maxwell_phase = repmat(maxwell_phase,[1,1,size(obj.data.sampB,3)]);
    maxwell_phase = exp(-1j*maxwell_phase(:));
    
else
    maxwell_phase  = maxxCorr2D(ifftshift(reshape(xhatGAMP_b.xhat,size(obj.data.sampB)))...
                        ,ifftshift(reshape(xhatGAMP_v.xhat,size(obj.data.sampB))));
    maxwell_phase = repmat(maxwell_phase,[1,1,size(obj.data.sampB,3)]);
    maxwell_phase  = fftshift(exp(1j*maxwell_phase));
    maxwell_phase = maxwell_phase(:);
    if strcmpi('gpu',obj.options.ReVEALOpts.compute)
        maxwell_phase = gpuArray(maxwell_phase);
    end
end

% =========================================================================
% Setup MRF 
% [xg,yg,zg] = meshgrid(1:size(obj.data.sampB,1),1:size(obj.data.sampB,2)...
%                       ,1:size(obj.data.sampB,3));
% coordinates = [xg(:),yg(:),zg(:)];
% beta = 0.8;
% betaz = beta;            % Horizontal inverse temperature
% betay = beta;            % Vertical inverse temperature
% betax = beta;
% alpha = 0;              % Sparsity parameter
% maxIter = 10;            % Max # of loopy BP iters per turbo iter
% mrf = markovField3D('coordinates',coordinates,'betax',betax,'betay'...
%     ,betay,'betaz',betaz,'alpha',alpha,'maxIter',maxIter);

% % =========================================================================
% First-order Markov chain properties
p01 = 0.05;     % Default active-to-inactive transition probability
learn_p01 = 'true';      % Learn p01 using EM alg. by default
dim = 'row';    % Each row forms a Markov chain
sparsity = 0.2;
sparsity = sparsity*ones(size(obj.data.sampB,1)*size(obj.data.sampB,2),size(obj.data.sampB,3));
mrc = MarkovChain1('p01',p01,'dim',dim);

tmp = gather(reshape(obj.data.x0,size(obj.data.sampB)));
level = graythresh(abs(tmp(:,:,1)));
mask = im2bw(abs(tmp(:,:,1)),level);
mask = repmat(mask,[1,1,size(obj.data.sampB,3)]);
clear tmp

gamma_in = gamma;
%==========================================================================
% Main Loop
for ind = 2:ReVEALOpts.nit;

    GAMPoptB.nit = 10;
    GAMPoptV.nit = 10;
    
    % Create input estimation class for xb
    EstimIn_b1 = ncCAwgnEstimIn(xhatGAMP_v.rhat,(xhatGAMP_v.rvar + sigmaSq));
    EstimIn_b2 = CAwgnEstimIn(bsxfun(@times,xhatGAMP_v.rhat,conj(maxwell_phase)),(xhatGAMP_v.rvar + sigmaSq));
    EstimIn_b = MixScaEstimIn(EstimIn_b1,gamma_in(:),EstimIn_b2);

    % Create input estimation class for xv
    EstimIn_v1 = ncCAwgnEstimIn(xhatGAMP_b.rhat,(xhatGAMP_b.rvar + sigmaSq));
    EstimIn_v2 = CAwgnEstimIn(bsxfun(@times,xhatGAMP_b.rhat,maxwell_phase),(xhatGAMP_b.rvar + sigmaSq));
    EstimIn_v = MixScaEstimIn(EstimIn_v1,gamma_in(:),EstimIn_v2);
    
    % EM update on Gaussian Mixture Components
    if ~L1
        % Learn weights for xb
        [ sigma1_b,sigma2_b,lambda_b] = CGMUpdate(8*xhatGAMP_b.zhat,64*xhatGAMP_b.zvar,size(obj.data.sampB),sigma1_b,sigma2_b,lambda_b);
        AnaEstimOut1 = CGaussMixEstimOut(zeros(8*numel(obj.data.sampB),1),sigma1_b,sigma2_b,lambda_b);
        EstimOut_b = EstimOutConcat({MeasEstimOut_b;AnaEstimOut1},[pMRI_b.M,W.M],...
            obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
        
        % Learn weights for xv
        [ sigma1_v,sigma2_v,lambda_v] = CGMUpdate(8*xhatGAMP_v.zhat,64*xhatGAMP_v.zvar,size(obj.data.sampB),sigma1_v,sigma2_v,lambda_v);
        AnaEstimOut1 = CGaussMixEstimOut(zeros(8*numel(obj.data.sampB),1),sigma1_v,sigma2_v,lambda_v);
        EstimOut_v = EstimOutConcat({MeasEstimOut_v;AnaEstimOut1},[pMRI_v.M,W.M],...
            obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);  
    end
   

    % Warm Start xb
    GAMPoptB.xhat0 = xhatGAMP_b.xhat;
    GAMPoptB.xvar0 = xhatGAMP_b.xvar;
    GAMPoptB.shat0 = xhatGAMP_b.shat;
    GAMPoptB.svar0 = xhatGAMP_b.svar;
    GAMPoptB.xhatPrev0 = xhatGAMP_b.xhatPrev;
    GAMPoptB.scaleFac = xhatGAMP_b.scaleFac;
    GAMPoptB.step = min(max(xhatGAMP_b.step,0.05),xhatGAMP_b.stepMax);
    
    % Warm Start xv
    GAMPoptV.xhat0 = xhatGAMP_v.xhat;
    GAMPoptV.xvar0 = xhatGAMP_v.xvar;
    GAMPoptV.shat0 = xhatGAMP_v.shat;
    GAMPoptV.svar0 = xhatGAMP_v.svar;
    GAMPoptV.xhatPrev0 = xhatGAMP_v.xhatPrev;
    GAMPoptV.scaleFac = xhatGAMP_v.scaleFac;
    GAMPoptV.step = min(max(xhatGAMP_v.step,0.05),xhatGAMP_v.stepMax);

    % Reconstruct x_b With GrAMPA
    xhatGAMP_b = gampEst(EstimIn_b,EstimOut_b,Op_b,GAMPoptB);
   
    % Reconstruct x_v With GrAMPA
    xhatGAMP_v = gampEst(EstimIn_v,EstimOut_v,Op_v,GAMPoptV);

    % Report Output if Desired
    change_b = norm(xhatGAMP_b.xhat-xhatb_prev)/norm(xhatGAMP_b.xhat);
    change_v = norm(xhatGAMP_v.xhat-xhatv_prev)/norm(xhatGAMP_v.xhat);
    fprintf('It = %d\t|dx_b|/|x_b| = %0.4f ----- |dx_v|/|x_v| = %0.4f\n', ind,change_b,change_v)
	
    % Save Previous Results
    xhatb_prev = xhatGAMP_b.xhat;
    xhatv_prev = xhatGAMP_v.xhat;
        
    em_gamma = 0;
    if em_gamma == 1 && ind > 3
        fprintf('Peforming LBP on MRF...')
        tic;
        v_1 = post_v(bsxfun(@times,xhatGAMP_b.rhat,maxwell_phase), xhatGAMP_v.rhat, xhatGAMP_b.rvar, xhatGAMP_v.rvar, sigmaSq, gamma);

%         % EM Update on Gamma
%         v_1 = post_v(xhatGAMP_b.rhat.*maxwell_phase, xhatGAMP_v.rhat, xhatGAMP_b.rvar, xhatGAMP_v.rvar, sigmaSq, gamma);
%         v_1 = reshape(v_1,size(obj.data.sampB));
% 
%         filt(1,1,:) = normpdf(-3:1:3,0,1)/sum(normpdf(-3:1:3,0,1));
%         gamma = convn(v_1,filt,'same');
%         gamma = circshift(gamma,[0,0,-3]);
%         gamma = gamma(:);
%         gamma(gamma>=1) = 0.99;

        % =================================================================
        % MRF
%         v_1 = mask(:).*v_1(:);
%         [gamma, ~, beta, alpha] = mrf.mrf_spd(v_1(:));
% %         gamma = max(reshape(gamma,size(obj.data.sampB)),[],3);
% %         gamma = repmat(gamma,[1,1,size(obj.data.sampB,3)]);
%         gamma = gamma(:);
%         fprintf(sprintf('Completed in %s s.\n',num2str(toc)));
%         if ind >3
%             mrf.betaz = beta;            % Horizontal inverse temperature
%             mrf.betay = beta;            % Vertical inverse temperature
%             mrf.betax = beta;
%             mrf.alpha = alpha;
%             fprintf(sprintf('EM Update alpha = %s, beta = %s\n',num2str(alpha),num2str(beta)));
%         end

        % =================================================================
        % MRC
        v_1 = reshape(v_1,[size(obj.data.sampB,1)*size(obj.data.sampB,2),size(obj.data.sampB,3)]);
        [gamma, S_POST, p01_upd] = mrc.binarymarkov_msgs(v_1, sparsity);
        sparsity = mean(S_POST(mask ==1));
        sparsity = sparsity*ones(size(obj.data.sampB,1)*size(obj.data.sampB,2),size(obj.data.sampB,3));
        if p01_upd<1e-6; p01_upd=1e-6 ;end
        mrc.p01 = p01_upd;
        gamma_in = max(reshape(gamma,size(obj.data.sampB)),[],3);
        gamma_in = repmat(gamma_in,[1,1,size(obj.data.sampB,3)]);
        gamma = gamma(:);
        fprintf(sprintf('\nEM Update p01 = %s, sparsity = %s\n',num2str(p01_upd),num2str(sparsity(1))));
    end
    
    % Check for Stoppping condition
    if mean([change_b,change_v])<obj.options.ReVEALOpts.tol;
%         break
    end
end
fprintf('ReVAMP Calculated in %s s\n', num2str(toc(reveal_start)))
% star_display('',0)

%==========================================================================
%Get Posterior on v and set outputs
v_1 = post_v(bsxfun(@times,xhatGAMP_b.rhat,maxwell_phase), xhatGAMP_v.rhat, xhatGAMP_b.rvar, xhatGAMP_v.rvar, sigmaSq, gamma);

% Calculate the Map Estimate theta
theta =  angle(bsxfun(@times,xhatGAMP_v.xhat,conj(maxwell_phase)).*conj(xhatGAMP_b.xhat));
theta_MAP = angle(bsxfun(@times,xhatGAMP_v.rhat,conj(maxwell_phase)).*conj(xhatGAMP_b.rhat));

% Reshape outputs
xhatb = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
xhatv = reshape(xhatGAMP_v.xhat,size(obj.data.sampB));
v_1 = reshape(v_1,size(obj.data.sampB));
theta = reshape(theta,size(obj.data.sampB));
theta_MAP = reshape(theta_MAP,size(obj.data.sampB));
% maxwell_phase = reshape(maxwell_phase,size(obj.data.sampB));

if size(gamma(:)) == size(v_1(:))
    gamma = reshape(gamma,size(obj.data.sampB));
end

% fftshift outputs
xhatb = ifftshift(ifftshift(xhatb,1),2);
xhatv = ifftshift(ifftshift(xhatv,1),2);
v_1 = ifftshift(ifftshift(v_1,1),2);
theta = ifftshift(ifftshift(theta,1),2);
theta_MAP = ifftshift(ifftshift(theta_MAP,1),2);
% maxwell_phase = ifftshift(ifftshift(maxwell_phase,1),2);

obj.outputs.xHatb = xhatb;
if strcmpi(dir,'x')
    obj.outputs.xHatx = xhatv;
    obj.outputs.vX = v_1;
    obj.outputs.thetaX = theta;
    obj.outputs.thetaMAPX = theta_MAP;
elseif strcmpi(dir,'y')
    obj.outputs.xHaty = xhatv;
    obj.outputs.vY = v_1;
    obj.outputs.thetaY = theta;
    obj.outputs.thetaMAPY = theta_MAP;
elseif strcmpi(dir,'z')
    obj.outputs.xHatz = xhatv;
    obj.outputs.vZ = v_1;
    obj.outputs.thetaZ = theta;
    obj.outputs.thetaMAPZ = theta_MAP;
end

if size(gamma(:)) == size(v_1(:))
    gamma= ifftshift(ifftshift(gamma,1),2);
end

%========================================================================
% Undo fftshifts
obj.data.Yb = ifftshift(ifftshift(obj.data.Yb,1),2);
obj.data.sampB = ifftshift(ifftshift(obj.data.sampB,1),2);
obj.data.sensMaps = ifftshift(ifftshift(obj.data.sensMaps,1),2);
obj.data.x0 =  ifftshift(ifftshift(obj.data.x0,1),2);
obj.options.ReVEALOpts.gamma = fftshift(fftshift(obj.options.ReVEALOpts.gamma,1),2);
obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);

% undo scaling
obj.data.sensMaps = obj.data.sensMaps/sqrt(R);
obj.data.Yb = obj.data.Yb/sqrt(R);
obj.options.ReVEALOpts.wvar = obj.options.ReVEALOpts.wvar/R;
obj.options.ReVEALOpts.sigmaSq = obj.options.ReVEALOpts.sigmaSq/R;

end

