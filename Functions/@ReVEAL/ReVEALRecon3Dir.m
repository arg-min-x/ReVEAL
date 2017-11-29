function ReVEALRecon3Dir(obj)
%REVEALRECON3DIR Summary of this function goes here
%   Detailed explanation goes here

%% fftshifts
obj.data.Yb = fftshift(fftshift(obj.data.Yb,1),2);
obj.data.Yx = fftshift(fftshift(obj.data.Yx,1),2);
obj.data.Yy = fftshift(fftshift(obj.data.Yy,1),2);
obj.data.Yz = fftshift(fftshift(obj.data.Yz,1),2);

obj.data.sampB = fftshift(fftshift(obj.data.sampB,1),2);
obj.data.sampX = fftshift(fftshift(obj.data.sampX,1),2);
obj.data.sampY = fftshift(fftshift(obj.data.sampY,1),2);
obj.data.sampZ = fftshift(fftshift(obj.data.sampZ,1),2);

obj.data.maxwellCorrX =  fftshift(fftshift(obj.data.maxwellCorrX,1),2);
obj.data.maxwellCorrY =  fftshift(fftshift(obj.data.maxwellCorrY,1),2);
obj.data.maxwellCorrZ =  fftshift(fftshift(obj.data.maxwellCorrZ,1),2);

obj.data.sensMaps = fftshift(fftshift(obj.data.sensMaps,1),2);
obj.data.x0 =  fftshift(fftshift(obj.data.x0,1),2);
obj.options.ReVEALOpts.gamma = fftshift(fftshift(obj.options.ReVEALOpts.gamma,1),2);
obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);

%% Setup Maxwell Corrections
maxwellCorrX = repmat(exp(1j*obj.data.maxwellCorrX),[1,1,size(obj.data.sampB,3)]);
maxwellCorrY = repmat(exp(1j*obj.data.maxwellCorrY),[1,1,size(obj.data.sampB,3)]);
maxwellCorrZ = repmat(exp(1j*obj.data.maxwellCorrZ),[1,1,size(obj.data.sampB,3)]);

maxwellCorrX = maxwellCorrX(:);
maxwellCorrY = maxwellCorrY(:);
maxwellCorrZ = maxwellCorrZ(:);

%% Normalize The columns of A to be unit norm
R = numel(obj.data.sampB)/length(find(obj.data.sampB ==1));
obj.data.sensMaps = obj.data.sensMaps*sqrt(R);
obj.data.Yb = obj.data.Yb*sqrt(R);
obj.data.Yx = obj.data.Yx*sqrt(R);
obj.data.Yy = obj.data.Yy*sqrt(R);
obj.data.Yz = obj.data.Yz*sqrt(R);
obj.options.ReVEALOpts.wvar = obj.options.ReVEALOpts.wvar*R;
obj.options.ReVEALOpts.sigmaSq = obj.options.ReVEALOpts.sigmaSq*R;
fprintf(sprintf('R = %s\n',num2str(R)))

% Set initializations
obj.options.GAMPOptB.xhat0 = obj.data.x0(:);
obj.options.GAMPOptX.xhat0 = obj.data.x0(:);
obj.options.GAMPOptY.xhat0 = obj.data.x0(:);
obj.options.GAMPOptZ.xhat0 = obj.data.x0(:);

% Downsample Data
dataB = downsample_data(obj.data.Yb,obj.data.sampB);
dataX = downsample_data(obj.data.Yx,obj.data.sampX);
dataY = downsample_data(obj.data.Yy,obj.data.sampY);
dataZ = downsample_data(obj.data.Yz,obj.data.sampZ);

% Create Input Estimation Class
pMRI_b = pMRI_Op_2D_t(obj.data.sensMaps,obj.data.sampB,'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);
pMRI_x = pMRI_Op_2D_t(obj.data.sensMaps,obj.data.sampX,'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);
pMRI_y = pMRI_Op_2D_t(obj.data.sensMaps,obj.data.sampY,'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);
pMRI_z = pMRI_Op_2D_t(obj.data.sensMaps,obj.data.sampZ,'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);

% Create nd-DWT Linear Transform
W = ndDWTLinTrans(obj.options.SparseTrans,size(obj.data.sampB),'uniform_var',...
    obj.options.ReVEALOpts.uniform_var,'compute',obj.options.ReVEALOpts.compute,...
    'precision',obj.options.ReVEALOpts.precision);

% Concatonate All The Linear Transform Operator Together
Op_b = LinTransConcat({pMRI_b;W},[1,1],obj.options.ReVEALOpts.precision,...
    obj.options.ReVEALOpts.compute);
Op_x = LinTransConcat({pMRI_x;W},[1,1],obj.options.ReVEALOpts.precision,...
    obj.options.ReVEALOpts.compute);
Op_y = LinTransConcat({pMRI_y;W},[1,1],obj.options.ReVEALOpts.precision,...
    obj.options.ReVEALOpts.compute);
Op_z = LinTransConcat({pMRI_z;W},[1,1],obj.options.ReVEALOpts.precision,...
    obj.options.ReVEALOpts.compute);

% Create Output Estimation classes
MeasEstimOut_b = CAwgnEstimOut(dataB,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
MeasEstimOut_x = CAwgnEstimOut(dataX,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
MeasEstimOut_y = CAwgnEstimOut(dataY,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);
MeasEstimOut_z = CAwgnEstimOut(dataZ,obj.options.ReVEALOpts.wvar,obj.options.ReVEALOpts.MAP);

% Set the Regularization per band
lambda = setLambda(size(obj.data.x0),obj.options.ReVEALOpts.lambda0);
AnaEstimOut1 = CplxLaplaceEstimOut(lambda);

% Create Output Estimation class
EstimOut_b = EstimOutConcat({MeasEstimOut_b;AnaEstimOut1},[pMRI_b.M,W.M],...
    obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
EstimOut_x = EstimOutConcat({MeasEstimOut_x;AnaEstimOut1},[pMRI_x.M,W.M],...
    obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
EstimOut_y = EstimOutConcat({MeasEstimOut_y;AnaEstimOut1},[pMRI_y.M,W.M],...
    obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);
EstimOut_z = EstimOutConcat({MeasEstimOut_z;AnaEstimOut1},[pMRI_z.M,W.M],...
    obj.options.ReVEALOpts.precision,obj.options.ReVEALOpts.compute);

% Create Input Estimation class
EstimIn = NullEstimIn(0,1);

% Clear some memory up
clear meanB meanX meanY meanZ varB varX varY varZ pvarb pvarx pvary pvarz...
      MeasEstimOut_b MeasEstimOut_x MeasEstimOut_y MeasEstimOut_z...
      AnaEstimOut1 pMRI_b pMRI_x pMRI_y pMRI_z W dataB dataX dataY dataZ ...
      lambda;

%% Initial Reconstruction 
tic
fprintf('\npB...')
xhatGAMP_b = gampEst(EstimIn,EstimOut_b,Op_b,obj.options.GAMPOptB);
xhatGAMP_b = clearOpts(xhatGAMP_b);
fprintf('pX...')
xhatGAMP_x = gampEst(EstimIn,EstimOut_x,Op_x,obj.options.GAMPOptX);
xhatGAMP_x = clearOpts(xhatGAMP_x);
fprintf('pY...')
xhatGAMP_y = gampEst(EstimIn,EstimOut_y,Op_y,obj.options.GAMPOptY);
xhatGAMP_y = clearOpts(xhatGAMP_y);
fprintf('pZ...\n')
xhatGAMP_z = gampEst(EstimIn,EstimOut_z,Op_z,obj.options.GAMPOptZ);
xhatGAMP_z = clearOpts(xhatGAMP_z);
fprintf('It = %d\t|dx_b|/|x_b| = %0.4f ----- |dx_v|/|x_v| = %0.4f\n', 1,1,1)

% % =========================================================================
% First-order Markov chain properties
p01 = 0.05;     % Default active-to-inactive transition probability
learn_p01 = 'true';      % Learn p01 using EM alg. by default
dim = 'row';    % Each row forms a Markov chain
sparsity = 0.2;
sparsity = sparsity*ones(size(obj.data.sampB,1)*size(obj.data.sampB,2),size(obj.data.sampB,3));
mrc = MarkovChain1('p01',p01,'dim',dim);

tmp = reshape(obj.data.x0,size(obj.data.sampB));
level = graythresh(abs(tmp(:,:,1)));
mask = im2bw(abs(tmp(:,:,1)),level);
mask = repmat(mask,[1,1,size(obj.data.sampB,3)]);
clear tmp

% Begin Main Loop
gamma = obj.options.ReVEALOpts.gamma;
%% Main Loop
for ind = 2:obj.options.ReVEALOpts.nit
    obj.options.GAMPOptB.nit = 10;
    obj.options.GAMPOptX.nit = 10;
    obj.options.GAMPOptY.nit = 10;
    obj.options.GAMPOptZ.nit = 10;
       
    %======================================================================
    % Create new input estimators
    pvarb = xhatGAMP_b.rvar + obj.options.ReVEALOpts.sigmaSq;
    pvarx = xhatGAMP_x.rvar + obj.options.ReVEALOpts.sigmaSq;
    pvary = xhatGAMP_y.rvar + obj.options.ReVEALOpts.sigmaSq;
    pvarz = xhatGAMP_z.rvar + obj.options.ReVEALOpts.sigmaSq;
    
    % Input Estimator for pb
    meanB = (abs(xhatGAMP_x.rhat).*pvary.*pvarz +...
             abs(xhatGAMP_y.rhat).*pvarx.*pvarz +...
             abs(xhatGAMP_z.rhat).*pvarx.*pvary)./...
             (pvarx.*pvary + pvarx.*pvarz + pvary.*pvarz) + eps;
    varB = (pvarx.*pvary.*pvarz)./...
             (pvarx.*pvary + pvarx.*pvarz + pvary.*pvarz) + eps;
    EstimIn_b1 = ncCAwgnEstimIn(meanB,varB);
    meanB = (xhatGAMP_x.rhat.*maxwellCorrX.*pvary.*pvarz +...
             xhatGAMP_y.rhat.*maxwellCorrY.*pvarx.*pvarz +...
             xhatGAMP_z.rhat.*maxwellCorrZ.*pvarx.*pvary)./...
             (pvarx.*pvary + pvarx.*pvarz + pvary.*pvarz) + eps;
    EstimIn_b2 = CAwgnEstimIn(meanB,varB);
    EstimIn_b = MixScaEstimIn(EstimIn_b1,obj.options.ReVEALOpts.gamma,EstimIn_b2);

    % Input Estimator for px
    meanX = (abs(xhatGAMP_b.rhat).*pvary.*pvarz +...
             abs(xhatGAMP_y.rhat).*pvarb.*pvarz +...
             abs(xhatGAMP_z.rhat).*pvarb.*pvary)./...
             (pvarb.*pvary + pvarb.*pvarz + pvary.*pvarz) + eps;
    varX = (pvarb.*pvary.*pvarz)./...
             (pvarb.*pvary + pvarb.*pvarz + pvary.*pvarz) + eps;     
    EstimIn_x1 = ncCAwgnEstimIn(meanX,varX);
    meanX = (xhatGAMP_b.rhat.*pvary.*pvarz +...
             xhatGAMP_y.rhat.*maxwellCorrY.*pvarb.*pvarz +...
             xhatGAMP_z.rhat.*maxwellCorrZ.*pvarb.*pvary)./...
             (pvarb.*pvary + pvarb.*pvarz + pvary.*pvarz) + eps;
    EstimIn_x2 = CAwgnEstimIn(meanX,varX);
    EstimIn_x = MixScaEstimIn(EstimIn_x1,obj.options.ReVEALOpts.gamma,EstimIn_x2);

    % Input Estimator for py
    meanY = (abs(xhatGAMP_b.rhat).*pvarx.*pvarz +...
             abs(xhatGAMP_x.rhat).*pvarb.*pvarz +...
             abs(xhatGAMP_z.rhat).*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvarz + pvarx.*pvarz) + eps;
    varY = (pvarb.*pvarx.*pvarz)./...
             (pvarb.*pvarx + pvarb.*pvarz + pvarx.*pvarz) + eps;     
    EstimIn_y1 = ncCAwgnEstimIn(meanY,varY);
    meanY = (xhatGAMP_b.rhat.*pvarx.*pvarz +...
             xhatGAMP_x.rhat.*maxwellCorrX.*pvarb.*pvarz +...
             xhatGAMP_z.rhat.*maxwellCorrZ.*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvarz + pvarx.*pvarz) + eps;
    EstimIn_y2 = CAwgnEstimIn(meanY,varY);
    EstimIn_y = MixScaEstimIn(EstimIn_y1,obj.options.ReVEALOpts.gamma,EstimIn_y2);

    % Input Estimator for py
    meanZ = (abs(xhatGAMP_b.rhat).*pvarx.*pvary +...
             abs(xhatGAMP_x.rhat).*maxwellCorrX.*pvarb.*pvary +...
             abs(xhatGAMP_y.rhat).*maxwellCorrY.*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvary + pvarx.*pvary) + eps;
    varZ = (pvarb.*pvarx.*pvary)./...
             (pvarb.*pvarx + pvarb.*pvary + pvarx.*pvary) + eps;
    EstimIn_z1 = ncCAwgnEstimIn(meanZ,varZ);
    meanZ = (xhatGAMP_b.rhat.*pvarx.*pvary +...
             xhatGAMP_x.rhat.*maxwellCorrX.*pvarb.*pvary +...
             xhatGAMP_y.rhat.*maxwellCorrY.*pvarb.*pvarx)./...
             (pvarb.*pvarx + pvarb.*pvary + pvarx.*pvary) + eps;
    EstimIn_z2 = CAwgnEstimIn(meanZ,varZ);
    EstimIn_z = MixScaEstimIn(EstimIn_z1,obj.options.ReVEALOpts.gamma,EstimIn_z2);
    
    % Clear some memory up
    clear meanB meanX meanY meanZ varB varX varY varZ pvarb pvarx pvary...
          pvarz EstimIn_b1 EstimIn_b2 EstimIn_x1 EstimIn_x2 EstimIn_y1...
          EstimIn_y2 EstimIn_z1 EstimIn_z2;
    
    %======================================================================
    % Warm Start xb
    obj.options.GAMPOptB.xhat0 = xhatGAMP_b.xhat;
    obj.options.GAMPOptB.xvar0 = xhatGAMP_b.xvar;
    obj.options.GAMPOptB.shat0 = xhatGAMP_b.shat;
    obj.options.GAMPOptB.svar0 = xhatGAMP_b.svar;
    obj.options.GAMPOptB.xhatPrev0 = xhatGAMP_b.xhatPrev;
    obj.options.GAMPOptB.scaleFac = xhatGAMP_b.scaleFac;
    obj.options.GAMPOptB.step = min(max(xhatGAMP_b.step,0.05),xhatGAMP_b.stepMax);

    % Warm Start xx
    obj.options.GAMPOptX.xhat0 = xhatGAMP_x.xhat;
    obj.options.GAMPOptX.xvar0 = xhatGAMP_x.xvar;
    obj.options.GAMPOptX.shat0 = xhatGAMP_x.shat;
    obj.options.GAMPOptX.svar0 = xhatGAMP_x.svar;
    obj.options.GAMPOptX.xhatPrev0 = xhatGAMP_x.xhatPrev;
    obj.options.GAMPOptX.scaleFac = xhatGAMP_x.scaleFac;
    obj.options.GAMPOptX.step = min(max(xhatGAMP_x.step,0.05),xhatGAMP_x.stepMax);

    % Warm Start xy
    obj.options.GAMPOptY.xhat0 = xhatGAMP_y.xhat;
    obj.options.GAMPOptY.xvar0 = xhatGAMP_y.xvar;
    obj.options.GAMPOptY.shat0 = xhatGAMP_y.shat;
    obj.options.GAMPOptY.svar0 = xhatGAMP_y.svar;
    obj.options.GAMPOptY.xhatPrev0 = xhatGAMP_y.xhatPrev;
    obj.options.GAMPOptY.scaleFac = xhatGAMP_y.scaleFac;
    obj.options.GAMPOptY.step = min(max(xhatGAMP_y.step,0.05),xhatGAMP_y.stepMax);

    % Warm Start xz
    obj.options.GAMPOptZ.xhat0 = xhatGAMP_z.xhat;
    obj.options.GAMPOptZ.xvar0 = xhatGAMP_z.xvar;
    obj.options.GAMPOptZ.shat0 = xhatGAMP_z.shat;
    obj.options.GAMPOptZ.svar0 = xhatGAMP_z.svar;
    obj.options.GAMPOptZ.xhatPrev0 = xhatGAMP_z.xhatPrev;
    obj.options.GAMPOptZ.scaleFac = xhatGAMP_z.scaleFac;
    obj.options.GAMPOptZ.step = min(max(xhatGAMP_z.step,0.05),xhatGAMP_z.stepMax);
    
    % Markov Model on velocity indicator
    if ind >4
        v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat.*maxwellCorrX, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
        v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat.*maxwellCorrY, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);
        v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat.*maxwellCorrZ, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, gamma);

        v = max([v_x(:),v_y(:),v_z(:)],[],2);
%         v = max(reshape(v,size(obj.data.sampB)),[],3);
%         v = repmat(v,[1,1,size(obj.data.sampB,3)]);
%         
%         v(v<1e-6) = 1e-6;
%         v(v>=1) = 0.999;
%         obj.options.ReVEALOpts.gamma = v(:) ;
        
        % =================================================================
        % MRC
        v(v<1e-6) = 1e-6;
        v(v>=1) = 0.999;
        v = reshape(v,[size(obj.data.sampB,1)*size(obj.data.sampB,2),size(obj.data.sampB,3)]);
        [v, S_POST, p01_upd] = mrc.binarymarkov_msgs(v, sparsity);
        sparsity = mean(S_POST(mask ==1));
        sparsity = sparsity*ones(size(obj.data.sampB,1)*size(obj.data.sampB,2),size(obj.data.sampB,3));
        mrc.p01 = p01_upd;
        if p01_upd<1e-6; p01_upd=1e-6 ;end
        obj.options.ReVEALOpts.gamma = max(reshape(v,size(obj.data.sampB)),[],3);
        obj.options.ReVEALOpts.gamma = repmat(obj.options.ReVEALOpts.gamma,[1,1,size(obj.data.sampB,3)]);
        obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);
        gamma = v(:);
        fprintf(sprintf('\nEM Update p01 = %s, sparsity = %s\n',num2str(p01_upd),num2str(sparsity(1))));
    end
    
    % Clear some memory up
    clear v xhatGAMP_b xhatGAMP_x xhatGAMP_y xhatGAMP_z;    
    
    %======================================================================
    % Reconstruct
    fprintf('pB...')
    xhatGAMP_b = gampEst(EstimIn_b,EstimOut_b,Op_b,obj.options.GAMPOptB);
	xhatGAMP_b = clearOpts(xhatGAMP_b);
    fprintf('pX...')
    xhatGAMP_x = gampEst(EstimIn_x,EstimOut_x,Op_x,obj.options.GAMPOptX);
	xhatGAMP_x = clearOpts(xhatGAMP_x);
    fprintf('pY...')
    xhatGAMP_y = gampEst(EstimIn_y,EstimOut_y,Op_y,obj.options.GAMPOptY);
	xhatGAMP_y = clearOpts(xhatGAMP_y);
    fprintf('pZ...\n')
    xhatGAMP_z = gampEst(EstimIn_z,EstimOut_z,Op_z,obj.options.GAMPOptZ);
	xhatGAMP_z = clearOpts(xhatGAMP_z);
    fprintf('It = %d\t|dx_b|/|x_b| = %0.4f ----- |dx_v|/|x_v| = %0.4f\n', ind,1,1)
    
end
fprintf('ReVAMP Calculated in %s s\n', num2str(toc))

v_x = post_v(xhatGAMP_b.rhat, xhatGAMP_x.rhat.*maxwellCorrX, xhatGAMP_b.rvar, xhatGAMP_x.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
v_y = post_v(xhatGAMP_b.rhat, xhatGAMP_y.rhat.*maxwellCorrY, xhatGAMP_b.rvar, xhatGAMP_y.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);
v_z = post_v(xhatGAMP_b.rhat, xhatGAMP_z.rhat.*maxwellCorrZ, xhatGAMP_b.rvar, xhatGAMP_z.rvar, obj.options.ReVEALOpts.sigmaSq, obj.options.ReVEALOpts.gamma);

%% Calculate velocity maps
obj.outputs.thetaX = angle(xhatGAMP_x.xhat.*conj(xhatGAMP_b.xhat).*maxwellCorrX);
obj.outputs.thetaY = angle(xhatGAMP_y.xhat.*conj(xhatGAMP_b.xhat).*maxwellCorrY);
obj.outputs.thetaZ = angle(xhatGAMP_z.xhat.*conj(xhatGAMP_b.xhat).*maxwellCorrZ);

% Reshape outputs
obj.outputs.xHatb = reshape(xhatGAMP_b.xhat,size(obj.data.sampB));
obj.outputs.xHatx = reshape(xhatGAMP_x.xhat,size(obj.data.sampX));
obj.outputs.xHaty = reshape(xhatGAMP_y.xhat,size(obj.data.sampY));
obj.outputs.xHatz = reshape(xhatGAMP_z.xhat,size(obj.data.sampZ));

obj.outputs.thetaX = reshape(obj.outputs.thetaX,size(obj.data.sampB));
obj.outputs.thetaY = reshape(obj.outputs.thetaY,size(obj.data.sampB));
obj.outputs.thetaZ = reshape(obj.outputs.thetaZ,size(obj.data.sampB));

obj.outputs.vX = reshape(v_x,size(obj.data.sampB));
obj.outputs.vY = reshape(v_y,size(obj.data.sampB));
obj.outputs.vZ = reshape(v_z,size(obj.data.sampB));

%% fftshifts
obj.data.Yb = fftshift(fftshift(obj.data.Yb,1),2);
obj.data.Yx = fftshift(fftshift(obj.data.Yx,1),2);
obj.data.Yy = fftshift(fftshift(obj.data.Yy,1),2);
obj.data.Yz = fftshift(fftshift(obj.data.Yz,1),2);

obj.data.sampB = fftshift(fftshift(obj.data.sampB,1),2);
obj.data.sampX = fftshift(fftshift(obj.data.sampX,1),2);
obj.data.sampY = fftshift(fftshift(obj.data.sampY,1),2);
obj.data.sampZ = fftshift(fftshift(obj.data.sampZ,1),2);

obj.data.maxwellCorrX =  fftshift(fftshift(obj.data.maxwellCorrX,1),2);
obj.data.maxwellCorrY =  fftshift(fftshift(obj.data.maxwellCorrY,1),2);
obj.data.maxwellCorrZ =  fftshift(fftshift(obj.data.maxwellCorrZ,1),2);

obj.data.sensMaps = fftshift(fftshift(obj.data.sensMaps,1),2);
obj.data.x0 =  fftshift(fftshift(obj.data.x0,1),2);
obj.options.ReVEALOpts.gamma = fftshift(fftshift(obj.options.ReVEALOpts.gamma,1),2);
obj.options.ReVEALOpts.gamma = obj.options.ReVEALOpts.gamma(:);

obj.outputs.xHatb = fftshift(fftshift(obj.outputs.xHatb,1),2);
obj.outputs.xHatx = fftshift(fftshift(obj.outputs.xHatx,1),2);
obj.outputs.xHaty = fftshift(fftshift(obj.outputs.xHaty,1),2);
obj.outputs.xHatz = fftshift(fftshift(obj.outputs.xHatz,1),2);

obj.outputs.thetaX = fftshift(fftshift(obj.outputs.thetaX,1),2);
obj.outputs.thetaY = fftshift(fftshift(obj.outputs.thetaY,1),2);
obj.outputs.thetaZ = fftshift(fftshift(obj.outputs.thetaZ,1),2);

obj.outputs.vX = fftshift(fftshift(obj.outputs.vX,1),2);
obj.outputs.vY = fftshift(fftshift(obj.outputs.vY,1),2);
obj.outputs.vZ = fftshift(fftshift(obj.outputs.vZ,1),2);

end


% Clear out unneeded GAMP outputs to save memory
function GAMPOpt = clearOpts(GAMPOpt)
    GAMPOpt.phat = [];
    GAMPOpt.pvar = [];
    GAMPOpt.zhat = [];
    GAMPOpt.zvar = [];
    GAMPOpt.Axhat = [];
    GAMPOpt.xhatPrev = [];
    GAMPOpt.xhatNext = [];
    GAMPOpt.xvarNext = [];
    GAMPOpt.xhatDamp = [];
    GAMPOpt.pvarOpt = [];
    GAMPOpt.rvarOpt = [];
    GAMPOpt.A2xvarOpt = [];
    GAMPOpt.shatNext = [];
    GAMPOpt.svarNext = [];
end
