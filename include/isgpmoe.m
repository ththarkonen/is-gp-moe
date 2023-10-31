function [resultObject] = isgpmoe( data, settings)

    J = settings.particleAmount;
    K = settings.numberOfExperts;
    N = settings.nInputDimensions;
        
    gatingMus = zeros( J, N * K);
    gatingSigmas = zeros( J, N * K);
    gatingPis = sampleGatingPi( J, [], K, settings);
        
    gpMeans = zeros( J, K);
    gpNoiseSigmas = zeros( J, K);
    gpLengthScales = zeros( J, K);
    gpSignalSigmas = zeros( J, K);
        
    counter = 1;
    
    for kk = 1:K
        
        mu_kk = sampleGatingMu( J, kk, K, settings);
        
        inds = counter:( counter + N - 1 );
        gatingMus( :, inds) = mu_kk;
        
        for n = 1:N
        
            sigma_kk = sampleGatingSigma( J, kk, K, settings);
            gatingSigmas( :, counter) = sigma_kk;
    
            counter = counter + 1;
        end
    end
    
    psiParticles = [ gatingMus, gatingSigmas, gatingPis];
    
    counter = 1;
    
    for kk = 1:K
    
        gpMu_kk = sampleExpertMean( J, kk, K, settings);
        gpSigma_k = sampleExpertNoiseSigma( J, kk, K, settings);
        gpSignalSigma_k = sampleExpertSignalSigma( J, kk, K, settings);
    
        gpMeans(:,kk) = gpMu_kk;
        gpNoiseSigmas(:,kk) = gpSigma_k;
        gpSignalSigmas(:,kk) = gpSignalSigma_k;
        
        for n = 1:N
            
            gpLengthScale_kk = sampleExpertLengthScale( J, kk, K, settings); 
            gpLengthScales( :, counter) = gpLengthScale_kk;
            
            counter = counter + 1;
        end 
    end
    
    thetaInitial = [ gpMeans, gpNoiseSigmas, gpLengthScales, gpSignalSigmas];
    
    partitions = cell( J, 1);
    gatingProbabilities = zeros( J, 1);
    
    theta = zeros( J, 4 * K);
    logLikelihoods = -Inf( J, 1);
    
    likelihoodEvaluations = 0;

    lfbgOptions = struct();
    lfbgOptions.GradObj = 'off';
    lfbgOptions.Display = 'off';
    lfbgOptions.LargeScale = 'off';

    lfbgOptions.HessUpdate = 'bfgs';
    lfbgOptions.InitialHessType = 'identity';
    lfbgOptions.GoalsExactAchieve = 1;
    lfbgOptions.GradConstr = 'false';

    tic
    parfor jj = 1:J
    
    
        gatingMus_jj = gatingMus( jj, :);
        gatingSigmas_jj = gatingSigmas( jj, :);
        gatingPis_jj = gatingPis( jj, :);
        
        theta_jj = thetaInitial( jj, :);
        [ partition_jj, gatingProbability_jj] = partition( data, gatingMus_jj, gatingSigmas_jj, gatingPis_jj);
    
        costFun = @(theta) -gpLogPosterior( partition_jj, theta, @thetaLogPrior, settings);
    
        try
            [thetaMAP, negLogPosterior, ~, output] = fminlbfgs( costFun, theta_jj, lfbgOptions);
            likelihoodEvaluations = likelihoodEvaluations + output.funccount;
        catch
            dims = size( theta_jj);
            thetaMAP = NaN( dims );
            negLogPosterior = Inf;
        end
    
        partitions{ jj } = partition_jj;
        gatingProbabilities( jj ) = gatingProbability_jj;
    
        thetaMAP = abs( thetaMAP );
        thetaMAP_mean = thetaMAP(1);
        
        validMean = ( 0 <= thetaMAP_mean );
        validMean = validMean && ( thetaMAP_mean <= settings.prior.theta.maxY );
        
        if( ~validMean )
            negLogPosterior = Inf;
        end
        
        theta( jj, :) = thetaMAP;
        logLikelihoods( jj, : ) = -negLogPosterior;
        
        jj
    end
    wallTime = toc;
    
    logLikelihoods = logLikelihoods - max( logLikelihoods );
    weights = exp( logLikelihoods );
    weights = weights / sum( weights );
    
    resampleInds = resampleResidual( weights );
    
    theta = theta( resampleInds, :);
    psiParticles = psiParticles( resampleInds, :);
    
    partitions = partitions( resampleInds );
    gatingProbabilities = gatingProbabilities( resampleInds );
    
    resultObject = struct();
    resultObject.psiParticles = psiParticles;
    resultObject.thetaParticles = theta;
    
    resultObject.partitions = partitions;
    resultObject.gatingProbabilities = gatingProbabilities;
    
    resultObject.J = J;
    resultObject.K = K;
    
    nX = 1024;
    nY = 1024;
    
    [ X, Y, yDataGrid, yGrid, posteriorMedian, posteriorDataQuantiles] = isMoePredictions( resultObject, data.x, data.y, nX, nY);        
    
    resultObject.X = X;
    resultObject.Y = Y;
    resultObject.data = data;
    
    resultObject.yGrid = yGrid;
    resultObject.yDataGrid = yDataGrid;
    resultObject.posteriorMedian = posteriorMedian;
    resultObject.posteriorDataQuantiles = posteriorDataQuantiles;
    
    resultObject.wallTime = wallTime;
    resultObject.likelihoodEvaluations = likelihoodEvaluations;
end