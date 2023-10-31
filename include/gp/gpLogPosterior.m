function [logPosterior] = gpLogPosterior( partition, theta, thetaLogPrior, param)

    K = length( partition.C );

    logLikelihood = gpExpertsLogLikelihood( partition, theta, 1);
    logPrior = 0;
        
    for kk = 1:K

        inds_kk = computeThetaInds( kk, K, 1);
        theta_kk = theta(inds_kk);
            
        logPrior = logPrior + thetaLogPrior( theta_kk, [], K, param);
    end
    
    logPosterior = logLikelihood + logPrior;
end