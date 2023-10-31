function [logP] = thetaLogPrior( theta_jj, kk, K, settings)
    
    nInputDimensions = settings.nInputDimensions;
    lengthScaleInds = 3:(3 + nInputDimensions - 1);
    
    noiseSigmaSigma = settings.prior.theta.noiseSigmaSigma;
    lengthScaleSigma = settings.prior.theta.lengthScaleSigma;
    signalSigmaSigma = settings.prior.theta.signalSigmaSigma;

    noiseSigma = theta_jj(2);
    lengthScales = theta_jj( lengthScaleInds );
    signalSigma = theta_jj(end);
    
    logP = 0;
    
    logP = logP - 0.5 * ( noiseSigma / noiseSigmaSigma ).^2;
    logP = logP - 0.5 * log( noiseSigmaSigma );
    logP = logP - 0.5 * log( 2*pi );
    
    logP = logP - 0.5 * ( signalSigma / signalSigmaSigma ).^2;
    logP = logP - 0.5 * log( signalSigmaSigma );
    logP = logP - 0.5 * log( 2*pi );
    
    for ii = 1:nInputDimensions
    
        lengthScale_ii = lengthScales(ii);
    
        logP = logP - 0.5 * ( lengthScale_ii / lengthScaleSigma ).^2;
        logP = logP - 0.5 * log( lengthScaleSigma );
        logP = logP - 0.5 * log( 2*pi );
    end
end