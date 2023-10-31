function [prior] = constructPrior( data, settings)

    minX = 0;
    maxX = 1;

    minY = 0;
    maxY = max( data.y );
    
    sdSigmaGP = 0.25 * (maxY - minY);
    sdLengthScalesGP = 0.125 * (maxX - minX);
    sdSignalSigmaGP = 0.25 * (maxY - minY);
    
    prior = {};
    prior.psi = {};
    prior.theta = {};

    prior.psi.samplers = {};
    prior.theta.samplers = {};

    prior.psi.samplers.mu = @sampleGatingMu;
    prior.psi.samplers.sigma = @sampleGatingSigma;
    prior.psi.samplers.pi = @sampleGatingPi;
        
    prior.psi.parameters.muSigma = settings.gatingMuSigma;
    prior.psi.parameters.sigmaSigma = settings.gatingSigmaSigma;
    prior.psi.parameters.concentrationParameter = settings.concentrationParameter;

    prior.theta.samplers.mu = @( N ) maxY * rand( N, 1);
    prior.theta.samplers.noiseSigma = @( N ) abs( sdSigmaGP * randn(N,1) );
    prior.theta.samplers.lengthScale = @( N ) abs( sdLengthScalesGP * randn(N,1) );
    prior.theta.samplers.signalSigma = @( N ) abs( sdSigmaGP * randn(N,1) );

    prior.theta.parameters.noiseSigmaSigma = sdSigmaGP;
    prior.theta.parameters.lengthScaleSigma = sdLengthScalesGP;
    prior.theta.parameters.signalSigmaSigma = sdSignalSigmaGP;

    prior.theta.logPrior = @thetaLogPrior;
end