function [ data, settings] = preprocess( data )

    nInputDimensions = size( data.x, 2);
    
    for ii = 1:nInputDimensions
        
        x_ii = data.x(:,ii);
        
        minX = min( x_ii );
        maxX = max( x_ii );
    
        data.x(:,ii) = ( x_ii - minX ) / ( maxX - minX );
    end
    
    data.R = computeDistanceMatrix( data.x );
    
    meanY = mean( data.y );

    data.y = data.y - meanY;
    data.y = data.y / std( data.y );
    data.y = data.y - min( data.y );

    maxY = max( data.y );
        
    settings = struct();
    settings.particleAmount = 100;
    settings.nInputDimensions = nInputDimensions;
    settings.numberOfExperts = 7;
        
    settings.prior.psi.muSigma = 0.25;
    settings.prior.psi.sigmaSigma = 0.25;
    settings.prior.psi.concentrationParameter = 0.5 * settings.numberOfExperts;

    settings.prior.theta.maxY = maxY;
    settings.prior.theta.noiseSigmaSigma = 0.25 * maxY;
    settings.prior.theta.lengthScaleSigma = 0.125;
    settings.prior.theta.signalSigmaSigma = 0.25 * maxY;
end