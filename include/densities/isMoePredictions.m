function [ X, Y, yDataGrid, yGrid, posteriorMedian, posteriorDataQuantiles] = isMoePredictions( inputObject, xData, yData, nX, nY)

    psiParticles = inputObject.psiParticles;
    thetaParticles = inputObject.thetaParticles;

    C = inputObject.partitions;
    gatingProbabilities = inputObject.gatingProbabilities;

    [ psiParticles, psiInds, countInds] = unique( psiParticles, "rows");
    thetaParticles = thetaParticles( psiInds, :);

    C = C( psiInds );
    gatingProbabilities = gatingProbabilities( psiInds );
    
    J = length( psiInds );
    K = inputObject.K;

    counts = accumarray( countInds, 1);

    xMin = min(xData);
    xMax = max(xData);
    
    yMin = min(yData);
    yMax = max(yData);
    
    yInterval = yMax - yMin;
    yPadding = 0.25 * yInterval;
    
    xStar = linspace( xMin, xMax, nX)';
    y0 = zeros( size(xStar) );
    
    yStar = linspace( yMin - yPadding, yMax + yPadding, nY)';
    
    yDataGrid = zeros( nY, nX);
    yGrid = zeros( nY, nX);

    for jj = 1:J
        
        psi_jj = psiParticles( jj, :);
        theta_jj = thetaParticles( jj, :);
        
        muInds = 1:K;
        sigmaInds = (K+1):2*K;
        piInds = (2*K + 1):3*K;
        
        mu_jj = psi_jj( muInds );
        sigma_jj = psi_jj( sigmaInds );
        pi_jj = psi_jj( piInds );
        
        C_jj = C{jj}.C;
        gatingProbability_jj = gatingProbabilities(jj);
        count_jj = counts(jj);

        data_jj = {};
        data_jj.x = xStar;
        data_jj.y = y0;
        data_jj.R = computeDistanceMatrix( xStar );

        [~, ~, classProbability, ~] = partition( data_jj, mu_jj, sigma_jj, pi_jj);
        
        for kk = 1:K

            C_jj_kk = C_jj{kk};
        
            if( isempty( C_jj_kk ) )
                continue;
            end
            
            x_kk = C_jj_kk(:,1);
            y_kk = C_jj_kk(:,2);
            
            classProbability_kk = classProbability(:,kk);
        
            thetaInds_kk = computeThetaInds( kk, K, 1);
            
            theta_jj_kk = theta_jj( thetaInds_kk );
            measurementSigma = theta_jj_kk(2);
            measurementVariance = measurementSigma.^2;
            
            [ yStar_m, yStarSigma_m] = createResultPredictions1D( x_kk, y_kk, xStar, theta_jj_kk);
        
            parfor ii = 1:nX

                mean_ii = yStar_m(ii);
                sigma_ii = yStarSigma_m(ii);
                dataPredictiveSigma_ii = sqrt( measurementVariance + sigma_ii.^2 );

                yPredictive_ii = normpdf( yStar, mean_ii, sigma_ii);
                yDataPredictive_ii = normpdf( yStar, mean_ii, dataPredictiveSigma_ii);
                
                yGrid(:,ii) = yGrid(:,ii) + count_jj * gatingProbability_jj * classProbability_kk(ii) .* yPredictive_ii;
                yDataGrid(:,ii) = yDataGrid(:,ii) + count_jj * gatingProbability_jj * classProbability_kk(ii) .* yDataPredictive_ii;
            end
        end
        
        progress = jj / J
    end
    
    yGridIntegral = sum( yGrid );
    yDataGridIntegral = sum( yDataGrid );
    
    for ii = 1:nX
        yGrid(:,ii) = yGrid(:,ii) / yGridIntegral(ii);
    end
    
    for ii = 1:nX
        yDataGrid(:,ii) = yDataGrid(:,ii) / yDataGridIntegral(ii);
    end

    [ X, Y] = meshgrid( xStar, yStar);
    
    posteriorMedian = computePosteriorMedian( yStar, yGrid );

    posteriorDataMedian = computePosteriorMedian( yStar, yDataGrid );
    posteriorDataLowerBound = computePosteriorQuantile( yStar, yDataGrid, 0.025);
    posteriorDataUpperBound = computePosteriorQuantile( yStar, yDataGrid, 0.975);

    posteriorDataQuantiles = {};
    posteriorDataQuantiles.median = posteriorDataMedian;
    posteriorDataQuantiles.lowerBound = posteriorDataLowerBound;
    posteriorDataQuantiles.upperBound = posteriorDataUpperBound;
end