function [fStar, fStarSigma] = createResultPredictions1D( x, f, xStar, theta)

    nData = length(f);
    nPredictions = length(xStar);
    
    R_xx = ones( nData, nData);
    R_xSx = ones( nPredictions, nData);
    R_xSxS = ones( nPredictions, nPredictions);
    
    for ii = 1:nData
        for jj = 1:nData

            x_ii = x(ii);
            x_jj = x(jj);

            R_xx( ii, jj) = abs( x_ii - x_jj );
        end
    end
    
    for ii = 1:nPredictions
        for jj = 1:nData

            x_ii = xStar(ii);
            x_jj = x(jj);

            R_xSx( ii, jj, 1) = abs( x_ii - x_jj );
        end
    end

    for ii = 1:nPredictions
        for jj = 1:nPredictions

            x_ii = xStar(ii);
            x_jj = xStar(jj);

            R_xSxS( ii, jj) = abs( x_ii - x_jj );
        end
    end
    
    mu = theta(1);
    s2 = theta(2).^2;
    ls = theta(3:end-1);
    s2F = theta(end).^2;

    K_xSxS = s2F * sqrExpCovMatrix( R_xSxS, ls);
    K_xSx = s2F * sqrExpCovMatrix( R_xSx, ls);
    K_xx = s2F * sqrExpCovMatrix( R_xx, ls) + s2 * eye(nData);

    [fStar, fStarSigma] = gpPrediction( f - mu, K_xSxS, K_xSx, K_xx );
    
    fStar = fStar + mu;
end