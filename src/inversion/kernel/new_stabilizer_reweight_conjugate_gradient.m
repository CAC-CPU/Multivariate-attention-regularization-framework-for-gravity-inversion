function inversion = new_stabilizer_reweight_conjugate_gradient(inversion,A)

m = inversion.mInitVector;
mapr = inversion.mPriorVector;
d = inversion.dataObs;
CN = inversion.maxIteration;
fmisfit = inversion.fmisfit;
misfit = inversion.misfit;
alpha = inversion.alpha;
gamma = inversion.gamma;
q = inversion.decayCoefficient;
e = inversion.focusFactor;
nx = inversion.xGridsNum;
nz = inversion.zGridsNum;
gridsNum = inversion.gridsNum;
reweight = inversion.reweightIndex;
theta = inversion.theta;
plotStablizersIndex = inversion.plotStablizersIndex;

reweightBegin = reweight(1);
reweightEnd = reweight(end);
F = A;
Wm = diag(F'*F).^0.25;
tic;
for  n=1:CN

    if n > 2 && norm(mapr)^2 - norm(m)^2 > 0.1*norm(m)^2
        theta(n) = norm(mapr./mw)/((norm(mapr)^2 - norm(m)^2)/norm(m)^2)^2;
        mapr = mapr/theta(n);
    end

    if n <= reweightBegin || reweightBegin == 0
        We = ones(gridsNum,1);

    elseif n <= reweightEnd + 1
        we = 1./sqrt(mPre.^2+e^2);
        spatiallyAdaptiveTerm = 1./(sqrt((mPre - mapr).^2) + e);
        temporalDecayTerm = 1;
        We = spatiallyAdaptiveTerm*temporalDecayTerm;
        fprintf('Reweight at iteration %3d \n',n)
    end

    if n == plotStablizersIndex
        minimumSupport = we;
        newStablizer = We;
    end

    W = Wm.*We;
    invW = 1./W;
    mw = W.*m;
    maprw = W.*mapr;
    snw = mw - maprw;
    Fw = F*spdiags(invW,0,length(invW),length(invW));
    Aw = Fw;
    rn = Aw*mw-d;

    if n == 2
        alpha(n) = norm(rn)^2/norm(snw)^2;
    end
    if n > 2
        gamma(n) = norm(snw)^2/norm(snwPre)^2;
        if gamma(n) > 1
             alpha(n) = alpha(n-1)/gamma(n);
        else
            alpha(n) = alpha(n-1)*q;
        end
    end
    
    misfit(n) = norm(rn)/norm(d);
    if misfit(n) <= fmisfit && n > reweightEnd + 1; break; end
    ln = Fw'*rn + alpha(n)*snw;
    if n == 1
        lln = ln;
    else
        beta = (norm(ln))^2/(norm(lnPre))^2;
        lln = ln + beta*llnPre;
    end
    kn = (lln'*ln)/(norm(Fw*lln)^2 + alpha(n)*norm(lln)^2);
    lnPre = ln;
    llnPre = lln;
    
    mPre = m;
    mw = mw - kn*lln;
    snwPre = mw - maprw;
    m = invW.*mw;
    fprintf('Iterration %d\tMisfit = %.8f\n',n,misfit(n));    
end
toc;
fprintf('Iterration %d\tMisfit = %.8f\n',n,misfit(n));
inversion.mPredictBuffuer = reshape(m,nz,nx);
inversion.mPredict = inversion.mPredictBuffuer(:,inversion.xIndexPlot);
inversion.dataPredict = A*m;
inversion.minimumSupport = reshape(minimumSupport,nz,nx);
inversion.newStablizer = reshape(newStablizer,nz,nx);
end