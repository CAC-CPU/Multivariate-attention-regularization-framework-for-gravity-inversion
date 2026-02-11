function inversion = trandition_reweight_conjugate_gradient_V1(inversion,A)

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
reweightGap = inversion.reweightGap;

F = A;
Wm = diag(F'*F).^0.5;
Wd = diag(F*F').^0.5;

for n = 1:CN

    if n == 3 || mod(n-3,reweightGap) == 0
        We = 1./sqrt((m+mapr).^2 + e^2);
    else
        We = ones(size(m));
    end

    W = Wm.*We;
    invW = 1./W;

    Fw = Wd.*F*spdiags(invW,0,length(invW),length(invW));
    Aw = Wd.*A*spdiags(invW,0,length(invW),length(invW));

    dw = Wd.*d;

    mw = W.*m;
    maprw = W.*mapr;

    rn = norm(Wd.*A*m - Wd.*d)^2;
    rnw = Aw*mw - dw;
    snw = mw - maprw;

    misfit(n) = norm(rnw)/norm(dw);

    if n == 2
        alpha(n) = norm(rnw)^2/norm(snw)^2;
    end
    if n == 3 || (n > 3 && mod(n-3,reweightGap) == 0)
        gamma(n) = norm(snw)^2/norm(snwPre)^2;
        if gamma(n) > 1
            alpha(n) = alpha(n - 1)/gamma(n);
        elseif norm(rnPre)^2 - norm(rn)^2 < 0.01*norm(rnPre)^2 && q < 1
            alpha(n) = alpha(n - 1)*q;
        end
    end

    if misfit(n) <= fmisfit; break; end

    ln = Fw'*rnw + alpha(n)*snw;

    if n == 1
        lln = ln;
    else
        beta = norm(ln)^2/norm(lnPre)^2;
        lln = ln + beta*llnPre;
    end

    kn = (lln'*ln)/(norm(Fw*lln)^2 + alpha(n)*norm(lln)^2);

    lnPre = ln;
    llnPre = lln;
    rnPre = Wd.*A*m - Wd.*d;
    snwPre = mw - maprw;

    mw = mw - kn*lln;
    m = invW.*mw;
    
    fprintf('Iterration %d\tMisfit = %.8f\n',n,misfit(n));
end
fprintf('Iterration %d\tMisfit = %.6f\n',n,misfit(n));
inversion.mPredictBuffuer = reshape(m,nz,nx);
inversion.mPredict = inversion.mPredictBuffuer(:,inversion.xIndexPlot);
inversion.dataPredict = A*m;
end