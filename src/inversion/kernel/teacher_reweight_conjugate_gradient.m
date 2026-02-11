function inversion = teacher_reweight_conjugate_gradient(inversion,A)

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

reweightBegin = reweight(1);
reweightEnd = reweight(end);
F = A;
Wm = diag(F'*F).^0.25;
tic;
for  n=1:CN
    if n <= reweightBegin || reweightBegin == 0
        We = ones(gridsNum,1);
    elseif n <= reweightEnd + 1
        We = 1./sqrt(mPre.^2+e^2);
        fprintf('Reweight at iteration %3d \n',n)
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
        alpha(n) = norm(Aw*mw - d)^2/norm(snw)^2;
    end

    if n>2
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
    kn = dot(ln,lln)/dot(lln,Fw'*Fw*lln + alpha(n)*lln);
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
end