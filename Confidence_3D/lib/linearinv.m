function showpar = linearinv(A,d,inverpar2)

% A: forward operator matrix which also server as the sensitivity matrix
% d: observed data for potential field
% inverpar2: inversion parameter
% showpar: output for inversion result

ao = inverpar2.InitialReg;
final_misfit = inverpar2.tol;
CN = inverpar2.maxiteration;
rew = inverpar2.reweight;
qq = inverpar2.DecayCoef;
e = inverpar2.focus;
num_cell = size(A,2);
m = inverpar2.mstart;
mapr = inverpar2.mapr;
nx = inverpar2.numx;
ny = inverpar2.numy;
nz = inverpar2.numz;
cut_block = inverpar2.cut_block;
istraditioninv = inverpar2.istraditioninv;


m_pre=zeros(num_cell,1);%previeous model;
l=zeros(num_cell,1);
ln1=zeros(num_cell,1);
ll=zeros(num_cell,1);
lln1=zeros(num_cell,1);


II=ones(num_cell,1);
misfit=zeros(1,CN);
MF=zeros(1,CN);
Pa=zeros(1,CN);%Parameteric functional
% apha=diag(zeros(1,CN));
apha = 0;
r=zeros(length(d),CN);
F=A;
Wm_diag=(diag(F'*F)).^0.25;
% Wm_diag=((sum(F.^2)).^0.25)';
mw = Wm_diag.*m;
maprw = Wm_diag.*mapr;

sm_fc=rew(1);
itw=1;
for  n=1:CN
    if n<=rew(1) || rew(1) == 0 %before the first re-weight point, do smooth inversion
        We_diag=ones(num_cell,1);
        %We=diag(We_diag)
    elseif n==sm_fc+1
        %         We_diag=(1./sqrt(diag((diag(m_pre))^2+diag(II)*e^2)));
        We_diag=(1./sqrt(m_pre.^2+II.*e^2));%Minimum support
        ep=0.01*e;
        %         We_diag=getw(m,nx,ny,nz,ep,[1 2 3]);%Minimum gradient support
        We_med_diag=We_diag;

        if itw<length(rew)
            itw=itw+1;
            sm_fc=rew(itw);
            fprintf('Reweight at iteration %3d \n',sm_fc)
        end

    else We_diag=We_med_diag;
    end
    W_diag=Wm_diag.*We_diag;
    %     inv_w_diag=diag(inv(diag(W_diag)));
    inv_w_diag=1./W_diag;
    mw=W_diag.*m;
    maprw=W_diag.*mapr;
    tmd2=norm(mw-maprw)^2;
    %     Fw=F*diag(inv_w_diag);
    Fw=F*spdiags(inv_w_diag,0,length(inv_w_diag),length(inv_w_diag));
    Aw=Fw;
    rn=Aw*mw-d;
    %====================Determine a
    %     if n==1
    %         apha(1)=0;
    %     end
    %     if n==2
    %         aow=norm(Aw*mw-d)^2/norm(mw-maprw)^2;
    %         apha(2)=aow;
    %     end
    %
    %     if n>2
    %         if (tmd2/tmd)>1
    %              apha(n)=apha(n-1)/(tmd2/tmd);
    % %             tmd2=tmd;
    %         else
    %             apha(n)=apha(n-1);
    %            %tmd=tmd2;
    %         end
    %        if (norm((Aw*(W_diag.*m_pre))-d).^2-norm(rn)^2)/ norm((Aw*(W_diag.*m_pre))-d).^2 < 0.001
    %            apha(n)=apha(n)*qq;
    %        end
    %     end
    
    if istraditioninv
        apha =  inverpar2.apha;
    else
        lambda=norm(Aw*mw-d)^2/norm(mw-maprw)^2;
        if n<=2
            [credibility_first,~,~,gamma] = generate_credibility(Aw,d,mw,maprw,n);
            apha = lambda .* diag(credibility_first)*gamma;
        else
            [credibility_history(n,:),~,~,gamma] = generate_credibility(Aw,d,m_pre,maprw,n);
            apha = lambda .* diag(credibility_history(n,:))*gamma;
        end
    end
    Apha(n) = max(diag(apha));
    %===============================
    misfit(n)=norm(rn)/norm(d);
    if misfit(n)<=final_misfit&&n>sm_fc+1 break; end  %when n<=sm_fc+1, even misfit is low, don't go out of iteration!
    MF(n)=(norm(rn))^2;
    Pa(n)=MF(n)+(mw-maprw)'*apha*(mw-maprw);
    l=Fw'*rn+apha*(mw-maprw);
    if n==1
        ll=l;
    else
        bet=(norm(l))^2/(norm(ln1))^2;
        ll=l+bet*lln1;
    end
    kn=dot(l,ll)/dot(ll,Fw'*Fw*ll+apha*ll);
    %store last l ll:
    ln1=l;
    lln1=ll;

    m_pre=m;
    mw=mw-kn*ll;
    tmd=norm(mw-maprw)^2;
    m=inv_w_diag.*mw;
    %     index=find(m<0);
    %     m(index)=0;
    %     save inversion_result m;
    fprintf('Iteration %3d    Misfit=%8.5f\n',n,misfit(n));
end

showpar.Model3d = reshape_back_cut(m,nx,ny,nz,cut_block);
showpar.ModelVector = showpar.Model3d(:);
% 添加先验模型到showpar结构体
showpar.PriorModel3d = reshape_back_cut(mapr,nx,ny,nz,cut_block);
showpar.PriorModelVector = showpar.PriorModel3d(:);
showpar.numx = nx;
showpar.numy = ny;
showpar.numz = nz;
showpar.x = inverpar2.x;
showpar.y = inverpar2.y;
showpar.z = inverpar2.z;
showpar.observed = d;
showpar.predicted = A*m;
showpar.num_iteration = n-1;
showpar.percentmisfit = misfit;
showpar.regpar = apha;
showpar.residual = MF;
showpar.objectfun = Pa;
showpar.component = inverpar2.component;
showpar.xobserved = inverpar2.xobserved;
showpar.yobserved = inverpar2.yobserved;
showpar.zobserved = inverpar2.zobserved;
showpar.Xtopography = inverpar2.Xtopography;
showpar.Ytopography = inverpar2.Ytopography;
showpar.Ztopography = inverpar2.Ztopography;
showpar.cutoff = inverpar2.cutoff;
showpar.Apha = Apha;