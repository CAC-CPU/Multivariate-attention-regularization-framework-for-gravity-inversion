function showpar = loginv(A,d,inverpar2)

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
upbound = inverpar2.upbound;
lowbound = inverpar2.lowbound;





m_pre=zeros(num_cell,1);%previous model for ml
m_plu=upbound*ones(num_cell,1); m_min=lowbound*ones(num_cell,1);
dm=m_plu-m_min;
mapr_l=log((mapr-(m_min))./((m_plu)-mapr));

l=zeros(num_cell,1);%l(n)
ln1=zeros(num_cell,1);%l(n-1)
ll=zeros(num_cell,1);%conjugate l(n)
lln1=zeros(num_cell,1);%conjugate l(n-1)
II=ones(num_cell,1);

misfit=zeros(1,CN);
MF=zeros(1,CN);
Pa=zeros(1,CN);%Parameteric functional
apha=zeros(1,CN);
F=A;
ml=log((m-(m_min))./((m_plu)-m));
sb=diag(exp(ml)./(1+exp(ml)).^2);
F=F*dm(1)*sb;
Wm_diag=(diag(F'*F)).^0.25;
mw = Wm_diag.*ml;
maprw = Wm_diag.*mapr_l;

sm_fc=rew(1);
itw=1;
for  n=1:CN
    if n<=rew(1)%before the first re-weight point, do smooth inversion
        We_diag=ones(num_cell,1);
    elseif n==sm_fc+1
%         We_diag=(1./sqrt(diag((diag(m_pre))^2+diag(II)*e^2)));
        We_diag=(1./sqrt(m_pre.^2+II.*e^2));%Minimum support
        We_med_diag=We_diag;
        
        if itw<length(rew)
            itw=itw+1;
            sm_fc=rew(itw);
            fprintf('Reweight at iteration %3d \n',sm_fc);
        end     
        
    else We_diag=We_med_diag;
    end
    W_diag=Wm_diag.*We_diag;
%     inv_w_diag=diag(inv(diag(W_diag)));
    inv_w_diag=1./W_diag;
    mw=W_diag.*ml;
    maprw=W_diag.*mapr_l;
    gamma2=norm(mw-maprw)^2;
%     Fw=F*diag(inv_w_diag);
    Fw=F*spdiags(inv_w_diag,0,length(inv_w_diag),length(inv_w_diag));
    Aw=Fw;
    rn=A*m-d;
    %====================Determine a
   
    if n==1
        apha(1)=0;
    elseif n==2 
        ao=norm(A*m-d)^2/norm(mw-maprw)^2;
        ao_med=ao;
    else
        ao=ao_med;
    end
    apha(2)=ao*qq;
    if n>2
        if (gamma2/gamma1)>1
            apha(n)=apha(n-1)/(gamma2/gamma1);
        else
            apha(n)=apha(n-1)*qq;
        end
    end
    %===============================    
    misfit(n)=norm(rn)/norm(d);
    if misfit(n)<=final_misfit&&n>sm_fc+1 break; end
    MF(n)=(norm(rn))^2;
    Pa(n)=MF(n)+apha(n)*(mw-maprw)'*(mw-maprw);
    l=Fw'*rn+apha(n)*(mw-maprw);
    if n==1
        ll=l;
    else
        bet=(norm(l))^2/(norm(ln1))^2;
        ll=l+bet*lln1;
    end
    kn=dot(l,ll)/dot(ll,Fw'*Fw*ll+apha(n)*ll);
    %store last l ll:
    ln1=l;
    lln1=ll;
    m_pre=ml;
    mw=mw-kn*ll;
    gamma1=norm(mw-maprw)^2;
    ml=inv_w_diag.*mw;
    m=(m_min(1)+m_plu(1)*exp(ml))./(1+exp(ml));
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