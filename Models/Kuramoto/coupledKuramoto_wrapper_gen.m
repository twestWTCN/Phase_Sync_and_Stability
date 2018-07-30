function [PLV dRPvar MsKappa LHat LVar RPvar rlxtime SRPeps] = coupledKuramoto_wrapper_gen(R,Klist,dvar,omvar,sigvar,dflag,surflag)
if nargin<7
    surflag=0;
end

dt = R.dt;
tend = R.tend;
tt = tend./dt;
fsamp = 1/dt;
burn = R.burn;
cfreq = R.cfreq;
period = R.period;
omega = randnbetween(cfreq,omvar,4,1);

A =[0 1 1 1;
    1 0 1 1
    1 1 0 1
    1 1 1 0];
A = rand(4).*A;
Nr = size(A,2);

if dflag == 1
    Da = floor(abs(randnbetween(0.25,dvar,sum(A(:)~=0))).*fsamp);
    D = [0      Da(1)   Da(2)   Da(3)
        Da(4)   0       Da(5)   Da(6)
        Da(7)   Da(8)   0       Da(9)
        Da(10)  Da(11)  Da(12)  0];
    D(D==0) = 1;
    y = [1*pi; 2*pi; -1*pi/2; pi/2];
    y = repmat(y,1,max(Da)+1);
else
    D = zeros(Nr);
    y = [1*pi; 2*pi; -1*pi/2; pi/2];
end


sppart = 5; L = 0;
for i=1:numel(Klist)
    [ystore{1} tvec{1} rlxtime(:,i)] = fx_Nnode_Kuramoto2O_gen(dt,tt,Nr,Klist(i),A,omega,sigvar,y,D,burn);
    if surflag ==1
        SRPeps(i) = 0.005; % phase_SurrStats(ystore{i},200,1)
    else
        SRPeps(i) = R.SRPeps;
    end
    
    a(1,:) = ystore{1}(1,:)-ystore{1}(2,:);
    a(2,:) = ystore{1}(2,:)-ystore{1}(3,:);
    a(3,:) = ystore{1}(3,:)-ystore{1}(4,:);
    a(4,:) = ystore{1}(4,:)-ystore{1}(1,:);
    a(5,:) = ystore{1}(4,:)-ystore{1}(2,:);
    a(6,:) = ystore{1}(3,:)-ystore{1}(1,:);
    for p = 1:size(a,1); PLV(p,i) = abs(mean(exp(1i*a(p,:)),2)); end
    for p = 1:size(a,1); dRPvar(p,i) = std(diff(a(p,:))); end
    for p = 1:size(a,1)
    [SRPlen{p} segRP{p}] =  SRP_Lengths(a(p,:),diff(a(p,:)),SRPeps(i),fsamp,period);
    end
    for p = 1:size(a,1); MsKappa(p,i) = nansum(SRPlen{p})./diff(tvec{1}([1 end])); end
    for p = 1:size(a,1); LHat(p,i) = nanmean(log(SRPlen{p})); end
    for p = 1:size(a,1); LVar(p,i) = nanstd(log(SRPlen{p}));
    
    end
    for p = 1:size(a,1); segRP{p} = segRP{p}(segRP{p}~=0);
        try RPvar(p,i) = sqrt(circ_var(segRP{p}'));  catch   RPvar(p,i) = NaN; disp('No Segments'); end;        end
    
    
    dRPvar = wghtcull(dRPvar);
    MsKappa = wghtcull(MsKappa);
    LHat = wghtcull(LHat);
    LVar = wghtcull(LVar);
    RPvar = wghtcull(RPvar);
    rlxtime = wghtcull(rlxtime);
    %         if rem(i,sppart) == 0
    %             L = L+1;
    %             figure(1)
    %             subplot(size(Klist,2)./sppart,2,(2*L)-1); plot(ystore{i}(1,:),ystore{i}(2,:));
    %             subplot(size(Klist,2)./sppart,2,(2*L)); plot(tvec{i},ystore{i}([1 2 3 4],:));
    %         end
end
a = 1;
function x = wghtcull(x)
colcull = sum(isnan(x))>(size(x,1)/2);
x(:,colcull) = NaN;

