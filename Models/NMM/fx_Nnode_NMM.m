function [ystore tvec rlxtime] = fx_Nnode_NMM(dt,tt,Nr,K,A,par,H,T,sigma,sigA,y,D,burn)
t = 0; hn = max(D);
% stvec  = sigA'.*(sqrt(dt).*normrnd(0.0,sigma,Nr,tt)); % noise input
stvec  = sigA'.*(sqrt(dt).*normrnd(0,0.5,Nr,tt)); % noise input
K = 200;
for i = size(y,3):tt
    dydt = dfx_NMM(y,A,K,Nr,D,par,H,T,stvec,i);
    y(:,:,i+1) = (y(:,:,i) +(dt*dydt) + [stvec(:,i) zeros(Nr,1)]); %wrapTo2Pi
    ystore(:,i) = (y(:,2,i));
    t = t+dt;
    tvec(i) = t;
end
close all

ystore = ystore(:,burn:end);
ystore = ystore-mean(ystore,2);
tvec = tvec(burn:end);
plot((repmat(linspace(0,size(ystore,2)*dt,size(ystore,2)),4,1)./1000)',ystore'); shg

R.frqz = 1:1:100;
close all
[F meannpd] = constructNPDMat(ystore,{'A','B','C','D'},{'A','B','C','D'},(1/dt)*1000,12,R)
npdplotter_110717({meannpd},{},F,R,1)
%% PERTUBATION(S)
% get end of sim (approx steady state)
y1 =y(:,end-(1/dt):end);
% Get possible combinations of channels
comblist = nchoosek(1:Nr,2);
eps = (pi/2); % pertubation size
epst = 0.25; % pertuation duration
simlength = 10;
for n = 1:Nr
    nrcmbind = find(sum(comblist==1,2));
    postind = size(y1,2):size(y1,2)+(simlength/dt);
    for i = postind
        dydt = dfx_kuramoto(y1,A,K,Nr,D,omega,i);

        y1(:,i+1) = (y1(:,i) +(dt*dydt') + (stvec(:,i)));
         if i<postind(1)+(epst/dt)
           y1(n,i+1) = y1(n,i+1) + eps;
        end
    end
    y1_post = y1(:,postind(1)+(epst/dt):end);
    phdiff_post = [];
    for i = nrcmbind'
        phdiff_post(i,:) = diff(unwrap(y1_post(comblist(i,:),:)));
    end
    rlx = diff(unwrap(phdiff_post),1,2);
    rlxt = [];
    for i = 1:3
        rlx(i,abs(rlx(i,:))>1.95*pi) = median(rlx(i,:));
        a = min(find(abs(rlx(i,:))<0.005));
        if isempty(a) || abs(median(rlx(i,1:10))-mean(rlx(i,:)))<0.005; a = NaN; end
        rlxt(i) = a;
    end
    rlxt = rlxt*dt;
    rlxtime(n) = nanmean(rlxt);
end

% rlxtime = nanmean(rlxtime);

% close all

