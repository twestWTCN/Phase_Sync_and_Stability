function [ystore tvec rlxtime] = fx_Nnode_Kuramoto2O_gen(dt,tt,Nr,K,A,omega,sigma,y,D,burn)
t = 0; hn = max(D);
stvec  = sqrt(dt).*normrnd(0,sigma,Nr,tt);
for i = size(y,2):tt
    dydt = dfx_kuramoto(y,A,K,Nr,D,omega,i);
    y(:,i+1) = (y(:,i) +(dt*dydt') + (stvec(:,i))); %wrapTo2Pi
    ystore(:,i) = (y(:,i));
    t = t+dt;
    tvec(i) = t;
end

ystore = ystore(:,burn:end);
tvec = tvec(burn:end);



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
    y1_post = y1(:,postind(1)+floor(epst/dt):end);
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

