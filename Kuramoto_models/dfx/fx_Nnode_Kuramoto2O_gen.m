function [ystore tvec rlxtime] = fx_Nnode_Kuramoto_gen(dt,tt,Nr,K,A,omega,sigma,y,D,burn)
t = 0; hn = max(D);
stvec  = normrnd(0,sigma,Nr,tt);
for i = size(y,2):tt
    dydt = kuramoto_d_fx(y,A,K,Nr,D,omega,i);
    y(:,i+1) = wrapTo2Pi(y(:,i) +(dt*dydt') + (stvec(:,i).*dt));
    ystore(:,i) = cos(y(:,i));
    t = t+dt;
    tvec(i) = t;
end

ystore = ystore(:,burn:end);
tvec = tvec(burn:end);

ypert =y(:,end-(1/dt):end);
comblist = nchoosek(1:Nr,2)
phdiff = [];
for i = 1:size(comblist,1)
    phdiff(i,:) = diff(unwrap(ypert(comblist(i,:),:)),1,1); % RP
    cmbname{i} = sprintf('%.0f : %.0f',comblist(i,:));
end
fpnt = mean(diff(phdiff,1,2),2); % Average dRP over channel pairs
subplot(2,1,1)
 plot(diff(phdiff,1,2)')
legend(cmbname)

% Kick it (do this individuallly!)
pertind = size(ypert,2)-floor(0.01/dt):size(ypert,2);
ypert(1,pertind) = repmat(ypert(1,end) + pi/2,size(pertind));
postind = size(ypert,2):size(ypert,2)+(5/dt);
for i = postind
    dydt = kuramoto_d_fx(ypert,A,K,Nr,D,omega,i);
    ypert(:,i+1) = wrapTo2Pi(ypert(:,i) +(dt*dydt') + (stvec(:,i).*dt));
end
ypert_post = ypert(:,postind);
phdiff_post = [];
for i = 1:size(comblist,1)
    phdiff_post(i,:) = diff(unwrap(ypert_post(comblist(i,:),:)));
end
rlx = sum(unwrap(diff(phdiff_post,1,2))-fpnt);
subplot(2,1,2)
 plot(unwrap(diff(phdiff_post,1,2))'); hold on; plot(rlx','k--','LineWidth',2); ylim([-0.25 0.25])
legend(cmbname)
if median(fpnt)<0.006
    rlxtime = min(find(rlx<0.001))*dt;
    if isempty(rlxtime)
        rlxtime = -inf;
    end
else
    rlxtime = inf;
end
close all

function dydt = kuramoto_d_fx(y,A,K,Nr,D,omega,i)
    for r = 1:Nr
        for d=1:Nr
            yp(d) = y(d,end-D(r,d));
        end
         I = (K/Nr).*sum( A(r,:).*(sin(yp-y(r,i)) + 0.25*sin(2*(yp-y(r,i)))) );
        dydt(r) = omega(r) + I;
    end