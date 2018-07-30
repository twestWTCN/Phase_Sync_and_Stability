function rlxtime = odepertubation(


ypert =y(:,end-(20/dt):end);
comblist = nchoosek(1:Nr,2);
phdiff = [];
for i = 1:size(comblist,1)
    phdiff(i,:) = diff(unwrap(ypert(comblist(i,:),:)),1,1); % RP
    cmbname{i} = sprintf('%.0f : %.0f',comblist(i,:));
end
fpnt =mean(diff(unwrap(phdiff),1,2),2); % Average dRP over channel pairs
% subplot(2,1,1)
% plot(repmat(tvec(end-(20/dt)+1:end),size(comblist,1),1)',diff(phdiff,1,2)')
% legend(cmbname)
% 
% Kick it (do this individuallly!)


rlxtime = mean(rlxtime);

pertind = size(ypert,2)-floor(0.01/dt):size(ypert,2);
for n = 1:Nr
    nrcmbind = find(sum(comblist==1,2));
    ypert(n,pertind) = repmat(ypert(1,end) + pi,size(pertind));
    postind = size(ypert,2):size(ypert,2)+(10/dt);
    for i = postind
        dydt = dfx_kuramoto(ypert,A,K,Nr,D,omega,i);
        ypert(:,i+1) = wrapTo2Pi(ypert(:,i) +(dt*dydt') + (stvec(:,i).*dt));
    end
    ypert_post = ypert(:,postind);
    ypert_post = ypert(:,floor(0.1/dt):end);
    phdiff_post = [];
    for i = nrcmbind'
        phdiff_post(i,:) = diff(unwrap(ypert_post(comblist(i,:),:)));
    end
    rlx = mean(unwrap(diff(phdiff_post,1,2))-fpnt(nrcmbind));
%     subplot(2,1,2)
%     plot(unwrap(diff(phdiff_post,1,2))'); hold on; plot(rlx','k--','LineWidth',2); ylim([-0.25 0.25])
%     legend(cmbname)
    rlx = rlx(floor(0.1/dt):end);
    if abs(median(fpnt))<0.01
        rlxt = (min(find(abs(rlx)<0.01))*dt)+0.1;
        if isempty(rlxt)
            rlxt = -inf;
        end
    else
        rlxt = inf;
    end
    rlxtime(n) = rlxt;
end
