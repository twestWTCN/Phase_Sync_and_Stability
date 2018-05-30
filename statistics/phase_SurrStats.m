function [SRPeps] = phase_SurrStats(phidata,NR,SRPeps_prctile)


    dphi_12_dt = wrapToPi(cumsum(phidata(1,:)-phidata(2,:)));

for N = 1:NR
    % Phase Shuffling
    Phishuff =[];
    Phishuff(1,:) = phidata(1,randperm(length(phidata(1,:))));
    Phishuff(2,:) = phidata(2,randperm(length(phidata(2,:))));
    dphi_12_dt = diff(unwrap(Phishuff(1,:)-Phishuff(2,:)));
    SRPbank(:,N) = mean(abs(dphi_12_dt));
end
% ppm.delete();
SRPeps = prctile(SRPbank(:),SRPeps_prctile);
