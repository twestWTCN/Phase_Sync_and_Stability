function [segL_ddt segRP] = SRP_Lengths(RP,dRP,SRPeps,fsamp,period)
segRP = NaN;
qstable = find(abs(dRP')<SRPeps);
consecSegs = SplitVec(qstable','consecutive');
% consecSegs = fixGap(consecSegs,minsep);
% lengths
segL_ddt = cellfun('length',consecSegs);
seglist = find(segL_ddt>(period));

seg_ddt = [consecSegs{segL_ddt>(period)}];
segL_ddt = segL_ddt(segL_ddt>(period))/fsamp;

for i = 1:numel(seglist)
    if segL_ddt(i)<(400) %%% || std(diff(dRP))>0.1
        segRP(i) = circ_mean(RP(consecSegs{seglist(i)})');
    end
end

if isempty(segL_ddt)
    segL_ddt = NaN; segRP = NaN;
end
