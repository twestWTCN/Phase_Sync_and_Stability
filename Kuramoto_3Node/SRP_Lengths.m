function [segL_ddt segRP] = SRP_Lengths(RP,dRP,SRPeps,fsamp,period)
qstable = find(abs(dRP')<SRPeps);
consecSegs = SplitVec(qstable','consecutive');
% consecSegs = fixGap(consecSegs,minsep);
% lengths
segL_ddt = cellfun('length',consecSegs);
seglist = find(segL_ddt>(period));

seg_ddt = [consecSegs{segL_ddt>(period)}];
segL_ddt = segL_ddt(segL_ddt>(period))/fsamp;

for i = 1:numel(consecSegs)
    segRP(i) = circ_mean(RP(consecSegs{i})');
end
