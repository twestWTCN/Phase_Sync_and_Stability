function coupleNMM_Generic(R,varnamez,varlist,dflag,Tag)
rng(R.seed)
rN = R.rN;
Klist = R.Klist;
% set defaults
omvar = 1;
sigvar =0;
dvar = 0.05;

for L = 1:3
    switch varnamez
        case 'omvar'
            omvar = varlist(L);
        case 'sigvar'
            sigvar = varlist(L);
        case 'dvar'
            dvar = varlist(L);
    end
        [d{1} d{2} d{3} d{4} d{5} d{6} d{7} SRPeps] = coupledNMM_wrapper_gen(R,Klist,dvar,omvar,sigvar,dflag,1);
    R.SRPeps = 0.005;
    parfor rp = 1:rN
        [PLV(:,:,rp) dRPvar(:,:,rp) MsKappa(:,:,rp) LHat(:,:,rp) LVar(:,:,rp) RPvar(:,:,rp) rlxtime(:,:,rp)] = coupledKuramoto_wrapper_gen(R,Klist,dvar,omvar,sigvar,dflag)
        disp(rp)
    end
    PLV_om{L} = PLV;
    dRPvar_om{L} = dRPvar;
    MsKappa_om{L} = MsKappa;
    LHat_om{L} = LHat;
    LVar_om{L} = LVar;
    RPvar_om{L} = RPvar;
    rlxtime_om{L} = rlxtime;
end
save([R.path '\Results\Kuramoto\' Tag '_simstats'],'PLV_om','dRPvar_om','MsKappa_om','LHat_om','LVar_om','RPvar_om','rlxtime_om')

% subplot(3,1,3); plot(ystore(1,:),ystore(2,:)); hold on
% plot(ystore(1,end-floor(500/dt):end),ystore(2,end-floor(500/dt):end),'r');