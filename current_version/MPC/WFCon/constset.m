function gp = constset(k,ap,gp)

yub = [];
ylb = [];
uub = [];
ulb = [];

% output constraints
% Pm <= P <= PM
% Fm <= F <= FM
for i=1:gp.Na
    v   = ap.SwMes{i};
    cf  = ap.cf(i,k);
    cp  = ap.cp(i,k);
    
    um  = 0;
    if k>2
        um=.1;
    end
    
    uM  = ap.uM;
    if k>ap.sd && i==3
        um = 0;
        uM = 0;
    end
    
    Pm  = um.*v.^3.*cp;
    PM  = uM.*v.^3.*cp;
    
    Fm  = -uM.*v.^2.*cf;
    FM  = -um.*v.^2.*cf;
    
    yub = [yub; [FM PM uM]' ];
    ylb = [ylb; [Fm Pm um]' ];
    
    uub = [uub; uM];
    ulb = [ulb; um];
    
end

% constraints wind farm error signal
yub = [yub; Inf];  
ylb = [ylb; -Inf];

% constaints for the complete horizon
gp.yub = repmat(yub, gp.Nh,1);
gp.ylb = repmat(ylb, gp.Nh,1);

% control constraints for the complete horizon
gp.uub = repmat(uub, gp.Nh,1);
gp.ulb = repmat(ulb, gp.Nh,1);

end

