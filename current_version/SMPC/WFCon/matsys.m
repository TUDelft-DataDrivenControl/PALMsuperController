
function gp = matsys(gp)


Am  = gp.a;
Bm  = gp.b;
Bmt = gp.bt;
Bmr = gp.br;
Cm  = gp.c;


Nh = gp.Nh;

A = Am;
Ai = A;
AA = Ai;
for ii = 2:Nh
    Ai = A*Ai;
    AA = [AA;Ai];
end
gp.A = AA;


AiB = Bm;
BB = kron(eye(Nh),AiB);
for ii = 1:Nh-1
    AiB = A*AiB;
    BB = BB+kron(diag(ones(Nh-ii,1),-ii),AiB);
end
gp.B = BB;

AiBt = Bmt;
q    = num2cell(AiBt,[1,2]);
BBt  = blkdiag(q{:});
for ii = 1:Nh-1
    AA = A;
    for jj = ii:Nh-1
       BBt(jj*nx+1:(jj+1)*nx , (ii-1)*nu+1:ii*nu) = AA*AiBt(:,:,ii);
       AA = AA*A;
    end   
end
gp.Bt = BBt;

AiB = Bmr;
BB = kron(eye(Nh),AiB);
for ii = 1:Nh-1
    AiB = A*AiB;
    BB = BB+kron(diag(ones(Nh-ii,1),-ii),AiB);
end
gp.Br = BB;


C = Cm;
Cc = kron(eye(Nh),C);
gp.C = Cc;

end

