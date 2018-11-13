
function gp = matsys(gp)


Am  = gp.a;
Bm  = gp.b;
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

