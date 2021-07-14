%------------------------------------------------
%defines the wedge fucntion used for holograms
%------------------------------------------------

function [ phi ] = wedge(x,y,angx,angy,apx,apy,lam)

phi=exp(1i*(x+apx)*(((tand(angx)*2*apx)/lam)/apx)).*exp(1i*(y+apy)*(((tand(angy)*2*apy)/lam)/apy));

end

