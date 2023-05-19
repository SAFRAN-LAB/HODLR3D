function [Tx,Ty,Tz] = gen_points(cx,cy,cz,L,kpx,kpy,kpz)

kxpanels = kpx;
kypanels = kpy;
kzpanels = kpz;
ax = cx - L;
bx = cx + L;
ay = cy - L;
by = cy + L;
az = cz - L;
bz = cz + L;
xh = linspace(ax,bx,kxpanels+2);
yh = linspace(ay,by,kypanels+2);
zh = linspace(az,bz,kzpanels+2);
xh(end) =[]; xh(1) =[];
yh(end) =[]; yh(1) =[];
zh(end) =[]; zh(1) =[];
Tx = zeros(kxpanels*kxpanels,1);
Ty = zeros(kypanels*kypanels,1);
Tz = zeros(kzpanels*kzpanels,1);

lk=1;

for i=1:length(xh)
        for j =1:length(yh)
            for k =1:length(zh)
                Tx(lk) = xh(i);
                Ty(lk) = yh(j);
                Tz(lk) = zh(k);
                lk=lk+1;
            end
        end
end