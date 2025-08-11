function [t2,y_ip,x_peaks,pidx,y_peaks,nidx]=slope_calc(t2,y_ip)

% Finding all the peaks % 
[x_peaks,pidx]=findpeaks(y_ip,t2);
[y_peaks,nidx]=findpeaks(-y_ip,t2);

z = 1;
% Slope calculation for interpolated data % (using Positive peaks)
for i = 1 : length(x_peaks)
    idx = pidx(i);
    if (find(t2==idx) > 2) && (find(t2==idx) < length(t2)-1)
        l1 = (find(t2 == idx)-1);
        l2 = (find(t2 == idx)-2);
        r1 = (find(t2 == idx)+1);
        r2 = (find(t2 == idx)+2);
        ml = (y_ip(l1)-y_ip(l2))/(t2(l1)-t2(l2));
        mr = (y_ip(r1)-y_ip(r2))/(t2(r1)-t2(r2));
        A = [-ml 1; -mr 1];
        B = [(y_ip(l1)) - (ml * t2(l1)) ; (y_ip(r1)) - (mr * t2(r1)) ];
        pos_points = A\B;
        px(z) = pos_points(1);
        py(z) = pos_points(2);
        tempp(z) = idx;
        z = z + 1;
    end
end

%%% new positive peaks replacement %%%
for i = 1 : length(tempp)
    idx = tempp(i);
    j = find(t2 == idx);
    t2(j) = px(i);
    y_ip(j) = py(i);
end

q = 1;
% Slope calculation for interpolated data % (using Negative peaks)
for i = 1 : length(y_peaks)
    idx = nidx(i);
    if (find(t2==idx) > 2) && (find(t2==idx) < length(t2)-1)
        l1 = (find(t2 == idx)-1);
        l2 = (find(t2 == idx)-2);
        r1 = (find(t2 == idx)+1);
        r2 = (find(t2 == idx)+2);
        ml = (y_ip(l1)-y_ip(l2))/(t2(l1)-t2(l2));
        mr = (y_ip(r1)-y_ip(r2))/(t2(r1)-t2(r2));
        A = [-ml 1; -mr 1];
        B = [(y_ip(l1)) - (ml * t2(l1)) ; (y_ip(r1)) - (mr * t2(r1)) ];
        pos_points = A\B;
        nx(q) = pos_points(1);
        ny(q) = pos_points(2);
        tempn(q) = idx;
        q = q + 1;
    end
end


%%% new negative peaks replacement %%%
for i = 1 : length(tempn)
    idx = tempn(i);
    j = find(t2 == idx);
    t2(j) = nx(i);
    y_ip(j) = ny(i);
end

