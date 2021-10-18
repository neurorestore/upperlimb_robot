%% Code to fix the trajectory and to reconstruct the step signal of the moving slide
% when the likelihood is too low

function pos_fix = Reconstruct_Trajectory(likelihood, pos,fS)

% find gaps that need to be reconstructed
index_low = find(likelihood<0.75);
index_out = find(isoutlier(pos)==1);
d_leng = length(pos);
% for outliers we can use fill missing, to have a flat reconstruction,
% because usually it is only one point
% pos_mis = pos;
% pos_mis(index_out)= NaN;
% pos_mis = fillmissing(pos_mis,'linear');
% pos(index_out) = pos_mis(index_out);

%
%index_low = [index_lik;index_out];
%index_low = sort(index_low);
differenza = diff(index_low);
n = 5; % number of points cut before st and after en 
ii_10 = find(differenza>=2*n); % join intervals too near
%ii = find(differenza>1.5);
if ~isempty(index_low) % it could be a file completly reconstructed
st = [index_low(1);index_low(ii_10+1)];
en = [index_low(ii_10);index_low(end)];

% considering also n points as NaN before and after likelihood intervals
st = st-3;
en = en+3;
% check border conditions
i_del = [];
for i_g = 1:length(st)
    if st(i_g)<25
        i_del = [i_del,i_g];
    end
    if en(i_g)>length(pos)-25
        i_del = [i_del,i_g];
    end
end
st(i_del) = [];
en(i_del) = [];

for i_g = 1:length(st)
    pos(st(i_g):en(i_g)) = NaN;
end
if sum(isnan(pos))~=0 % if there is only one part at the extremity of the recording with low likelihood we don't consider it
figure
hold on
plot(pos-nanmean(pos),'LineWidth',1) %plot around zero

% if the gap is in a part of flat signal, a fillgaps reconstruction is used
dev = nanstd (pos);

speed = derivative(pos,1/fS);
th_s = nanmean(speed(isoutlier(speed,'movmedian',100)==0))+nanstd(speed(isoutlier(speed,'movmedian',100)==0));

% find index start and end of intervals for the speed because they are
% bigger
index = find(isnan(speed)==1);
ii_s = find(diff(index)>1.5);
st_s = [index(1);index(ii_s+1)];
en_s = [index(ii_s);index(end)];

% scroll all the gaps
pos_fix = pos;
for i_g = 1:length(st)
    % find the nearest start and end for the speed
    [~,i_m] = min(abs(st(i_g)-st_s));
    
    dist = abs(pos(en(i_g)+1)-pos(st(i_g)-1));
    if dist<3  %dev/8 % signal enough flat in the gap
        int = fillgaps(pos,5);
        pos_fix(st(i_g)-1:en(i_g)+1) = int(st(i_g)-1:en(i_g)+1);
    else %the signal is missing in a step
        if abs(speed(en_s(i_m)+1))<th_s && abs(speed(st_s(i_m)-1))<th_s 
            % step in the middle between two flat parts
            p_st = polyfit([st_s(i_m)-5:st_s(i_m)-1]',[pos(st_s(i_m)-5:st_s(i_m)-1)],1);
            p_en = polyfit([en_s(i_m)+1:en_s(i_m)+5]',[pos(en_s(i_m)+1:en_s(i_m)+5)],1);
            x_mid = (en_s(i_m)-st_s(i_m))/2+st_s(i_m);
            y_mid = (pos(st_s(i_m))-pos(en_s(i_m)))/2+pos(en_s(i_m));
            b = y_mid-(-3)*x_mid;
            p_mid = [-3,b];
            %calculate intersection
            x1_intersect = round(fzero(@(x) polyval(p_st-p_mid,x),0));
            y1_intersect = round(polyval(p_st,x1_intersect));
            x2_intersect = round(fzero(@(x) polyval(p_mid-p_en,x),0));
            y2_intersect = round(polyval(p_en,x2_intersect));
            pos(x1_intersect) = y1_intersect;
            pos(x2_intersect) = y2_intersect;
            % divide filling gap from start to intersection1
            int = fillgaps(pos(1:x1_intersect),5);
            pos_fix(st(i_g)-1:x1_intersect) = int(st(i_g)-1:x1_intersect);
            % filling gap from intersection1 to intersection2
            int = fillgaps(pos);
            pos_fix(x1_intersect:x2_intersect) = int(x1_intersect:x2_intersect);
            % filling gap from intersection2 to end
            int = fillgaps(pos(x2_intersect:end),5);
            d = en(i_g)+2-x2_intersect;
            pos_fix(x2_intersect:en(i_g)+1) = int(1:d);
            
        elseif abs(speed(en_s(i_m)+1))>=th_s && abs(speed(st_s(i_m)-1))<th_s
            % step during a pulling movement
            p_st = polyfit([st_s(i_m)-5:st_s(i_m)-1]',[pos(st_s(i_m)-5:st_s(i_m)-1)],1);
            p_en = polyfit([en_s(i_m)+1:en_s(i_m)+5]',[pos(en_s(i_m)+1:en_s(i_m)+5)],1);
            %calculate intersection
            x_intersect = round(fzero(@(x) polyval(p_st-p_en,x),0));
            y_intersect = round(polyval(p_st,x_intersect));          
            pos(x_intersect) = y_intersect;
            % divide filling gap from start to intersection
            int = fillgaps(pos(1:x_intersect),5);
            pos_fix(st(i_g)-1:x_intersect) = int(st(i_g)-1:x_intersect);
            % filling gap from intersection to end
            int = fillgaps(pos);
            pos_fix(x_intersect:en(i_g)+1) = int(x_intersect:en(i_g)+1);
        elseif abs(speed(en_s(i_m)+1))<th_s && abs(speed(st_s(i_m)-1))>=th_s
            % step during a pulling movement, it is missing the start of
            % the movement
            p_st = polyfit([st_s(i_m)-5:st_s(i_m)-1]',[pos(st_s(i_m)-5:st_s(i_m)-1)],1);
            p_en = polyfit([en_s(i_m)+1:en_s(i_m)+5]',[pos(en_s(i_m)+1:en_s(i_m)+5)],1);
            %calculate intersection
            x_intersect = round(fzero(@(x) polyval(p_st-p_en,x),0));
            y_intersect = round(polyval(p_st,x_intersect));          
            pos(x_intersect) = y_intersect;
            % divide filling gap from start to intersection
            int = fillgaps(pos(1:x_intersect));
            pos_fix(st(i_g)-1:x_intersect) = int(st(i_g)-1:x_intersect);
            int = fillgaps(pos(x_intersect:end),5);
            d = en(i_g)+2-x_intersect;
            pos_fix(x_intersect:en(i_g)+1) = int(1:d);            
        elseif abs(speed(en_s(i_m)+1))>=th_s && abs(speed(st_s(i_m)-1))>=th_s
            if speed(en_s(i_m)+1)>= th_s && speed(st_s(i_m)-1)<-th_s
                % situation when it is missing a big part of the signal and DLC find back the signal when it is already started the pushing phase
                p_st = polyfit([st_s(i_m)-5:st_s(i_m)-1]',[pos(st_s(i_m)-5:st_s(i_m)-1)],1);
                p_en = [0,pos(en_s(i_m)+1)];
                %calculate intersection
                x_intersect = round(fzero(@(x) polyval(p_st-p_en,x),0));
                y_intersect = round(polyval(p_st,x_intersect));          
                pos(x_intersect) = y_intersect;
                % divide filling gap from start to intersection
                int = fillgaps(pos(1:x_intersect));
                pos_fix(st(i_g)-1:x_intersect) = int(st(i_g)-1:x_intersect);
                int = fillgaps(pos(x_intersect:end),5);
                d = en(i_g)+2-x_intersect;
                pos_fix(x_intersect:en(i_g)+1) = int(1:d);            
            else % both derivative are negative, it is missing only a little part of the signal along the same track
                int = fillgaps(pos);
                pos_fix(st(i_g)-1:en(i_g)+1) = int(st(i_g)-1:en(i_g)+1);
            end
        end        
    end
end
else 
    pos_fix = pos;
end
else
    pos_fix = pos;
end

pos_fix = sgolayfilt(pos_fix(1:d_leng),3,21);
plot(pos_fix-nanmean(pos_fix));
hold on
plot(derivative(pos_fix,1/fS)/3);

end