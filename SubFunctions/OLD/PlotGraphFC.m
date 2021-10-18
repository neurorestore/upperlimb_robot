
function PlotGraphFC (Matrix,cell_type,info_days,Mforce,days_force,MENO,PIU,fScalcium,fSrobot)
ind_ty = find(cell_type == 1); % excitatory
t_tot = MENO+PIU;
st = MENO/2;
X = [st:1/fScalcium:t_tot];
Xf = [st:1/fSrobot:t_tot];
for i_c = 1:length(ind_ty)% cycle for all the units
    M = Matrix{ind_ty(i_c)};
    % verify in which day the cell has been recorded and consider only
    % those peaks of force
    day = unique(info_days{ind_ty(i_c)});
    in_force = find(ismember(days_force,day));
    if size(M,1)==length(in_force)
        Area = NaN (size(M,1),2); % 1 col -> area force peak, 2 col -> area calcium peak
        MM = -Mforce(in_force,:);
        for i_ind = 1:size (M,1) %cycle for all the index
            % move the interval for the calcium peaks around the real peak
            % of the calcium activity
            [~,mas] = max(M(i_ind,:));
            mas = mas-MENO*fScalcium;
            stc = mas-MENO/2*fScalcium;
            enc = mas+PIU*fScalcium;
            if enc > size(M,2)
                enc = size(M,2);
                stc = enc-t_tot*fScalcium-1;
            elseif enc < length(X)
                stc = 1;
                enc = length(X);
            end
            Area(i_ind,2) = trapz(X,M(i_ind,stc:enc));
            %bring the force peak at zero level, to calculate only the
            %integration of the area inside the force peak
            line = MM(i_ind,:);%-min(MM(i_ind,:));
            Area(i_ind,1) = trapz(Xf,line(st*fSrobot:t_tot*fSrobot));
        end
    else
        disp('Error - not same number of index')
    end
    figure
    scatter (Area(:,1),Area(:,2))
    xlabel('Integration of force peaks')
    ylabel('Integration of calcium peaks')
end

end