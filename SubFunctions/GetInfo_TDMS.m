%% Function to obatin information by the _info.tdms files %%

function info = GetInfo_TDMS(Info)

info = struct();

if Info.Info_task.Forelimb_sx.data
    info.Paw = 'left';
else
    info.Paw = 'right';
end

if Info.Info_task.Stand_up.data
    info.Pos = 'bipedal';
else
    info.Pos = 'all4';
end

if Info.Info_task.Status.data ==0
    info.Status = 'active';
elseif Info.Info_task.Status.data ==1
    info.Status = 'Push_active';
elseif Info.Info_task.Status.data ==2
    info.Status = 'Pull_active';
elseif Info.Info_task.Status.data ==3
    info.Status = 'Passive';
else
    disp ('error Status task')
end

end