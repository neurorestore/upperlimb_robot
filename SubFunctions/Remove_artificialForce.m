% Function to remove the artificial force that the motors that move the
% sled create during the movement

function [ff] = Remove_artificialForce(force,status,good_trials,fS)

start1 = round(good_trials(:,2)*fS);
end2 = round(good_trials(:,3)*fS);
%find the start and the end of each trial
% index1 = find (status ==1);
% index2 = find (status ==2);
% 
% st1 = find(diff(index1)>1);
% st2 = find(diff(index2)>1);
% 
% start1 = [index1(1),index1(st1+1)];
% end2 = [index2(st2),index2(end)];
% 
% if length(start1) ~= length(end2)
%     start1(end) =[];
% elseif index2(end)==length(force)
%     start1(end) = [];
%     end2(end) = [];
% end

% resample to have the same number of frames for every cicle
% (in theory all cicles should be exactly the same but sometimes it can go
% slower for various reasons)
cicles = end2-start1+1;
max_length = max(cicles);

matrix_force = NaN (length (start1), max_length);
for i_c = 1:length(start1)
    f_points = resample(force(start1(i_c):end2(i_c)),max_length,cicles(i_c));
    matrix_force (i_c,:) = f_points;
end

valueToSub = sgolayfilt(median(matrix_force),3,21);
figure ; plot(valueToSub)

ff = force;
for i_c = 1:length(start1)
    valueToSubSampled = resample(valueToSub,cicles(i_c),max_length);
    ff(start1(i_c):end2(i_c)) = force(start1(i_c):end2(i_c))-valueToSubSampled;
end

end
    