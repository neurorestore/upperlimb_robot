% eliminate peaks speed in the trials that are not good and add information
% about the number of trial and the status of the trial

function [Analysis_post] = InsideGoodTrialsSpeed (Analysis,cycles,good_trials,status,fSrobot,fSKIN)

good_trials = round(good_trials(:,1));
Analysis_post = Analysis;
Analysis_post.name_pks = [Analysis_post.name_pks,'num_cicle','task_phase'];

for n_point = 1:size(Analysis.pks,1)
    index = [];
    cycle_column = [];
    status_column = [];
    peaks = Analysis.pks{n_point,2}; % index of the peaks
    for i_peak = 1: size(peaks,2)
        onset = round(peaks(i_peak)*fSrobot);
        cycle = cycles(onset);
        stat = status (onset);
        if sum(ismember(good_trials,cycle))>0.5
            index = [index , i_peak];
            cycle_column = [cycle_column; cycle];
            status_column = [status_column; stat];
        end
    end
    Analysis_post.pks{n_point,1}= Analysis.pks{n_point,1}(index);
    Analysis_post.pks{n_point,2}= Analysis.pks{n_point,2}(index);
    Analysis_post.pks{n_point,3}= Analysis.pks{n_point,3}(index);
    Analysis_post.pks{n_point,4}= cycle_column';
    Analysis_post.pks{n_point,5}= status_column';
end

end