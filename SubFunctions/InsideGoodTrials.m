% eliminate peaks force in the trials that are not good and add information
% about the number of trial and the status of the trial

function [Analysis_post] = InsideGoodTrials (Analysis,cicles,good_trials,status)

names = fieldnames(Analysis);
good_trials = round(good_trials(:,1));
Analysis_post = Analysis;
Analysis_post.name_peaks_matrix = [Analysis_post.name_peaks_matrix,'num_cicle','task_phase'];

for n = 1 :size (names,1)
    if ~contains (names{n},'name') && ~isempty(Analysis.(names{n}))
        index = [];
        cicle_column = [];
        status_column = [];
        peaks = Analysis.(names{n});
        for i_peak = 1: size(peaks,1)
            onset = peaks(i_peak,2);
            cicle = cicles(onset);
            stat = status (onset);
            if sum(ismember(good_trials,cicle))>0.5
                index = [index , i_peak];
                cicle_column = [cicle_column; cicle];
                status_column = [status_column; stat];
            end
        end
        Analysis_post.(names{n})= [peaks(index,:),cicle_column,status_column];
    end
end

end
