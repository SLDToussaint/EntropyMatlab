function Res=EntropyWithBootstrap(event_list, N_bootstraps)
% EntropyWithBootstrap  Compute entropy on generated bootstraped samples of events.
%   Res = EntropyWithBootstrap(event_list, N_bootstraps)
%   Arguments:
%   - event_list: raw list of events observed (with observed repetitions).
%     Can be an array or a cell array.
%   - N_bootstraps (optional): Number of boostraps to perform. Defaults to
%     1 (which is still resampled so not equal to the original event_list).
%     if N_bootstraps=0, then an entropy computation is done without
%     bootstrapping
%
%
%   Output:
%   - Res: Array of entropy values for each bootstrap

if nargin<2
    N_bootstraps=1;
end


[~, ~, sampling_vector] = unique(event_list);
N_events = length(sampling_vector);


% we generate the bootstraped lists of events and compute their entropy
if N_bootstraps>0
    Res = nan(N_bootstraps,1);
    for boot=1:N_bootstraps
        
        % resample with replacement to make a bootstrapped sample
        resampled_v = sampling_vector(randi(N_events,N_events,1));
        
        % compute frequency of events
        bootstraped_event_counts = accumarray(resampled_v,ones(length(resampled_v),1));
        freqs = bootstraped_event_counts(bootstraped_event_counts>0); % events not appearing are not counted
        freqs = freqs/sum(freqs);
        
        % compute entropy
        Res(boot)=-sum(freqs.*log2(freqs));
    end
    
else
    % simply compute entropy without bootstrap
    event_counts = accumarray(sampling_vector,ones(length(sampling_vector),1));
    freqs=event_counts(event_counts>0); % events not appearing are nto counted
    freqs = freqs/sum(freqs);
    
    % compute entropy
    Res = -sum(freqs.*log2(freqs));
end

end
