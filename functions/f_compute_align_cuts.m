function [params] = f_compute_align_cuts(params)
if ~isfield(params, 'align_pulse_crop_method'); params.align_pulse_crop_method = 0; end % default
if ~isfield(params, 'align_pulse_buff'); params.align_pulse_buff = 10; end % in frames
if ~isfield(params, 'align_pulse_min_frag_size'); params.align_pulse_min_frag_size = 10; end % in frames
if ~isfield(params, 'align_pulse_thresh'); params.align_pulse_thresh = 0.5; end % from normalized ca ave trace
if ~isfield(params, 'align_pulse_plot_out'); params.align_pulse_plot_out = true; end

mean_trace = mean(params.ave_trace);
max_trace = max(params.ave_trace);
norm_ave_trace = (params.ave_trace - mean_trace)/(max_trace-mean_trace);

T = numel(norm_ave_trace);
vid_cuts_trace = false(T,1);

if params.align_pulse_crop_method
    if params.align_pulse_crop_method == 1 % auto find 
        
        pulse_trace = norm_ave_trace > params.align_pulse_thresh;
        
        d_pulse_trace = diff(pulse_trace);
    
        pulse_on = find(d_pulse_trace>0)+1;
        pulse_off = find(d_pulse_trace<0)+1;
    
        % quality check if light turns off in beginning 
        if and(pulse_trace(1) == 1, pulse_on(1) > pulse_off(1))
            pulse_on = [1; pulse_on];  % removing first off ?
        end
        
        pulse_on_buff = max(pulse_on - params.align_pulse_buff, 1);
        pulse_off_buff = min(pulse_off + params.align_pulse_buff, T);
        
        num_frag = round(numel(pulse_on)+1);
        vid_cuts = ones(num_frag,2);
        vid_cuts(end,2) = T;
    
        for n_frag = 1:(num_frag-1)
            vid_cuts(n_frag,2) = pulse_on_buff(n_frag);
            vid_cuts(n_frag+1,1) = pulse_off_buff(n_frag);
        end
        
        throw_idx = diff(vid_cuts,1,2) < params.align_pulse_min_frag_size;
        vid_cuts(throw_idx,:) = [];
        num_frag = size(vid_cuts,1);
    
        for n_frag = 1:num_frag
            vid_cuts_trace(vid_cuts(n_frag,1):vid_cuts(n_frag,2)) = 1;
        end
        
        title_tag = 'Automatic';

    elseif params.align_pulse_crop_method == 2 % manual
        f1 = figure;
        plot(norm_ave_trace);
        axis tight;
        title('how many fragments?');
        num_frag = input('how many fragments? (int):');
        vid_cuts = zeros(num_frag,2);
        for n_frag = 1:num_frag
            title(sprintf('Select fragment %d/%d (2 clicks)', n_frag,num_frag));
            [temp_cuts, ~] = ginput(2);
            vid_cuts(n_frag,:) = round(temp_cuts);
            if vid_cuts(n_frag,1) < 1
                vid_cuts(n_frag,1) = 1;
            end
            if vid_cuts(n_frag,2) > numel(ave_trace)
                vid_cuts(n_frag,2) = numel(ave_trace);
            end
            vid_cuts_trace(vid_cuts(n_frag,1):vid_cuts(n_frag,2)) = 1;
            plot(norm_ave_trace);
            hold on;
            plot(vid_cuts_trace);
            hold off;
            axis tight;
        end
        close(f1)
        
        title_tag = 'Manual';
        
    end
    
    if params.align_pulse_plot_out
        figure;plot(norm_ave_trace);
        hold on; plot(vid_cuts_trace);
        axis tight;
        title(sprintf('%s %d frag selected %s', title_tag, num_frag, params.title_tag), 'interpreter', 'none');
    end
else    % no pulses
    vid_cuts = [1, T];
    vid_cuts_trace(vid_cuts(1):vid_cuts(2)) = 1;
end

params.vid_cuts = vid_cuts;
params.vid_cuts_trace = vid_cuts_trace;
end
