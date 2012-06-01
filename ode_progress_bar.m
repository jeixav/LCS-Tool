function status = ode_progress_bar(t,~,flag,progressBar,initialTime)

if nargin < 3 || isempty(flag)
    progressBar.setValue(abs(t(end)-initialTime))
else
    switch(flag)
        case 'init'
            progressBar.setElapsedTimeVisible(1);
            progressBar.setRemainedTimeVisible(1);
            progressBar.start
        case 'done'
            progressBar.stop
            fprintf('\n')
    end
end

status = 0;
