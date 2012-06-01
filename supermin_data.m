function [supermin_data_,varargout] = supermin_data(show_plot)
%SUPERMIN_DATA   Generate data for superminimization.

nposition = 100;
position = linspace(-1.5,1.5,nposition);

% True value of function to minimize
value = sin(2*pi*position + 3*pi/2) + position.^2 + .6*position;

% Value with noise
noise_amplitude = .5;
noisy_value = value + .25*sin(2*pi*position*4) ...
    + noise_amplitude*(rand(size(position)) - .5);

% Sampled values from values with noise
sampling_proportion = .75;
samples = round(sampling_proportion*nposition);
position_sample = randsample(position,samples);
position_sample = sort(position_sample);
noisy_value_sample = interp1(position,noisy_value,position_sample);

supermin_data_.position = position_sample;
supermin_data_.value = noisy_value_sample;

if show_plot
    figure
    a1 = axes('box','on','nextplot','add','xgrid','on','ygrid','on');
    varargout{1} = a1;
    plot(a1,position,value,'color','k')
    plot(a1,position,noisy_value)
    plot(a1,position_sample,noisy_value_sample,'linestyle','none',...
        'marker','o','markerfacecolor','b')
end
