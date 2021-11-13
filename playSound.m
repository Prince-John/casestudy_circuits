%% function playSound(y, Fs)
% scales the audio signal in y to range from -1.0 to 1.0 (if y is too
% large), then plays the entire sequence. Control is returned once the
% sequence is finished playing.
%
% y - audio time sequence 
% Fs - sampling rate for y (Hz). The sampling rate Fs is related to the
% sampling interval by the relation Fs = 1/h.
%
% Matthew Lew 10/25/2018
% updated 11/27/2020

function playSound(y, Fs)
if max(abs(y))>1
    y=y/max(abs(y));
end
obj = audioplayer(y, Fs);
playblocking(obj);
end