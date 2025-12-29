recObj = audiorecorder(22050, 8, 1); % 22.05kHz, 8 bits per sample, 1 channel
recObj.StartFcn = 'disp(''Start speaking.'')';
recObj.StopFcn = 'disp(''End of recording.'')';
recordblocking(recObj, 5); % Record for 5 seconds (do not change)
play(recObj); % Play recorded audio
audioRaw = getaudiodata(recObj); % Extract audio to vector
save('u2227951_lab_Audio.mat', 'audioRaw'); 