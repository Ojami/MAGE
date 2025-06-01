function randStr = getRandomName
% Generates random string

charSet = ['a':'z',upper('a':'z'),'0':'9'];
randStr = string(charSet(randi(numel(charSet),1, 10)));