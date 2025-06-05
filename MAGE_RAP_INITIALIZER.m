function MAGE_RAP_INITIALIZER
% a helper function for install MATLAB dependencies inside docker container
% within a ttyd instance

system("wget https://mathworks.com/mpm/glnxa64/mpm");
system("chmod +x mpm");
system("sudo chmod +x /usr/local/MATLAB/R2024b");

% by default on ttyd: /home/matlab/Documents/MATLAB/MAGE/mpm_input_r2024b.txt
system("./mpm install --inputfile=/home/matlab/Documents/MATLAB/MAGE/mpm_input_r2024b.txt");

end % END