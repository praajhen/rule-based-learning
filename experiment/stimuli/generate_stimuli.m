% ============================================================
% Stimulus Generation Script for Rule-Based Learning Experiment
%
% Author: Praghajieeth Raajhen Santhana Gopalan
%% Institution: University of Jyväskylä, Finland
%
% Description:
% This script generates 3 shapes (circle, square, triangle)
% in 3 colors (red, green, blue) on a black background.
% Total: 9 stimuli used in rule-based learning experiment.
%
% Output:
% PNG images (400 x 400 pixels)
%
% Created: 2026
% ============================================================

clear
close all
clc

img_size = 400;
shape_size = 200;

% RGB colors
colors = [255 0 0;     % red
          0 180 0;     % green
          0 0 255];    % blue

color_names = {'red','green','blue'};
shapes = {'circle','square','triangle'};

% output folder (GitHub-friendly)
outdir = 'stimuli';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

for s = 1:3
    for c = 1:3
        
        % black background
        img = zeros(img_size,img_size,3,'uint8');
        
        center = img_size/2;
        half = shape_size/2;
        
        if s == 1   % circle
            
            [x,y] = meshgrid(1:img_size,1:img_size);
            mask = (x-center).^2 + (y-center).^2 <= half^2;
            
        elseif s == 2   % square
            
            mask = false(img_size);
            mask(center-half:center+half,center-half:center+half) = true;
            
        else   % triangle
            
            mask = poly2mask([center center-half center+half],...
                             [center-half center+half center+half],...
                             img_size,img_size);
        end
        
        % apply color
        for ch = 1:3
            temp = img(:,:,ch);
            temp(mask) = colors(c,ch);
            img(:,:,ch) = temp;
        end
        
        % save image
        filename = fullfile(outdir, ...
            sprintf('%s_%s.png',shapes{s},color_names{c}));
        
        imwrite(img,filename);
        
    end
end

disp('Stimuli generated successfully');