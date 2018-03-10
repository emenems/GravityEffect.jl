function shp = shp2txt(filein,fileout,varargin)
%SHP2TXT convert shapefile to txt binary file 
% Use this function to convert shapefiles to either 
% ascii or Matlab's v7 binary file that can be loaded 
% to Julia using MAT package. 
% When using Octave, IO and MAPPING packages are required.
% When using Matlab, Mapping toolbox is required
%
% Input:
%   filein  ... input file name (full)
%   fileout ... output file name (full). Set to []
%               reading only.
%   varargin{1} : optional output precision (by default '%.4f');
% 
% Output:
%   shp     ... loaded/witten shapefile
%
% Example:
%   shp2txt('exclusion_polygon.shp','exclusion_polygon.txt','%.2f');
%

%% Code
    % Check if matlab or octave is used
    v = version;
    if ~strcmp(v(end),')')
        pkg load io;
        pkg load mapping;
    end
    % load all 
    shp = shaperead(filein);
    % write
    if ~isempty(fileout)
		switch fileout(end-2:end)
			case 'mat'
				if ~strcmp(v(end),')')
					save(fileout,'shp','-mat7-binary');
				else
					save(fileout,'shp');
				end
			case 'txt'                             
                if nargin > 2 
                    outprec = varargin{1}
                else
                    outprec = '%.4f';
                end
                outprec = [outprec,' ',outprec,'\n'];
                fid = fopen(fileout,'w');
                fprintf(fid,'%%X Y\n');
                for i = 1:length(shp)
                    tempx = shp(i).X;
                    tempy = shp(i).Y;
                    for j = 1:length(tempx)
                        fprintf(fid,outprec,tempx(j),tempy(j));
                    end
                    if ~isnan(tempx(end)) || ~isnan(tempy(end))
                        fprintf(fid,outprec,NaN,NaN);
                    end
                end
                fclose(fid);
		end 
    end
end