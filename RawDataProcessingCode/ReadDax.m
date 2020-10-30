function [movie, infoFile, infoFileRoi] = ReadDax(varargin)
%--------------------------------------------------------------------------
% [movie, infoFiles] = ReadDax(fileName, varargin)
% This function loads a STORM movies from the dax file associated with the
% provided .inf file
%--------------------------------------------------------------------------
% Outputs:
% movies/LxMxN array: A 3D array containing the specified movie
% infoFile: infoFile structure for the specified daxfile
% infoFileRoi: modified infoFile corresponding to the daxfile
%--------------------------------------------------------------------------
% Input:
% fileName/string: a path to the dax or inf file 
%
%--------------------------------------------------------------------------
% Variable Inputs:
%
% 'file'/string ([]): A path to the associated .inf file
%
% 'path'/string ([]): Default path to look for .inf files
%
% 'startFrame'/double  (1): first of movie to load.  
%
% 'endFrame'/double ([]): last frame of the movie to load.  If empty will
% be max.   
%
% 'infoFile'/info file structure ([]): An info file for
% the files to be loaded.  
%
% 'imageDimensions'/3x1 integer array ([]): The size of the movies to be
% loaded.  
%
% 'verbose'/boolean (true): Display or hide function progress
%
% 'orientation'/string ('normal'): Control the relative rotation of the data
%   structure
%--------------------------------------------------------------------------
% Siyuan (Steven) Wang
% sywang1984@gmail.com
% October 30, 2020
%
% Version 2
%-------------------Updates:
% 10/30/2020: SW
% Combined the ReadInfoFile function and the other auxiliary functions into
% this function, so that the ReadDax function can be run alone without
% setting path to the matlab-storm folder and the co-existence of the
% ReadInfoFile.m
%
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% September 7, 2012
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
orientationValues = {'normal','nd2'};
flags = {'file', 'infoFile', 'startFrame','endFrame', 'verbose', ...
    'orientation', 'path'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
dataPath = defaultDataPath;
startFrame = 1;
endFrame = []; 
fileName = [];
infoFile = [];
verbose = true;
orientation = 'normal';

maxMemory = 10E9; % 1 Gb


%--------------------------------------------------------------------------
% Parse Required Input
%--------------------------------------------------------------------------
if nargin >= 1
    if ~ismember(varargin{1}, flags)
        fileName =varargin{1};
        varargin = varargin(2:end);
    end
end

%--------------------------------------------------------------------------
% Parse Variable Input
%--------------------------------------------------------------------------
if (mod(length(varargin), 2) ~= 0 )
    error(['Parameter names and parameters not in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch parameterName
        case 'file'
            fileName = parameterValue;
        case 'startFrame'
            startFrame = parameterValue;
        case 'endFrame'
            endFrame = parameterValue;
        case 'infoFile'
            infoFile = parameterValue;
        case 'maxMemory'
            maxMemory = parameterValue;
        case 'verbose'
            verbose = parameterValue;
        case 'orientation'
            orientation = parameterValue;
        case 'path'
            dataPath = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end

%--------------------------------------------------------------------------
% Check parameter consistency
%--------------------------------------------------------------------------
if ~isempty(infoFile) && ~isempty(fileName)
    error('You cannot specify both info file and file name.');
end

%--------------------------------------------------------------------------
% Load info files if needed
%--------------------------------------------------------------------------
if isempty(infoFile)
    if isempty(fileName)
        infoFile = ReadInfoFile('path', dataPath, 'verbose', verbose);
    else
        infoFile = ReadInfoFile(fileName, 'verbose', verbose);
    end
    
    if isempty(infoFile)
        display('Canceled');
        movie = {};
        return;
    end 
end

%--------------------------------------------------------------------------
% Load Dax Files
%--------------------------------------------------------------------------

% --------- Determine number of frames to load
framesInDax = infoFile.number_of_frames;

% parse startFrame and endFrame
if isempty(endFrame)
    endFrame = framesInDax;
end
if endFrame > framesInDax;
    if verbose
        warning('input endFrame greater than total frames in dax_file.  Using all available frames after startFrame');
    end
    endFrame = framesInDax; 
end
framesToLoad = endFrame - startFrame + 1;
frameDim = [infoFile.frame_dimensions(1)/infoFile.binning(1),...
            infoFile.frame_dimensions(2)/infoFile.binning(2)];
frameSize = frameDim(1)*frameDim(2);

% Check memory requirements.  Ask for confirmation if > maxMemory.  
memoryRequired = frameSize*framesToLoad*16/8;
DoThis = 1; 
if memoryRequired > maxMemory
    warning([fileName,'  ', 'is > ',num2str(maxMemory/1E6,3),' Mbs.']);
    DoThis = input(['Requested file requires ',...  
        num2str(memoryRequired/1E6,3),' Mbs. ',...
        'Are you sure you want to load it? ',... 
        '(Filling memory may crash the computer) ',...
        '0 = abort, 1 = continue, n = new end frame  ']);
    if DoThis > 1 
        endFrame = DoThis; % read in the new last frame
        framesToLoad = endFrame - startFrame + 1;
        DoThis = true;
    end
end
  
%--------------------------------------------------------
% Proceed to load specified poriton of daxfile
%--------------------------------------------------------
if DoThis
    fileName = [infoFile.localName(1:(end-4)) '.dax'];
    if verbose
        display(['Loading ' infoFile.localPath fileName ]);
    end

    if ~isempty( strfind(infoFile.data_type,'little endian') );
        binaryFormat = 'l';
    else
        binaryFormat = 'b';
    end
    
    %----------------------------------------------------------------- 
    % Read all pixels from selected frames
    %----------------------------------------------------------------- 
    fid = fopen([infoFile.localPath fileName]);
    if fid < 0
        error(['Invalid file: ' infoFile.localPath fileName]);
    end
    
    fseek(fid,(frameSize*(startFrame - 1))*16/8,'bof'); % bits/(bytes per bit)
    dataSize = frameSize*framesToLoad;
    movie = fread(fid, dataSize, '*uint16', binaryFormat);
    fclose(fid);
    
    try % Catch corrupt files
        if framesToLoad == 1
            movie = reshape(movie, frameDim)';
        else
            switch orientation % Change orientation
                case 'normal'
                    movie = permute(reshape(movie, [frameDim framesToLoad]), [2 1 3]);
                case 'nd2'
                    movie = permute(reshape(movie, [fliplr(frameDim) framesToLoad]), [2 1 3]);
                otherwise
            end
        end
    catch
        display('Serious error somewhere here...check file for corruption');
        movie = zeros(frameDim);
    end
    
    if verbose
        display(['Loaded ' infoFile.localPath fileName ]);
        display([num2str(framesToLoad) ' ' num2str(frameDim(1)) ' x ' num2str(frameDim(2)) ...
            ' frames loaded']);
    end
else
    error('User aborted load dax due to memory considerations '); 
end

function infoFile = ReadInfoFile(varargin)
%--------------------------------------------------------------------------
% infoFile = ReadInfoFile(fileName, varargin)
% This function returns a structure, info, containing the elements of an
% .inf file.
%--------------------------------------------------------------------------
% Outputs: 
% info/struct: A structure array containing the elements of the info file
%
%--------------------------------------------------------------------------
% Inputs:
% fileName/string or cell array of strings ([]): A path to a .dax or .ini
%   file
% Field: explanation
%                  localName: inf filename matlab found / should save as 
%                  localPath: where matlab found / should save this file
%                   uniqueID: ?
%                       file: full path to daxfile. 
%               machine_name: e.g. 'storm2'
%            parameters_file: Full pathname of pars file used in Hal
%              shutters_file: e.g. 'shutters_default.xml'
%                   CCD_mode: e.g. 'frame-transfer'
%                  data_type: '16 bit integers (binary, big endian)'
%           frame_dimensions: [256 256]
%                    binning: [1 1]
%                 frame_size: 262144
%     horizontal_shift_speed: 10
%       vertical_shift_speed: 3.3000
%                 EMCCD_Gain: 20
%                Preamp_Gain: 5
%              Exposure_Time: 0.1000
%          Frames_Per_Second: 9.8280
%         camera_temperature: -70
%           number_of_frames: 10
%                camera_head: 'DU897_BV'
%                     hstart: 1
%                       hend: 256  
%                     vstart: 1
%                       vend: 256
%                  ADChannel: 0
%                    Stage_X: 0
%                    Stage_Y: 5
%                    Stage_Z: 0
%                Lock_Target: 0
%                   scalemax: 4038
%                   scalemin: 0
%                      notes: ''
%--------------------------------------------------------------------------
% Variable Inputs:
% 'file'/string or cell array: The file name(s) for the .ini file(s) to load
%   Path must be included. 
%
% 'verbose'/boolean(true): Determines if the function hides progress
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% September 5, 2012
% jeffmoffitt@gmail.com
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------
quiet = 1;
flags = {'file', 'verbose', 'path'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Define default parameters
%--------------------------------------------------------------------------
infFileName = [];
dataPath = defaultDataPath;
verbose = false;

%--------------------------------------------------------------------------
% Parse Required Input
%--------------------------------------------------------------------------
if nargin >= 1
    if ~ismember(varargin{1}, flags)
        infFileName =varargin{1};
        varargin = varargin(2:end);
    end
end

%--------------------------------------------------------------------------
% Parse Variable Input Arguments
%--------------------------------------------------------------------------
if (mod(length(varargin), 2) ~= 0 )
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch parameterName
        case 'file'
            infFileName = parameterValue; 
        case 'path'
            dataPath = parameterValue; 
        case 'verbose'
            verbose = parameterValue; 
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end

%--------------------------------------------------------------------------
% Get file if needed
%--------------------------------------------------------------------------
if isempty(infFileName)
    [infFileName infFilePath] = uigetfile([dataPath '*.inf']);
    if isempty(infFileName)
        display('Loading canceled');
        infoFile = [];
        return;
    end
    if infFilePath(end) ~= '\'
        infFilePath = [infFilePath '\']; % All paths must end in '\'
    end
    infFileName = [infFilePath infFileName];
end

%--------------------------------------------------------------------------
% Open Inf File
%--------------------------------------------------------------------------
if strcmp(infFileName((end-3):end), '.dax')
    infFileName = [infFileName(1:(end-4)) '.inf'];
end

% Open file
fid = fopen(infFileName);
if fid == -1
    error([infFileName ' is not a valid .inf file']);
end

%--------------------------------------------------------------------------
% Read Inf File
%--------------------------------------------------------------------------
count = 1;
text = {};
while ~feof(fid)
    text{count} = fgetl(fid);
    count = count + 1;
end
fclose(fid);

%--------------------------------------------------------------------------
% Create Info File
%--------------------------------------------------------------------------
infoFile.localName = '';
infoFile.localPath = '';
infoFile.uniqueID = now;
infoFile.file = '';
infoFile.machine_name = '';
infoFile.parameters_file = '';
infoFile.shutters_file = '';
infoFile.CCD_mode = '';
infoFile.data_type = '16 bit integers (binary, big endian)';
infoFile.frame_dimensions = [0 0];
infoFile.binning = [1 1];
infoFile.frame_size = 0;
infoFile.horizontal_shift_speed = 0;
infoFile.vertical_shift_speed = 0;
infoFile.EMCCD_Gain = 1;
infoFile.Preamp_Gain = 1;
infoFile.Exposure_Time = 1;
infoFile.Frames_Per_Second = 1;
infoFile.camera_temperature = 1;
infoFile.number_of_frames = 1;
infoFile.camera_head = '';
infoFile.hstart = 1;
infoFile.hend = 256;
infoFile.vstart = 1;
infoFile.vend = 256;
infoFile.ADChannel = 0;
infoFile.Stage_X = 0;
infoFile.Stage_Y = 0;
infoFile.Stage_Z = 0;
infoFile.Lock_Target = 0;
infoFile.scalemax = 0;
infoFile.scalemin = 0;

[infFilePath, name, extension] = fileparts(infFileName);
infoFile.localName = [name extension];
infoFile.localPath = [infFilePath, filesep];

%--------------------------------------------------------------------------
% Parse each line and build ini structure
%--------------------------------------------------------------------------
for j=1:length(text)

    % Does the line contain a definition
    posEqual = strfind(text{j}, '=');
    if ~isempty(posEqual)

        %Parse value
        value = strtrim(text{j}((posEqual+1):end)); % Read value
        posX = strfind(value, ' x '); % Find a potential 'X'--a flag of a 2 element entry
        posColon = strfind(value, ':'); % Is there a colon?
        if ~isempty(posX) && isempty(posColon)% Parse both elements
            value1 = value(1:(posX-1));
            value2 = value((posX+2):end);
        end

        %Prepare field name
        fieldName = CoerceFieldName(text{j}(1:(posEqual-1)));

        %Prepare value
        if ~isempty(posX) && isempty(posColon)
            infoFile.(fieldName) = [str2num(value1) str2num(value2)];
        else
            fieldValue = str2num(value);
            if isempty(fieldValue) %The value is a string
                fieldValue = value;
            end
            infoFile.(fieldName) = fieldValue;
        end

    elseif strcmp(text{j}, 'information file for')  %If the line does not contain a definition, then the next line is the file name
        infoFile.file = text{j+1};
    end
end

if verbose
    display(['Loaded ' infFilePath, filesep, name, '.inf']);
end

%--------------------------------------------------------------------------
% Check frame dimensions
%--------------------------------------------------------------------------
if any(infoFile.frame_dimensions == 0)
    warning('matlabSTORM:corruptedInfoFile', 'Unexpected frame dimensions');
    infoFile.frame_dimensions = [infoFile.hend - infoFile.hstart + 1, ...
        infoFile.vend - infoFile.vstart + 1];
end

function fieldName = CoerceFieldName(inputName)
%--------------------------------------------------------------------------
% fieldName = CoerceFieldName(inputName) 
% This removes illegal characters from an inputName string forming a
% string that can be used as the name of a structure field
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% September 5, 2012
% jeffmoffitt@gmail.com
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Remove text between paranthesis, brackets, or any of the like
%--------------------------------------------------------------------------
paraSignals = {'[', ']', '{', '}', '(', ')'};
ind = [];
for i=1:length(paraSignals)
    ind = union(ind, strfind(inputName, paraSignals{i}));
end
if length(ind>1)
    inputName = inputName(setdiff(1:length(inputName), ind(1):ind(end)));
elseif length(ind) == 1
    inputName = inputName(setdiff(1:length(inputName), ind));
end

%--------------------------------------------------------------------------
% Trim whitespace and replace internal whitespace with underscores
%--------------------------------------------------------------------------
inputName = strtrim(inputName);
inputName(isspace(inputName)) = '_';

%--------------------------------------------------------------------------
% Remove all remaining non-alphabetic letters
%--------------------------------------------------------------------------
fieldName = inputName(isletter(inputName)|(inputName == '_'));