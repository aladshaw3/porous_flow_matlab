%% @package create_gif
%
%   @brief create_gif is a function to create a *.gif file from a set of
%   user provided commands
%
%   @details Function to build a gif animation from a set of user defined
%   commands for any type of plot in Matlab. The set of commands MUST be
%   provided as a character array, where each individual command must be
%   separated by the `;` character. Within the scope of the commands, there
%   should be an `INDEX` that represents how each frame of the `*.gif` is
%   created and changed. For instance, if you are ploting an image that
%   updates at each timestep, then the `INDEX` in your command set should
%   define a time index in the command set to plot each individual frame
%   on. 
%
%   @author Austin Ladshaw
%
%   @date 03/29/2023
%
%   @copyright This software was designed and built by Austin Ladshaw,
%   2023. It is openly available through the MIT License. 

%% Function to create gif animates from user commands and variables
%
%       Calling this function will create and save a `*.gif` animation file
%       based on user provided commands (which define how a frame of the
%       gif is to be created) and a variable set (which defines what the
%       variables are in the command set). The command set must have an
%       `INDEX` item somewhere within it to define what in the plotting
%       will change in each frame. For instance, the most basic thing to
%       update is a time index in a variable or data set to define which
%       data belongs to the current time index. 
%
%       EXAMPLE:
%
%           % Create list of amplitudes and x values
%           a = linspace(1,-1,30);
%           x = linspace(0,2*pi,30);
%           idx_limit = 30;
%
%           % Create matrix y where y(i,:) is the sinewave of the i-th frame
%           y = a'*sin(x);
%
%           % Define a command_set for how each frame is to be plotted
%           command_set = ['plot(x,y(INDEX,:));','axis([0 2*pi -1 1]);'];
%
%           % Define a variable_set for what each variable is in the commands
%           variable_set{1,1} = 'x'; variable_set{1,2} = x;
%           variable_set{2,1} = 'y'; variable_set{2,2} = y;
%
%           % Call the 'create_gif' function
%           output = create_gif(command_set,variable_set,idx_limit);
%
%       END EXAMPLE:
%
%   @param command_set A character array of commands separated by ';'
%   @param variable_set An Nx2 cell array where column 1 is a variable name
%                       and column 2 is the data for that variable
%   @param INDEX_LIMIT A scalar for how many frames there are in the data
%                       set. NOTE: If you provide a number that goes
%                       outside the bounds of your data, this will cause
%                       an error.
%   @param INDEX_START What frame-index to start from in the data
%                       (default = 1)
%   @param filename The name you wish to save the gif animation as 
%                   NOTE: The gif will be saved in
%                   'output/Gifs/{filename}.gif'.
%                   NOTE 2: default filename is 'create_gif_result'
function output_file = create_gif(command_set, variable_set, INDEX_LIMIT, INDEX_START, filename)
    %Validate user inputs
    arguments
        command_set (1,:) {mustBeText}
        variable_set (:,2) { iscell(variable_set) }
        INDEX_LIMIT (1,1) {mustBeNonnegative}
        INDEX_START (1,1) {mustBeNonnegative} = 1;
        filename (1,:) {mustBeText} = 'create_gif_result';
    end

    loc = 'output/Gifs';
    [status, msg, msgID] = mkdir(loc);
    if (status == 0)
        causeException = MException(msgID,msg);
        throw(causeException);
    end

    % Create the output file path and name
    output_file = join([ loc,'/', filename,'.gif'], '');

    % Setup all workspace variables 
    for i=1:size(variable_set,1)
        eval([variable_set{i,1},'=[];'])
        localupdate(variable_set{i,1},variable_set{i,2});
    end

    % Initialize the index and evaluate commands
    INDEX=INDEX_START;
    eval(command_set);
    gif(output_file,'overwrite',true);

    for INDEX=INDEX_START+1:INDEX_LIMIT
        eval(command_set);
        gif;
    end

    % Automatically close the open figure before exiting 
    close gcf

end

%% This helper function is used to assign variable values in the workspace of the above function
function localupdate(name,value)
    assignin('caller',name,value);
end