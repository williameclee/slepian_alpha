function defval(nameOfVar, valueOfVar)
    % DEFVAL(name,value)
    %
    % Assigns a default value to the named variable
    %
    % INPUT:
    %
    % name    A string, enclosed in single quotes, with a variable name
    % value   The value, whatever it is, that you want the variable to have
    %
    % OUTPUT:
    %
    %      None. The variables appear as if by magic into your workspace or
    %      will be available inside your function.
    %
    % NOTE:
    %
    % This won't work for an unassigned structure variable.
    %
    % Last modified by ebrevdo-at-alumni-princeton.edu, 05/28/2011
    % Last modified by fjsimons-at-alum.mit.edu, 12/20/2012
    %
    % It appears that defval('bla',functioncall) evaluates the function call
    % regardless of whether or not 'bla' has been assigned.

    % Note: basically this is an 'assign if DNE' function

    p = inputParser;
    addRequired(p, 'name', @ischar);
    addRequired(p, 'value');
    parse(p, nameOfVar, valueOfVar);
    nameOfVar = p.Results.name;
    valueOfVar = p.Results.value;

    existsAndNotEmpty = evalin('caller', ['exist(''' nameOfVar ''',''var'')']) && ... % exists and ...
        ~evalin('caller', ['isempty(' nameOfVar ')']); % ... is not empty

    if ~existsAndNotEmpty
        assignin('caller', nameOfVar, valueOfVar);
    end

end
