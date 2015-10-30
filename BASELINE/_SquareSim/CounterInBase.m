function cb = CounterInBase(inbase)
% COUNTERINBASE Counter in the base specified as the argument.
% 
% Description:
%   c = CounterInBase(3)
%      Creates a new counter with base 3 and resets its state to 0.
%
%   c.next(), or
%   c.next(1)
%      Returns the state of the counter and increments it by 1.
%
%   c.next(0)
%      Returns the state of the counter and does not change its state.
%
%   c.next(20)
%      Returns a matrix of the next 20 numbers and advances the state of
%      the counter by 20.
%
%   c.save('filename')
%      Saves the counter state to disk. If no filename is given, it uses
%      the default name 'CounterInBaseState'.
%
%   c.load('filename')
%      Loads the previously stored counter state from disk. If no filename
%      is given, it loads from file with the default name
%      'CounterInBaseState'. 
%      NOTE: To be able to call the c.load() method, the counter must first
%      exist: c.load() must be preceded by the call to c=CounterInBase(x),
%      i.e. the c counter may be in any state with any base.
%
% Example:
%
%   >> c = CounterInBase(3);
%   >> c.next()
%   ans =
%        0
%   >> c.next()
%   ans =
%        1
%   >> c.next()
%   ans =
%        2
%   >> c.next()
%   ans =
%        1     0
%   >> c.next(3)
%   ans =
%        1     1
%        1     2
%        2     0
%   >> c.next(4)
%   ans =
%        0     2     1
%        0     2     2
%        1     0     0
%        1     0     1
%   >> c.next(0)
%   ans =
%        1     0     2
%   >> c.next(0)
%   ans =
%        1     0     2

% Filename     : CounterInBase.m
% Author       : Petr Posik (posik#labe.felk.cvut.cz, replace # with @ to get email)
% Created      : 04-Feb-2010 16:25:19
% Modified     : $Date: 2010-02-05 19:34:28 +0100 (p?, 05 II 2010) $
% Revision     : $Revision: 4750 $
% Developed in : 7.5.0.342 (R2007b)
% $Id: CounterInBase.m 4750 2010-02-05 18:34:28Z posik $

    %% Public methods
    cb.next = @next;
    cb.set = @set;
    cb.save = @mysave;
    cb.load = @myload;
    cb.getBase = @getBase;
    %% Public data
    %% Private data
    base = 10;
    numberbase = 0;
    numberdec = 0;
    defaultfname = 'CounterInBaseState';
%% Initialization
    if inbase < 2, 
        error('CounterInBase:WrongArgument', ...
            'The base must be greater or equal to 2.'); 
    end
    base = inbase;
%% Public methods implementation
    function b = getBase()
        b = base;
    end
    function n = next(nextn)
        % Return the 'nextn' states of the counter, advance the counter by
        % 'nextn' states.
        if nargin < 1 || isempty(nextn), nextn = 1; end
        if nextn < 0, 
            error('CounterInBase:WrongArgument',...
                'When calling next(n), n must be nonnegative.'); 
        end
        if nextn == 0, n = numberbase; return; end        
        % Advance the counter by the nextn
        n = zeros(nextn, numel(numberbase));
        for inext = 1:nextn,
            % Update the size of the output matrix if needed
            if numel(numberbase) > size(n,2),
                n = [zeros(nextn,1) n];
            end
            % Assign the counter state to the appropriate output slot 
            n(inext,:) = numberbase;
            % Add 1 to decimal state
            numberdec = numberdec + 1;
            % Add 1 to the state in the base-b
            for i = numel(numberbase):-1:1,
                numberbase(i) = mod(numberbase(i) + 1, base);
                if numberbase(i) ~= 0, break; end
                if i == 1,
                    numberbase = [1 numberbase];
                end
            end
        end
    end
    function set(n)
        % Set the state of the counter
        numberdec = n;
        numberbase = dec2base(n);
    end
    function mysave(fname)
        if nargin < 1 || isempty(fname),
            fname = defaultfname;
        end
        save(fname, 'base', 'numberdec', 'numberbase');
    end
    function myload(fname)
        if nargin < 1 || isempty(fname),
            fname = defaultfname;
        end
        S = load(fname, 'base', 'numberdec', 'numberbase');
        base = S.base;
        numberdec = S.numberdec;
        numberbase = S.numberbase;
    end
%% Private methods implementation
    function nb = dec2base(nd)
        % This is a code snippet from the MATLAB dec2base function,
        % however, it is not its precise copy.
        nd = nd(:);
        ndigits = max(1,round(log(max(nd)+1)/log(base)));
        while any(base.^ndigits <= nd)
            ndigits = ndigits + 1;
        end
        nb = zeros( numel(nd), ndigits );
        for i = ndigits:-1:1,
            nb(:,i) = mod(nd,base);
            nd = floor(nd/base);
        end
    end
end