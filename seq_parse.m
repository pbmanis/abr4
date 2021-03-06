function seq = seq_parse(list)
% parse the list of the format:
% 12;23/10 etc... like nxtrec in datac
% now also parses matlab functions and array formats, using eval
%
% first arg is starting number for output array
% second arg is final number
% / indicates the skip arg type
% basic: /n means skip n : e.g., 1;10/2 = 1,3,5,7,9
% special: /##r means randomize order (/##rn means use seed n for randomization)
% special: /##l means spacing of elements is logarithmic
% special: /##s means spacing is logarithmic, and order is randomized. (/##sn means use seed n for randomization)
% special: /##o means octave (or fractions) between the first and last
% point. 1o is 1-octave spacing, based on the first value
% e.g., 2;16/o should yield 2,4,8,16 
% 2;16/0.5o should yield 2, 2^0.5*2, 2^1*2, 2^1.5*2, etc. ... 16. 
% 2;16/0.25o would be quarter octave spacing, 2, 2^0.25*2, 2^0.5*2,
% 2^0.75*2, etc. 
%
% multiple sequences are returned in a cell array

% 3 ways for list to be structured:
% 1. standard datac record parses. List is enclosed inbetween single quotes
% 2. matlab : (array) operator expressions. [0:10:100], for example
% 3. matlab functions (not enclosed in quotes). Each function generates a new list
% note that matlab functions and matrices are treated identically
%
% Updated 9/07/2000, 11/13/2000
% Paul B. Manis, Ph.D.
% pmanis@med.unc.edu

seq={};

if(nargout < 1)
    QueMessage('seq_parse: insufficient output arguments', 1);
    return;
end;
if(isnumeric(list))
    QueMessage('seq_parse: input must be string, not numeric', 1);
    return;
end;

if(strcmp(unblank(list), '') | isempty(list))
    QueMessage('seq_parse: input string appears to be empty', 1);
    return;
end;

% obtain tokens from input.
% There are two forms available.
% if there are spaces in the individual commands, then tokens must be separated by &
% Each token is pulled, and deblanked prior to use.
% if there is no &, then tokens are assumed to be separated by spaces.

if(find(list == '&')) % if list is &'d, then remove all spaces
    list=deblank(list); % first, remove extraneous stuff...
    token = '&'; % and let the & act as the separator
else
    token = ' '; % lists with no & use space delimiter...
end;

ndim = 1;
while(any(list))   % run through the whole list
    [alist, list] = strtok(list, token); % space separated lists are returned in sequential output variables
    % Now we must determine present list type as follows:
    % if list contains only members consisting of '0123456789.-+;/nrls' then list is datac
    % if list contains [ and ] then list is matlab array
    % is list contains neither, attempt to evaluate it as a matlab function,
    %   which must return a single dimensional matrix (array)
    %
    listtype = 0; % default is a matlab command - note we can mix types!

    k = find(any(ismember(alist, '[]')) == 1);
    if(rem(length(k), 2) == 0)
        listtype = 1; % Possibly matlab array - check details later
    end;

    if(all(ismember(unblank(alist), '0123456789.+-/;,nrlso!'))) % datac command
        listtype = 2; % datac format list
    end;

    switch(listtype)

        case {0, 1} % operate on matlab arrays or functions in the list
            try
                rp = eval([alist]);
            catch
                QueMessage(sprintf('seq_parse: Matlab was unable to evaluate %s\n', alist),1);
                return;
            end
            if(size(rp, 1) ~= 1)
                QueMessage(sprintf('seq_parse: Result must be array, not matrix: %s\n', alist),1);
                return;
            end;
            xo{ndim} = rp;
            ndim = ndim + 1;

        case 2 % handle lists with datac-format

            rp = [];
            while(any(alist))
                [plist, alist] = strtok(alist,',');
                [rpn, rperr] = recparse(unblank(plist));
                if(isempty(rpn) | rperr ~= 0) % catch errors from recparse
                    return;
                end;

                rp = [rp rpn];
            end
            xo{ndim} = rp;
            ndim = ndim + 1;

        otherwise
            QueMessage(sprintf('seq_parse: Failed to identify list type for %s\n', alist),1);
            return;
    end;
end;

% now place output in cell array out

% case for single output is special
l=length(xo);
if(l == 1)
    seq{1} = xo{1};
    return;
end;

% otherwise, setup
ins='';
outs='';
for i = 1:l
    in=sprintf('xo{%d}', i);
    if(i < l)
        in = [in ','];
    end;
    ins = [ins in];

    out=sprintf('s{%d}', i);
    if (i < l)
        out = [out ','];
    end;
    outs = [outs out];
end;
s = cell(1,l);
seq = cell(1,l);
cmd = sprintf('[%s] = ndgrid(%s);', outs, ins);
try
    eval(cmd, 'disp(''error'')');
catch
    QueMessage(sprintf('seq_parse: Failed to evaluate %s \n  -- ???? FATAL ???? --\n', cmd),1);
    seq={};
    return;
end;

for i = 1:l
    seq{i}=reshape(s{i}, 1, prod(size(s{i})));
end;
return;


function [recs, err] = recparse(list)
% function to parse basic word unit of the list - a;b/c or the like
%
recs=[];
err = 1; % assume an error
if(isempty(list))
    return;
end

fn=[]; ln=[]; sn=[]; arg=[]; seed=0;
[fn, rest] = strtok(list, ';');
[ln, rest]=strtok(rest, ';/');
[sn, rest] = strtok(rest,'/!');
[an, rest] = strtok(rest, '!');
u=isletter(sn); % see if skip argument has a type in its syntax
if(any(u)) % must have another kind of skip argument - l, n, r, s
    arg=sn(u==1);
    if(length(arg) > 1)
        QueMessage(sprintf('seq_parse: incorrect skip argument: %s', arg), 1);
        return;
    end;
    fu = find(u > 0);
    wsn = sn(fu(1)+1:end); % save whole skip text
    sn=sn(1:fu(1)-1);	% get the first number
    if(length(wsn) >= 1) % check for a seed in the argument (e.g: /10l5 would use a seed of 5)
        ps = str2num(wsn);
        if(~isempty(ps))
            seed = ps;
        else
            QueMessage(sprintf('seq_parse: Bad seed argument: %s', ps), 1);
            return;
        end;
    end;
end;

if(~isempty(rest))
    QueMessage(sprintf('seq_parse: Bad Record list to parse: %s', list),1)
    return;
end
if(~isempty(fn) & isempty(str2num(fn)))
    QueMessage(sprintf('seq_parse: Bad first arg: "%s" in "%s"', fn, list),1);
    return;
end
if(~isempty(ln) & isempty(str2num(ln)))
    QueMessage(sprintf('seq_parse: Bad last arg: "%s" in "%s"', ln, list),1);
    return;
end
if(~isempty(sn) & isempty(str2num(sn)))
    QueMessage(sprintf('seq_parse: Bad skip arg: "%s" in "%s"', sn, list),1);
    return;
end
if(~isempty(an) & isempty(str2num(an)))
    QueMessage(sprintf('seq_parse: Bad alternation arg: "%s" in "%s"', an, list),1);
    return;
end


if(isempty(sn) & ~isempty(fn) & ~isempty(ln))
    recs=[str2num(char(fn)):str2num(char(ln))];
    recs = alt(recs, an);
    err = 0;
    return;
end
if(~isempty(fn) & ~isempty(ln) & ~isempty(sn) & isempty(arg))
    recs=[str2num(char(fn)):str2num(char(sn)):str2num(char(ln))];
    recs = alt(recs, an);
    err = 0;
    return;
end
if(~isempty(fn) & ~isempty(ln) & ~isempty(sn) & ~isempty(arg))
    a=str2num(char(fn));
    b=str2num(char(ln));
    n=str2num(char(sn));
    switch(arg)
        case 'n' % just n steps between a and b
            sk=(b-a)/(n-1);
            recs=a:sk:b;
            recs = alt(recs, an);
        case 'r' % steps spaced at n between a and b, but randomize order..
            sk=n; % (b-a)/(n-1);
            recs=a:sk:b;
            rand('state', seed); % force systematic state
            [ignore,v] = sort(rand(1,length(recs)));
            recs=recs(v); % reorder result!
            recs = alt(recs, an);
        case 'l' % make log steps between the elements: the skip arg is the length of the resulting vector.
            if(a <= 0 | b <= 0)
                QueMessage(sprintf('seq_parse: arguments to log are negative or 0: %6.1f %6.1f', a, b));
                return;
            end;
            la = log(a);
            lb = log(b);
            if(n < 2)
                QueMessage(sprintf('seq_parse: not enough steps in log seq: %d', n));
                return;
            end;
            ls = (lb-la)/(n-1);
            lrecs=la:ls:lb;
            recs=exp(lrecs);
            recs = alt(recs, an);
        case 's' % make log steps, but then randomize them
            if(a <= 0 | b <= 0)
                QueMessage(sprintf('seq_parse: arguments to log are negative or 0: %6.1f %6.1f', a, b));
                return;
            end;
            la = log(a);
            lb = log(b);
            if(n < 2)
                QueMessage(sprintf('seq_parse: not enough steps in log seq: %d', n));
                return;
            end;
            ls = (lb-la)/(n-1);
            lrecs=la:ls:lb;
            recs=exp(lrecs);
            rand('state', seed); % force systematic state
            [ignore,v] = sort(rand(1,length(recs)));
            recs=recs(v); % reorder result!
            recs = alt(recs, an);
        case 'o' % octave spacing/ fractional octave spacing
            if(a <= 0 | b <= 0)
                QueMessage(sprintf('seq_parse: arguments to octave generator are negative or 0: %6.1f %6.1f', a, b));
                return;
            end;
            fo = a;
            f = fo;
            i = 1;
            recs=[];
            while(f <= b)
                recs(i) = f;
                i = i + 1;
                f = fo* 2^((i-1)*n);
            end;

                
    otherwise
            QueMessage(sprintf('seq_parse: unrecognized spacing/order: %s\n', sn));
            return;
    end;
    err = 0; % the otherwise handled the error: after the switch we should be ok
    return;
end

if(isempty(sn) & isempty(ln) & ~isempty(fn))
    recs=str2num(char(fn));
    err = 0; % indicate success
    return;
end

QueMessage(sprintf('seq_parse: Failed to parse: %s', list),1);
return;

function reco = alt(recs, an)
% put an in in alternate positions in recs (alternation with a single value)
%
if(isempty(an) | ~ischar(an))
    reco = recs;
    return;
end;
index = [1:length(recs)];
index1 = index * 2 - 1;
index2 = index * 2;
reco(index1) = recs;
reco(index2) = str2num(an);
return;

function [y]=unblank(x)
y=deblank(fliplr(deblank(fliplr(x))));
return;

function QueMessage(str, arg)

fprintf(1, '%s\n', str);
return;
