function [rmap, freqs, levels] = ABRA_1(arg)
%
% dispatch function that does all work for ABRA.m
% Modified 6/6/07 P. Manis for new ABR3 output - tone pips, spls, etc.
%

if(nargin == 0)
    return;
end
sr=1/100000; % 4/97635;
maxt = 8;
persistent abr_path
persistent filename
persistent t
persistent d

%if(isempty(abr_path))
switch(computer)
    case {'PC', 'WIN'}
        abr_path = 'c:\mat_datac\abr\';
    case {'MACI'}
        abr_path = '/users/pbmanis/documents/CurrentProjects/AR/Alex-ABR/';
    case 'MAC'
        abr_path = '/users/pmanis/desktop/abr/';
end
%end;

switch(arg)

    case 'open'
        [filename, pathname ] = uigetfile([abr_path '*.txt'], 'Open ABR text file');
        abr_path = pathname;
        if(filename==0)
            return;
        end
        [p, f] = fileparts(filename);
        %        fb = f(1:(length(f)-2));
        fb = f(1:14);
        fe = f(16:end);
        fn = [pathname fb 'n' fe '.txt'];
        fp = [pathname fb 'p' fe '.txt']; % make filenames for both polarities, regardless of which one was picked.
        fs = [pathname fb 'SPL.txt']; % sound pressure file
        ff = [pathname fb 'kHz.txt'];
        try
            dn = load(fn, '-ascii');
            nf = 1;
        catch exception
            nf = 0;
        %    dn = 0 * dp 
        end
        dp = load(fp, '-ascii');
        ds = load(fs, '-ascii');
        d = (dn+dp)/2;

        n = length(d);
        t = (0:sr:(n-1)*sr)*1000;
        [t0, i0] = find(t > 1);
        i0 = min(i0);
        [t1, i1] = find(t >= 6);
        i1 = min(i1);
        [t2, i2] = find(t > 6);
        i2 = min(i2);
        [t3, i3] = find(t > 10);
        i3 = min(i3);
        for i = 1:size(d,2)
            d(:,i) = d(:,i) - mean(d(1:i0,i));
        end

        htxt = findobj('tag', 'abra_info');
        [pfp, ffp, efp] = fileparts(fp);
        [pfn, ffn, efn] = fileparts(fn);
        if(~isempty(htxt) && ishandle(htxt))
            set(htxt, 'String', sprintf('Dir:%s\nFiles:\n%s\n%s', ...
                abr_path, ffn, ffp));
        end
        % make plot with data scaled all the same
        %
        mx = max(max(d));
        mn = min(min(d));
        if(mx > mn)
            mx1 = mx;
        else
            mx1 = mn;
        end
        m = size(d,2);
        bl=zeros(n,m); % matching baseline
        fullscale = (m+2)*mx1; % max plot...
        for i = 1:m
            %offset(i) = fullscale - (i)*mx1;
            offset(i) = (m-i+1)*mx1;
        end
        hf = findobj('tag', 'ABRA_g2'); % get the graphic window
        if(isempty(hf))
            subplot('Position', [0.75, 0.1, 0.2, 0.25]);
            set(gca, 'Tag', 'ABRA_g2');
            hf = findobj('tag', 'ABRA_g2'); % get the graphic window
        end
        axes(hf);
        hold off
        for i = 1:m
            plot(t, d(:,i)+offset(i));
            hold on
            plot(t, bl(:,i)+offset(i), '--');
            text(t(1), bl(1,i)+offset(i), sprintf('%6.1f dB', ds(m-i+1)), ...
                'VerticalAlignment', 'top');

        end
        if(fullscale < 0)
            ylims = [fullscale 0];
        else
            ylims = [0 fullscale];
        end
        set(gca, 'YLim', ylims);
        set(gca, 'xlim', [0 maxt]);
        set(gca, 'tag', 'ABRA_g2');
        xlabel('ms')
        calX = [0 0.0];
        calY = [0 2.0e-6];
        plot(calX, calY, 'k-', 'linewidth', 2);
        hf = findobj('tag', 'ABRA_g1'); % get the graphic window
        if(isempty(hf))
            hf = subplot('Position', [0.1 0.1 0.5 0.75]);
            set(hf, 'tag', 'ABRA_g1');
        end
        axes(hf);
        hold off
        for i = 1:m
            ks = mx1/max(abs(d(:,i)));
            plot(t, ks*d(:,i)+offset(i));
            hold on
            plot(t, bl(:,i)+offset(i), '--');
            text(t(1), bl(1,i)+offset(i), sprintf('%6.1f dB', ds(m-i+1)), ...
                'VerticalAlignment', 'top');
        end
        set(gca, 'YLim', ylims);
        set(gca, 'xlim', [0 maxt]);
        set(gca, 'tag', 'ABRA_g1');

        subplot('Position', [0.75, 0.5, 0.2, 0.25]);
        set(gca, 'Tag', 'ABRA_g3');
        hf = findobj('tag', 'ABRA_g3'); % get the graphic window
        axes(hf);
        hold off
        n = size(d, 2);
        pord = n:-1:1;
        pmax = max(d((i0:i1),:));
        pmax(pord) = pmax;
        nmax = abs(min(d((i0:i1),:)));
        nmax(pord) = nmax;
        bln = std(d((1:i0),:));
        bln(pord) = bln;
        plot(ds, pmax, 'go-', 'markersize', 2.0);
        hold on
        plot(ds, nmax, 'rs-', 'markersize', 2.0);
        plot(ds, bln, 'bx-', 'markersize', 1.8);

        hamax=findobj('tag', 'abra_attnmax');
        hamin=findobj('tag', 'abra_attnmin');
        hastep = findobj('tag', 'abra_attnstep');
        %         if(~isempty(hamax) & ~isempty(hamin))
        %             attmax = str2num(get(hamax, 'String'));
        %         attmin = str2num(get(hamin, 'String'));
        %         attstep = str2num(get(hastep, 'String'));
        %     else
        %         attmax = 80;
        %         attmin = 0;
        %         attstep = 5;
        %     end;
        %     if(attmax > 0) attmax = -attmax; end;
        %         if(attmin > 0) attmin = -attmin; end;
        %         if(attmax < -80) atm = attmax;
        %         else
        %             atm = min(at);
        %         end;
        set(gca, 'XLim', [min(ds) max(ds)]);
        yl = get(gca, 'Ylim');
        set(gca, 'Ylim', [0 yl(2)]);
        set(gca, 'tag', 'ABRA_g3');


    case 'Export' % write files for Igor Pro
        [p, f] = fileparts(filename);
        filen = [f(1:13) '.dat'];
        export2(t, d, filen, filename);




    case 'openmap' % get mulitple files and read the response map
        %         [filename, pathname ] = uigetfile([abr_path '*.txt'], ...
        % 'Open ABR text file', 'multiselect', 'on');
        [filename, pathname ] = uigetfile([abr_path '*.txt'], ...
            'Open ABR text file');
        if isequal(filename,0) || isequal(pathname,0)
            return;
        end;
        abr_path = pathname;
        %         if(~iscell(filename))
        %             ABRA_1('open');
        %             return;
        %         end;
        [p, f, e] = fileparts(filename);
        fb = f(1:14);
        fe = f(16:end);
        fn = [pathname fb 'n' fe '.txt'];
        fp = [pathname fb 'p' fe '.txt']; % make filenames for both polarities, regardless of which one was picked.
        fs = [pathname fb 'SPL.txt']; % sound pressure file
        ff = [pathname fb 'kHz.txt']; % frequency list file
        try
            dn = load(fn, '-ascii');
            nf = 1;
        catch exception
            dp = load(fp, '-ascii');
            dn = 0 * dp;
            nf = 0;
        end
        dsize = size(dn);
        ds = load(fs, '-ascii');
        if(exist(ff, 'file'))
            df = load(ff, '-ascii');
        else
            disp 'missing frequency file list'
            df=[];
        end;
        rmapstd = zeros(length(df), length(ds), 1); % 2 d map with response amplitudes
        rmapp2p = rmapstd;
        rmapnoise = rmapstd;
        rmapspec = rmapstd;
        hs = spectrum.yulear(1024);
        disp('df: ')
        df
        for i = 1:length(df) % for each frequency
            % construct the expected file for this frequency
            %        fb = f(1:(length(f)-2));
            fe = sprintf('-%.3f', df(i));
            fn = [pathname fb 'n' fe '.txt'];
            fp = [pathname fb 'p' fe '.txt']; % make filenames for both polarities, regardless of which one was picked.

            dp = load(fp, '-ascii');
            if nf == 1
                dn = load(fn, '-ascii');
            else
                dn = 0*dp;
            end;
            d = (dn+dp)/2;
            c=find(isnan(d));
            [ci, cj] = ind2sub(size(d), c);
            d(ci, cj) = d(ci-1, cj); % remove NaNs...
            n = length(d);
            t = [0:sr:(n-1)*sr]*1000;
            [t0, i0] = find(t > 1);
            i0 = min(i0);
            [t1, i1] = find(t >= maxt);
            i1 = min(i1);
            [t2, i2] = find(t > maxt);
            i2 = min(i2);
            [t3, i3] = find(t > 10);
            i3 = min(i3);
            ttwo = find(t >= 2, 1, 'first');
            tfour = find(t <= 4, 1, 'last');
            %      d = detrend(d);
            for j = 1:size(d,2)
                d(:,j) = detrend(d(:,j)) - mean(d(1:i0,j));
            end;
            if(i == 1)
                disp 'ddmap making'
                ddmap = zeros(length(df), size(d, 2), size(d, 1));
            else
                disp 'did not make ddmap'
                fprintf(1, 'i = %d', i);
            end;
            ddmap(i,:,:) = d'; % all traces into a big matrix now
            rmapstd(i, :) = std(d); % rms signal
            rmapp2p(i, :) = max(d(ttwo:tfour, :)) - min(d(ttwo:tfour, :)); % peak - peak signal
            rmapnoise(i, :) = std(d(1:i0, :)); % baseline rms

            for j = 1:size(d, 2)
                abrspec = psd(hs, d(:,j), 'Fs', 100000);
                if(j == 1)
                    if0 = find(abrspec.Frequencies > 800, 1, 'first');
                    if0 = min(if0);
                    if1 = find(abrspec.Frequencies < 1200, 1, 'last');
                    if1 = max(if1);
                end;
                rmapspec(i,j) = sum(abrspec.Data(if0:if1))/(if1-if0);
            end;

            htxt = findobj('tag', 'abra_info');
            [pfp, ffp, efp] = fileparts(fp);
            [pfn, ffn, efn] = fileparts(fn);
            if(~isempty(htxt) && ishandle(htxt))
                set(htxt, 'String', sprintf('Dir:%s\nFiles:\n%s\n%s', ...
                    pathname, ffn, ffp));
            end;
        end;

        % make plot with data scaled all the same
        %
        mx = max(max(max(ddmap)));
        mn = min(min(min(ddmap)));
        if(mx > mn)
            mx1 = mx;
        else
            mx1 = mn;
        end;
        % ddmap dim 1 is freq
        % ddmap dim 2 is intensity

        bl=zeros(size(ddmap, 3), 1); % matching baseline
        fullscale = (size(ddmap, 2) + 2)*mx1; % max plot...
        for i = 1:size(ddmap, 2)
            %offset(i) = fullscale - (i)*mx1;
            yoffset(i) = i*mx1;
        end;
        for i = 1:size(ddmap, 1)
            xoffset(i) = (i-1) * (maxt+2);
        end;
        hf = findobj('tag', 'ABRA_g2'); % get the graphic window
        if(isempty(hf))
            hf = subplot('Position', [0.05 0.1 0.4 0.7]);
            set(hf, 'tag', 'ABRA_g2');
        else
            set(hf, 'position', [0.05 0.1 0.4 0.7]);
        end;
        axes(hf);
        hold off
        for j = 1:size(ddmap, 1)
            for i = 1:size(ddmap, 2)
                %  dm=ddmap(j, i, 1:i1);
                plot(t(1:i1)+xoffset(j), shiftdim(ddmap(j, i, 1:i1))+yoffset(i));
                hold on
                plot(t(1:i1)+xoffset(j), bl(1:i1)+yoffset(i), '--');
            end;
        end;
        if(fullscale < 0)
            ylims = [fullscale 0];
        else
            ylims = [0 fullscale];
        end;
        set(gca, 'YLim', ylims);
        set(gca, 'xlim', [0 (maxt+2)*(size(ddmap, 1))]);
        set(gca, 'tag', 'ABRA_g2');

        hf = findobj('tag', 'ABRA_g3'); % get the graphic window
        if(~isempty(hf))
            delete(hf);
        end;
        hf = findobj('tag', 'ABRA_g1'); % get the graphic window
        if(~isempty(hf))
            delete(hf);
        end;
        %         subplot('Position', [0.7, 0.5, 0.25, 0.25]);
        %         set(gca, 'Tag', 'ABRA_g3');
        %         hf = findobj('tag', 'ABRA_g3'); % get the graphic window
        %         axes(hf);
        %         hold off
        %
        %         pmax = max(d((i0:i1),:));
        %         nmax = abs(min(d((i0:i1),:)));
        %         bln = std(d((1:i0),:));
        %         at = ds; % -[0:5:5*(length(pmax)-1)];
        %         plot(at, pmax, 'go-', 'markersize', 2.0);
        %         hold on
        %         plot(at, nmax, 'rs-', 'markersize', 2.0);
        %         plot(at, bln, 'bx-', 'markersize', 1.8);
        %
        %         hamax=findobj('tag', 'abra_attnmax');
        %         hamin=findobj('tag', 'abra_attnmin');
        %         hastep = findobj('tag', 'abra_attnstep');
        %         set(gca, 'XLim', [min(ds) max(ds)]);
        %         yl = get(gca, 'Ylim');
        %         set(gca, 'Ylim', [0 yl(2)]);
        %         set(gca, 'tag', 'ABRA_g3');


        subplot('Position', [0.55, 0.1, 0.4, 0.7]);
        set(gca, 'Tag', 'ABRA_g4');
        hf = findobj('tag', 'ABRA_g4'); % get the graphic window
        axes(hf);
        hold off

        dfk = df/1000;
        %     quiver(dfk, ds, zeros(size(rmapstd')), rmapstd');
        quiver(dfk, ds, zeros(size(rmapspec')), rmapspec');
        hold on
        quiver(dfk, ds, rmapp2p', rmapp2p');
        quiver(dfk, ds, rmapnoise', zeros(size(rmapnoise')));
        %cmap = [meshgrid(dfk, ds), rmapstd];
        contour(dfk, ds, rmapspec');
        set(gca, 'xlim', [0.8*min(dfk) 1.2*max(dfk)]);
        set(gca, 'xscale', 'log');
        drawnow
        % threshold detection
        thr = 1.25; % standard deviations from lowest stimulus level
        wt_p2p = 0.3;
        wt_pspec = (1-wt_p2p);
        noise = sqrt(wt_pspec*mean(rmapspec(:, 1))^2 + wt_p2p*mean(rmapp2p(:, 1))^2);
        tc = NaN * zeros(length(dfk), 1); % threshold curve points
        for j = 1: length(dfk)
            s = find(sqrt(wt_pspec*rmapspec(j, :).^2+wt_p2p*rmapp2p(j, :).^2) > thr*noise); % indices into those with response
            slr = fliplr(s);
            u  = find (diff(slr) == -1, 1, 'last');

            if(isempty(u))
                tc(j) = NaN; % set to a high value
            else
                tc(j) = ds(slr(u)-1);
            end;
        end;
        plot(dfk, tc, 'k--', 'LineWidth', 2); % dashed thick black line
        fnout = [pathname fb '-thr.txt'];

        hout = fopen(fnout, 'w');
        for j = 1: length(dfk)
            fprintf(hout, '%.3f, %.3f\n', dfk(j), tc(j));
        end;
        fclose(hout);
        if(nargout > 0)
            rmap = rmapspec;
            freqs = dfk;
            levels = ds;
        end;

    otherwise
        return;

end;

% export the data to a text file

function export2(t, d, fn, filename)
    %
    % function to export data traces to a binary file to be read by Igor Pro.
    %

    rec = str2num(filename(15:18));
    hw=fopen(fn, 'w');
    if(isempty(hw) || hw <= 0)
        fprintf(1, 'Unable to open output file %s', fn);
        return;
    end;
    prefix = 'd';
    fprintf(hw, '%% File: %s %s', fn, filename);
    for k = 1:length(d)
        fprintf(hw, 'd_%d ', k);
    end;
    fprintf(hw, '\r\n');
    fprintf(hw, '"t%d"', rec);
    for k = 1:size(d,2)
            fprintf(hw, '\t"%sR%dI%d"', prefix, rec, k); % titles of variables for igor
    end;
    fprintf(hw, '\r\n');
    for j = 1:length(t) % now we write the data array
        fprintf(hw, '%f', t(j));
        for k = 1:size(d,2)
            fprintf(hw, '\t%f', 1.0e6*d(j,k));
        end;
        fprintf(hw, '\r\n');
    end;
    fclose(hw);



