%% Dynamic time warping
dista = zeros(length(pTot));
for i = 1:length(pTot)
    for j = 1:length(pTot)
        dista(i,j) = (pTot(i)-M1Tot(j))^2;
    end
end

% Calculate Accumulative dists (DTW)
AccCos = zeros(length(pTot));
for i = 2:length(pTot)
    AccCos(1,i) = dista(1,i) + AccCos(1, i-1);
end
for i = 2:length(pTot)
    AccCos(i,1) = dista(i,1) + AccCos(i-1, 1);
end
for i = 2:length(pTot)
    for j = 2:length(pTot)
        AccCos(i,j) = min([AccCos(i-1,j-1), AccCos(i-1,j), AccCos(i,j-1)]) + dista(i,j);
        if (i>=j)
            AccCos(i,j) = AccCos(i,j)*100;
        end
        if (j>i+3)
            AccCos(i,j) = AccCos(i,j)*(j-i);
        end
    end
end

% Find path (DTW)
path = [];
i = length(pTot)-1;
j = length(M1Tot)-1;
while i>1 && j>1
    if i==1
        j = j - 1;
    elseif j==1
        i = i - 1;
    else
        if AccCos(i-1,j) == min([AccCos(i-1,j-1), AccCos(i-1,j), AccCos(i,j-1)])
            i = i - 1;
        elseif AccCos(i,j-1) == min([AccCos(i-1,j-1), AccCos(i-1,j), AccCos(i,j-1)])
            j = j-1;
        else
            i = i - 1;
            j= j- 1;
        end
    end
    path = [path; j i];
end
path = [path; 1 1];
path = flip(path);
warp_p53 = pTot(path(:,2));
warp_mdm2 = M1Tot(path(:,1));

%% Plot warped signals
figure('Position',[959   694   958   300])
nexttile(1)
plot(warp_p53,'LineWidth',1.5);
xlim([1 max([length(warp_p53),length(warp_mdm2)])]);

hold on
plot(warp_mdm2,'LineWidth',1.5);
ylabel('Level')

%% Calculate velocity
vp=[NaN;diff(warp_p53)];
vm=[NaN;diff(warp_mdm2)];

%% Skip detection

skip_p53=NaN(size(M1Tot));
skip_mdm2=NaN(size(M1Tot));
skip_start=NaN(size(M1Tot));
skip_p53_dtw=NaN(size(warp_mdm2));
skip_mdm2_dtw=NaN(size(warp_mdm2));
skip_start_dtw=NaN(size(warp_mdm2));
skip_duration=[];
skip_duration_dtw=[];
detection_point=cell(2,1);
fft_range=[];
riseloc=[];
skip_per_p53=[];
skip_ratio=[];


for i=2:length(plocs)-1
    
    % Find Mdm2 peak between 2 p53 peaks
    idx=find(and(M1locs>=plocs(i)...
        ,M1locs<=plocs(i+1)));
    
    % If no Mdm2 peak between 2 p53 peaks
    if isempty(idx)
        
        x2=plocs(i);    % p53 peaks of interests
        
        x3=plocs(i-1);  % Previous p53 peak
        x4=plocs(i+1);  % Next p53 peak
        
        x5=ptrlocs(and(ptrlocs>x3,ptrlocs<x2)); % Previous p53 valley
        x6=ptrlocs(and(ptrlocs>x2,ptrlocs<x4)); % Next p53 valley
        
        % Convert to warped coordinates
        xx2=find(path(:,2)==x2,1,'first');

        xx3=find(path(:,2)==x3,1,'last');

        xx4=find(path(:,2)==x4,1,'first');
        if isempty(xx4)
            continue
        end
        
        xx5=find(path(:,2)==x5,1,'last');

        xx6=find(path(:,2)==x6,1,'first');

        % Filter with relative velocity
        vp1=vp(xx5:xx2);
        vp_pos=sum(vp1(vp1>0));
        vm1=vm(xx5:xx2);
        vm_pos=sum(vm1(vm1>0));
        riv=vm_pos/vp_pos;

        if riv>10^(-1.2)
            continue
        end
        
        % FFT filter
        L=x4-x3+1;
        Fs=2;
        X=pTot(x3:x4);
        X=detrend(X,2);
        Y=fft(X);
        P2 = abs(Y/L);
        if mod(L,2)==0
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(L/2))/L;
        else
            P1 = P2(1:(L-1)/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(L-1)/2)/L;
        end
        
        [maxpw,idx]=max(P1(f>0));
        maxfreq=f(P1==maxpw);
        main_per_p53=1/maxfreq;
        
        if main_per_p53>=7.375|| main_per_p53<2.5
            continue
        end
        
        skip_per_p53=vertcat(skip_per_p53,main_per_p53);
        skip_ratio=vertcat(skip_ratio,riv);
        fft_range=vertcat(fft_range,[path(x3,2), path(x4,2)]);
        
        skip_p53(path(xx5,2):path(xx6,2))=pTot(path(xx5,2):path(xx6,2));
        skip_mdm2((path(xx5,1)):(path(xx6,1)))=M1Tot((path(xx5,1)):(path(xx6,1)));
        
        skip_p53_dtw(xx5:xx6)=warp_p53(xx5:xx6);
        skip_mdm2_dtw(xx5:xx6)=warp_mdm2(xx5:xx6);

        skip_start(path(xx2,2))=50;
        
        detection_point{1}=vertcat(detection_point{1},plocs(i));
        
    end
    
end

%% Add markers to figure
figure(fg)
nexttile(1);
hold on
plot(t,skip_p53,'b','LineWidth',2.5);
plot(t,skip_start-20,'LineStyle','none','Marker','v' ...
    ,'MarkerFaceColor',c(2,:),'MarkerEdgeColor','k'...
    ,'MarkerSize',6);

nexttile(2);
hold on
plot(t,skip_mdm2,'r','LineWidth',2.5);
plot(t,skip_start,'LineStyle','none','Marker','v'...
    ,'MarkerFaceColor',c(2,:),'MarkerEdgeColor','k'...
    ,'MarkerSize',6);