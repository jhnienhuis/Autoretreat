%changes in delta width?
%predefined initial bathymetry

clr

%Grid setup
noNodes = 200; % number of modeling nodes
dx = 20000; % x increment, unit is meter
dxp = dx; % previous dx for sed partition calculation
dxn = 1/noNodes; % distretized interval in moving boundary coordinates.
xn=0:dxn:1; %dimensionless grid

%simulation time-stepping
tf = 0.5; % ratio to adjust number of time steps and sea level rate
noTimeSteps = 3500/tf; % modeling time period
ti = 50/tf;  % sample interval for ploting figure
dt = 60*60*24*365*tf; % time increment of 1 year X tf
nt = (0:noTimeSteps-1)*tf;
%sampling simulation results
noSkipSteps = ti; %grab system information every noSkipSteps*dt time interval

%system thresholds (for autoretreat calculation)
rS_threshold = 0.1;
rS = ones(1,noNodes+1); %initial ratio of Sfrix/Sx

%channel system slope
Sfi = 2e-5; % initial slope of fluvial reach%system slopes
Sbb = 5e-5; % Slope of subaerial bedrock reach
Ssb = 0; % slope of subaqueous basement
Sfore = 2e-4; % slope of foreset.

% model settings (like the Mississippi)
%seaLevel = cumsum([0 repmat(10e-3*tf,1,2000/tf) repmat(2e-3*tf,1,2000/tf) repmat(0.5e-3*tf,1,(8000/tf)-1)]); %sea level over the modelling period (in meter)
seaLevel = cumsum([0 repmat(5e-3*tf,1,noTimeSteps)]);
River_Discharge = 100; %m3/s (bankfull, following Lamb et al., JGR 2012)
River_SusLoad = 10; %kg/s
River_Width = 10; %m
rB = 60; % ratio between channel and flood plain width
If = 0.5; %intermittancy
Cf = 1/22.5^2; % friction coefficient
kc = 75e-3; % Roughness height
ar = 8.1; % Coefficient in Manning-Strickler resistance relation
g = 9.81; % gravitational acceleration
D50 = 120e-6; % medium grain size
R = 1.65; % Submerged specific gravity
au = 1; % explicit vs. implicit
mud_fraction = rB; % fraction of mud deposited in flood plain
Sinu = 1.7; % Sinuosity
porosity = 0.4; % bed porosity

%dependent variables
qw = River_Discharge/River_Width;
un = qw^0.4*ar^0.6*g^0.3*Sfi^0.3*kc^-0.1; % uniform flow velocity m/s 
Qw = ones(1,noNodes+1).*qw;
sp = zeros(1,noNodes); % Flow surface profile
qs = zeros(1,noNodes); % sediment flux (bed load, sand)
Hn = qw/un; % normal flow depth

% fitted initial channel bed elevation and x node location
x = 0:dx:dx*noNodes;
x = x./1000;
xp = x; % x coordinates for previous time step
be = zeros(1,noNodes+1); %preallocate bed elevation
be(end) = -Hn; %initial bed elevation (simply normal flow depth?
be(1) = dx*noNodes*Sfi+be(end);%bed elevation
be = be(1):((be(end)-be(1))/noNodes):be(end);%bed elevation
bep = be; % bed elevation for previous time step
bei = be; % set the initial bed elevation for comparison

%inital basin configuration
bebasei = -20; % initial bed elevation at foreset-basement break
sbai = -dx*noNodes; % initial x coordinate at bed rock/alluvial contact
stfi = 0; % initial x coordinate at shoreline break
ssbi = stfi+(be(end)-bebasei)/Sfore; % initial x coodinate of foreset bottom set boundary
sba = sbai; % x coordinate at bed rock/alluvial contact
stf = stfi; % x coordinate at shoreline break
ssb = ssbi; % x coordinate at foreset/bottom set contact.
bebase = be(end)-Sfore*(ssb-stf); % elevation of foreset/basement break, the value is updated each step, the value here is not important.
bebasep = bebase; %elevation of basement break for previous step, foreset increment calculation
autosb = nan(1,noTimeSteps); % index of shoreline break during shoreline break
stfauto = []; % location of shoreline location during autoretreat
    
%basin rate variables
rsba = 0; % rate of migration of bedrock/alluvial contact 
rstf = 0; % rate of migration of shoreline break
rssb = 0; % rate of migration of foreset/bottom set contact

%storage variables
STF = zeros(1,noTimeSteps); % x coordinate at shoreline break throughout the modeling time frame
RSFT = zeros(1,noTimeSteps);% rate of migration of shoreline break for all the timed interval
STFauto = zeros(1,noTimeSteps); % collection of shoreline location during autoretreat

%preallocation of derivatives
dndt = zeros(1,noNodes+1); %elevation time der.
dqdxn = zeros(1,noNodes+1);  %sed flux time der.
Sx = zeros(1,noNodes+1); % slope over x space
dSxdx = zeros(1,noNodes+1); % slope change over x space
Sfrix = zeros(1,noNodes+1); % friction slope over x space
dSfdx = zeros(1,noNodes+1); % friction slope change over x space

x = xn.*(stf-sba)+sba;%initial node location



% flow depth calculation preallocation
H(1:noNodes) = Hn; % flow depth.
H(noNodes+1) = seaLevel(1)-be(end); % flow depth this sets downstream boundary value as well
kh1 = zeros(1,noNodes+1);%4th oder runge-kutta method
kh2 = zeros(1,noNodes+1);
kh3 = zeros(1,noNodes+1);
kh4 = zeros(1,noNodes+1);
dHdx = zeros(1,noNodes+1);

% pre-allocate model results output
%these fields are nxmxl; n = simulation number; m = grid nodes; 1 = time-step
noObservations = noTimeSteps/noSkipSteps;
tmp = nan(noObservations,noNodes+1);

flowDepth = tmp;
channelSlope = tmp;
aggradationRate = tmp;
grainSize = tmp;
sedimentFlux = tmp;
channelElevation = tmp;
nodeLocation = tmp;
backwaterLength = nan(1,noObservations);
backwaterLengthxIndex = nan(1,noObservations);
shorelineLocation = nan(1,noObservations);
shorelineLocation2 = nan(1,noObservations);
shorelineMigRate = nan(1,noObservations);
sedPartition =  nan(1,noObservations); 
basementPosition = nan(2,noObservations); 
maxAggradationX = nan(1,noObservations); 
shorelineIdx = nan(1,noObservations);

v = VideoWriter('delta2.avi');
v.FrameRate = 10;
v.Quality = 95;
open(v); f = figure(1);
f.Units = "inches";
f.Position = [2 2 11 6];


% backwater calculation
for tt=1:noTimeSteps
    idxu = (1:noNodes+1); % node index upstream of auto shoreline
    idxd = (1:noNodes+1); % node index downstream of auto shoreline
    idxu = idxu(rS>rS_threshold);
    idxd = idxd(rS<=rS_threshold);
    [colum, idxul]= size(rS(rS>rS_threshold));
    idx = (1:noNodes);

% solving the backwater equation using 4th-order runge-kutta method
    kh1(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./H(idx+1).^(10/3)))./(1-qw^2/g./H(idx+1).^3)*dxn*(stf-sba);      
    kh2(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./(H(idx+1)-kh1(idx)/2).^(10/3)))./(1-qw^2/g./(H(idx+1)-kh1(idx)/2).^3)*dxn*(stf-sba);        
    kh3(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./(H(idx+1)-kh2(idx)/2).^(10/3)))./(1-qw^2/g./(H(idx+1)-kh2(idx)/2).^3)*dxn*(stf-sba);  
    kh4(idx) = ((-be(idx+1)+be(idx))/dxn/(stf-sba)...
                -(kc^(1/3)*qw^2/ar^2/g./(H(idx+1)-kh3(idx)).^(10/3)))./(1-qw^2/g./(H(idx+1)-kh3(idx)).^3)*dxn*(stf-sba);  
    H(idx) = H(idx+1)-kh1(idx)/6-kh2(idx)/3-kh3(idx)/3-kh4(idx)/6;    
    H(H<=Hn) = Hn; % upstream flow depth might get below Hn due to model instability        
    qs = (R*g*D50.^3).^0.5./(R*g*D50).^3*0.0355/ar^4.*(kc./H).^(2/3).*(qw./H).^6; % calculate qs
    
    autosb(tt) = idxul;%storing location index of shoreline
    %{
    if idxu(end)~=noNodes+1
        H_auto(tt) = interp1(rS(autosb(tt):end),H(autosb(tt):end),rS_threshold,'spline');
        
    else
        H_auto(tt) = H(idxu(end));%storing flow depth at shoreline break
    end
    %}
    Qw = qw;
       
    % end of grain size calculation        
    stfauto = x(autosb(tt));
    STFauto(tt) = stfauto;

    if tt==1
        dqdxn(end) = (qs(end)-qs(end-1))/dxn;
        dqdxn(1) = (qs(2)-qs(1))/dxn;
        dndt(end) = 1e-5*(-Sfi)-dqdxn(end)*If*(1+mud_fraction)*Sinu/rB/(1-porosity)/(stf-sba);
        dndt(1) = 1e-5*(-Sfi)-dqdxn(1)*If*(1+mud_fraction)*Sinu/rB/(1-porosity)/(stf-sba);
    end

    rstf = 1/Sfore*(If*Sinu*(1+mud_fraction)*qs(end)/rB/(1-porosity)/(ssb-stf)-dndt(end)); % rate of migration of shoreline break
    rssb = 1/(Sfore-Ssb)*(Sfore*rstf+dndt(end)); % rate of rate of migration of foreset/bottom set contact
    rsba = -1/Sbb*dndt(1); % rate of migration of bedrock/alluvial contact

    if qs(end)==0
        rstf = 0;
        rssb = 0;
    end
    
    RSFT(tt) = rstf;
    Sxp = Sx(end); % Sx at downstream end from previous iteration, for sediment partition calculation 
    bet = be(end); % temporary be(end) for mass balance calculation
    Sxendi = Sx(end); % downstream slope at previous time step

    %calculate spatial derivatives of sediment flux, grain size fraction,
    %flow depth and active layer thickness
    %downstream point
    dqdxn(noNodes+1) = (qs(noNodes+1)-qs(noNodes))/dxn;
    dndt(noNodes+1) = (xn(noNodes+1)*rstf+(1-xn(noNodes+1))*rsba)/(stf-sba)*(be(noNodes+1)-be(noNodes))/dxn...
    -dqdxn(noNodes+1)*If*(1+mud_fraction)*Sinu/rB/(1-porosity)/(stf-sba);
    dHdx(noNodes+1) = (H(noNodes+1)-H(noNodes))/(x(noNodes+1)-x(noNodes));

    % upstream point
    dqdxn(1) = (qs(2)-qs(1))/dxn;
    dndt(1) = (xn(1)*rstf+(1-xn(1))*rsba)/(stf-sba)*(be(2)-be(1))/dxn...
    -dqdxn(1)*If*(1+mud_fraction)*Sinu/rB/(1-porosity)/(stf-sba);
    dHdx(1) = (H(2)-H(1))/(x(2)-x(1));

    % cases 2 to N points
    dqdxn(2:noNodes) = au*(qs(2:noNodes)-qs(1:noNodes-1))/dxn+(1-au)*(qs(3:noNodes+1)-qs(1:noNodes-1))/2/dxn;
    dndt(2:noNodes) = (xn(2:noNodes)*rstf+(1-xn(2:noNodes))*rsba)/(stf-sba).*(be(2:noNodes)-be(1:noNodes-1))/dxn...
    -dqdxn(2:noNodes)*If*(1+mud_fraction)*Sinu/rB/(1-porosity)/(stf-sba);
     dHdx(2:noNodes) = (H(2:noNodes)-H(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    
    %update bed elevation and flow surface profile
    be = be+dndt*dt;
    sp = be+H;

    Sx(noNodes+1) = -(be(noNodes+1)-be(noNodes))/(x(noNodes+1)-x(noNodes));
    Sx(1) = -(be(2)-be(1))/(x(2)-x(1));
    Sx(2:noNodes) = -(be(2:noNodes)-be(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    Sfrix = kc^(1/3)/ar^2.*Qw.^2./g./H.^(10/3);

    dSfdx(noNodes+1) = (Sfrix(noNodes+1)-Sfrix(noNodes))/(x(noNodes+1)-x(noNodes));
    dSfdx(1) = (Sfrix(2)-Sfrix(1))/(x(2)-x(1));
    dSfdx(2:noNodes) = (Sfrix(2:noNodes)-Sfrix(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    dSxdx(noNodes+1) = (Sx(noNodes+1)-Sx(noNodes))/(x(noNodes+1)-x(noNodes));
    dSxdx(1) = (Sx(2)-Sx(1))/(x(2)-x(1));
    dSxdx(2:noNodes) = (Sx(2:noNodes)-Sx(1:noNodes-1))./(x(2:noNodes)-x(1:noNodes-1));
    rS = Sfrix./Sx;
    
    H(end)=seaLevel(tt)-be(end);% update river mouth flow depth

    x = xn.*(stf-sba)+sba; % convert normalized x distance back to dimensional
    
    ssbp = ssb; % temporary ssb for mass balance calculation
    ssbt = ssb; % temporary ssb for mass balance calculation
    stft = stf; % temporary stf for mass balance calculation
    bebaset = bebase; % temporary bebase for mass balance calculation

    sba = sba+rsba*dt; % update sba
    ssb = ssb+rssb*dt; % update ssb
    stf = stf+rstf*dt; % update stf   
    STF(tt) = stf;
    ssby = be(end)-Sfore*(ssb-stf);%elevation of basement foreset break

    % calculate bed elevation increment
    if rstf>0
        BEInt = interp1(x,be,xp,'spline');  
        BEI = BEInt-bep;
    else
        BEInt = interp1(xp,bep,x,'spline');
        BEI = be - BEInt;
    end
    [depoFront,idx] = max(BEI);
    depoFrontX = x(idx);
    bep = be; % update BEP for next round of iteration

    % Calculate location of BW transition
    idx = (100:noNodes+1);
    idx = idx(dHdx(100:noNodes+1)<5e-6);
    bwl = x(autosb(tt))-x(idx(end));
    bwlxIdx = idx(end);
    %bexc(tt,:) = interp1(x,be,x(autosb(tt))-200000:10000:x(autosb(tt)),'spline');

    
    % Mass Balance calculation
    if rstf>0 %progradation
        fa1 = sum(BEI) * dxp - 0.5*BEI(end)*dxp;
        if be(end)<bep(end)
            fa2 = 0.5 * Sx(end) * (x(end)-xp(end))^2 - 0.5 * (bep(end)-be(end))*(x(end)-xp(end));
            da1 = 0.5 * (x(end)-xp(end)-(bep(end)-be(end))/Sfore) * (bep(end)-be(end));
            da2 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (bep(end)-bebasep);
        else
            fa2 = 0.5 * Sx(end) * (x(end)-xp(end))^2 + 0.5 * (be(end)-bep(end)) * (x(end)-xp(end));
            da1 = 0.5 * (x(end)-xp(end)+(bep(end)-be(end))/Sfore) * (be(end)-bep(end));
            da2 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (bep(end)-bebasep);
        end
        else  %rstf<0 shoreline retreat
            fa1 = sum(BEI) * dx - 0.5*BEI(end)*dx;
            fa2 = 0.5 * (xp(end)-x(end)) * (be(end)-bep(end)) - 0.5 * Sxp * (xp(end)-x(end))^2;
            da1 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (be(end)-bep(end));
            da2 = 0.5 * (ssb-ssbp - (bebasep-bebase)/Sfore) * (bep(end)-bebasep);  
    end
    da3 = 0.5 * (bebasep-bebase) * ((ssb-ssbp) - (bebasep-bebase)/Sfore);
    fa = fa1 + fa2;
    da = da1 + da2 + da3;
    fp = fa / (fa+da);
    sedPar = fp;
    
    xp = x; % update xp so for next iteration of calculation, this is placed after ploting because plotting needs previous x location.
    bep = be; % update BEP so for next iteration of calculation
    bebasep = bebase; % temporary bebase for mass balance calculation

   


    
    
    %store results and plot
    if mod(tt,noSkipSteps) == 0
        ic = floor(tt./noSkipSteps)+1
        flowDepth(ic,:) = H;
        channelSlope(ic,:) = Sx;
        aggradationRate(ic,:) = BEI;
        grainSize(ic,:) = D50;
        sedimentFlux(ic,:) = qs;
        channelElevation(ic,:) = be;
        nodeLocation(ic,:) = x;
        
        backwaterLength(ic) = bwl;
        backwaterLengthxIndex(ic) = bwlxIdx;
        shorelineLocation(ic) = stf;
        shorelineLocation2(ic) = stfauto;
        shorelineIdx(ic) = autosb(tt);
        shorelineMigRate(ic) = rstf;
        sedPartition(ic) = sedPar;
        basementPosition(1,ic) = ssb; %x coordinate of basement foreset break
        basementPosition(2,ic) = ssby; %y coordinate of basement foreset break
        maxAggradationX(ic) = depoFrontX; % location of maximum deposition
        
        
    % plotting

    xx = Hn;

    ax = gca;


    hold on,box on
    set(ax,'PlotBoxAspectRatio',[2 1 1],'FontSize',15,'FontName','Myriad Pro')

    blub(1) = area(ax,[-1000,1000],[seaLevel(tt) seaLevel(tt)]-H(autosb(tt))+Hn,-5,'FaceColor',[91 155 213]/255);
    

    %aa=plot(x./1e3,be+xx,'k'); %plot bed elevation

    aa = area(ax,x./1e3,be+Hn,-5,'FaceColor',[112 173 71]./255);
    aa(2) = area(ax,nodeLocation(2,:)./1e3,channelElevation(2,:)+Hn,-5,'FaceColor',[92 153 51]./255);


    bb = plot(ax,x(1:autosb(tt))./1e3,be(1:autosb(tt))+Hn,'Color',[34 63 89]./255,'Linewidth',4);
    cc = scatter(ax,x(autosb(tt))./1e3,be(autosb(tt))+Hn,40,'filled','MarkerFaceColor',[34 63 89]./255,'MarkerEdgeColor',[34 63 89]./255);



    %p1 = plot(stf/1e3,be(end),'bo');
    title(ax,['Time = ' num2str(tt/10) ' yrs'])
    if rssb>0
    cc=area(ax,[x(end)/1e3,ssb/1e3],[be(end),ssby]+xx,-5,'FaceColor',[112 173 71]./255); % plot delta foreset
    cc(2)=area(ax,[nodeLocation(2,end)/1e3,basementPosition(1,2)/1e3],[channelElevation(2,end),basementPosition(2,2)]+xx,-5,'FaceColor',[92 153 51]./255); % plot delta foreset

    
    end
    xlabel('distance (km)');
    ylabel('elevation (m)');
    p2 = [];

    %if rS(end)<rS_threshold
    %p2 = plot(stfauto/1e3,be(autosb(tt))+xx,'bo');
    %legend([p1,p2],{'shoreline before autobreak','shoreline after autobreak'});
    %end
    xlim(ax,[-500 500])
    ylim(ax,[-5 20])

    %blub(1) = plot(x./1e3,be+H,'-b');
    %blub(1) = plot([stfauto/1e3,1000],[be(autosb(tt)),be(autosb(tt))]+xx,'-b');
    

    %blub(1) = plot([stfauto/1e3,1000],[seaLevel(tt), seaLevel(tt)],'-b');

        %H(end)=seaLevel(tt)-be(end);% update river mouth flow depth




    frame = getframe(gcf);
    writeVideo(v,frame);    

    %cla(blub(1))
    delete(p2)
    delete(blub)
    delete(cc)
    delete(aa)
    delete(bb)

    %bb=plot(ax,x./1e3,be+xx,'Color',[0.5 0.5 0.5]);
    %if rssb>0
    %dd=plot(ax,[x(end)/1e3,ssb/1e3],[be(end),ssby]+xx,'Color',[0.5 0.5 0.5]); % plot delta foreset
    %end









        
    end

end
close(v);

%return values in the form of a structure
results.flowDepth = flowDepth;
results.channelSlope = channelSlope;
results.aggradationRate = aggradationRate;
results.maxAggradationX = maxAggradationX;
results.grainSize = grainSize;
results.sedimentFlux = sedimentFlux;
results.channelElevation = channelElevation;
results.nodeLocation = nodeLocation;

results.backwaterLength =  backwaterLength;
results.backwaterLengthIndex =  backwaterLengthxIndex;
results.shorelineLocation =  shorelineLocation;
results.shorelineLocation2 =  shorelineLocation2;

results.sedPartition = sedPartition;
results.basementPosition = basementPosition;
results.shorelineIdx = shorelineIdx;

%calculate autostratigraphic time and length scales
L_auto = qs(1)/10e-3*60*60*24*365/60/10*2*1.7/1e3/0.6;%autostratiraphic length scale
T_auto = Sfi * qs(1)/(10e-3)^2/60/10*2*1.7/0.6*(60*60*24*365);%autostratigraphic time scale



% plotting results
figure(1)
hold on
plot(x./1e3,be+H,'-b')
plot([stf/1e3,ssb/1e3],[seaLevel(end),seaLevel(end)],'-b')

figure(2)
yyaxis left
plot((0:noSkipSteps:noTimeSteps).*tf,(shorelineLocation2-shorelineLocation2(end)+backwaterLength(end)))
ylabel('shoreline trajectory')
yyaxis right
plot((1:noTimeSteps).*tf,diff(seaLevel)./tf)
ylabel('sea-level rise (m/yr)')
%}
%close all
