%% clusters

clear all; 

bases={'fr2020';'al2020';'it2020';'es2020';'en2020'};


for len = 1:length(bases)
    language = bases{len};
    load(['./' language '/TIMELINE'])

    minword = 9;

    desde = find(double(TIMELINE.words(1).years)==1750);
    hasta = find(double(TIMELINE.words(1).years)==2000);

    TREND = [TIMELINE.words.trend]';
    OSC = [TIMELINE.words.smoothed]';
    
    OSCshort = OSC(:,desde:hasta);
    TRENDshort = TREND(:,desde:hasta);
    for indp=1:length(TIMELINE.words)
        OSCshort(indp,:) = (OSCshort(indp,:)-mean(OSCshort(indp,:)))/max(OSCshort(indp,:));
    end
    
    %trend clusters
    cutoff = 0.045;
    PD = pdist(TRENDshort,'correlation');
    Z = linkage(PD,'average');
    T = cluster(Z,'Cutoff',cutoff,'Criterion','distance');
    if shuffle==0
        size_com = zeros(max(T),1);
        for ind = 1:max(T)
            size_com(ind) = length(find(T==ind));
        end
        megustan = find(size_com>minword & size_com<maxword);
        [~,index] = sort(size_com(megustan),'descend');
        megustan = megustan(index);
        
        save(sprintf('../%s/clusters_trend',language)...
            ,'T','megustan','desde','hasta','cutoff')
        
        Ncom = length(megustan);
        ParOrden = nan(Ncom,length(desde:hasta));
        Ang = nan(Ncom,length(desde:hasta));
        for comunidad = 1:Ncom
            index = T == megustan(comunidad);
            osc = OSCshort(index,:);
            [par,ang] = PolarFun.ParOrden(osc);
            ParOrden(comunidad,:) = par;
            Ang(comunidad,:) = ang;
        end
        save(sprintf('./%s/stat_true',language),'ParOrden','Ang');
        
        %shuffle clusters
        T = T(randperm(length(T)));
        for ind = 1:max(T)
            size_com(ind) = length(find(T==ind));
        end
        
        megustan = find(size_com>minword);
        [~,index] = sort(size_com(megustan),'descend');
        megustan = megustan(index);
        
        Ncom = length(megustan);
        ParOrden = nan(Ncom,length(desde:hasta));
        Ang = nan(Ncom,length(desde:hasta));
        for comunidad = 1:Ncom
            index = T==megustan(comunidad);
            osc = OSCshort(index,:);
            [par,ang] = PolarFun.ParOrden(osc);
            ParOrden(comunidad,:) = par;
            Ang(comunidad,:) = ang;
        end
        
        
        save(sprintf('./%s/stat_shuffle',language),'ParOrden','Ang');
        
        %osc clusters
        cutoff = 0.5;
        PD = pdist(OSCshort,'correlation');
        Z = linkage(PD,'average');
        T = cluster(Z,'Cutoff',cutoff,'Criterion','distance');
        
        % Comunidades con distancia menor que cutoff y más de minword palabras
        for ind = 1:max(T)
            size_com(ind) = length(find(T==ind));
        end
        megustan = find(size_com>minword);
        [~,index] = sort(size_com(megustan),'descend');
        megustan = megustan(index);
        save(sprintf('./%s/clusters_osc',language),'T','megustan','desde','hasta')
    end
end



%% figura 1

clear all
close all
figure(1);clf
set(gcf,'color','w')
set(gcf,'position',[680 298 1366 768])
bordeizq = 0.055;
ancho = 0.29;
bordeinf = 0.07;
alto = 0.41;
sepv = 0.08;
seph = 0.03;

clear handles
handles(1)=axes('position',[bordeizq                bordeinf+alto+sepv ancho alto     ]);
handles(2)=axes('position',[bordeizq                bordeinf           ancho alto     ]);
handles(3)=axes('position',[bordeizq+ancho+seph     bordeinf+alto      ancho alto+sepv]);
handles(4)=axes('position',[bordeizq+ancho+seph     bordeinf           ancho alto     ]);
handles(5)=axes('position',[bordeizq+2*ancho+2*seph bordeinf+alto      ancho alto+sepv]);
handles(6)=axes('position',[bordeizq+2*ancho+2*seph bordeinf           ancho alto     ]);
letsize = 12; %font size

% Panel A. Top
language = 'en2020';
load(['./' language '/TIMELINE'])

gpalabras = { 'time' 'work' 'god'}; 
index = arrayfun(@(x) find(ismember({TIMELINE.words.word},x)),gpalabras);

set(gcf,'currentaxes',handles(1)); hold all
hasta = 2000;
line([TIMELINE.desde hasta],[0 0],'color',[.5 .5 .5])
line([1750 1750], [-1.5e-4 2.5e-4],'color','k');
colores  = [75 156 211; 15 82 186; 0 35 102]/256;
for indp = 1:length(gpalabras)
    set(gcf,'currentaxes',handles(1))
    hold all
    
    years = double(TIMELINE.words(index(indp)).years);
    freqrel = TIMELINE.words(index(indp)).freqrel;
    trend = TIMELINE.words(index(indp)).trend;
    osc = TIMELINE.words(index(indp)).smoothed;
    total = trend+osc;    
    plot(years,TIMELINE.words(index(indp)).freqrel,'.','Color',colores(indp,:))    
    plot(years,trend,'linewidth',.5,'color','k');
    plot(years,total,'linewidth',.5,'color',colores(indp,:),'linewidth',1.5)
    ylabel('Word frequency','Fontsize',letsize)
    xlabel('Year','Fontsize',letsize)
    xlim([TIMELINE.desde+50 hasta])
    ylim([0 2e-3])
    set(gca,'FontSize',letsize)
    indexpos = length(years)-8;
    if indp == 3
        text(years(indexpos)-10,trend(indexpos)+0.0002,gpalabras{indp},...
            'horizontalalignment','right','fontsize',15)
    elseif indp == 2
        text(years(indexpos)-8,trend(indexpos)+.00015,gpalabras{indp},...
            'horizontalalignment','right','fontsize',15)
    else
        text(years(indexpos)-8,trend(indexpos)+.0001,gpalabras{indp},...
            'horizontalalignment','right','fontsize',15)        
    end
   
end
%

% PAnel A. bottom
load('./Wavelets')

ancho = 50;
step = 0.25;
Scales = 1:step:ancho;

bases = {'en2020';'es2020';'fr2020';'al2020';'it2020'};

set(gcf,'currentaxes',handles(2))
hold all

for len=1:length(bases)
    
    PerMax = ((8*pi^2/5)^(1/2))*M(len).ScalesMax;
    Periodos = ((8*pi^2/5)^(1/2))*Scales;
    [N,edges] = histcounts(PerMax,Periodos);
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    scatter(edges,N,20,'filled');
    [mmax,indmax] = max(N(2:end-1));
    fprintf('%s %2.1f\n',bases{len},edges(indmax))
    box on 
    set(gca, 'YGrid', 'off', 'XGrid', 'on')
end

legend({'English','Spanish','French','German','Italian'})
set(gca,'ytick',0:2000:20000,'fontsize',letsize)
xlim([5.5 80])
ylim([0 9999])
xlabel('Periods (years)','fontsize',letsize)
ylabel('Count','fontsize',letsize)


% Panel B

CuantasCloud = 15;
load(['./' language '/clusters_trend'])
load(['./' language '/TIMELINE'])
elegidatr = 42;

comunidad = elegidatr; 
index = find(T==megustan(comunidad));
palabras = {TIMELINE.words(index).word};
index_com = find(ismember({TIMELINE.words.word},palabras));
word = 'copper';
indcolor = find(ismember({TIMELINE.words(index_com).word},word));

wordcolor = zeros(length(index),3);
wordcolor(indcolor,:) = [1 0 0];
PAL = struct('palabra',palabras,'count',num2cell([TIMELINE.words(index_com).tot]));
wordcloud({PAL.palabra},[PAL.count],'position',[.38 .44 .28 .16],...
    'MaxDisplayWords',CuantasCloud,'Color',wordcolor);

set(gcf,'currentaxes',handles(3));
totales_com = [TIMELINE.words(index_com).tot];
[~,indsorted] = sort(totales_com,'descend');
for indp = 1:CuantasCloud
    yplot = TIMELINE.words(index_com(indsorted(indp))).trend(desde:hasta);
    yplot = yplot/max(yplot);
    plot(TIMELINE.words(1).years(desde:hasta),yplot); hold on;
end

box off
set(gca,'xcolor','none','fontsize',letsize)
set(gca,'ytick',[],'fontsize',letsize)
ylabel('Trends','fontsize',letsize)
ylim([-.238 1])
yticks([0 1])
set(gca,'ytick',[])
set(gca,'xtick',[])

% Panel B. Bottom
ZZ = [TIMELINE.words(index_com).smoothed]';
ZZ = ZZ(:,desde:hasta);
set(gcf,'currentaxes',handles(4));
hold all
for indp = 1:CuantasCloud
    plot(TIMELINE.words(1).years(desde:hasta),ZZ(indsorted(indp),:)/max(abs(ZZ(indp,:))));
end

box off
ylim([-1.1 1.3])
ylabel('Oscilations','fontsize',letsize)
yticks(0)
set(gca,'ytick',[],'fontsize',letsize)
linkaxes(handles([3 4]),'x')
xlim([1750 1999])
xlabel('Year','fontsize',letsize)


% Panel C

CuantasCloud = 15;
load(['./' language '/clusters_osc'])
load(['./' language '/TIMELINE'])
word = 'copper';
indword = find(ismember({TIMELINE.words.word},word));
index = find(T==T(indword));

load('./en2020/cluster_osc_subcom')
palabras = {TIMELINE.words(index_com).word};

indcolor = find(ismember({TIMELINE.words(index_com).word},word));
wordcolor = zeros(length(index_com),3);
wordcolor(indcolor,:) = [1 0 0];
PAL = struct('palabra',palabras,'count',num2cell([TIMELINE.words(index_com).tot]));
wordcloud({PAL.palabra},[PAL.count],'position',[.70 .42 .28 .16],...
    'MaxDisplayWords',CuantasCloud,'Color',wordcolor);

set(gcf,'currentaxes',handles(5));
for indp = 1:length(index_com)
    yplot = TIMELINE.words(index_com(indp)).trend(desde:hasta);
    yplot = yplot/max(yplot);
    plot(TIMELINE.words(1).years(desde:hasta),yplot); hold on;
end
box off
set(gca,'xcolor','none','fontsize',letsize)
set(gca,'ytick',[])
ylabel('Trends','fontsize',letsize)
ylim([-.07 1])
yticks([0 1])
set(gca,'ytick',[])
set(gca,'xtick',[])

%Panel C. Bottom
ZZ = [TIMELINE.words(index_com).smoothed]';
ZZ = ZZ(:,desde:hasta);
set(gcf,'currentaxes',handles(6));
hold all

for indp = 1:length(index_com)
    plot(TIMELINE.words(1).years(desde:hasta),ZZ(indp,:)/max(abs(ZZ(indp,:)))); hold on;
end

box off
ylim([-1.1 1.33])
ylabel('Oscilations')
yticks(0)
set(gca,'ytick',[],'fontsize',letsize)
linkaxes(handles([5 6]),'x')
xlim([1750 1999])
xlabel('Year','fontsize',letsize)


letras='abc';
for indp=1:length(letras)
    set(gcf,'currentaxes',handles(2*indp-1))
    if indp==4
        posx = min(xlim)+range(xlim)*(.05-.04);
        posy = max(ylim)-range(ylim)*(.05-.03);
    else
        posx = min(xlim)-range(xlim)*.08;
        posy = max(ylim)-range(ylim)*.07;
    end
    text(posx,posy,letras(indp),'fontsize',20)   
end


%% figura 2

clear all
close all
figure(2);clf
set(gcf,'color','w')
set(gcf,'position',[680 298 1366 768])
bordeizq = 0.05;
ancho = 0.45;
bordeinf = 0.07;
alto = 0.9;
sepv = 0.08;
seph = 0.025;
letsize = 12;

clear handles
handles(1)=axes('position',[bordeizq            bordeinf+.05 ancho alto-.1]);
handles(2)=axes('position',[bordeizq+seph+ancho bordeinf ancho alto]);

language='en2020';
load(['./' language '/AjustesIndo13'])
medidas={'corrtotal';'corrosc';'disttotal';'distosc';'SumaCODT'};
load(['./' language '/TIMELINE'])

% PANEL A
set(gcf,'currentaxes',handles(1));
hold all


ylabel('\tau (years)','fontsize',letsize)
xlabel('R (years^{-1})','fontsize',letsize)
xx=.2:.01:1;
plot(xx,4./xx,'k','linewidth',4)
plot(xx,8./(27*xx),'k','linewidth',4)
xlim([.21 0.98])
ylim([0.1 11.5])
xticks(0.3:0.1:1)
posmaximos = [Todo.posmaximos];
Razon = [Todo.razon];
maximos = [Todo.maximos];

indp = 5;
x = rs([posmaximos.(medidas{indp})]);
x = x+(rand(size(x))-1/2).*0.01;
y = taus([posmaximos.(medidas{indp})]);
C = [Razon.(medidas{indp})];
y = y+(rand(size(y))-1/2).*0.1;
indice = [Razon.(medidas{indp})]>0.75 & [Razon.(medidas{indp})]<1.25;
Indices = find(indice);
set(gca,'fontsize',letsize)

h = scatter(x(indice),y(indice),10,'k','o','filled');
h.CData = [.5 .5 .5];
colormap bone

imagedir='./Images/';


axes('position',[.13 .25 .09 .19]);
I=importdata([imagedir 'Imagen1.png']);
h=image(I.cdata);
set(gca,'xcolor','none','ycolor','none','color','none')
set(h, 'AlphaData', I.alpha);

%HOPF
axes('position',[.28 .7 .09 .19]);
I=importdata([imagedir 'Imagen3.png']);
h=image(I.cdata);
set(gca,'xcolor','none','ycolor','none','color','none')
set(h, 'AlphaData', I.alpha);

%REALES
axes('position',[.028 .023 .09 .19]);
I=importdata([imagedir 'Imagen2.png']);
h=image(I.cdata);
set(gca,'xcolor','none','ycolor','none','color','none')
set(h, 'AlphaData', I.alpha);

set(gcf,'currentaxes',handles(2));
I = importdata([imagedir 'FigCoordPol.png']);
h = imshow(I.cdata);
set(gca,'xcolor','none','ycolor','none','color','none')
set(h, 'AlphaData', I.alpha);

letras='ab';
for indp = 1:length(letras)
    set(gcf,'currentaxes',handles(indp))
    if indp==1
        posx = min(xlim)-range(xlim)*.08;
        posy = max(ylim)-range(ylim)*.07;
    else
        posx = min(xlim)-range(xlim)*.01;
        posy = max(ylim)-range(ylim)*0.92;
    end
    text(posx,posy,letras(indp),'fontsize',20)   
end



%% figura 3



clear all
close all
figure(3);clf
set(gcf,'color','w')
set(gcf,'position',[680 298 1366 768])
bordeizq = 0.05;
ancho = 0.2;
bordeinf = 0.08;
alto = 0.91;
sepv = 0;
seph = 0.05;
seph2 = 0.1;
letsize = 12;
clear handles
handles(1)=axes('position',[bordeizq               bordeinf+alto/2 ancho   alto/2]);
handles(2)=axes('position',[bordeizq               bordeinf        ancho   alto/2]);
handles(3)=axes('position',[bordeizq+seph+ancho    bordeinf        2*ancho alto  ]);
handles(4)=axes('position',[bordeizq+seph2+3*ancho alto/5          ancho   3*alto/4  ]);

fprintf('------\n')
pp = annotation('rectangle',[bordeizq,bordeinf+alto/2,3*ancho+seph+0.01,...
    alto/2],'color','b');
set(pp,'FaceColor','b','FaceAlpha',0.1,'linestyle','none');


% PANEL A top
set(gcf,'currentaxes',handles(1));
hold all
bases = {'en2020','es2020','fr2020','al2020','it2020'};
Lenguajes = {'English','Spanish','French','German','Italian'};
for len = 1:length(bases)
    language = bases{len};
    load(['./' language '/clusters_osc'])
    cuantascom = length(megustan);
    tamanios = zeros(1,cuantascom);
    for indcom = 1:cuantascom
        tamanios(indcom) = sum(T==megustan(indcom));
    end
    tamsort = sort(tamanios,'descend');
    plot(1:cuantascom,tamsort)
end
etiqx = [10 100];
xticks(etiqx)
etiqy = [10 100 1000];
yticks(etiqy)
ylim([10 2000])
yticklabels({'10^1' '10^2' '10^3'})
ylabel('Size of Topic Communities','fontsize',letsize)
set(gca,'xscale','log','fontsize',letsize)
set(gca,'yscale','log','fontsize',letsize)

% PANEL A bottom
set(gcf,'currentaxes',handles(2));
hold all

for len = 1:length(bases)
    language = bases{len};
    load(['./' language '/clusters_trend'])
    cuantascom = length(megustan);
    tamanios = zeros(1,cuantascom);
    for indcom = 1:cuantascom
        tamanios(indcom) = sum(T==megustan(indcom));
    end
    tamsort = sort(tamanios,'descend');
    plot(1:cuantascom,tamsort)
end
etiqx = [1 10 100];
xticks(etiqx)
xticklabels({'10^0' '10^1' '10^2'})
etiqy = [10 100 1000];
ylabel('Size of Keyword Communities','fontsize',letsize)
yticks(etiqy)
yticklabels({'10^1' '10^2' '10^3'})
ylim([10 2000])
linkaxes(handles([1 2]),'x')
xlim([1 150])
xlabel('Rank')
set(gca,'yscale','log','fontsize',letsize)
set(gca,'xscale','log','fontsize',letsize)
legend(Lenguajes)

% Panel B and C
set(gcf,'currentaxes',handles(3));hold all

for len = 1:length(bases)
    language = bases{len};
    load(sprintf('./%s/clusters_trend',language))
    load(sprintf('./%s/Stats/stat_true',language))
    ParExp = ParOrden;
    load(sprintf('./%s/Stats/stat_shuffle',language))
    ParShuf = ParOrden;
    
    load(sprintf('./%s/TIMELINE',language))
    clusters = T;
    carpetasave = sprintf('./%s/Stats/',language);
    
    desde = 50;
    hasta = 300;
    mincom = 9;
    maxcom = 701;
    
    tamanios = zeros(length(megustan),1);
    for i = 1:length(megustan)
        tamanios(i) = length(find(clusters==megustan(i)));
    end
    indcoms = find(tamanios>mincom & tamanios<maxcom);
    ParExp = ParExp(indcoms,:);
    ParShuf = ParShuf(indcoms,:);
    
    mediaexp = mean(ParExp,2)';
    mediashuf = mean(ParShuf,2)';
    
    load(sprintf('%sStat_Acople',carpetasave))

    indm = [M.Ncom]>mincom & [M.Ncom]<maxcom;
    parpol = [M(indm).ParSinTr];
    parpol = parpol(11:end-20,:);
    mediastr = mean(parpol);
    
    parpol = [M(indm).ParConTr];
    parpol = parpol(11:end-20,:);
    mediactr = mean(parpol);


    parpol = [M(indm).ParAc];
    parpol = parpol(11:end-20,:);
    mediaac = mean(parpol);
    
    medias = [mean(mediaexp) nanmean(mediashuf) mean(mediastr) mean(mediactr) mean(mediaac)];
    err = [nansem(mediaexp) nansem(mediashuf) nansem(mediastr) nansem(mediactr) nansem(mediaac)];
    
    errorbar(0.7+0.1*len:4.7+0.1*len,medias,err...
        ,'linestyle','none','linewidth',3,'capsize',0)
    
    [~,p] = ttest(mediaac,mediaexp);
    if p>0.05
    else
        fprintf('%s Trend ConAc/Exp NO\n',Lenguajes{len})
    end

end

row1 = {'Exp','Shuffle','  No Trend','    Trend','  Trend'};
row2 = {'','','No Coupling','No Coupling','Coupling'};
xlim([0.7 5.2])
xpos = 1:length(row1);
xticks(xpos)
labelArray = [row1; row2];
tickLabels = (sprintf('%s\\newline%s\n', labelArray{:}));
set(gca,'xtick',xpos,'xticklabel',tickLabels,'fontsize',letsize)


colores = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125;...
    0.494 0.184 0.556; 0.466 0.674 0.188];
fprintf('------\n')

% Shaded Area in B and C
for len = 1:length(bases)
    language = bases{len};
    load(sprintf('./%s/clusters_osc',language))
    load(sprintf('./%s/TIMELINE',language))
    clusters = T;
    carpetasave = sprintf('./%s/Stats/',language);
    Osc = [TIMELINE.words.smoothed]';
    
    mincom = 9;
    maxcom = 10001;
    
    tamanios = zeros(length(megustan),1);
    for i = 1:length(megustan)
        tamanios(i) = length(find(clusters==megustan(i)));
    end
    indcoms = find(tamanios>mincom & tamanios<maxcom);
    
    ParExp = zeros(length(indcoms),length(desde:hasta));
    for i = 1:length(indcoms)
        comexp = find(clusters==megustan(indcoms(i)));
        oscy = Osc(comexp,desde:hasta);
        ParExp(i,:) = PolarFun.ParOrden(oscy);
    end
    mediaexp = mean(ParExp,2)';
    
    
    load(sprintf('%sStat_Osc',carpetasave))
    indm = [M.Ncom]>mincom;
    parpol = [M(indm).ParAc];
    parpol = parpol(11:end-20,:);
    mediaac = mean(parpol);

    medias = [mean(mediaexp) mean(mediaac)];
    err = [nansem(mediaexp) nansem(mediaac)];
    
    errorbar([0.7+0.1*len 4.7+0.1*len],medias,err,'color',colores(len,:)...
        ,'linestyle','none','linewidth',3,'capsize',0)
    
    [~,p] = ttest(mediaac,mediaexp);
    if p>0.05
    else
        fprintf('%s Osc ConAc/Exp NO\n',Lenguajes{len})
    end
end
ylim([0.15 .62])
grid on
line([2.5 2.5],ylim,'color','k')
ylabel('Phase Coherence \rho','fontsize',letsize)

 
set(gcf,'currentaxes',handles(4));hold all
colores  = [75 156 211;50 120 200; 15 82 186; 0 35 102]/256;
xlabel('\lambda_w','fontsize',letsize)
ylabel('\lambda_s','fontsize',letsize)

M = struct([]);
M(1).legend = 'N';
M(2).legend = 'N/2';
M(3).legend = 'N/4';
M(4).legend = 'N/8';

M(1).lw = 0.0675;
M(2).lw = 0.07;
M(3).lw = 0.075;
M(4).lw = 0.075;

M(1).ls = 0.1125;
M(2).ls = 0.1175;
M(3).ls = 0.125;
M(4).ls = 0.1325;

M(1).elw = 0.0075;
M(2).elw = 0.005;
M(3).elw = 0.005;
M(4).elw = 0.005;

M(1).els = 0.0025;
M(2).els = 0.0025;
M(3).els = 0.0025;
M(4).els = 0.0025;

pp = [];
for i = 1:4
    pp(end+1) = scatter(M(i).lw,M(i).ls,[],colores(5-i,:),'filled');
    pp(end+1) = errorbar(M(i).lw,M(i).ls,M(i).els,M(i).els,M(i).elw,M(i).elw,...
        'color',colores(5-i,:));
end
legend(pp(1:2:7),{M.legend},'location','northwest')
set(gca,'fontsize',letsize)
yticks(0.09:0.005:0.13)
xlim([0.064 0.082])
ylim([0.088 0.135])


letras='abcd';
for indp = 1:length(letras)
    if indp==1
        set(gcf,'currentaxes',handles(indp))
        posx = min(xlim)*0.7;
        posy = max(ylim)*.75;
    elseif indp == 2
        set(gcf,'currentaxes',handles(3))
        posx = 0.28;
        posy = max(ylim)*.98;
    elseif indp == 3
        set(gcf,'currentaxes',handles(indp))
        posx = 2.3;
        posy = max(ylim)*.98;  
    else
        set(gcf,'currentaxes',handles(4))
        posx = 0.0625;
        posy = max(ylim)*.99;
    end
    text(posx,posy,letras(indp),'fontsize',20)   
end





