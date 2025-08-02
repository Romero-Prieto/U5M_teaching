clear
pATh                  = "/Users/lshjr3/Documents/SARMAAN/";
lISt                  = {'bASe','IGMEnigeria','DHSnigeria','MICSnigeria'};
for i = 1:numel(lISt)
    options    = detectImportOptions(char(pATh + string(lISt{i}) + ".csv"));
    for j = 1:numel(options.VariableTypes)
        if isequal(options.VariableTypes{j},'char')
            options.VariableTypes{j} = 'categorical';
        end
        
        if isequal(options.VariableTypes{j},'datetime') && ~isequal(options.VariableNames{j},'sTArt') && ~isequal(options.VariableNames{j},'eNd') && ~isequal(options.VariableNames{j},'sUBmiSSion')
            options.VariableOptions(1,j).InputFormat = 'dd/MM/yyyy';
        end
    end

    sET        = readtable(char(pATh + string(lISt{i}) + ".csv"),options);
    assignin('base',lISt{i},sET);
    clear options j sET
end

clear lISt i
bASe.interview                      = datetime(year(bASe.sUBmiSSion),month(bASe.sUBmiSSion),day(bASe.sUBmiSSion));
bASe.B_min                          = datetime(year(bASe.B),month(bASe.B),1);
bASe.B_max                          = datetime(year(bASe.B),month(bASe.B) + 1,1) - 1;
sEL                                 = ~isnat(bASe.B_min);
bASe.B_max(sEL)                     = max(min(bASe.B_max(sEL),bASe.interview(sEL)),bASe.B_min(sEL));
bASe.DOB_min                        = datetime(year(bASe.DOB),month(bASe.DOB),1);
bASe.DOB_max                        = datetime(year(bASe.DOB),month(bASe.DOB) + 1,1) - 1;
bASe.age                            = max(floor((years(bASe.interview - bASe.DOB_min) + years(bASe.interview - bASe.DOB_max))/2),10);
bASe.GO                             = 1 + min(floor((bASe.age - 10)/5),8);
bASe.cluster(bASe.cluster == "C-1") = "Kano";
bASe.cluster(bASe.cluster == "C-2") = "Kano";
bASe.cluster(bASe.cluster == "C-3") = "Kaduna";
save(char(pATh + "Results/bASe.mat"),'bASe','IGMEnigeria','DHSnigeria','MICSnigeria');


clear
pATh                  = "/Users/lshjr3/Documents/SARMAAN/";
RESolUTioN            = 300;
load(char(pATh + "Results/bASe.mat"),'bASe');
date                  = max(bASe.interview);
Ts                    = {datetime([year(date) - 5 year(date)]',month(date),day(date)),datetime([year(date) - 4 year(date)]' - 6,month(date),day(date)),datetime([year(date) - 2 year(date)]' - 4,month(date),day(date)),datetime([year(date) - 2 year(date)]' - 2,month(date),day(date)),datetime([year(date) - 1 year(date)]' - 1,month(date),day(date)),datetime([year(date) - 1 year(date)]',month(date),day(date))};
x{1}                  = [(0:7:28)/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
x{2}                  = [string(0);string([(7:7:28)';(2:1:11)';(12:3:24)';(36:12:60)']) + char([kron('d',ones(4,1));kron('m',ones(18,1))])];
n                     = x{1}(2:end) - x{1}(1:end - 1);
models                = {'$\textit{15-49}$','$\textit{10-55}$'};
cLUsTEr               = {'Kano','Kaduna'};
MO                    = {bASe.age >= 15 & bASe.age < 50,bASe.age >= 0};
R                     = 250;
aGEs                  = [5 16 23];

data                  = bASe(bASe.k == 1,{'cluster','moTHer','K'});
W                     = NaN(size(data,1),R + 1);
rng(0);    
for j = 1:numel(cLUsTEr)
    sEL             = (data.cluster == cLUsTEr{j});
    S               = [(1:sum(sEL))',unidrnd(sum(sEL),sum(sEL),R)];
    for r = 1:R + 1
        temp     = tabulate([S(:,1);S(:,r)]);
        W(sEL,r) = temp(:,2) - 1;
        clear temp
        clc;
        r/(R + 1)
    end
    clear sEL S r  
end
W                     = W(repelem((1:size(data,1))',data.K),:);
clear data j
for j = 1:numel(cLUsTEr)
    sEL             = (bASe.survival == 'dead' | bASe.survival == 'alive') & bASe.cluster == cLUsTEr{j};
    pOP{j}          = char("$\textrm{" + cLUsTEr{j} + ", B = " + sum(sEL) + "}$");
    sEL             = bASe.k == 1 & bASe.cluster == cLUsTEr{j};
    pOPw{j}         = char("$\textrm{" + cLUsTEr{j} + ", W = " + sum(sEL) + "}$");
    sEL             = bASe.survival == 'dead' & bASe.cluster == cLUsTEr{j};
    pOPd{j}         = char("$\textrm{" + cLUsTEr{j} + ", D = " + sum(sEL) + "}$");
end

for i = 1:2
    for j = 1:numel(models)
        h          = j + (i - 1)*numel(models);
        dATa{h}    = {bASe.cluster == cLUsTEr{i} & MO{j},W};
        for k = 1:numel(Ts)
            date      = eXAcTTime(Ts{k});
            date      = char("$\textrm{" + cLUsTEr{i} + ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2))) + "}$");
            pOPs{k,h} = {date;''};
            clear date
        end
    end
end
save(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','Ts','-v7.3','-nocompression');


sex              = [2 1];
for i = 1:numel(dATa)
    pOPs{numel(Ts) + 1,i}    = pOPs{1,i};
    pOPs{numel(Ts) + 1,i}{1} = char(pOPs{1,i}{1} + " (female)");
    pOPs{numel(Ts) + 2,i}    = pOPs{1,i};
    pOPs{numel(Ts) + 2,i}{1} = char(pOPs{1,i}{1} + " (male)");
end

p                     = rand(size(bASe,1),R + 1);
d                     = (p.*bASe.D_min + (1 - p).*bASe.D_max)/365.25;
p                     = rand(size(bASe,1),R + 1);
bR                    = (bASe.birth ~= 'livebirth');
B                     = p.*eXAcTTime(bASe.B_min) + (1 - p).*eXAcTTime(bASe.B_max);
B(bR,:)               = NaN;
bR                    = (bASe.birth ~= 'stillbirth');
sB                    = p.*eXAcTTime(bASe.B_min) + (1 - p).*eXAcTTime(bASe.B_max);
sB(bR,:)              = NaN;

interview             = eXAcTTime(bASe.interview);
D                     = B + d;
O                     = min(D,interview);
p                     = rand(size(bASe,1),R + 1);
ageS                  = p.*years(bASe.interview - bASe.DOB_min) + (1 - p).*years(bASe.interview - bASe.DOB_max);
ageB                  = ageS - (interview - B);
ageSB                 = ageS - (interview - sB);
sEL                   = (bASe.survival == 'dead' | bASe.survival == 'alive');
        
for i = 1:numel(dATa)
    s = dATa{i}{1};
    w = dATa{i}{2};
    for h = 1:size(pOPs,1)
        if h <= numel(Ts)
            TS            = Ts{h};
            T             = eXAcTTime(TS);
            sEX           = ones(size(bASe,1),1) == 1;
        else
            TS            = Ts{1};
            T             = eXAcTTime(TS);
            sEX           = (bASe.sex == sex(h - numel(Ts)));
        end
        sH              = dATa{i}{1} & sEX;

        for j = 1:numel(x{1}) - 1
            a             = max(min(B(sH,:) + x{1}(j),min(O(sH,:),T(2))),T(1));
            o             = max(min(B(sH,:) + x{1}(j + 1),min(O(sH,:),T(2))),T(1)); 
            exposure(j,:) = sum((o - a).*(w(sH,:).*sEL(sH)));
            events(j,:)   = sum((d(sH,:) >= x{1}(j) & d(sH,:) < x{1}(j + 1) & D(sH,:) >= T(1) & D(sH,:) < T(2)).*(w(sH,:).*sEL(sH,:)));
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end
        
        stillbirths     = sum((bASe.FlagBI(s) ~= 1 & bASe.gestation(s) >= 28 & sB(s,:) >= T(1) & sB(s,:) < T(2)).*w(s,:));
        births          = sum((B(s,:) >= T(1) & B(s,:) < T(2)).*w(s,:));
        stillbirths     = stillbirths./(births + stillbirths);        
        SRB             = sum((B(s,:) >= T(1) & B(s,:) < T(2) & bASe.sex(s) == 1).*w(s,:))./sum((B(s,:) >= T(1) & B(s,:) < T(2) & bASe.sex(s) == 2).*w(s,:));
        m               = events./exposure;
        q               = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        
        TaBle.m{h,i}    = m;
        TaBle.q{h,i}    = q;
        TaBle.s{h,i}    = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBle.SNR{h,i}  = [TaBle.s{h,i}(1,:);TaBle.q{h,i}(2,:);TaBle.q{h,i}(5,:);TaBle.s{h,i}(1,:)./TaBle.q{h,i}(5,:)/1000;TaBle.q{h,i}(2,:)./TaBle.q{h,i}(5,:)/1000];
        TaBle.SRB{h,i}  = SRB;
        TaBle.tAU{h,i}  = TS;
        TaBle.BTHW{h,i} = births./sum(w(s,:));
        clear T TS sEX events exposure q m births stillbirths SRB sH
    end
    
    A                    = [max(min(bASe.age(dATa{i}{1})),15),min(max(bASe.age(dATa{i}{1})),49) + 1];
    Ax                   = (min(bASe.age(dATa{i}{1})):max(bASe.age(dATa{i}{1})))';
    Ax                   = 10:55;
    for j = 1:numel(Ax)
        a               = max(min(ageS(s,:) - 3,Ax(j) + 1),Ax(j));
        o               = max(min(ageS(s,:),Ax(j) + 1),Ax(j));
        exposure        = sum((bASe.k(s) == 1).*(o - a).*w(s,:));
        events          = sum((interview(s) - B(s,:) <= 3 & interview(s) - B(s,:) > 0).*(ageB(s,:) >= Ax(j) & ageB(s,:) < Ax(j) + 1).*w(s,:));
        TFR(j,:)        = events./max(exposure,eps);
        clear events exposure a o
        clc;
        j/numel(Ax)
    end
    
    for j = 1:numel(Ax)
        a               = max(min(ageS(s,:) - 3,Ax(j) + 1),Ax(j));
        o               = max(min(ageS(s,:),Ax(j) + 1),Ax(j));
        exposure        = sum((bASe.k(s) == 1).*(o - a).*w(s,:));
        events          = sum((interview(s) - sB(s,:) <= 3 & interview(s) - sB(s,:) > 0).*(ageSB(s,:) >= Ax(j) & ageSB(s,:) < Ax(j) + 1).*w(s,:));
        TSR(j,:)        = events./max(exposure,eps);
        clear events exposure a o
        clc;
        j/numel(Ax)
    end
    
    TaBle.F{1,i}         = TFR;
    TaBle.TFR{1,i}       = sum(TFR(Ax >= A(1) & Ax < A(2),:),1);
    TaBle.TSR{1,i}       = sum(TSR(Ax >= A(1) & Ax < A(2),:),1);
    TaBle.sRB{1,i}       = sum((bASe.birth(s) == 'livebirth' & bASe.sex(s) == 1).*w(s,:))./sum((bASe.birth(s) == 'livebirth' & bASe.sex(s) == 2).*w(s,:)); 
    TaBle.parity{1,i}    = sum((bASe.birth(s) == 'livebirth').*w(s,:))./sum((bASe.k(s) == 1).*w(s,:));
    TaBle.childless{1,i} = sum((bASe.mother(s) ~= 1 & bASe.k(s) == 1).*w(s,:))./sum((bASe.k(s) == 1).*w(s,:))*100;
    clear TFR TSR w s A
end
clear b ageB ageSB ageS bR h i j k r

for i = 1:numel(dATa)
    TaBle.sEX{1,i} = log(1 - TaBle.q{numel(Ts) + 2,i})./log(1 - TaBle.q{numel(Ts) + 1,i});
end

lABelS                = {'Sex','Place of residence','Education','Electricity'};
lABelSd               = {{'female' 'male'} {'urban' 'rural'} {'less than complete secondary' 'complete secondary or more'} {'access' 'no access'}}; 
pOPh{1}               = {'All women'};
for i = 1:numel(lABelS)
    pOPh{end + 1,1} = {char(lABelS{i} + ": " + lABelSd{i}{1})};
    for j = 2:numel(lABelSd{i})
        pOPh{end + 1,1} = {char(lABelSd{i}{j})};
    end
end

setH                  = ones(size(bASe.sex));
setH(:,end + 1)       = (bASe.sex == 2);
setH(:,end + 1)       = (bASe.sex == 1);
setH(:,end + 1)       = (bASe.UR == 1);
setH(:,end + 1)       = (bASe.UR == 2);
setH(:,end + 1)       = (bASe.Education == 1 | bASe.Education == 2);
setH(:,end + 1)       = (bASe.Education == 3);
setH(:,end + 1)       = (bASe.Electricity == 2);
setH(:,end + 1)       = (bASe.Electricity == 1);
T                     = eXAcTTime(Ts{1});

for i = 1:numel(dATa)
    for h = 1:size(setH,2)
        s               = dATa{i}{1} & setH(:,h); 
        wH              = dATa{i}{2}(s,:);
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(s,:) + x{1}(j),min(O(s,:),T(2))),T(1));
            o             = max(min(B(s,:) + x{1}(j + 1),min(O(s,:),T(2))),T(1)); 
            exposure(j,:) = sum((o - a).*wH.*sEL(s));
            events(j,:)   = sum((d(s,:) >= x{1}(j) & d(s,:) < x{1}(j + 1) & D(s,:) >= T(1) & D(s,:) < T(2)).*wH.*sEL(s));
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end
        
        stillbirths     = sum((bASe.FlagBI(s) ~= 1 & bASe.gestation(s) >= 28 & sB(s,:) >= T(1) & sB(s,:) < T(2)).*wH);
        births          = sum((B(s,:) >= T(1) & B(s,:) < T(2)).*wH);
        stillbirths     = stillbirths./(births + stillbirths);        
        SRB             = sum((B(s,:) >= T(1) & B(s,:) < T(2) & bASe.sex(s) == 1).*wH)./sum((B(s,:) >= T(1) & B(s,:) < T(2) & bASe.sex(s) == 2).*wH);
        m               = events./exposure;
        q               = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        
        TaBle.m{h,i}    = m;
        TaBlE.q{h,i}    = q;
        TaBlE.s{h,i}    = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBlE.tAU{h,i}  = T;
        clear events exposure q m births stillbirths SRB wH s
    end
end
clear T i sH interview sH sEL p d D B O j h sB



BrassTimeD            = datetime([kron((2015:2024)',ones(12,1)),kron(ones(2024 - 2015 + 1,1),(1:12)'),15*ones((2024 - 2015 + 1)*12,1)]);
BrassTime             = (2015.0 + 1/24:1/12:2025.0 - 1/24)';
BrassP                = [BrassTime < 2021.0,BrassTime >= 2021.0];
BrassP                = BrassP./sum(BrassP);
CD                    = {'North','South','East','West'};
family                = 1;

for i = 1:2
    H             = i + numel(dATa);
    dATaB{i}      = {bASe.cluster == cLUsTEr{i} & bASe.k == 1 & MO{1},W};
    for k = 1:size(BrassP,2)
        date       = string(sprintf('%0.1f',BrassP(:,k)'*BrassTime - 2.5)) + "-" + string(sprintf('%0.1f',BrassP(:,k)'*BrassTime + 2.5));
        pOPs{k,H}  = {char("$\textrm{" + cLUsTEr{i} + " SBH (Brass: " + string(CD{family}) + "), " + date + "}$");''};
        clear date
    end
    pOPs{k + 1,H} = {char("$\textrm{" + cLUsTEr{i} + " SBH (Brass: " + string(CD{family}) + ")}$");''};
    clear h H
end

for i = 1:numel(dATaB)
    s                       = dATaB{i}{1};
    w                       = dATaB{i}{2}(s,:);
        
    T                       = eXAcTTime(bASe.interview(s));
    T                       = (T'*w)./sum(w);
    xx                      = bASe.age(s) >= 15:5:45 & bASe.age(s) < 20:5:50;
    W                       = xx'*w;
    B                       = xx'*((bASe.sons(s) + bASe.daughters(s)).*w);
    D                       = xx'*((bASe.sonsD(s) + bASe.daughtersD(s)).*w);
    
    h                       = i + numel(dATa);
    TaBle.sRB{1,h}          = (bASe.sons(s)'*w)./(bASe.daughters(s)'*w);
    children                = bASe.sons(s) + bASe.daughters(s);
    TaBle.parity{1,h}       = sum(children.*w)./sum(w);        
    TaBle.childless{1,h}    = ((children == 0)'*w)./((children >= 0)'*w)*100;
    
    temp                    = BraSsTruSseLl(W,B,D,T);
    for k = 1:numel(CD)
        for j = 1:2
            BraSs_s{k,i}{j} = LinInterPol(flip(temp{k,end}),flip(temp{k,j + 1}),BrassTime);
            BraSs_r{k,i}{j} = flip(temp{k,j + 1});
        end
        BraSs_r{k,i}{j + 1} = flip(temp{k,end});
    end
    
    temp                    = temp(family,:);
    for j = 1:numel(temp) - 1
        q                   = LinInterPol(flip(temp{end}),flip(temp{j}),BrassTime);
        for k = 1:size(BrassP,2)
            if isequal(j,1)
                TaBle.q{k,h} = NaN(max(aGEs),R + 1);
            end
            TaBle.q{k,h}(aGEs(j),:) = BrassP(:,k)'*q;
            TaBle.tAU{k,h}          = BrassP(:,k)'*BrassTime;
        end
        TaBle.q{k + 1,h}{j} = q; 
    end
    TaBle.q{k + 1,h}{j + 1} = BrassTimeD;
    TaBle.tAU{k + 1,h}      = BrassTimeD;
    clear s w W T xx W B D h children temp q
end




load(char(pATh + "Results/bASe.mat"),'DHSnigeria');
mIn                   = min(DHSnigeria.interview);
mAx                   = max(DHSnigeria.interview);
mAx                   = mean([mIn mAx]);
mIn                   = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
dATe{1}               = [mIn mAx];
dATe{2}               = eXAcTTime(dATe{1});
Tdhs                  = dATe{2};
models                = {'$\textit{National}$','$\textit{North West}$','$\textit{Kano}$','$\textit{Kaduna}$','$\textit{Sokoto}$'};
dATaDHS{1}            = {DHSnigeria.age >= 0,DHSnigeria.W};
dATaDHS{2}            = {DHSnigeria.Region == 3,DHSnigeria.W};
dATaDHS{3}            = {DHSnigeria.State == 31 | DHSnigeria.State == 32,DHSnigeria.W};
dATaDHS{4}            = {DHSnigeria.State == 29 | DHSnigeria.State == 30,DHSnigeria.W};
dATaDHS{5}            = {DHSnigeria.State == 37 | DHSnigeria.State == 38,DHSnigeria.W};

for i = 1:numel(dATaDHS)
    data      = "$\textrm{DHS VII, " + string(sprintf('%0.1f',dATe{2}(1))) + "-" + string(sprintf('%0.1f',dATe{2}(2))) + "}$";
    h         = i + numel(dATa) + numel(dATaB);
    pOPs{1,h} = {char(data);models{i}};
    pOPs{2,h} = {char(data{1}(1:end - 2) + " (female)}$");models{i}};
    pOPs{3,h} = {char(data{1}(1:end - 2) + " (male)}$");models{i}};
end
clear mIn mAx models data

rng(0);
s                     = DHSnigeria(DHSnigeria.k == 1, {'cluster','K','woman','Women','iNDeX'});
w                     = rand(size(s,1),R);
w                     = [(1:size(s,1))',ceil(s.Women.*w) + s.iNDeX];
Wdhs                  = NaN(size(DHSnigeria,1),R + 1);
for r = 1:R + 1
    S         = tabulate([w(:,1);w(:,r)]);
    Wdhs(:,r) = repelem(S(:,2) - 1,s.K);
    clc;
    r/(R + 1)
end
clear r s S wDHS;
DHSnigeria.births      = zeros(size(DHSnigeria,1),1);  
DHSnigeria.stillbirths = zeros(size(DHSnigeria.births,1),1);
sET                    = find(~isnan(DHSnigeria.row));
for j = 1:numel(sET)
    cal = char(DHSnigeria.CAL(sET(j)));
    if ismember('T',cal) || ismember('B',cal)
        cal = cal(DHSnigeria.row(sET(j)) + 1:end);
        sT  = find(ismember(cal,'T'));
        for h = 1:numel(sT)
            dT = datetime(year(DHSnigeria.interview(sET(j))),month((DHSnigeria.interview(sET(j)))) - sT(h) + 1,day(DHSnigeria.interview(sET(j))),'Format','dd/MM/yyyy');
            gA = cal(sT(h):min(sT(h) + 6,end));
            if isequal(gA,'TPPPPPP') && dT >= dATe{1}(1) && dT  < dATe{1}(2)
               DHSnigeria.stillbirths(sET(j)) = DHSnigeria.stillbirths(sET(j)) + 1;
            end
            clear dT gA
        end
        sB  = find(ismember(cal,'B'));
        for h = 1:numel(sB)
            dB = datetime(year(DHSnigeria.interview(sET(j))),month((DHSnigeria.interview(sET(j)))) - sB(h) + 1,day(DHSnigeria.interview(sET(j))),'Format','dd/MM/yyyy');
            gA = cal(sB(h):min(sB(h) + 6,end));
            if dB >= dATe{1}(1) && dB  < dATe{1}(2)
               DHSnigeria.births(sET(j))      = DHSnigeria.births(sET(j)) + 1;
            end
            clear dB gA
        end
        clear sT sB
    end
    clear cal
    j/numel(sET)
end


for i = 1:numel(dATaDHS)
    dATaDHS{i}{2} = Wdhs.*dATaDHS{i}{2};
end
clear Wdhs
save(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','dATaDHS','Ts','-v7.3','-nocompression');


p                     = rand(size(DHSnigeria,1),R + 1);
d                     = (p.*DHSnigeria.D_min + (1 - p).*DHSnigeria.D_max)/365.25;
p                     = rand(size(DHSnigeria,1),R + 1);
B                     = datetime(datenum(DHSnigeria.B_min) + datenum(DHSnigeria.B_max - DHSnigeria.B_min).*p,'ConvertFrom','datenum');
B                     = eXAcTTime(B);
D                     = B + d;
O                     = min(D,Tdhs(2));

p                     = rand(size(DHSnigeria,1),R + 1);
dob                   = p.*eXAcTTime(DHSnigeria.DOB) + (1 - p).*eXAcTTime(datetime(year(DHSnigeria.DOB),month(DHSnigeria.DOB) + 1,day(DHSnigeria.DOB)));
ageS                  = dATe{2}(2) - dob;
ageB                  = B - dob;

for i = 1:numel(dATaDHS)
    h                    = i + numel(dATa) + numel(dATaB);
    s                    = dATaDHS{i}{1};
    w                    = dATaDHS{i}{2};
    for k = 1:3
        if k == 1
            sEX           = ones(size(DHSnigeria,1),1);
        else
            sEX           = (DHSnigeria.sex == sex(k - 1));
        end
        sH                   = dATaDHS{i}{1} & sEX;
        
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(sH,:) + x{1}(j),min(O(sH,:),Tdhs(2))),Tdhs(1));
            o             = max(min(B(sH,:) + x{1}(j + 1),min(O(sH,:),Tdhs(2))),Tdhs(1));
            exposure(j,:) = sum((o - a).*w(sH,:));
            events(j,:)   = sum((d(sH,:) >= x{1}(j) & d(sH,:) < x{1}(j + 1) & D(sH,:) >= Tdhs(1) & D(sH,:) < Tdhs(2)).*w(sH,:));
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end

        births               = sum(DHSnigeria.births(s).*w(s,:));
        stillbirths          = sum(DHSnigeria.stillbirths(s).*w(s,:));
        stillbirths          = stillbirths./(stillbirths + births);
        SRB                  = ((DHSnigeria.sex(s) == 1 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*w(s,:))./((DHSnigeria.sex(s) == 2 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*w(s,:));
        m                    = events./exposure;
        q                    = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];

        TaBle.m{k,h}         = m;
        TaBle.q{k,h}         = q;    
        TaBle.s{k,h}         = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBle.SNR{k,h}       = [TaBle.s{k,h}(1,:);TaBle.q{k,h}(2,:);TaBle.q{k,h}(5,:);TaBle.s{k,h}(1,:)./TaBle.q{k,h}(5,:)/1000;TaBle.q{k,h}(2,:)./TaBle.q{k,h}(5,:)/1000];
        TaBle.SRB{k,h}       = SRB;
        TaBle.tAU{k,h}       = dATe{1}';
        clear sEX events exposure births stillbirths m q SRB
    end
    
    A                    = [max(min(DHSnigeria.age(dATaDHS{i}{1})),15),min(max(DHSnigeria.age(dATaDHS{i}{1})),49)];
    Ax                   = (min(DHSnigeria.age(dATaDHS{i}{1})):max(DHSnigeria.age(dATaDHS{i}{1})))';
    Ax                   = 10:55;
    for j = 1:numel(Ax)
        a               = max(min(ageS(s,:) - 3,Ax(j) + 1),Ax(j));
        o               = max(min(ageS(s,:),Ax(j) + 1),Ax(j));
        exposure        = sum((DHSnigeria.k(s) == 1).*(o - a).*w(s,:));
        events          = sum((dATe{2}(2) - B(s,:) <= 3 & dATe{2}(2) - B(s,:) > 0).*(ageB(s,:) >= Ax(j) & ageB(s,:) < Ax(j) + 1).*w(s,:));
        TFR(j,:)        = events./max(exposure,eps);
        clear events exposure a o
        clc;
        j/numel(Ax)
    end

    TaBle.F{1,h}         = TFR;
    TaBle.TFR{1,h}       = sum(TFR(Ax >= A(1) & Ax < A(2),:),1);
    TaBle.sRB{1,h}       = sum((DHSnigeria.sex(s) == 1).*w(s,:))./sum((DHSnigeria.sex(s) == 2).*w(s,:));
    TaBle.parity{1,h}    = sum(~isnan(DHSnigeria.bidx(s)).*w(s,:))./sum((DHSnigeria.k(s) == 1).*w(s,:));
    TaBle.childless{1,h} = sum((DHSnigeria.mother(s) ~= 1 & DHSnigeria.k(s) == 1).*w(s,:))./sum((DHSnigeria.k(s) == 1).*w(s,:))*100;
    clear h j TFR w A
end

for i = 1:numel(dATaDHS)
    h              = i + numel(dATa) + numel(dATaB);
    TaBle.sEX{1,h} = log(1 - TaBle.q{3,h})./log(1 - TaBle.q{2,h});
end

setH             = ones(size(DHSnigeria.W));
setH(:,end + 1)  = (DHSnigeria.sex == 2);
setH(:,end + 1)  = (DHSnigeria.sex == 1);
setH(:,end + 1)  = (DHSnigeria.UR == 1);
setH(:,end + 1)  = (DHSnigeria.UR == 2);
setH(:,end + 1)  = (DHSnigeria.Education == 1 | DHSnigeria.Education == 2);
setH(:,end + 1)  = (DHSnigeria.Education == 3);
setH(:,end + 1)  = (DHSnigeria.Electricity == 2);
setH(:,end + 1)  = (DHSnigeria.Electricity == 1);

for i = 1:numel(dATaDHS)
    k                = i + numel(dATa) + numel(dATaB);
    for h = 1:size(setH,2)
        s               = dATaDHS{i}{1} & setH(:,h); 
        wH              = dATaDHS{i}{2}(s,:);
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(s,:) + x{1}(j),min(O(s,:),Tdhs(2))),Tdhs(1));
            o             = max(min(B(s,:) + x{1}(j + 1),min(O(s,:),Tdhs(2))),Tdhs(1));
            exposure(j,:) = sum((o - a).*wH);
            events(j,:)   = sum((d(s,:) >= x{1}(j) & d(s,:) < x{1}(j + 1) & D(s,:) >= Tdhs(1) & D(s,:) < Tdhs(2)).*wH);
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end
        
        births          = sum(DHSnigeria.births(s).*wH);
        stillbirths     = sum(DHSnigeria.stillbirths(s).*wH);
        stillbirths     = stillbirths./(births + stillbirths);        
        SRB             = ((DHSnigeria.sex(s) == 1 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*wH)./((DHSnigeria.sex(s) == 2 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*wH);
        m               = events./exposure;
        q               = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        
        TaBle.m{h,k}    = m;
        TaBlE.q{h,k}    = q;
        TaBlE.s{h,k}    = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBlE.tAU{h,k}  = Tdhs;
        clear events exposure q m births stillbirths SRB s wH
    end
    clear h j k
end
clear p d b B D O dob ageS ageB sET i



load(char(pATh + "Results/bASe.mat"),'MICSnigeria');
mIn                   = min(MICSnigeria.interview);
mAx                   = max(MICSnigeria.interview);
mAx                   = mean([mIn mAx]);
mIn                   = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
dATe{1}               = [mIn mAx];
dATe{2}               = eXAcTTime(dATe{1});
Tmics                 = dATe{2};
models                = {'$\textit{National}$','$\textit{North West}$','$\textit{Kano}$','$\textit{Kaduna}$','$\textit{Sokoto}$'};
dATaMICS{1}           = {MICSnigeria.age >= 0,MICSnigeria.W};
dATaMICS{2}           = {MICSnigeria.Region == 1,MICSnigeria.W};
dATaMICS{3}           = {MICSnigeria.State == 19,MICSnigeria.W};
dATaMICS{4}           = {MICSnigeria.State == 18,MICSnigeria.W};
dATaMICS{5}           = {MICSnigeria.State == 33,MICSnigeria.W};

for i = 1:numel(dATaMICS)
    data      = "$\textrm{MICS 6, " + string(sprintf('%0.1f',dATe{2}(1))) + "-" + string(sprintf('%0.1f',dATe{2}(2))) + "}$";
    h         = i + numel(dATa) + numel(dATaB) + numel(dATaDHS);
    pOPs{1,h} = {char(data);models{i}};
    pOPs{2,h} = {char(data{1}(1:end - 2) + " (female)}$");models{i}};
    pOPs{3,h} = {char(data{1}(1:end - 2) + " (male)");models{i}};
end
clear mIn mAx models data

rng(0);
s                     = MICSnigeria(MICSnigeria.k == 1, {'cluster','K','woman','Women','iNDeX'});
w                     = rand(size(s,1),R);
w                     = [(1:size(s,1))',ceil(s.Women.*w) + s.iNDeX];
Wmics                 = NaN(size(MICSnigeria,1),R + 1);
for r = 1:R + 1
    S          = tabulate([w(:,1);w(:,r)]);
    Wmics(:,r) = repelem(S(:,2) - 1,s.K);
    clc;
    r/(R + 1)
end
clear r s S w

for i = 1:numel(dATaMICS)
    dATaMICS{i}{2} = Wmics.*dATaMICS{i}{2};
end
clear Wdhs
load(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','dATaDHS','Ts'); 
save(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','dATaDHS','dATaMICS','Ts','-v7.3','-nocompression');


p                     = rand(size(MICSnigeria,1),R + 1);
d                     = (p.*MICSnigeria.D_min + (1 - p).*MICSnigeria.D_max)/365.25;
p                     = rand(size(MICSnigeria,1),R + 1);
B                     = datetime(datenum(MICSnigeria.B_min) + datenum(MICSnigeria.B_max - MICSnigeria.B_min).*p,'ConvertFrom','datenum');
B                     = eXAcTTime(B);
D                     = B + d;
O                     = min(D,Tmics(2));

p                     = rand(size(MICSnigeria,1),R + 1);
dob                   = p.*eXAcTTime(MICSnigeria.DOB_min) + (1 - p).*eXAcTTime(MICSnigeria.DOB_max);
ageS                  = dATe{2}(2) - dob;
ageB                  = B - dob;

for i = 1:numel(dATaMICS)
    h                    = i + numel(dATa) + numel(dATaB) + numel(dATaDHS);
    s                    = dATaMICS{i}{1};
    w                    = dATaMICS{i}{2};
    for k = 1:3
        if k == 1
            sEX           = ones(size(MICSnigeria,1),1);
        else
            sEX           = (MICSnigeria.sex == sex(k - 1));
        end
        sH                   = dATaMICS{i}{1} & sEX;
        
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(sH,:) + x{1}(j),min(O(sH,:),Tmics(2))),Tmics(1));
            o             = max(min(B(sH,:) + x{1}(j + 1),min(O(sH,:),Tmics(2))),Tmics(1));
            exposure(j,:) = sum((o - a).*w(sH,:));
            events(j,:)   = sum((d(sH,:) >= x{1}(j) & d(sH,:) < x{1}(j + 1) & D(sH,:) >= Tmics(1) & D(sH,:) < Tmics(2)).*w(sH,:));
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end

        SRB                  = ((MICSnigeria.sex(s) == 1 & B(s,:) >= Tmics(1) & B(s,:) < Tmics(2))'*w(s,:))./((MICSnigeria.sex(s) == 2 & B(s,:) >= Tmics(1) & B(s,:) < Tmics(2))'*w(s,:));
        m                    = events./exposure;
        q                    = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];

        TaBle.m{k,h}         = m;
        TaBle.q{k,h}         = q;    
        TaBle.SRB{k,h}       = SRB;
        TaBle.tAU{k,h}       = dATe{1}';
        clear sEX events exposure m q SRB
    end
    
    A                    = [max(min(MICSnigeria.age(dATaMICS{i}{1})),15),min(max(MICSnigeria.age(dATaMICS{i}{1})),49)];
    Ax                   = (min(MICSnigeria.age(dATaMICS{i}{1})):max(MICSnigeria.age(dATaMICS{i}{1})))';
    Ax                   = 10:55;
    for j = 1:numel(Ax)
        a               = max(min(ageS(s,:) - 3,Ax(j) + 1),Ax(j));
        o               = max(min(ageS(s,:),Ax(j) + 1),Ax(j));
        exposure        = sum((MICSnigeria.k(s) == 1).*(o - a).*w(s,:));
        events          = sum((dATe{2}(2) - B(s,:) <= 3 & dATe{2}(2) - B(s,:) > 0).*(ageB(s,:) >= Ax(j) & ageB(s,:) < Ax(j) + 1).*w(s,:));
        TFR(j,:)        = events./max(exposure,eps);
        clear events exposure a o
        clc;
        j/numel(Ax)
    end
    
    TaBle.F{1,h}         = TFR;
    TaBle.TFR{1,h}       = sum(TFR(Ax >= A(1) & Ax < A(2),:),1);
    TaBle.sRB{1,h}       = sum((MICSnigeria.sex(s) == 1).*w(s,:))./sum((MICSnigeria.sex(s) == 2).*w(s,:));
    TaBle.parity{1,h}    = sum(~isnan(MICSnigeria.bidx(s)).*w(s,:))./sum((MICSnigeria.k(s) == 1).*w(s,:));
    TaBle.childless{1,h} = sum((MICSnigeria.mother(s) ~= 1 & MICSnigeria.k(s) == 1).*w(s,:))./sum((MICSnigeria.k(s) == 1).*w(s,:))*100;
    clear h j TFR w A
end

for i = 1:numel(dATaMICS)
    h              = i + numel(dATa) + numel(dATaB) + numel(dATaDHS);
    TaBle.sEX{1,h} = log(1 - TaBle.q{3,h})./log(1 - TaBle.q{2,h});
end

setH             = ones(size(MICSnigeria.W));
setH(:,end + 1)  = (MICSnigeria.sex == 2);
setH(:,end + 1)  = (MICSnigeria.sex == 1);
setH(:,end + 1)  = (MICSnigeria.UR == 1);
setH(:,end + 1)  = (MICSnigeria.UR == 2);
setH(:,end + 1)  = (MICSnigeria.Education == 1 | MICSnigeria.Education == 2);
setH(:,end + 1)  = (MICSnigeria.Education == 3);
setH(:,end + 1)  = (MICSnigeria.Electricity == 2);
setH(:,end + 1)  = (MICSnigeria.Electricity == 1);

for i = 1:numel(dATaMICS)
    k                = i + numel(dATa) + numel(dATaB) + numel(dATaDHS);
    for h = 1:size(setH,2)
        s               = dATaMICS{i}{1} & setH(:,h); 
        wH              = dATaMICS{i}{2}(s,:);
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(s,:) + x{1}(j),min(O(s,:),Tmics(2))),Tmics(1));
            o             = max(min(B(s,:) + x{1}(j + 1),min(O(s,:),Tmics(2))),Tmics(1));
            exposure(j,:) = sum((o - a).*wH);
            events(j,:)   = sum((d(s,:) >= x{1}(j) & d(s,:) < x{1}(j + 1) & D(s,:) >= Tdhs(1) & D(s,:) < Tdhs(2)).*wH);
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end
        
        SRB             = ((MICSnigeria.sex(s) == 1 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*wH)./((MICSnigeria.sex(s) == 2 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*wH);
        m               = events./exposure;
        q               = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        
        TaBle.m{h,k}    = m;
        TaBlE.q{h,k}    = q;
        TaBlE.tAU{h,k}  = Tmics;
        clear events exposure q m births stillbirths SRB s wH
    end
    clear h j k
end
clear p d b B D O dob ageS ageB sET i



load(char(pATh + "Results/bASe.mat"),'IGMEnigeria');
IGME             = mat2cell(table2array(IGMEnigeria),size(IGMEnigeria,1),[1 ones(1,(size(IGMEnigeria,2) - 1)/3)*3]);
for i = 2:numel(IGME)
    IGME{i} = IGME{i}/1000;
end

A{1}             = log(IGME{2}(1:22,:));
A{2}             = log(IGME{5}(1:22,:));
for i = 1:2
    for j = 1:size(A{i},1)
        delta = 0.005;
        mu    = (A{i}(j,2) + A{i}(j,3))/2;
        sigma = 1;
        E     = 10;
        while E^2 > eps
            dE    = (normcdf(A{i}(j,2),mu,sigma*exp(delta)) - normcdf(A{i}(j,2),mu,sigma*exp(-0.005)))/(2*delta);
            sigma = sigma*exp(-dE*E);
            E     = normcdf(A{i}(j,2),mu,sigma) - 0.05;
        end
        A{2 + i}(j,:) = exp(normrnd(mu,sigma,1,10000));
        clear E sigma mu delta
    end
end
A{4}             = A{4}./(1 + A{4});
A{3}             = A{4} + (1 - A{4}).*A{3};
IGME{7}          = prctile(A{3}',[50 5 95])';
IGME{5}          = IGME{5}(1:22,:)./(1 + IGME{5}(1:22,:));
IGME{6}          = NaN(size(IGME{5}));

h                = 1 + numel(dATa) + numel(dATaB) + numel(dATaDHS) + numel(dATaMICS);
k                = 4;
for i = 1:k
    q               = NaN(numel(x{1}),3);
    q(aGEs,:)       = [IGME{2}(k - i + 1,:);IGME{3}(k - i + 1,:);IGME{4}(k - i + 1,:)];
    TaBle.q{i,h}    = q;
    s               = [IGME{5}(k - i + 1,:);IGME{6}(k - i + 1,:);IGME{7}(k - i + 1,:)];
    TaBle.s{i,h}    = s;
    pOPs{i,h}       = {char("UN IGME, " + string(IGME{1}(k - i + 1)));''};
    TaBle.tAU{i,h}  = datetime(IGME{1}(k - i + 1) - .5,7,1);
    clear q s 
end
TaBle.q{i + 1,h}   = {IGME{2},IGME{3},IGME{4},datetime(IGME{1} - .5,7,1)};
TaBle.s{i + 1,h}   = {IGME{5},IGME{6},IGME{7},datetime(IGME{1} - .5,7,1)};
pOPs{i + 1,h}      = {'UN IGME';''};
TaBle.tAU{i + 1,h} = datetime(IGME{1} - .5,7,1);



mAX                = 2;
tHtiles            = [50 2.5 97.5];
selection          = {[1 1 2],[1 3 2],[1 5 1],[1 6 1],[1 7 1],[1 8 1],[1 9 1],[1 10 1],[1 12 1],[1 13 1],[1 14 1],[1 15 1]};
lISt               = {TaBle.childless,TaBle.parity,TaBle.TFR,TaBle.sRB};
for i = 1:numel(selection)
    temp       = pOPs{selection{i}(1),selection{i}(2)};
    temp       = {char(temp{1});temp{2}};
    if numel(strfind(temp{1},"(")) > 0
        temp{1}    = char(string(temp{1}(1:strfind(temp{1},"(") - 2)) + "}$");
    end
    if numel(strfind(temp{1},",")) > 0
        temp{1}    = char(string(temp{1}(1:strfind(temp{1},",") - 1)) + "}$");
    end
    label{i,1} = temp;
    clear temp
    for j = 1:numel(lISt)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end

        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = lISt{j}{selection{i}(1),selection{i}(2) + k - 1};
            if size(temp,2) == R + 1 
                table{i,h} = prctile(temp(:,2:end)',tHtiles);
            elseif numel(temp) > 0
                table{i,h} = temp;
            end
        end
    end
end

sEt          = {'$\textrm{Childless Women \%}$','$\textrm{Average Parity}$','$\textrm{TFR}$','$\textrm{SRB}$'};
models       = {'$\textit{15-49}$','$\textit{10-55}$'};
vARs         = {models models models models};
foRMaT       = {'%0.2f','%0.2f','%0.2f'};
nOTe         = {'$\textrm{Sample}$','$\textit{Bootstrapping}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$'};
lABs         = {{1} {2} {3} {4} {5 6 7 8} {9 10 11 12}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.120,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_1.png"),'Resolution',RESolUTioN);
clear table label temp

selection  = {[2 1 2],[3 1 2],[4 1 2],[5 1 2],[6 1 2],[2 3 2],[3 3 2],[4 3 2],[5 3 2],[6 3 2],[1 17 1],[2 17 1],[3 17 1],[4 17 1]};
mAX        = 2;
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:numel(aGEs)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.q{selection{i}(1),selection{i}(2) + k - 1}(aGEs(j),:);
            if size(temp,2) > 3
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

sEt          = {'NMR - $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','IMR - $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','U5MR - $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};
vARs         = {models models models};
lABs         = {{1 2 3 4 5} {6 7 8 9 10} {11 12 13 14}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.140,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_2.png"),'Resolution',RESolUTioN);
clear table label temp

selection  = {[1 1 2],[1 3 2],[1 5 1],[2 5 1],[1 6 1],[2 6 1],[1 7 1],[1 8 1],[1 9 1],[1 10 1],[1 12 1],[1 13 1],[1 14 1],[1 15 1],[1 17 1],[2 17 1],[3 17 1],[4 17 1]};
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:numel(aGEs)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.q{selection{i}(1),selection{i}(2) + k - 1}(aGEs(j),:);
            if size(temp,2) == R + 1
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

sEt          = {'NMR - $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','IMR - $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','U5MR - $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};;
vARs         = {models models models};
lABs         = {{1} {2} {3 4} {5 6} {7 8 9 10} {11 12 13 14} {15 16 17 18}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.200,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_3.png"),'Resolution',RESolUTioN);
clear table label temp

tAU{1}       = ", " + string(sprintf('%0.1f',TaBlE.tAU{1,1}(1))) + "-" + string(sprintf('%0.1f',TaBlE.tAU{1,1}(2)));
tAU{2}       = ", " + string(sprintf('%0.1f',TaBlE.tAU{1,3}(1))) + "-" + string(sprintf('%0.1f',TaBlE.tAU{1,3}(2)));
tAU{3}       = ", " + string(sprintf('%0.1f',TaBlE.tAU{1,7}(1))) + "-" + string(sprintf('%0.1f',TaBlE.tAU{1,7}(2)));
tAU{4}       = ", " + string(sprintf('%0.1f',TaBlE.tAU{1,12}(1))) + "-" + string(sprintf('%0.1f',TaBlE.tAU{1,12}(2)));
sEt          = {"$\textrm{" + cLUsTEr{1} + tAU{1} + "}$","$\textrm{" + cLUsTEr{2} + tAU{2} + "}$","$\textrm{DHS VII" + tAU{3} + "}$","$\textrm{MICS 6" + tAU{4} + "}$"};
DHS          = {'$\textit{National}$','$\textit{North West}$','$\textit{Kano}$','$\textit{Kaduna}$'};
vARs         = {models models DHS DHS};
selection    = [1 2 3 4 7 8 9 10 11 12 13 14];
for i = 1:size(TaBlE.q,1)
    for j = 1:numel(selection)
        bOx{i,j}  = prctile(TaBlE.q{i,selection(j)}(5,2:end),[50 2.5 97.5])*1000;
        bOx2{i,j} = prctile(TaBlE.q{i,selection(j)}(end,2:end),[50 2.5 97.5])*1000;
    end
end
lABs         = {{1} {4 5} {6 7} {8 9}};
nOTe         = {'$\textrm{Attributes}$/$\textrm{Sample}$','$\textit{Bootstrapping}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$'};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOPh,cell2mat(bOx),0.175,0.070,[]);
exportgraphics(gcf,char(pATh + "Results/Table_4_NMR.png"),'Resolution',RESolUTioN);

tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOPh,cell2mat(bOx2),0.175,0.070,[]);
exportgraphics(gcf,char(pATh + "Results/Table_5_U5MR.png"),'Resolution',RESolUTioN);



P                        = cell(0);
selection                = {[1 1],[1 3],[2 5],[2 6],[1 7],[1 12],[1 17]};
selectionT               = {[1 1],[1 3],[3 5],[3 6],[1 7],[1 12],[5 17]};
LAB                      = {'Neonatal mortality rate $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','Infant mortality rate $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','Under-five mortality rate $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};
load(char(pATh + "Results/paleTTe.mat"),'paleTTe');
coloR                    = paleTTe([1 2 3 5 7 4 8]);

mPIX                     = 538756;
pix                      = 1/37.7952755906;
z                        = min(sqrt(mPIX/((10*3)*(10*2)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 3 2]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for i = 1:6
    nexttile(i)
    if i > 3
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.XScale                = 'log';
        ax{i}.XAxis.TickLabelFormat = '%.1f';
        ax{i}.YAxis.TickLabelFormat = '%.2f';
        ax{i}.XTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.XAxis.MinorTickValues = 10:10:200;
        xlim([6.25 200])
        xlabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$ (log scale)','Interpreter','latex','FontSize',11*z);
        ylabel('$\mathit{kernel}$ $\mathit{density}$','Interpreter','latex','FontSize',11*z);
        title(char(string(char(96 + i)) + ". " + string(LAB{i - 3})),'Interpreter','latex');
    else
        d                           = scatter(datetime('01-Jan-2021'),0);
        delete(d)
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.YAxis.TickLabelFormat = '%.1f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickLabelFormat = 'yyyy';
        ax{i}.XTick                 = datetime(2012:2:2026,1,1);
        ax{i}.XAxis.MinorTickValues = datetime(2012:2026,1,1);
        ax{i}.YScale                = 'log';
        ax{i}.YAxis.TickLabelFormat = '%.1f';
        ax{i}.YTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.YAxis.MinorTickValues = 10:10:200;
        ylim([6.25 200])
        xlabel('$\mathit{year}$','Interpreter','latex','FontSize',11*z);
        ylabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$ (log scale)','Interpreter','latex','FontSize',11*z);
        xlim([datetime(2012,1,1) datetime(2026,1,1)]);
        title(char(string(char(96 + i)) + ". " + string(LAB{i})),'Interpreter','latex','FontSize',12*z);
    end
    grid on;
    box on;
    hold on;
end

for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) == R + 1
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.00*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00*z,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
        end
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_1.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P

P                        = cell(0);
selection                = {[7 1],[8 1],[2 9],[3 9],[2 14],[3 14]};
selectionT               = {[7 1],[8 1],[2 9],[3 9],[2 14],[3 14]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) == R + 1
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.00*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00*z,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
        end
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_2.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P

P                        = cell(0);
selection                = {[7 3],[8 3],[2 10],[3 10],[2 15],[3 15]};
selectionT               = {[7 3],[8 3],[2 10],[3 10],[2 15],[3 15]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) == R + 1
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.00*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00*z,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
        end
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_3.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P

coloR                    = paleTTe([1 1 1 1 1 2 7 4 8]);
P                        = cell(0);
selection                = {[2 1],[3 1],[4 1],[5 1],[6 1],[2 5],[1 9],[1 14],[1 17]};
selectionT               = {[2 1],[3 1],[4 1],[5 1],[6 1],[3 5],[1 9],[1 14],[5 17]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) == R + 1
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.00*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00*z,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
        end
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_4.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P


P                        = cell(0);
selection                = {[2 3],[3 3],[4 3],[5 3],[6 3],[2 6],[1 10],[1 15],[1 17]};
selectionT               = {[2 3],[3 3],[4 3],[5 3],[6 3],[3 6],[1 10],[1 15],[5 17]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) == R + 1
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.00*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00*z,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
        end
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_5.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P


LAB                      = {'Age-Specific Fertility Rates','Age-Specific Mortality Rates','Cumulative Probability of Dying'};
z                        = min(sqrt(mPIX/((10*3)*(10*2)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 3 2]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for i = 1:6
    nexttile(i)
    if isequal(mod(i,3),1)
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.XAxis.TickLabelFormat = '%.0f';
        ax{i}.YAxis.TickLabelFormat = '%.2f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickValues      = 10.00:5.00:55.00;
        ax{i}.XAxis.MinorTickValues = 10.00:2.50:55.00;
        xlim([10 55])
        ylim([0 0.35])
        xlabel('$\textit{age}$ $\textit{(in years)}$','Interpreter','latex','FontSize',11*z);
        ylabel('$\mathit{_nf_x}$','Interpreter','latex','FontSize',11*z);
        title(char("$\textrm{" + string(char(96 + i)) + ". " + string(LAB{1}) + ", " + cLUsTEr{floor((i - 1)/3) + 1} + "}$"),'Interpreter','latex');
    elseif isequal(mod(i,3),2)
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.YScale                = 'log';
        ax{i}.XAxis.TickLabelFormat = '%.1f';
        ax{i}.YAxis.TickLabelFormat = '%.2f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickValues      = 0.00:1.00:5.00;
        ax{i}.XAxis.MinorTickValues = x{1};
        xlim([-0.1 5])
        xlabel('$\textit{age}$ $\textit{(in years)}$','Interpreter','latex','FontSize',11*z);
        ylabel('$\mathit{_nM_x}$','Interpreter','latex','FontSize',11*z);
        title(char("$\textrm{" + string(char(96 + i)) + ". " + string(LAB{2}) + ", " + cLUsTEr{floor((i - 1)/3) + 1} + "}$"),'Interpreter','latex');
    else
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.XAxis.TickLabelFormat = '%.1f';
        ax{i}.XAxis.TickLabelFormat = '%.1f';
        ax{i}.YAxis.TickLabelFormat = '%.2f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickValues      = 0.00:1.00:5.00;
        ax{i}.XAxis.MinorTickValues = x{1};
        xlim([-0.1 5])
        xlabel('$\textit{age}$ $\textit{(in years)}$','Interpreter','latex','FontSize',11*z);
        ylabel('$\mathit{q}\mathrm{(}\mathit{x}\mathrm{)}$','Interpreter','latex','FontSize',11*z);
        title(char("$\textrm{" + string(char(96 + i)) + ". " + string(LAB{3}) + ", " + cLUsTEr{floor((i - 1)/3) + 1} + "}$"),'Interpreter','latex');
    end
    grid on;
    box on;
    hold on;
end

for i = 1:6
    nexttile(i)
    if isequal(floor((i - 1)/3) + 1,1)
        sEL    = [1 9 14];
        leGend = {'$\textrm{SARMAAN, Kano}$','$\textrm{DHS VII, Kano}$','$\textrm{MICS 6, Kano}$'};
        coloR  = paleTTe([2 7 4]);
    else
        sEL    = [3 10 15];
        leGend = {'$\textrm{SARMAAN, Kaduna}$','$\textrm{DHS VII, Kaduna}$','$\textrm{MICS 6, Kaduna}$'};
        coloR  = paleTTe([3 7 4]);
    end

    if isequal(mod(i,3),1)
        Y   = TaBle.F(1,sEL);
        X   = 10:55;
    elseif isequal(mod(i,3),2)
        Y   = TaBle.m(1,sEL);
        X   = (x{1}(1:end - 1) + x{1}(2:end))/2;
    else
        Y   = TaBle.q(1,sEL);
        X   = x{1};
    end

    for j = 1:numel(Y)
        plot(NaN,NaN,'color',coloR{j},'LineWidth',1.00*z);
    end
    for j = numel(Y):-1:1
        plot(X,Y{j},'color',[coloR{j} 0.025],'LineWidth',0.50*z);
        plot(X,prctile(Y{j}(:,2:end)',50)','color',coloR{j},'LineWidth',1.50*z);
    end
    if isequal(mod(i,3),2)
        legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_6.png"),'Resolution',RESolUTioN);



options                  = detectImportOptions(char(pATh + "GPSlimits.csv")); 
options.VariableTypes{1} = 'string';
lIMiTs                   = readtable(char(pATh + "GPSlimits.csv"),options);
country                  = {'NG'};

leGend                   = {};
coloR                    = paleTTe([2 3]);
z                        = min(sqrt(mPIX/((20)*(16)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*[0 0 20 16]/pix,'Theme','light');
axes1                    = geoaxes('Parent',fi,'Position',[0.05 0.05 0.95 0.95]);
hold(axes1,'on');

geobasemap grayland
ax                       = gca;
ax.FontName              = 'Times New Roman';
ax.FontSize              = 8*z;
limits                   = lIMiTs(lIMiTs.country == 'NG',:);
geolimits([limits.minLAT limits.maxLAT],[limits.minLONG limits.maxLONG])

for j = 1:numel(cLUsTEr)
    leGend{j}  = char("$\textrm{" + cLUsTEr{j} + "}");
    sEL        = (bASe.k == 1 & bASe.cluster == cLUsTEr{j});
    S          = geoscatter(bASe.latitude(sEL),bASe.longitude(sEL),2.5,'filled','MarkerEdgeColor',coloR{j},'MarkerFaceColor',coloR{j},'MarkerFaceAlpha',.65);
end
legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
exportgraphics(gcf,char(pATh + "Results/mAP_1.png"),'Resolution',RESolUTioN);

geobasemap satellite
exportgraphics(gcf,char(pATh + "Results/mAP_2.png"),'Resolution',RESolUTioN);
exportgraphics(gcf,char(pATh + "Results/mAP_3.png"),'Resolution',RESolUTioN);



clear
pATh         = "/Users/lshjr3/Documents/SARMAAN/";
RESolUTioN   = 300;
load(char(pATh + "Results/bASe.mat"),'bASe','DHSnigeria','MICSnigeria');
load(char(pATh + "Results/ReSAmPLiNG.mat"));
dATaDHS      = dATaDHS(end - 4:end - 1);
dATaMICS     = dATaMICS(end - 4:end - 1);

cLUsTEr      = {'Kano','Kaduna'};
list         = {char("$\textrm{" + cLUsTEr{1} + "}$"),char("$\textrm{" + cLUsTEr{2} + "}$")};
models       = {'$\textit{15-49}$','$\textit{10-55}$'};
listDHS      = {'$\textrm{DHS VII}$'};
listMICS     = {'$\textrm{MICS 6}$'};
modelsDHS    = {'$\textit{National}$','$\textit{North West}$','$\textit{Kano}$','$\textit{Kaduna}$'};
modelsMICS   = modelsDHS;

for i = 1:numel(list)
    for j = 1:numel(models)
        I          = numel(models)*(i - 1) + j;
        date       = max(bASe.interview);
        ts{I,1}    = datetime([year(date) - 5 year(date)]',[1 month(date)]',[1 day(date)]');
        date       = eXAcTTime(ts{I,1});
        date       = string(list{i}(1:end - 2)) + ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2))) + "}$";
        pOPsD{I,1} = {char(date);models{j}};
        pOPsS{I,1} = list{i};
        clear data I
    end
end

for i = 1:numel(modelsDHS)
    mIn              = min(DHSnigeria.interview);
    mAx              = max(DHSnigeria.interview);
    mAx              = mean([mIn mAx]);
    mIn              = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
    date             = [mIn mAx]';
    ts{end + 1,1}    = date;
    date             = eXAcTTime(date);    
    date             = ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
    date             = char("DHS VII" + date);
    pOPsD{end + 1,1} = {date;modelsDHS{i}};
    pOPsS{end + 1,1} = listDHS{1};
    clear deta mIn mAx
end

for i = 1:numel(modelsMICS)
    mIn              = min(MICSnigeria.interview);
    mAx              = max(MICSnigeria.interview);
    mAx              = mean([mIn mAx]);
    mIn              = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
    date             = [mIn mAx]';
    ts{end + 1,1}    = date;
    date             = eXAcTTime(date);    
    date             = ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
    date             = char("MICS 6" + date);
    pOPsD{end + 1,1} = {date;modelsMICS{i}};
    pOPsS{end + 1,1} = listMICS{1};
    clear deta mIn mAx
end


pACk         = {{bASe,dATa},{DHSnigeria,dATaDHS},{MICSnigeria,dATaMICS}};
lABelS       = {'Place of residence','Age','Education','Electricity','Drinking water'};
lABelSd      = {{'urban' 'rural'} {'10-14' '15-19' '20-24' '25-29' '30-34' '35-39' '40-44' '45-49' '50+'} {'less than complete primary' 'incomplete secondary' 'complete secondary or more'} {'access' 'no access'} {'safe source' 'other source'}}; 
outcomes     = {[1 2],[1 2 3 4 5 6 7 8 9],[1 2 3],[2 1],[2 1]};

H            = 0;
for h = 1:numel(pACk)
    d    = pACk{h}{1};
    data = [d.UR,d.GO,d.Education,d.Electricity,d.Water];
    for i = 1:numel(pACk{h}{2})
        H  = H + 1;
        s  = pACk{h}{2}{i}{1} & (d.k == 1);
        sW = pACk{h}{2}{i}{2}(s,:);
        I  = 0;
        for j = 1:numel(outcomes)
            for k = 1:numel(outcomes{j})
                I        = I + 1;
                BST      = ((data(s,j) == outcomes{j}(k))'*sW)./sum(sW,1);
                bOx{I,H} = prctile(100*BST(2:end),[50 2.5 97.5]);
                
                if isequal(h,1) & isequal(i,1)
                    if k == 1
                        pOPsd{I,1} = {char(string(lABelS{j}) + ": " + string(lABelSd{j}{k}))};
                    else
                        pOPsd{I,1} = {lABelSd{j}{k}};
                    end
                end
            end
        end
        
        if isequal(h,1) & isequal(i,1)
            bOx{end + 1,H}   = [sum(s) NaN NaN];
            pOPsd{end + 1,1} = {'Observations'};            
        else
            bOx{end,H}       = [sum(s) NaN NaN];
        end
    end
end

selection    = 1 + [0 cumsum([ones(1,2)*numel(models) numel(modelsDHS)])];
for i = 1:numel(selection)
    sEt{i} = {pOPsS{selection(i)}};
end


lABs         = {{1} {3 4 5 6 7 8 9 10 11} {12 13 14} {15} {17} {19}};
vARs         = {models models modelsMICS modelsDHS};
foRMaT       = {'%0.2f','%0.2f','%0.2f'};

nOTe         = {'$\textrm{Attributes}$/$\textrm{Instrument}$','$\mathit{Bootstrapping}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$'};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOPsd,cell2mat(bOx),0.190,0.065,[]);
exportgraphics(gcf,char(pATh + "Results/Table_6.png"),'Resolution',RESolUTioN);



clear
pATh                  = "/Users/lshjr3/Documents/SARMAAN/";
load(char(pATh + "Results/bASe.mat"),'bASe','DHSnigeria','MICSnigeria');
RESolUTioN            = 300;
coloR                 = {[0.05 0.05 0.05],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.00 0.00 0.75],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65]};
sAMplE                = {'Children','Women'};
vAR                   = {'year','month of the year','day of the month'};
sOUrCe                = {'$\textrm{SARMAAN}$','$\textrm{DHS VII}$','$\textrm{MICS 6}$'};
aLPha                 = [.2 .1 .1];

pix                   = 1/37.7952755906*0.75;
fi                    = figure('Color',[1 1 1]);
fi.Theme              = 'light';
fi.Position           = [0 0 numel(sAMplE) numel(vAR)]*7/pix;
axes1                 = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                    = tiledlayout(numel(vAR),numel(sAMplE),'Padding','compact','TileSpacing','compact');

for i = 1:numel(sAMplE)*numel(vAR)
    nexttile(i)
    ax{i}                          = gca;
    ax{i}.FontName                 = 'Times New Roman';
    ax{i}.FontSize                 = 10;
    ax{i}.XAxis.TickLabelFormat    = '%.0f';
    ax{i}.YAxis.TickLabelFormat    = '%.2f';
    ax{i}.XAxis.MinorTick          = 'off';
    ax{i}.YAxis.MinorTick          = 'off';
    ax{i}.LabelFontSizeMultiplier  = 1;
    
    if i <= numel(sAMplE)
        title(char("$\textrm{DOB: " + sAMplE{2 - mod(i,2)} + "}$"),'Interpreter','latex','FontSize',11);
    end
    xlabel(char("$\textit{" + vAR{ceil(i/2)} + "}$"),'Interpreter','latex','FontSize',10);
    if isequal(mod(i,numel(sAMplE)),1)
        ylabel('$\textit{probability density function}$','Interpreter','latex','FontSize',10);
    end
    grid on;
    box on;
    hold on;
    
    if isequal(mod(i,numel(sAMplE)),1)
        base = datevec(bASe.B(bASe.birth == 'livebirth'));
        base = base(:,1:3);
        DHS  = datevec(DHSnigeria.B_min(~isnat(DHSnigeria.B_min)));
        DHS  = DHS(:,1:3);
        MICS = datevec(MICSnigeria.B_min(~isnat(MICSnigeria.B_min)));
        MICS = MICS(:,1:3);
    else
        base = datevec(bASe.DOB(bASe.k == 1));
        base = base(:,1:3);
        DHS  = datevec(DHSnigeria.DOB(DHSnigeria.k == 1));
        DHS  = [DHS(:,1:2) NaN(size(DHS,1),1)];
        MICS = datevec(MICSnigeria.DOB_min(MICSnigeria.k == 1));
        MICS = [MICS(:,1:2) NaN(size(MICS,1),1)];
    end
    j    = ceil(i/numel(sAMplE));    
    dATa = {base(:,j),DHS(:,j),MICS(:,j)};
    
    if isequal(j,2)
        xlim([.5 12.5])
    elseif isequal(j,3)
        xlim([.5 31.5])
    end
    
    for j = 1:numel(dATa)
        histogram(dATa{j},'Normalization','pdf','FaceColor',coloR{j},'FaceAlpha',aLPha(j),'EdgeColor',coloR{j},'EdgeAlpha',min(aLPha(j)*1.75,1));
    end
    if isequal(i,numel(sAMplE)*numel(vAR) - 1)
        legend(sOUrCe,'Interpreter','latex','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',numel(sOUrCe),'Box','off');
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_8.png"),'Resolution',RESolUTioN);





clear
pATh                  = "/Users/lshjr3/Documents/SARMAAN/";
load(char(pATh + "Results/bASe.mat"),'bASe','DHSnigeria','MICSnigeria');
RESolUTioN            = 300;

bASe                  = bASe(bASe.age >= 15 & bASe.age <  50,:);
bASe.id               = cumsum(bASe.k == 1);
bASe                  = bASe(bASe.birth == 'livebirth',{'B_min','B_max','DOB_min','DOB_max','id'});
bASe.t                = zeros(size(bASe,1),1);
bASe.t(2:end)         = bASe.id(1:end - 1) == bASe.id(2:end) & bASe.B_min(1:end - 1) == bASe.B_min(2:end);
bASe                  = bASe(~bASe.t,:);
rng(0);
p                     = rand(size(bASe,1),1);
bASe.B                = eXAcTTime(bASe.B_min).*p + eXAcTTime(bASe.B_max).*(1 - p);
bASe                  = sortrows(bASe,{'id','B'});
for i = 1:bASe.id(end)
    sEL                = find(bASe.id == i);
    K                  = numel(sEL);
    bASe.k(sEL)        = (1:K)';
    bASe.K(sEL)        = K;
end

sEL                   = (bASe.k == 1);
p                     = rand(sum(sEL),1);
bASe.DOB              = repelem(eXAcTTime(bASe.DOB_min(sEL)).*p + eXAcTTime(bASe.DOB_max(sEL)).*(1 - p),bASe.K(sEL));
bASe.id               = cumsum(sEL);
dATa{1}               = bASe(:,{'B','DOB','k','K','id'});
clear p ID bASe sEL K i

DHSnigeria            = DHSnigeria(~isnat(DHSnigeria.B_min),:);
sEL                   = (DHSnigeria.k == 1);
DHSnigeria.id         = cumsum(sEL);
p                     = rand(sum(sEL),1);
DHSnigeria.DOB        = repelem(eXAcTTime(DHSnigeria.DOB(sEL)).*p + eXAcTTime(datetime(year(DHSnigeria.DOB(sEL)),month(DHSnigeria.DOB(sEL)) + 1,1) - 1).*(1 - p),DHSnigeria.K(sEL));
DHSnigeria.t          = zeros(size(DHSnigeria,1),1);
DHSnigeria.t(2:end)   = DHSnigeria.id(1:end - 1) == DHSnigeria.id(2:end) & DHSnigeria.B_min(1:end - 1) == DHSnigeria.B_min(2:end);
DHSnigeria            = DHSnigeria(~DHSnigeria.t,:);
DHSnigeria.B          = eXAcTTime(DHSnigeria.B_min);
DHSnigeria            = sortrows(DHSnigeria(:,{'B','DOB','id'}),{'id','B'});
for i = 1:DHSnigeria.id(end)
    sEL                = find(DHSnigeria.id == i);
    K                  = numel(sEL);
    DHSnigeria.k(sEL)  = (1:K)';
    DHSnigeria.K(sEL)  = K;
end
dATa{2}               = DHSnigeria(:,{'B','DOB','k','K','id'});
clear p DHSnigeria sEL

MICSnigeria           = MICSnigeria(~isnat(MICSnigeria.B_min),:);
sEL                   = (MICSnigeria.k == 1);
MICSnigeria.id        = cumsum(sEL);
p                     = rand(sum(sEL),1);
MICSnigeria.DOB       = repelem(eXAcTTime(MICSnigeria.DOB_min(sEL)).*p + eXAcTTime(MICSnigeria.DOB_max(sEL)).*(1 - p),MICSnigeria.K(sEL));
MICSnigeria.t         = zeros(size(MICSnigeria,1),1);
MICSnigeria.t(2:end)  = MICSnigeria.id(1:end - 1) == MICSnigeria.id(2:end) & MICSnigeria.B_min(1:end - 1) == MICSnigeria.B_min(2:end);
MICSnigeria           = MICSnigeria(~MICSnigeria.t,:);
p                     = rand(size(MICSnigeria,1),1);
MICSnigeria.B         = eXAcTTime(MICSnigeria.B_min).*p + eXAcTTime(MICSnigeria.B_max).*(1 - p);
MICSnigeria           = sortrows(MICSnigeria(:,{'B','DOB','id'}),{'id','B'});
for i = 1:MICSnigeria.id(end)
    sEL                = find(MICSnigeria.id == i);
    K                  = numel(sEL);
    MICSnigeria.k(sEL) = (1:K)';
    MICSnigeria.K(sEL) = K;
end
dATa{3}               = MICSnigeria(:,{'B','DOB','k','K','id'});
clear p MICSnigeria sEL

for i = 1:numel(dATa)
    mATrIx      = NaN(dATa{i}.id(end),max(dATa{i}.K) + 2);
    mATrIx(:,1) = dATa{i}.DOB(dATa{i}.k == 1);
    for j = 1:max(dATa{i}.K)
        sEL               = dATa{i}.id(dATa{i}.k == j);
        mATrIx(sEL,1 + j) = dATa{i}.B(dATa{i}.k == j);
    end
    mATrIx      = mATrIx(:,2:end) - mATrIx(:,1:end - 1);
    dATa{i}     = mATrIx;
end

coloR                 = {[0.05 0.05 0.05],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.00 0.00 0.75],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65]};
sOUrCe                = {'$\textrm{SARMAAN}$','$\textrm{DHS VII}$','$\textrm{MICS 6}$','$\textit{unlikely}$'};
aLPha                 = [.2 .1 .1];
pix                   = 1/37.79527559068*0.75;
fi                    = figure('Color',[1 1 1]);
fi.Position           = [0 0 3 2]*7/pix;
fi.Theme              = 'light';
axes1                 = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                    = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for i = 1:6
    nexttile(i)
    ax{i}                          = gca;
    ax{i}.FontName                 = 'Times New Roman';
    ax{i}.FontSize                 = 10;
    ax{i}.XAxis.TickLabelFormat    = '%.0f';
    ax{i}.YAxis.TickLabelFormat    = '%.2f';
    ax{i}.XTickLabelRotation       = 0;
    ax{i}.YTickLabelRotation       = 0;
    ax{i}.LabelFontSizeMultiplier  = 1;

    if isequal(i,1)
        f                              = 10;
        title('$\textrm{Maternal debut}$','Interpreter','latex','FontSize',11);
        xlabel('$\textit{age}$','Interpreter','latex','FontSize',10);
        ax{i}.XAxis.TickValues         = 5:5:50;
        xlim([5 35]);
    else
        f                              = 26*7/365.25;
        title(char("$\textrm{Birth interval " + (i - 1) + "}$"),'Interpreter','latex','FontSize',11);
        xlabel('$\textit{years}$','Interpreter','latex','FontSize',10);
        ax{i}.XAxis.TickValues         = 0:1:10;
        xlim([0 7]);
    end

    if isequal(mod(i,3),1)
        ylabel('$\textit{probability density function}$','Interpreter','latex','FontSize',9);
    end
    grid on;
    box on;
    hold on;
        
    for j = 1:numel(dATa)
        histogram(dATa{j}(:,i),'Normalization','pdf','FaceColor',coloR{j},'FaceAlpha',aLPha(j),'EdgeColor',coloR{j},'EdgeAlpha',min(aLPha(j)*1.75,1));
    end

    F                              = ax{i}.YLim;
    fill([0 f f 0],kron(F,ones(1,2)),coloR{7},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineWidth',0.25,'LineStyle','-','EdgeColor',coloR{7});
    ylim(F);
    if isequal(i,5)
        legend(sOUrCe,'Interpreter','latex','FontSize',9,'Location','southoutside','NumColumns',numel(sOUrCe),'Box','off');
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_9.png"),'Resolution',RESolUTioN);


clear
pATh                  = "/Users/lshjr3/Documents/SARMAAN/";
load(char(pATh + "Results/bASe.mat"),'bASe');
load(char(pATh + "Results/paleTTe.mat"),'paleTTe');
RESolUTioN            = 300;

coloR                 = paleTTe([2 3 7 4 5 1]);
bASe                  = bASe(bASe.k == 1,{'sTArt','eNd','sUBmiSSion','enumeratorid','cluster'});
cLUsTEr               = {'Kano','Kaduna'};

for i = 1:numel(cLUsTEr)
    tEMp{i}            = bASe(bASe.cluster == cLUsTEr{i},:);
    tEMp{i}.S          = datetime(year(tEMp{i}.sUBmiSSion),month(tEMp{i}.sUBmiSSion),day(tEMp{i}.sUBmiSSion));
    tEMp{i}.duration   = minutes(tEMp{i}.sUBmiSSion - tEMp{i}.sTArt);
    T                  = tabulate(days(tEMp{i}.S - min(tEMp{i}.S)));
    S{i}               = table(T(:,1) + min(tEMp{i}.S),T(:,2),'VariableNames',{'date','women'});
    S{i}.cumulative    = cumsum(S{i}.women);

    T                  = tabulate(tEMp{i}.enumeratorid);
    T                  = T(T(:,2) > 0,1);
    for j = 1:numel(T)
        s         = tEMp{i}.duration(tEMp{i}.enumeratorid == T(j));
        s         = [T(j) numel(s) prctile(s,[50 2.5 5 10 25 75 90 95 97.5])];
        P{i}(j,:) = table(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10),s(11),'VariableNames',{'enumerator','women','p50','p2','p5','p10','p25','p75','p90','p95','p97'});
    end
    P{i}               = flip(sortrows(P{i},{'p50'}));
end


mPIX                     = 538756;
pix                      = 1/37.7952755906;
LAB                      = {'Data collection schedule','Time to completion','Time to completion (enumerator)'};
z                        = min(sqrt(mPIX/((10*3)*(10*2)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 3 2]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for i = 1:6
    nexttile(i)
    if isequal(mod(i,3),1)
        d                           = scatter(datetime('01-Jan-2021'),0);
        delete(d)
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.YAxis.TickLabelFormat = '%.1f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickLabelFormat = 'dd.MM.yyyy';
        xlim(datetime(2025,[5 8],1))
        xlabel('$\textit{date of submission}$','Interpreter','latex','FontSize',11*z);
        ylabel('$\textit{women interviewed (000)}$','Interpreter','latex','FontSize',11*z);
        title(char("$\textrm{" + string(char(96 + i)) + ". " + string(LAB{1}) + ", " + cLUsTEr{floor((i - 1)/3) + 1} + ", W = " + sum(S{floor((i - 1)/3) + 1}.women) + "}$"),'Interpreter','latex');
    elseif isequal(mod(i,3),2)
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.XAxis.TickLabelFormat = '%.0f';
        ax{i}.YAxis.TickLabelFormat = '%.1f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickValues      = 0:100:1200;
        ax{i}.XAxis.MinorTickValues = 0:25:1200;
        xlim([-100 800])
        xlabel('$\textit{time (in munutes)}$','Interpreter','latex','FontSize',11*z);
        ylabel('$\textit{probability density function}$','Interpreter','latex','FontSize',11*z);
        title(char("$\textrm{" + string(char(96 + i)) + ". " + string(LAB{2}) + ", " + cLUsTEr{floor((i - 1)/3) + 1} + "}$"),'Interpreter','latex');
    else
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10*z;
        ax{i}.XAxis.TickLabelFormat = '%.0f';
        ax{i}.YAxis.TickLabelFormat = '%.0f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickValues      = 0:100:1200;
        ax{i}.XAxis.MinorTickValues = 0:25:1200;
        xlim([-100 800])
        xlabel('$\textit{time (in munutes)}$','Interpreter','latex','FontSize',11*z);
        ylabel('$\textit{enumerator}$','Interpreter','latex','FontSize',11*z);
        title(char("$\textrm{" + string(char(96 + i) + ". " + string(LAB{3})) + ", " + cLUsTEr{floor((i - 1)/3) + 1} + "}$"),'Interpreter','latex');
    end
    grid on;
    box on;
    hold on;
end

for i = 1:6
    nexttile(i)
    h = floor((i - 1)/3) + 1;

    if isequal(mod(i,3),1)
        plot(S{h}.date,S{h}.women/1000,'color',coloR{h},'LineWidth',1.25*z);
    elseif isequal(mod(i,3),2)
        sEL = tEMp{h}.duration(tEMp{h}.duration > -100 & tEMp{h}.duration < 800);
        histogram(sEL,200,'Normalization','pdf','FaceColor',coloR{h},'FaceAlpha',.15,'EdgeColor',coloR{h},'EdgeAlpha',.35);
        F   = ax{i}.YLim;
        fill([-100 0 0 -100],kron(F,ones(1,2)),coloR{3},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineWidth',0.25,'LineStyle','-','EdgeColor',coloR{3});
        ylim(F);
        legend({char("$\textrm{ W = " + sum(S{h}.women) + ", E = " + size(P{h},1) + "}$"),'$\textit{unlikely duration}'},'Interpreter','latex','FontSize',9,'Location','southoutside','NumColumns',1,'Box','off');
    else
        plot(NaN,NaN,'color','k','LineWidth',2.0);
        F   = [0 size(P{h},1) + 1];
        fill([-100 0 0 -100],kron(F,ones(1,2)),coloR{h},'FaceAlpha',.10,'EdgeAlpha',0.25,'LineWidth',0.25,'LineStyle','-','EdgeColor',coloR{3});
        ylim(F);

        for j = 1:size(P{h},1)
            A = table2array(P{h}(j,4:end));
            A = [A(1:numel(A)/2)' flip(A(numel(A)/2 + 1:end))'];
            for k = 1:size(A,1)
                fill([A(k,1) A(k,2) A(k,2) A(k,1)],kron(j - [1 0],ones(1,2)),coloR{k + 2},'FaceAlpha',.15,'EdgeAlpha',0.00,'LineWidth',0.01,'LineStyle','-','EdgeColor',coloR{k + 2});
            end
            plot([1 1]*P{h}.p50(j),j - [1 0],'color','k','LineWidth',2.0);
        end
        legend({'$\textrm{median duration}$','$\textit{unlikely duration}'},'Interpreter','latex','FontSize',9,'Location','southoutside','NumColumns',1,'Box','off');
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_7.png"),'Resolution',RESolUTioN);



clear
pATh                     = "/Users/lshjr3/Documents/SARMAAN/";
load(char(pATh + "Results/bASe.mat"),'bASe');
load(char(pATh + "Results/paleTTe.mat"),'paleTTe');
RESolUTioN               = 300;
coloR                    = paleTTe([2 3 7 4 5 1]);
cLUsTEr                  = {'Kano','Kaduna'};

mPIX                     = 538756;
pix                      = 1/37.7952755906;
z                        = min(sqrt(mPIX/((15)*(15*2)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*15*[0 0 2 1]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

for i = 1:2
    nexttile(i)
    d                           = scatter(datetime('01-Jan-2021'),0);
    delete(d)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10*z;
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.YMinorGrid            = 'on';
    ax{i}.XMinorGrid            = 'on';
    ax{i}.XAxis.TickLabelFormat = 'yyyy';
    ax{i}.XTick                 = datetime(2016:2:2026,1,1);
    ax{i}.XAxis.MinorTickValues = datetime(2016:2026,1,1);
    ax{i}.YTick                 = 0:1:5;
    ax{i}.YAxis.MinorTickValues = [(0:7:28)/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
    ylabel('$\textit{age (in years)}$','Interpreter','latex','FontSize',11*z);
    xlabel('$\textit{year}$','Interpreter','latex','FontSize',11*z);
    xlim([datetime(2016,1,1) datetime(2026,1,1)]);
    ylim([0 8])
    title(char("$\textrm{" + string(char(96 + i)) + ". " + cLUsTEr{i} + "}$"),'Interpreter','latex','FontSize',12*z);
    clear d
    grid on;
    box on;
    hold on;
end

bASe                     = bASe(:,{'B_min','B_max','D_min','D_max','interview','cluster'});
sEL                      = ~isnan(bASe.D_min);
bASe.D_max(sEL)          = max(min(max(days(bASe.interview(sEL) - bASe.B_max(sEL)),0),bASe.D_max(sEL)),bASe.D_min(sEL));

r                        = rand(size(bASe,1),2);
bASe.B                   = bASe.B_min + r(:,1).*days(min(bASe.B_max,bASe.interview) - bASe.B_min);
s                        = max(days(bASe.interview - bASe.B),0);
bASe.xd                  = bASe.D_min + r(:,2).*min(s,(bASe.D_max - bASe.D_min));
bASe.D                   = bASe.B + bASe.xd;

for i = 1:2
    nexttile(i);
    sEL = (bASe.cluster == cLUsTEr{i} & ~isnat(bASe.D));
    scatter(bASe.D(sEL),bASe.xd(sEL)/365.25,2.5*z,'filled','MarkerFaceColor',coloR{3},'MarkerFaceAlpha',.25,'MarkerEdgeColor',coloR{3},'MarkerEdgeAlpha',.25);
    
    A   = [min(bASe.interview(bASe.cluster == cLUsTEr{i})) max(bASe.interview(bASe.cluster == cLUsTEr{i})) datetime(2026,1,1)];
    fill([A(1) A(2) A(2) A(1)],kron([0 8],ones(1,2)),coloR{1},'FaceAlpha',.25,'EdgeAlpha',0.00,'LineWidth',0.01,'LineStyle','-','EdgeColor',coloR{1});
    fill([A(2) A(3) A(3) A(2)],kron([0 8],ones(1,2)),coloR{2},'FaceAlpha',.25,'EdgeAlpha',0.00,'LineWidth',0.01,'LineStyle','-','EdgeColor',coloR{2});
    legend({char("$\textrm{deaths = " + sum(sEL) + "}$"),'$\textrm{data collection}',char("$\textrm{unlikely reports, D = " + sum(bASe.B_min(sEL) + bASe.D_min(sEL) > bASe.interview(sEL)) + ", B = " + sum(bASe.B(sEL) > bASe.interview(sEL)) + "}$")},'Interpreter','latex','FontSize',9,'Location','southoutside','NumColumns',1,'Box','off');
end
exportgraphics(gcf,char(pATh + "Results/Figure_10.png"),'Resolution',RESolUTioN);



