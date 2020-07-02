%{
This pipeline generates scatter plots and 
significance between different groups.
INSTRUCTIONS: 
1) Go to your folder of interest 
2) Define the variable location, which could either be 'Cellular' or
'Nuclear'
3) Double check the Filenames of the CellProfiler output (e.g.
'MyExpt.csv')
%}

clear all; close all; 
location = 'Cellular'; 
folder = cd; 
nuclei = readtable('MyExpt100nM_Nuclei.csv');
cells = readtable('MyExpt100nM_Cells.csv'); 
[rr, cc] = size(cells);
colors = {[0.4, 0.4, 0.4], [0.6, 0.1, 0.0], [0.8, 0.4, 0.7], [0.2, 0.5, 0.4], [0.4, 0.5, 0.2], [0.9, 0.6, 0.3], [0.2, 0.4, 0.6]};
drugs = {'CTRL', 'DAS', 'IMA', 'NIL', 'VAN', 'ERL', 'BOS','CTRL2'};
drug_label = {'CTRL', 'DAS', 'IMA', 'NIL', 'VAN', 'ERL', 'BOS'};
row = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'}; 
col = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'};
if location == 'Cellular'
    %FOR CELLULAR MORPHOMETRY 
    compartment = cells; 
    var = {'AreaShape_Area', 'AreaShape_Eccentricity', 'AreaShape_FormFactor', 'AreaShape_Extent', 'AreaShape_Solidity',...
    'Intensity_IntegratedIntensity_Actin'}; 
    ax_tit = {'Cell Spreading Area (um^2)', 'Cell Eccentricity', 'Cell Form Factor', 'Cell Extent', 'Cell Solidity', 'Mean actin intensity/cell (a.u.)'};
    comp = 'Cell'; 
elseif location == 'Nuclear'
    %FOR NUCLEAR MORPHOMETRY
    compartment = nuclear;     
    var = {'AreaShape_Area', 'AreaShape_FormFactor', 'AreaShape_Extent', 'AreaShape_Eccentricity', 'Math_nucYAP_cellYAP'};
    ax_tit = {'Nuclear Spreading Area (um^2)', 'Nuclear Form Factor', 'Nuclear Extent', 'Nuclear Eccentricity', 'Math_nucYAP_cellYAP'}; 
    comp = 'Nuclear';
end

numArrays = 12; 

for vv = 1:length(var)
    matrix = nan(rr,8); 
    total = cell(numArrays, 1);
    for i = 1:length(col)
        for j = 1:length(row)   
            well = strcat(row(j), col(i));
            ro = strcmp(compartment.Metadata_Well, well) == 1; 
            vect = table2array(compartment(ro, var(vv)));
            matrix(1:length(vect), j) = vect; 
            means(1,j) = nanmean(vect);
        end
        total{i} = matrix; 
    end
 
    total_copy= total;
 
    for n = 1:length(total)
        if mod(n,2) == 0 
            total_copy{n} = fliplr(total{n});
        end
    end
 
 
    tp_48 = [total_copy{1}; total_copy{2}; total_copy{3}]; % Combining the three trials
    tp_48(tp_48 == 0) = NaN;
    tp_24 = [total_copy{4}; total_copy{5}; total_copy{6}]; 
    tp_24(tp_24 == 0) = NaN;
    tp_1h = [total_copy{7}; total_copy{8}; total_copy{9}];
    tp_1h(tp_1h == 0) = NaN;
    tp_30m = [total_copy{10}; total_copy{11}; total_copy{12}];
    tp_30m(tp_30m == 0) = NaN;
 
 
    [tp_48_c, ns_tp_48] = condense(tp_48,0); 
    [tp_24_c, ns_tp_24] = condense(tp_24,0);
    [tp_1h_c, ns_tp_1h] = condense(tp_1h,0);
    [tp_30m_c, ns_tp_30] = condense(tp_30m,0); 
 
    all = {tp_48_c, tp_24_c, tp_1h_c, tp_30m_c};
    nss = {ns_tp_48, ns_tp_24, ns_tp_1h, ns_tp_30};
    title ={'48 h', '24 h', '1 hr', '30 min'};
    facecolors = [[0 0 1];[1 0 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 0 1]];
    %violin plots
    for kk = 1:length(all)
        figure; 
        plotSpread(all{kk}, 'distributionColors', colors); 
        boxplot(all{kk}, 'Colors', 'k', 'plotStyle', 'compact', 'Labels', drug_label, 'LabelOrientation', 'horizontal'); 
        %violin(all{kk}, 'facecolor', facecolors, 'medc', 'k', 'mc', []);
        set(gca, 'XTickLabel', drugs);
        ylabel(ax_tit(vv));
        suptitle(title{kk});
        legend off  
        [p, tbl, stats] = kruskalwallis(all{kk}, [], 'off'); 
        c = multcompare(stats, 'displayopt', 'off'); 
        [vect, ps] = analyze_c(c); 
        sigstar(vect, ps); 
        figuresize(5.5,4,'inches'); 
        print(strcat(folder,'\',comp, var{vv}, '_', title{kk}), '-dpdf'); 
        csvwrite(strcat(folder,'\',comp, var{vv}, '_', title{kk}, '.csv'), nss{kk});
    end
end


function [all_counts,ns] = condense(numbers,second_c) 
% Condense function organizes all measurements into a cell array, each
% array signify 1 drug treatment. 
    CTRL1 = numbers(:,1);
    CTRL1 = CTRL1(isnan(CTRL1)~=1); 
    DAS = numbers(:,2); 
    DAS = DAS(isnan(DAS)~=1);
    IMA = numbers(:,3); 
    IMA = IMA(isnan(IMA)~=1);
    NIL = numbers(:,4);
    NIL = NIL(isnan(NIL)~=1);
    VAN = numbers(:,5); 
    VAN = VAN(isnan(VAN)~=1);
    ERL = numbers(:,6); 
    ERL = ERL(isnan(ERL)~=1);
    BOS = numbers(:,7); 
    BOS = BOS(isnan(BOS)~=1);
    CTRL2 = numbers(:,8); 
    CTRL2 = CTRL2(isnan(CTRL2)~=1);
    everything = {CTRL1, DAS, IMA, NIL, VAN, ERL, BOS, CTRL2};
    for i = 1:length(everything)
        length_e(i) = size(everything{i},1);
    end
    max_e = max(length_e);
    for k = 1:length(everything)
        everything{k} = padnan(everything{k}, max_e);
    end
    everything = cell2mat(everything); 
    if second_c == 1 %if you want to include the second control (ROW H)
        all_counts = everything;
        ns = [length(CTRL1), length(DAS), length(IMA), length(NIL), length(VAN), length(ERL), length(BOS),...
            length(CTRL2)];
    else
        all_counts = everything(:,1:7);
        ns = [length(CTRL1), length(DAS), length(IMA), length(NIL), length(VAN), length(ERL), length(BOS)];
    end
end

function [vect, ps]  = analyze_c(c)
    ind = find(c(:,1) == 1 & c(:,6) < 0.01);
    pair = c(ind, 2);  
    pair2 = cell(1, length(pair));
    for kk = 1:length(pair)
        pair2{kk} = [1, pair(kk)]; 
    end 
    vect = pair2;
    ps = c(ind, 6); 
end

