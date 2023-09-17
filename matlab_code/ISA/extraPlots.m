% Start by reading in the necessary data
%Xbar = readtable("../../data/isaMetadata/surrogateModelWithFixedSampleMetadataRandomOrder.txt");
folder = 'sampleFeaturesCorrTol0.001';

colormap('parula');
scriptfcn;

X = readtable("./" + folder + "/coordinates.csv");
Yvals = readtable("./" + folder + "/algorithm_raw.csv");
YbinSVM = readtable("./" + folder + "/algorithm_svm.csv");
YSVMprediction = readtable("./" + folder + "/portfolio_svm.csv", 'Delimiter',',');
featVals = readtable("./" + folder + "/feature_process.csv");
XallFeat = readtable("./" + folder + "/metadataAllFeatures.csv");
undecided = readtable("./" + folder + "/undecidedInstances.csv");
rulesSelector = readtable("./" + folder + "/rulesSelector.csv");


% Plot rules based selector
clf;
set(gcf,'position',[0,0,500,500])
drawPortfolioSelectionsCustom(table2array(X(:, 2:3)), rulesSelector.("Best_Algorithm"), ["Kriging", "Co-Kriging"], 'Rules-based prediction');
print(gcf,'-dpng',["./" + folder + "/customSelectorRulesBased.png"]);

% Plot sources giving priority to certain sources
% Reorder data to do this
Zsolar = X(strcmp(XallFeat.("Source"), "SOLAR"), :);
Zlit = X(strcmp(XallFeat.("Source"), "Literature"), :);
Zdist = X(strcmp(XallFeat.("Source"), "Disturbance-based"), :);
Zsources = [Zdist array2table(repelem({'Disturbance'}, size(Zdist,1),1));
            Zlit array2table(repelem({'Literature'}, size(Zlit,1),1));
            Zsolar array2table(repelem({'SOLAR'}, size(Zsolar,1),1))];
S = categorical(table2array(Zsources(:, 4)));
clf;
drawSourcesNico(table2array(Zsources(:, 2:3)), S);
print(gcf,'-dpng',["./" + folder + "/customDistributionSources.png"]);


% Plot desired features
clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_highFiBudgetRatio"), "B^r_h");
print(gcf,'-dpng',["./" + folder + "/customExtraFeature_highFiBudgetRatio.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_CC"), "CC");
print(gcf,'-dpng',["./" + folder + "/customExtraFeature_CC.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_budgetRatio"), "B^r");
print(gcf,'-dpng',["./" + folder + "/customFeature_budgetRatio.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_high_ela_level_mmce_lda_50"), "f_h, MMCE^{0.5}_{lda}");
print(gcf,'-dpng',["./" + folder + "/customFeature_mmce.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_high_ic_m0"), "f_h, M_0");
print(gcf,'-dpng',["./" + folder + "/customFeature_ic.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_LCC_sd"), "LCC^{0.2}_{sd}");
print(gcf,'-dpng',["./" + folder + "/customFeature_lcc_sd.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_LCCrel_0_4"), "LCC^{0.2^{1/d}}_{0.4}");
print(gcf,'-dpng',["./" + folder + "/customFeature_lcc_rel_4.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_LCCrel_0_95"), "LCC^{0.2^{1/d}}_{0.95}");
print(gcf,'-dpng',["./" + folder + "/customFeature_lcc_rel_95.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_mid_ela_meta_lin_simple_adj_r2"), "f_h - f_l, R^2_L");
print(gcf,'-dpng',["./" + folder + "/customFeature_R_linear_simple.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_mid_ela_meta_lin_w_interact_adj_r2"), "f_h - f_l, R^2_{LI}");
print(gcf,'-dpng',["./" + folder + "/customFeature_R_linear_interact.png"]);

clf;
drawScatterNico(table2array(X(:, 2:3)), XallFeat.("feature_sample_RRMSE"), "RRMSE");
print(gcf,'-dpng',["./" + folder + "/customFeature_RRMSE.png"]);


% Plot predictions
clf;
drawPortfolioSelectionsCustom(table2array(X(:, 2:3)), table2array(YSVMprediction(:, 2)), ['Kriging', "Co-Kriging"], "Predicted best algorithm")
print(gcf,'-dpng',["./" + folder + "/customSelector.png"]);

% Plot algorithm performance 
clf;
drawScatterNico(table2array(X(:, 2:3)), table2array(Yvals(:, 2)), "Kriging");
print(gcf,'-dpng',["./" + folder + "/customKriging.png"]);
clf;
drawScatterNico(table2array(X(:, 2:3)), table2array(Yvals(:, 3)), "Co-Kriging");
print(gcf,'-dpng',["./" + folder + "/customCoKriging.png"]);

clf;
drawBinaryPerformanceNico(table2array(X(:, 2:3)), table2array(Yvals(:, 2)) >= 0.5, "Kriging");
print(gcf,'-dpng',["./" + folder + "/customKrigingBinary.png"]);
clf;
drawBinaryPerformanceNico(table2array(X(:, 2:3)), table2array(Yvals(:, 3)) >= 0.5, "Co-Kriging");
print(gcf,'-dpng',["./" + folder + "/customCoKrigingBinary.png"]);

clf;
drawBinaryPerformanceNico(table2array(X(:, 2:3)), table2array(YbinSVM(:, 2)) >= 0.5, "Kriging");
print(gcf,'-dpng',["./" + folder + "/customKrigingBinarySVM.png"]);
clf;
drawBinaryPerformanceNico(table2array(X(:, 2:3)), table2array(YbinSVM(:, 3)) >= 0.5, "Co-Kriging");
print(gcf,'-dpng',["./" + folder + "/customCoKrigingBinarySVM.png"]);






% Plot critical instances
labels = readtable("../../data/isaMetadata/instanceFiltering/instanceFilteringSampleFeaturesCorrTol0.001eps0.3.txt");
labels = table2array(labels);
% First remove the rows of instances which were precluded
labels = labels(strcmp(labels(:,1), 'FALSE'), :);
% Now find visa and dissimlar instances
Xvisa = X(strcmp(labels(:,2), 'TRUE'), :);
Xdiss = X(strcmp(labels(:,3), 'TRUE'), :);
Xvisa = [Xvisa array2table(repelem({'ViSA'}, size(Xvisa,1),1))];
Xdiss = [Xdiss array2table(repelem({'Dissimilar'}, size(Xdiss,1),1))];
Xbar = [Xdiss; Xvisa];
S = categorical(table2array(Xbar(:,4)));
clf;
drawCritical(table2array(Xbar(:, 2:3)), S);
print(gcf,'-dpng',["./" + folder + "/customCriticalInstances.png"]);












function handle = drawCritical(Z, S)
    ubound = ceil(max(Z));
    lbound = floor(min(Z));
    sourcelabels = cellstr(unique(S));
    sourcelabels = sourcelabels([2, 1], :);
    nsources = length(sourcelabels);
    clrs = [[1.0 0.6471 0.0]; [0.0 0.0 1.0]];
    handle = zeros(nsources,1);
    for i=nsources:-1:1
    line(Z(S==sourcelabels{i},1), ...
         Z(S==sourcelabels{i},2), ...
         'LineStyle', 'none', ...
         'Marker', '.', ...
         'Color', clrs(i,:), ...
         'MarkerFaceColor', clrs(i,:), ...
         'MarkerSize', 4);
    handle(i) = patch([0 0],[0 0], clrs(i,:), 'EdgeColor','none');
    end
    xlabel('z_{1}'); ylabel('z_{2}'); title('Instance type');
    legend(handle, sourcelabels, 'Location', 'SouthEast');
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end


function handle = drawScatterNico(Z, X, titlelabel)
    X = (X-min(X))./range(X);
    ubound = ceil(max(Z));
    lbound = floor(min(Z));
    handle = scatter(Z(:,1), Z(:,2), 6, X, 'filled');
    caxis([0,1])
    xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);
    colorbar('EastOutside');
end

function handle = drawSourcesNico(Z, S)

ubound = ceil(max(Z));
lbound = floor(min(Z));
sourcelabels = cellstr(unique(S));
nsources = length(sourcelabels);
%clrs = flipud(lines(nsources));
clrs = flipud(parula(nsources));
handle = zeros(nsources,1);
for i=1:nsources
    line(Z(S==sourcelabels{i},1), ...
         Z(S==sourcelabels{i},2), ...
         'LineStyle', 'none', ...
         'Marker', '.', ...
         'Color', clrs(i,:), ...
         'MarkerFaceColor', clrs(i,:), ...
         'MarkerSize', 6);
    handle(i) = patch([0 0],[0 0], clrs(i,:), 'EdgeColor','none');
end
xlabel('z_{1}'); ylabel('z_{2}'); title('Sources');
legend(handle, sourcelabels, 'Location', 'SouthEast');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end

function h = drawBinaryPerformanceNico(Z, Ybin, titlelabel)
ubound = ceil(max(Z));
lbound = floor(min(Z));
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if any(~Ybin)
    h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
end
if any(Ybin)
    h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
end
for i=1:size(Z,1)
    if Ybin(i)
        line(Z(i,1), Z(i,2), 'LineStyle', 'none', ...
                               'Marker', '.', ...
                               'Color', blue, ...
                               'MarkerFaceColor', blue, ...
                               'MarkerSize', 6);
    else
        line(Z(i,1), Z(i,2), 'LineStyle', 'none', ...
                                 'Marker', '.', ...
                                 'Color', orange, ...
                                 'MarkerFaceColor', orange, ...
                                 'MarkerSize', 6);
    end
end

xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
legend(h(h~=0), lbls(h~=0), 'Location', 'SouthEast');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end


function drawPortfolioSelectionsCustom(Z, P, algolabels, titlelabel)

ubound = ceil(max(Z));
lbound = floor(min(Z));
nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
h = zeros(1,nalgos+1);
isworthy = sum(bsxfun(@eq, P, 0:nalgos))~=0;
clr = flipud(lines(nalgos+1));

for i=1:size(Z,1)
    line(Z(i,1), Z(i,2), 'LineStyle', 'none', ...
                               'Marker', '.', ...
                               'Color', clr(P(i)+1,:), ...
                               'MarkerFaceColor', clr(P(i)+1,:), ...
                               'MarkerSize', 6);
end



for i=0:nalgos
    if ~isworthy(i+1)
        continue;
    end
%     line(Z(P==i,1), Z(P==i,2), 'LineStyle', 'none', ...
%                                'Marker', '.', ...
%                                'Color', clr(i+1,:), ...
%                                'MarkerFaceColor', clr(i+1,:), ...
%                                'MarkerSize', 6);
    h(i+1) = patch([0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
    if i==0
        algolbls{i+1} = 'None';
    else
        algolbls{i+1} = strrep(algolabels{i},'_',' ');
    end
end
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
legend(h(isworthy), algolbls(isworthy), 'Location', 'SouthEast');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis([lbound(1)-1 ubound(1)+1 lbound(2)-1 ubound(2)+1]);

end

