%%%%%% QUPATH PIPELINE %%%%%%
% Zachary Pranske
% Paradis Lab
<<<<<<< HEAD
% Rev. 06/07/2023 

%% Description
%
=======
% Rev. 11/2022  
>>>>>>> 024c743ab0dabfa782115650a250f4909c39a1db

%% Combine files
defaultpath = 'C:\Users\Zachary_Pranske\Desktop\Datasets\2023-04-21 analysis rnascope\output';
uiwait(msgbox('Open analysis folder containing output files ending with "Detections.txt"'))
d = dir(uigetdir(defaultpath));
% while(~max(ismember({d(1:end).name},"Detections")))
%     error('Folder must be parent folder containing "measure" and "summary" folders')
%     d = dir(uigetdir(defaultpath));
% end; 

% %DELETE THIS LINE LATER
%  d = dir('C:\Users\Zachary_Pranske\Desktop\pipeline_inputs');

 disp('Found file folders...')
 
 d_measure = d(3:end);
 disp(['Combining measure files for ' d(5).folder '...'])
 warning('off','all') 
 
 T = [];
 for i = 1:length(d_measure)
    temp_T = readtable([d_measure(i).folder '\' d_measure(i).name]);
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_IN_Area")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_IN_Area");
    end
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_IN_MeanChannelIntensity")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_IN_MeanChannelIntensity");
    end
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_SemaPlexin_Area")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_SemaPlexin_Area");
    end
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_SemaPlexin_MeanChannelIntensity")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_SemaPlexin_MeanChannelIntensity");
    end
    %if~(temp_M.Properties.VariableNames == measure_format)
     %   error(['Header format for file ' measure_files(i).name ' is not correct!'])
    %end;
    T = [T; temp_T];
 end 

 %ONE WAY of identifying PV+ cells (but both methods below converge on
    %exact same answer!!)
 %T_PVpos = T((T.Subcellular_Pvalb_NumClusters>0 & T.Nucleus_PvalbODSum>100),:);
 
%STAIN 2 = CELL TYPE MARKER, RED    STAIN 3 = SEMA/PLEXIN, BLUE
 T_allcells = T(T.ROI=="Polygon",:);
 
 T_allcells = [T_allcells table(repmat(0,height(T_allcells),1))];
 T_allcells = renamevars(T_allcells,'Var1','Qual_Score2');
 T_allcells = [T_allcells table(repmat(0,height(T_allcells),1))];
 T_allcells = renamevars(T_allcells,'Var1','Qual_Score3');

 for i=1:height(T_allcells)
     T_allcells(i,:).Qual_Score2 = getscore(T_allcells(i,:).Subcellular_IN_NumSpotsEstimated);
     T_allcells(i,:).Qual_Score3 = getscore(T_allcells(i,:).Subcellular_SemaPlexin_NumSpotsEstimated);
 end
 
%  T_allcells = [T_allcells table(repmat(0,height(T_allcells),1))];
%  T_allcells = renamevars(T_allcells,'Var1','Region');
%  for i=i:height(T_allcells)
%      T_allcells(i,:).Region = getregion(T_allcells(i,:).Image{1});
%  end
 
 %T_allcells.Region = getregion(T_allcells(:,:).Image{:})
 
 %To count as IN+ a cell must have at least 1 cluster (not just spot) and a
 %sufficiently high OD sum to nuclear area ratio
 T_celltypemarker_pos = T_allcells((T_allcells.Subcellular_IN_NumClusters>0 & T_allcells.Nucleus_INODSum./T_allcells.Nucleus_Area>1),:);
 T_celltypemarker_neg = T_allcells((T_allcells.Subcellular_IN_NumClusters == 0 | T_allcells.Nucleus_INODSum./T_allcells.Nucleus_Area<=1),:);
 T_semaplexin = T_allcells((T_allcells.Subcellular_SemaPlexin_NumClusters>0 | T_allcells.Subcellular_SemaPlexin_NumSingleSpots>4),:);
 T_highsemaplexin = T_allcells((T_allcells.Subcellular_SemaPlexin_NumClusters>0 & T_allcells.Nucleus_SemaPlexinODSum./T_allcells.Nucleus_Area>2),:);

 % Table of all IN+ cells with at least one Sema/Plexin cluster
 T_semaplexin_celltypemarker_pos = T_celltypemarker_pos(T_celltypemarker_pos.Subcellular_SemaPlexin_NumSpotsEstimated>=4,:);
 T_semaplexin_celltypemarker_neg = T_celltypemarker_neg(T_celltypemarker_neg.Subcellular_SemaPlexin_NumSpotsEstimated>=4,:);

  % Table of all "high transcript level" Sema/Plexin cells with that meet
  % IN+ criteria (gives info about what proportion of high transcript cells
  % are INs, i.e. do INs constitute a disproportionate share of high transcript cells?)
 T_coloc2 = T_highsemaplexin((T_highsemaplexin.Subcellular_IN_NumClusters>0 & T_highsemaplexin.Nucleus_INODSum./T_highsemaplexin.Nucleus_Area>1) ,:);

 n_cells = height(T_allcells)
 n_celltypemarker_pos = height(T_celltypemarker_pos)
 n_celltypemarker_neg = height(T_celltypemarker_neg)
 n_semaplexin_celltypemarker_pos = height(T_semaplexin_celltypemarker_pos)
 n_semaplexin_celltypemarker_neg = height(T_semaplexin_celltypemarker_neg)
 
 n_highsemaplexin = height(T_highsemaplexin)
 n_coloc2 = height(T_coloc2)
 
 mean_od_semaplexin = mean(T_allcells.Nucleus_SemaPlexinODMean)
 mean_od_semaplexin_celltypemarker_pos = mean(T_celltypemarker_pos.Nucleus_SemaPlexinODMean)
 mean_od_semaplexin_celltypemarker_neg = mean(T_celltypemarker_neg.Nucleus_SemaPlexinODMean)

 mean_puncta_semaplexin = mean(T_allcells.Subcellular_SemaPlexin_NumSpotsEstimated)
 mean_puncta_semaplexin_celltypemarker_pos = mean(T_celltypemarker_pos.Subcellular_SemaPlexin_NumSpotsEstimated)
 mean_puncta_semaplexin_celltypemarker_neg = mean(T_celltypemarker_neg.Subcellular_SemaPlexin_NumSpotsEstimated)
 
 mean_qualscore_semaplexin = mean(T_allcells.Qual_Score3)
 mean_qualscore_semaplexin_celltypemarker_pos = mean(T_celltypemarker_pos.Qual_Score3)
 mean_qualscore_semaplexin_celltypemarker_neg = mean(T_celltypemarker_neg.Qual_Score3)
 
 sd_qualscore_semaplexin_celltypemarker_pos = std(T_celltypemarker_pos.Qual_Score3)
 sd_qualscore_semaplexin_celltypemarker_neg = std(T_celltypemarker_neg.Qual_Score3)
 
 [x,p] = ttest2(T_celltypemarker_pos.Qual_Score3,T_celltypemarker_neg.Qual_Score3)
 
 vacell = {n_cells n_celltypemarker_pos n_celltypemarker_neg n_semaplexin_celltypemarker_pos n_semaplexin_celltypemarker_neg ...
     n_highsemaplexin n_coloc2 mean_od_semaplexin_celltypemarker_pos mean_od_semaplexin_celltypemarker_neg ...
     mean_puncta_semaplexin_celltypemarker_pos mean_puncta_semaplexin_celltypemarker_neg ...
     mean_qualscore_semaplexin_celltypemarker_pos mean_qualscore_semaplexin_celltypemarker_neg ...
     sd_qualscore_semaplexin_celltypemarker_pos sd_qualscore_semaplexin_celltypemarker_neg}

%% GRAPHING

 %plot(T_allcells.Qual_Score2, T_allcells.Qual_Score3, '.k')
 %xlim([0 max(T_allcells.Qual_Score2)]); ylim([0 max(T_allcells.Qual_Score3)])
 %lsline
 
%% STATS

% runstats2(T, "NormCount", "Condition", ["GFP" "Sema4D"]);
% runstats2(T, "AverageSize", "Condition", ["GFP" "Sema4D"]);
% runstats2(T, "PercentArea", "Condition", ["GFP" "Sema4D"]);
% runstats2(T, "mean_Mean_P", "Condition", ["GFP" "Sema4D"]);

 function [sc] = getscore(est_spots)
     if est_spots < 1 || est_spots == NaN
         sc = 0;
     elseif est_spots >= 1 & est_spots < 4
         sc = 1;
     elseif est_spots >= 4 & est_spots < 10
         sc = 2;
     elseif est_spots >= 10 & est_spots < 15
         sc = 3;
     else
         sc = 4;
    end
 end

    
 
