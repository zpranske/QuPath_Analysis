%%%%%% PIPECAT TM %%%%%%
% Zachary Pranske
% Paradis Lab
% Rev. 11/2022  

%% Combine files
defaultpath = 'D:\RNASCOPE\output';
uiwait(msgbox('Open analysis folder containing output files ending with "Detections.txt"'))
d = dir(uigetdir('D:\RNASCOPE\output\plxnb1 pvalb'));
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
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_Stain2_Area")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_Stain2_Area");
    end
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_Stain2_MeanChannelIntensity")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_Stain2_MeanChannelIntensity");
    end
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_Stain3_Area")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_Stain3_Area");
    end
    if max(temp_T.Properties.VariableNames=="SubcellularCluster_Stain3_MeanChannelIntensity")==0
        temp_T = [temp_T table(string(repmat(0,height(temp_T),1)))];
        temp_T = renamevars(temp_T,'Var1',"SubcellularCluster_Stain3_MeanChannelIntensity");
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
     T_allcells(i,:).Qual_Score2 = getscore(T_allcells(i,:).Subcellular_Stain2_NumSpotsEstimated);
     T_allcells(i,:).Qual_Score3 = getscore(T_allcells(i,:).Subcellular_Stain3_NumSpotsEstimated);
 end
 
%  T_allcells = [T_allcells table(repmat(0,height(T_allcells),1))];
%  T_allcells = renamevars(T_allcells,'Var1','Region');
%  for i=i:height(T_allcells)
%      T_allcells(i,:).Region = getregion(T_allcells(i,:).Image{1});
%  end
 
 T_allcells.Region = getregion(T_allcells(:,:).Image{:})
 T_celltypemarker = T_allcells((T_allcells.Subcellular_Stain2_NumClusters>0 & T_allcells.Nucleus_Stain2ODSum./T_allcells.Nucleus_Area>1),:);
 T_semaplexin = T_allcells((T_allcells.Subcellular_Stain3_NumClusters>0 | T_allcells.Subcellular_Stain3_NumSingleSpots>4),:);
 T_highsemaplexin = T_allcells((T_allcells.Subcellular_Stain3_NumClusters>0 & T_allcells.Nucleus_Stain3ODSum./T_allcells.Nucleus_Area>2),:);

 T_coloc = T_celltypemarker(T_celltypemarker.Subcellular_Stain3_NumClusters>0 ,:);
 T_coloc2 = T_highsemaplexin((T_highsemaplexin.Subcellular_Stain2_NumClusters>0 & T_highsemaplexin.Nucleus_Stain2ODSum./T_highsemaplexin.Nucleus_Area>1) ,:);


 n_cells = height(T_allcells)
 n_celltypemarker = height(T_celltypemarker)
 n_semaplexin = height(T_semaplexin)
 n_coloc = height(T_coloc)
 
 n_highsemaplexin = height(T_highsemaplexin)
 n_coloc2 = height(T_coloc2)
 
 mean_od_semaplexin = mean(T_allcells.Nucleus_Stain3ODMean)
 mean_od_semaplexin_celltypemarker = mean(T_celltypemarker.Nucleus_Stain3ODMean)
 
 mean_puncta_semaplexin = mean(T_allcells.Subcellular_Stain2_NumSpotsEstimated)
 mean_puncta_semaplexin_celltypemarker = mean(T_celltypemarker.Subcellular_Stain2_NumSpotsEstimated)
 
 mean_qualscore_semaplexin = mean(T_allcells.Qual_Score3)
 mean_puncta_semaplexin_celltypemarker = mean(T_celltypemarker.Qual_Score3)
   
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
     if est_spots < 1
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

    
 
