
sigma=1;
fileFolder=fullfile('D:\Share\0401_gold163_all_soma_sort_strict_swc_only_6_29');
dirOutput=dir(fullfile(fileFolder,'*'));
LengthFileFolder=length(dirOutput);
tempFileFolder=[];
for m=3:LengthFileFolder
    tempFileFolder=cat(2,tempFileFolder,str2num(dirOutput(m).name));
end
[NoFileFolder in]=sort(tempFileFolder);

for m=1:length(NoFileFolder)
    m
    num2str(NoFileFolder(m))
    goldFiles = dir(fullfile(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\*.swc')));
    swc1=struct([]);
    for r=1:length(goldFiles)
        if ~isempty(strfind(goldFiles(r).name,'consensus'))
            swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',goldFiles(r).name));
             Constructfile=goldFiles(r).name;
        else
            swc1=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',goldFiles(r).name));
            Goldfile=goldFiles(r).name;
        end
    end
    if isempty(swc1)
        pause;
    end
    data1=[];
    data1=swc1.data;
    
    Files = dir(fullfile(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\processed\*.swc')));
    LengthFiles = length(Files);
    
    for n= mod(length(goldFiles),2):LengthFiles
        if n>=1
            n;
        swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\processed\',Files(n).name));
        Constructfile=Files(n).name;
        end
        if iscell(swc2)
            data2=[];
            for k=2:size(swc2,1)
                data2=cat(1,data2,str2num(swc2{k,1}));
            end
        else
            data2=swc2.data;
        end
        sumn=0;
        if size(data2,1)>10 % If the row number is less than 10, let sumn be 0
            for i=1:size(data1,1)
                i;
                point1=data1(i,3:5);
                if data1(i,7)==-1
                    parent1=data1(find(data1(:,7)==1,1),3:5);
                else
                    parent1=data1(find(data1(:,1)==data1(i,7)),3:5);
                end
                tangent1=(parent1-point1)/norm(parent1-point1);
                tempmat=repmat(point1,size(data2,1),1)-data2(:,3:5);
                [dist,ind]=min(sum(tempmat.*tempmat,2));
                di=sqrt(dist);
                point2=data2(ind,3:5);
                if data2(ind(1),7)==-1
                    parent2=data2(find(data2(:,7)==1,1),3:5);
                else
                    parent2=data2(find(data2(:,1)==data2(ind,7)),3:5);
                end
                tangent2=(parent2-point2)/norm(parent2-point2);
                sumn=sumn+sqrt(abs(sum(tangent1.*tangent2))*exp(-di*di/(2*sigma*sigma)));
            end
        end
        if m==1 && n==0
            fp = fopen('.\nblast_score_0401_6_29.txt','wt');
        else
            fp = fopen('.\nblast_score_0401_6_29.txt','at');
        end
        fprintf(fp, '%d  %s  %s   %5.3f \n', NoFileFolder(m),Constructfile,Goldfile, sumn);
        fclose(fp);
    end
end


