
sigma=1;
fileFolder=fullfile('D:\Yangjian\Deskdop\BigNeuron\gold_163_all_swcs_only');
dirOutput=dir(fullfile(fileFolder,'*'));
LengthFileFolder=length(dirOutput);
tempFileFolder=[];
for m=4:LengthFileFolder
    tempFileFolder=cat(2,tempFileFolder,str2num(dirOutput(m).name));
end
[NoFileFolder in]=sort(tempFileFolder);

for m=1:length(NoFileFolder)
    m
    num2str(NoFileFolder(m))
    Files = dir(fullfile(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\*.swc')));
    LengthFiles = length(Files);
    swc1=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',Files(1).name));
    data1=swc1.data;
    for n= 2:LengthFiles
        n;
        if  ~isempty(strfind(Files(n).name,'app1.swc'))
            swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',Files(n).name),' ',17);
        elseif ~isempty(strfind(Files(n).name,'app2new1.swc')) ||~isempty(strfind(Files(n).name,'app2new2.swc')) ...
                || ~isempty(strfind(Files(n).name,'smartTracing.swc'))
            swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',Files(n).name), ' ', 25);
        elseif ~isempty(strfind(Files(n).name,'LCMboost_updated.swc')) ||~isempty(strfind(Files(n).name,'LCMboost.swc'))||~isempty(strfind(Files(n).name,'app2.swc'))...
                ||~isempty(strfind(Files(n).name,'LCMboost_3.swc')) ||~isempty(strfind(Files(n).name,'app2new3.swc'))
            swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',Files(n).name), ' ', 24);
        elseif ~isempty(strfind(Files(n).name,'nctuTW.swc'))
            swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',Files(n).name), ' ', 2);
        else
            swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',Files(n).name));
        end
        if ~isstruct(swc2)
            data2=swc2;
        else
            data2=swc2.data;
        end
        if data2(1,1)~=1 && ~isempty(strfind(Files(n).name,'smartTracing.swc'))
            swc2=importdata(strcat(fileFolder,'\',num2str(NoFileFolder(m)),'\',Files(n).name), ' ', 24);
            data2=swc2.data;
        end
        
        if ~isempty(strfind(Files(n).name,'Rivulet.swc'))
            %             disp(['For image ', num2str(NoFileFolder(m)), ', ', Files(n).name, ', the number of node with parent -1 is ' num2str(length(find(data2(:,7)==-1)))]);
            %             disp(['For image ', num2str(NoFileFolder(m)), ', ', Files(n).name, ', the number of node with parent -2 is ' num2str(length(find(data2(:,7)==-2)))]);
            data2(data2(:,7)==-2,7)=-1;
        end
        if data2(1,1)~=1 && (isempty(strfind(Files(n).name,'simple.swc'))  && isempty(strfind(Files(n).name,'autotrace.swc'))  && isempty(strfind(Files(n).name,'MST_Tracing.swc')))
            disp(['For image ', num2str(NoFileFolder(m)), ', ', num2str(m), ', ', Files(n).name, ' there is an error.']);
            pause;
        end
        sumn=0;
        if size(data2,1)>10
            
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
                [dist,ind]=sort(sum(tempmat.*tempmat,2));
                q=1;
                while 1
                    di=sqrt(dist(q));
                    point2=data2(ind(q),3:5);
                    if data2(ind(q),7)==-1
                        parent2=data2(find(data2(:,7)==data2(ind(q),1),1),3:5);
                        if ~isempty(parent2)
                            break;
                        else
                            q=q+1;
                        end
                    else
                        temparent2=data2(find(data2(:,1)==data2(ind(q),7)),[3:5 7]);
                        while ~isempty(temparent2) && norm(temparent2(1:3)-point2)==0
                            temparent2=data2(find(data2(:,1)==temparent2(4)),[3:5 7]);
                        end
                        if ~isempty(temparent2)
                            parent2=temparent2(1:3);
                            break;
                        else
                            q=q+1;
                        end
                        
                    end
                end
                
                tangent2=(parent2-point2)/norm(parent2-point2);
                sumn=sumn+sqrt(abs(sum(tangent1.*tangent2))*exp(-di*di/(2*sigma*sigma)));
            end
        end
        if m==1 && n==2
            fp = fopen('.\nblast_nstrict_score.txt','wt');
        else
            fp = fopen('.\nblast_nstrict_score.txt','at');
        end
        fprintf(fp, '%d  %s  %s   %5.3f \n', NoFileFolder(m),Files(n).name,Files(1).name, sumn);
        fclose(fp);
    end
end


