clear;
L = 5.8277997971*6;
minD = sqrt(2)*5.8277997971/2;
color = ['k' 'r' 'g' 'y' 'b' 'm' 'c' ];
nparallel = 1;
nvertical = 0;
px = [];
sAve = [];
s=[];
S=[];
count =1;
data = zeros(9,14);

% for ii = 1:2
% pdirs = dir(['Np_216_',num2str(ii),'_p_0.05-1/p_*']);
pdirs = dir(['Np_216_',num2str(2),'_p_0.05-1/p_*']);
h = figure;
o = 1;
calculatedIndex = [];
for pidx = 1:1:length(pdirs)
    if pidx ~=1 && pidx ~= 5 && pidx ~= 9 && pidx ~= 11 && pidx ~= 13 && pidx ~= 17 && pidx ~= 21
        continue
    end
    
%     for ppidx = 1:5
        
    table1 = csvread([pdirs(pidx).folder, '/', pdirs(pidx).name, '/',num2str(5),'/positionP.csv']);
   
    table1 = table1(:,2:4);

    table2 = csvread([pdirs(pidx).folder, '/', pdirs(pidx).name, '/',num2str(5),'/positionPt.csv']);
    
    table2 = table2(:,2:3);
    vector =  zeros(length(table1),2);
    vectorPostion =  zeros(length(table1),2);
    for i = 1:length(table1)
%         if ismember(i,calculatedIndex)
%             continue
%         end
        vector(i,1)=table1(i,1) - table1(table1(i,3),1);
%         calculatedIndex = [calculatedIndex table1(i,3)];
        vectorPostion(i,1) = (table1(i,1) + table1(table1(i,3),1))/2;

        if abs(vector(i,1)) > L/2
            vectorPostion(i,1) = (L+ table1(i,1) + table1(table1(i,3),1))/2;
            if vector(i,1) > 0
                vector(i,1) = vector(i,1)-L;
            else
                vector(i,1) = vector(i,1)+L;
            end
        end

        
        vector(i,2)=table1(i,2)  - table1(table1(i,3),2);
        vectorPostion(i,2) = (table1(i,2) + table1(table1(i,3),2))/2;
        if abs(vector(i,2)) > L/2
            vectorPostion(i,2) = (L+table1(i,2) + table1(table1(i,3),2))/2;
            if vector(i,2) > 0
                vector(i,2) = vector(i,2)-L;
            else
                vector(i,2) = vector(i,2)+L;
            end
        end 
    end
    k = 1;
%     v = nonzeros(vector');
%     newVector = reshape(v,2,72)';
    innerProduct = zeros(length(table1)*(length(table1)-1)/2 , 1);
    distance = zeros(length(table1)*(length(table1)-1)/2 , 1);
    for i = 1:length(table1)
        for j = i+1 :length(table1)
             tempx=abs(vectorPostion(i,1)-vectorPostion(j,1)); 
             tempy=abs(vectorPostion(i,2)-vectorPostion(j,2));
%              if tempx>=L/2            %For boundary conditions
%                  tempx=L-tempx;
%               end
%               if tempy>=L/2
%                  tempy=L-tempy;
%               end
             distance(k) = sqrt(tempx*tempx+tempy*tempy);
             modi = sqrt(vector(i,1)^2+vector(i,2)^2);
             modj = sqrt(vector(j,1)^2+vector(j,2)^2);
            innerProduct(k) = (vector(i,1)*vector(j,1)+vector(i,2)*vector(j,2))/(modi*modj);
            innerProduct(k) = abs(innerProduct(k));
            k = k+1;
        end
    end
    
%     for j = 2:length(newVector)
%         modi = sqrt(newVector(1,1)^2+newVector(1,2)^2);
%         modj = sqrt(newVector(j,1)^2+newVector(j,2)^2);
%         tempInner = (newVector(1,1)*newVector(j,1)+newVector(1,2)*newVector(j,2))/(modi*modj);
%         tempInner = abs(tempInner);     
%         if tempInner >0.5
%            nparallel = nparallel+1;
%         else
%            nvertical = nvertical+1;
%         end
%     end
    
%     s = [s nparallel-nvertical];
%     if ppidx == 5
%         sAve = [sAve mean(s)];
%         s = [];
%   
%     end
%     nparallel = 1;
%     nvertical = 0;
%     end
%     px(o) = (pidx-1)*0.05;
% 
%   
%     o = o+1;

    m=ceil(((max(distance))-minD)/minD)+1;
    innerProduct_R=zeros(m,1); %The number of particles at R
    x=zeros(1,m);

    g2_R=zeros(1,m); %Create a matrix to save the value of g2

    s=1;
    number = 0;
    for R=minD:minD:max(distance) 
       delta_R=0.1*R;
       for i=1:length(distance)
          if distance(i)<=(R+(delta_R/2)) && distance(i)>=(R-(delta_R/2))
             innerProduct_R(s) = innerProduct_R(s) + innerProduct(i);
             number = number +1;
          end
       end
       if number ~= 0
         g2_R(s)=innerProduct_R(s)/number;   %Calculate g2
       end
       x(s)=R;
       s=s+1;
       number = 0;
    end
    
  data(:,count)=x(1:9);
    data(:,count+1)=g2_R(1:9);
    count = count + 2;
    
    figure(h)
     ax = gca;
     set(gcf, 'Position',  [0, 0, 1000, 1000])
     set(ax,'FontSize',50);
     axis square;
     box on;
     set(ax,'linewidth',2);
    axis([0,40,0,1.2])
    p =plot(x(1:9),g2_R(1:9),'o--','MarkerFaceColor',color(mod(pidx,7)+1),'MarkerEdgeColor',color(mod(pidx,7)+1),'Markersize',10,'LineWidth',1.5);
     
    p.Color = color(mod(pidx,7)+1);
    p.DisplayName = pdirs(pidx).name;
    title('Np_216_random','Interpreter','none')
    xlabel('R');
    ylabel('C(R)');
    hold on;

end
% S(ii,:) = sAve;
% sAve = [];
% end
% S = S/72;
%figure(h);
% ax = gca;
%set(gcf, 'Position',  [0, 0, 500, 500])
% set(ax,'FontSize',30);
% axis square;
% box on;
% set(ax,'linewidth',2);
% hold on;
% p =plot(px,S(1,:),'o--k','MarkerFaceColor','k','MarkerEdgeColor','k','Markersize',10,'LineWidth',1.5);
% p =plot(px,S(2,:),'o--r','MarkerFaceColor','r','MarkerEdgeColor','r','Markersize',10,'LineWidth',1.5);
%p.Color = 'r';
%title('Np_216_ordered','Interpreter','none')
% xlabel('p');
% ylabel('S','Rotation',360,'Position',[-0.25,0.4]);
hold off;
% legend({'P5-R456','P5-C46'},'FontSize',20,'TextColor','black')
legend('Interpreter','none','FontSize',30);
%saveas(gcf,'Np_216_random_Cr.png')
