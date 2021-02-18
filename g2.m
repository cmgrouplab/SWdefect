clear;   
N=216; %Total number of particles
L=5.8277997971*6; %The length of box
A=L*L; %The area of box
Rho=N/A; %Density of number of particles 
color = ['k' 'r' 'g' 'b' 'y' 'm' 'c' ];
count = 1;
data = zeros(165,12);

pdirs = dir('Np_216_1_p_0.05-1/p_*');
h = figure;
for pidx = 1:4:length(pdirs)
 
%Reading data from files
    %rng(4);
    %c = rand(length(pdirs),3);       %随机生成颜色。RGB随机。
    table1 = csvread([pdirs(pidx).folder, '/', pdirs(pidx).name, '/5/positionP.csv']); %The coordinates of particle 1
    %table1 = csvread('Np_216/p_0.3/5/positionP.csv');
    table1 = table1(:,2:3);

    
    table2 = csvread([pdirs(pidx).folder, '/', pdirs(pidx).name, '/5/positionPt.csv']); %The coordinates of particle 1
    %table2 = csvread('Np_216/p_0.3/5/positionPt.csv');
    table2 = table2(:,2:3);

%Create a new table for saving coordinates of all particles
    table=[table1;table2]; 

    N_D=(N*(N-1))/2; %The number of distances (How many different distances)
    D=zeros(N_D,1); %Create a matrix to record distances between any two particles

    k=1;
    for i=1:1:N   %Calculate the distance 
        for j=i+1:1:N
            tempx=abs(table(i,1)-table(j,1)); 
            tempy=abs(table(i,2)-table(j,2));
            if tempx>=L/2            %For boundary conditions
                tempx=L-tempx;
            end
            if tempy>=L/2
                tempy=L-tempy;
            end
            D(k)=sqrt(tempx*tempx+tempy*tempy); %Calculate distance between two particles
            k=k+1;
        end
    end

    m=ceil(((L/2)-1)/0.1)+1;
    N_R=zeros(m,1); %The number of particles at R
    x=zeros(1,m);

    g2_R=zeros(1,m); %Create a matrix to save the value of g2

    s=1;
    for R=1:0.1:(L/2)  
        delta_R=0.1*R;
        for i=1:1:N_D
            if D(i)<=(R+(delta_R/2)) && D(i)>=(R-(delta_R/2))
                N_R(s)=N_R(s)+2;
            end
        end
        g2_R(s)=N_R(s)/(N*Rho*2*pi*R*delta_R);   %Calculate g2
        x(s)=R;
        s=s+1;
    end
    
    data(:,count)=x(1:165);
    data(:,count+1)=g2_R(1:165);
    count = count + 2;
   
    figure(h)
    ax = gca;
    set(gcf, 'Position',  [0, 0, 1000, 1000])
    set(ax,'FontSize',50);
    axis square;
    box on;
    set(ax,'linewidth',2);
    axis([0,18,0,6])
    %p = semilogx(x(1:165),log(g2_R(1:165)),'o-','MarkerFaceColor',color(mod(pidx,7)+1),'MarkerEdgeColor',color(mod(pidx,7)+1),'Markersize',10,'LineWidth',1.5)
    p = plot(x(1:165),g2_R(1:165),'o-','MarkerFaceColor',color(mod(pidx,7)+1),'MarkerEdgeColor',color(mod(pidx,7)+1),'Markersize',10,'LineWidth',1.5);
    p.Color = color(mod(pidx,7)+1);
     p.DisplayName = pdirs(pidx).name;
    title('Np_216_random','Interpreter','none')
    xlabel('R');
    ylabel('g2');
    hold on;
end
hold off;
legend('Interpreter','none','FontSize',50);
%saveas(gcf,'Np_216_random_g2.png')

    
    
    