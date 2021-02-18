clear;
Length=5.8277997971*6; %length of box
Number=216; %number of particles;
color = ['k' 'r' 'g' 'b' 'y' 'm' 'c' ];
count = 1;
data = zeros(195,12);
pdirs = dir('Np_216_2_p_0.05-1/p_*');
h = figure;
for pidx = 1:4:length(pdirs)
    clear kx ky;
    c=1;r=1;
    n=1:20;
    kx(1,:)=n*((2*pi)/Length);
    [kx,ky] = ndgrid(kx,kx);
    k = [kx(:),ky(:)];
    %[kx,ky,kz] = ndgrid(kx,kx,kx);
    %k=[kx(:),ky(:),kz(:)];

    %table1 = csvread('Np_216/p_0.3/3/positionP.csv'); %The coordinates of particle 1
    table1 = csvread([pdirs(pidx).folder, '/', pdirs(pidx).name, '/5/positionP.csv']);
    table1 = table1(:,2:3);

    %table2 = csvread('Np_216/p_0.3/3/positionPt.csv'); %The coordinates of particle 1
    table2 = csvread([pdirs(pidx).folder, '/', pdirs(pidx).name, '/5/positionPt.csv']); %The coordinates of particle 1
    table2 = table2(:,2:3);

    A=[table1;table2]; 

    positions = A;
    positions=positions';

    product=k*positions;
    exp1=exp(1i*product);

    sum1=sum(exp1,2);
    kmod=abs(sum1);
    sk=(kmod.^2)/Number;

    k1=k.^2;    
    k2=sum(k1,2);
    ksqrt=k2.^(1/2);

    func=[ksqrt,sk];
    func=sortrows(func,1);
    B=tabulate(func(:,1));
    [row,col]=size(func);
    [rowb,colb]=size(B);
    final=zeros(rowb,2);
    while c<row+1
        d=B(r,2);
        final(r,:)=(sum(func(c:(c+d-1),:),1))/d;
        c=c+d;
        r=r+1;
    end
    finals=final';
 
    x(1,:)=finals(1,:);
    y(1,:)=finals(2,:);
    

    figure(h);
    ax = gca;
    set(gcf, 'Position',  [0, 0, 1000, 1000])
    set(ax,'FontSize',50);
    axis square;
    box on;
    set(ax,'linewidth',2);
    hold on;
    %axis([0,6,-10,10])
    
    data(:,count)=x;
    data(:,count+1)=y;
    count = count + 2;
    
    p = plot(x,y,'o--','MarkerFaceColor',color(mod(pidx,7)+1),'MarkerEdgeColor',color(mod(pidx,7)+1),'Markersize',10,'LineWidth',1.5);
    p.Color = color(mod(pidx,7)+1);
    p.DisplayName = pdirs(pidx).name;
    title('Np_216_ordered','Interpreter','none');
    xlabel('k');
    ylabel('S(k)');
    hold on;
    %saveas(gcf,'~/Desktop/Sk.JamT1.png');
end
hold off;
legend('Interpreter','none','FontSize',40);
%saveas(gcf,'Np_216_ordered_sk.png')