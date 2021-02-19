load carsmall
x1 = Weight;
x2 = Horsepower;    % Contains NaN data
y = MPG;

X = [ones(size(x1)) x1 x2 x1.*x2];
b = regress(y,X)    % Removes NaN data

scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):100:max(x1);
x2fit = min(x2):10:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Weight')
ylabel('Horsepower')
zlabel('MPG')
view(50,10)
hold off

%% projection
% https://www.mathworks.com/help/stats/fitting-an-orthogonal-regression-using-principal-components-analysis.html
% fitting a plane to 3D data
rng(5,'twister');
X = mvnrnd([0 0 0], [1 .2 .7; .2 1 0; .7 0 1],50);
plot3(X(:,1),X(:,2),X(:,3),'bo');
grid on;
maxlim = max(abs(X(:)))*1.1;
axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
axis square
view(-9,12);

[coeff,score,roots] = pca(X);
basis = coeff(:,1:2);
normal = coeff(:,3);

pctExplained = roots' ./ sum(roots);

[n,p] = size(X);
meanX = mean(X,1);
Xfit = repmat(meanX,n,1) + score(:,1:2)*coeff(:,1:2)';
residuals = X - Xfit;

error = abs((X - repmat(meanX,n,1))*normal);
sse = sum(error.^2);

[xgrid,ygrid] = meshgrid(linspace(min(X(:,1)),max(X(:,1)),5), ...
                         linspace(min(X(:,2)),max(X(:,2)),5));
zgrid = (1/normal(3)) .* (meanX*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);

hold on
above = (X-repmat(meanX,n,1))*normal < 0;
below = ~above;
nabove = sum(above);
X1 = [X(above,1) Xfit(above,1) nan*ones(nabove,1)];
X2 = [X(above,2) Xfit(above,2) nan*ones(nabove,1)];
X3 = [X(above,3) Xfit(above,3) nan*ones(nabove,1)];
plot3(X1',X2',X3','-', X(above,1),X(above,2),X(above,3),'o', 'Color',[0 .7 0]);
nbelow = sum(below);
X1 = [X(below,1) Xfit(below,1) nan*ones(nbelow,1)];
X2 = [X(below,2) Xfit(below,2) nan*ones(nbelow,1)];
X3 = [X(below,3) Xfit(below,3) nan*ones(nbelow,1)];
plot3(X1',X2',X3','-', X(below,1),X(below,2),X(below,3),'o', 'Color',[1 0 0]);

hold off
maxlim = max(abs(X(:)))*1.1;
axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
axis square
view(-9,12);


% fitting a line to 3D data
dirVect = coeff(:,1)
Xfit1 = repmat(meanX,n,1) + score(:,1)*coeff(:,1)';
t = [min(score(:,1))-.2, max(score(:,1))+.2];
endpts = [meanX + t(1)*dirVect'; meanX + t(2)*dirVect'];
plot3(endpts(:,1),endpts(:,2),endpts(:,3),'k-');

X1 = [X(:,1) Xfit1(:,1) nan*ones(n,1)];
X2 = [X(:,2) Xfit1(:,2) nan*ones(n,1)];
X3 = [X(:,3) Xfit1(:,3) nan*ones(n,1)];
hold on
plot3(X1',X2',X3','b-', X(:,1),X(:,2),X(:,3),'bo');
hold off
maxlim = max(abs(X(:)))*1.1;
axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
axis square
view(-9,12);
grid on
